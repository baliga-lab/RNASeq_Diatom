# Make tidyverse package imports explicit so we have
# a cleaner package
library(tidyr)
library(dplyr)
library(tibble)
library(purrr)
library(readr)

library(fs)

# Input parameters
#organism1 = 'Past'
#organism2 = 'Smic'
#organisms = paste(organism1, organism2, sep='_')
#cat("ORGANISMS: ", organisms, "\n", sep='')

# derived constants
#KALLISTO_ABUNDANCE_PATH = paste("/results_Kallisto_", organisms, "-reefGenomics/abundance.tsv", sep="")
#BWA_SALMON_QUANT_SF_PATH = paste("/results_bwa_Salmon_", organisms, "-reefGenomics/salmon_quant/quant.sf", sep="")
#BWA_SALMON_REGEXP = paste('*', BWA_SALMON_QUANT_SF_PATH, sep='')

## Function to extract TPM from each file
quant_tpm_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)
  # create and return a data frame like this:
  return(tibble::tibble(data$Name, data$TPM))
}
## Function to extract NumReads from each file
quant_numreads_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)
  # create and return a data frame like this:
  return(tibble::tibble(data$Name, data$NumReads))
}


# Remove the path components from the column names of the tibble
# x: a vector of column names
remove_quant_path = function(x) {
   return (unlist(lapply(x, function(s) {
     if (is.character(s) && endsWith(s, "quant.sf")) {
       return (unlist(strsplit(s, '/'))[1])
     }
     return(s)
   })));
}

extract_salmon_quant_tpms = function(salmon_files, analysis_dir) {
  # get the file list and pipe it into our extractor function
  salmon_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(quant_tpm_extractor, .id = "file.ID") %>%
    tidyr::pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    # strip the input directory part out of the column name
    dplyr::rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    # remove the last path part from the name
    dplyr::rename_with(~ remove_quant_path(.)) %>%
    dplyr::rename(gene_id=`data$Name`)
}

#' Extract number of reads from salmon quant.sf files
#'
#' `extract_salmon_quant_numreads()` extract number of reads
#'
#' @importFrom tidyr %>%
#' @param salmon_files the set of quant.sf files to process
#' @param analysis_dir path to the analysis directory
#' @return data frame with the reads
extract_salmon_quant_numreads = function(salmon_files, analysis_dir) {
  # get the file list and pipe it into our extractor function
  salmon_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(quant_numreads_extractor, .id = "file.ID") %>%
    tidyr::pivot_wider(names_from = `file.ID`, values_from=`data$NumReads`) %>%
    # strip the input directory part out of the column name
    dplyr::rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    # remove the last path part from the name
    dplyr::rename_with(~ remove_quant_path(.)) %>%
    dplyr::rename(gene_id=`data$Name`)
}

write_out_tables <- function(merged_df, outdir, algo, typename, organism1, organism2) {
  organisms = paste(organism1, organism2, sep='_')
  df_org1 <- merged_df %>% dplyr::filter(grepl(organism1, gene_id))
  df_org2 <- merged_df %>% dplyr::filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_Merged.csv", sep='')
  org1_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_", organism2, ".csv", sep='')
  readr::write_csv(merged_df, file=merged_file)
  readr::write_csv(df_org1, file=org1_file)
  readr::write_csv(df_org2, file=org2_file)
}

# default regular expression for extract_salmon_quants()
STAR_SALMON_REGEXP = "*/*salmon_quant/quant.sf"

#' Extract TPMS and number of reads from salmon generated quant.sf files
#'
#' `extract_salmon_quants()` extracts all the salmon_quant.sf files
#' in the given analysis directory and writes 2 output files, 1 for
#' the TPMs and one for the number of reads
#'
#' @param organism1 the first organism of the coral
#' @param organism2 the second organism of the coral
#' @param analysis_dir the analysis directory
#' @param outdir the output directory
#' @param regexp the regular expression to match files. Default
#'     value: "*/*salmon_quant/quant.sf"
#' @export
extract_salmon_quants <- function(organism1, organism2, analysis_dir, outdir,
		                  regexp=STAR_SALMON_REGEXP) {
  # Step 1. Sanitize the analysis directory path, otherwise we
  # get weird effects
  analysis_dir = normalizePath(analysis_dir)
  # input dir path *must* end with the path slash
  if (!endsWith(analysis_dir, '/')) {
     analysis_dir <- paste(analysis_dir, '/', sep='');
  }
  # Step 2. Create the output directory if it does not exist
  dir.create(outdir, showWarnings=FALSE)

  message("\nExtract Salmon calculated quants")
  message(paste("dir: [", analysis_dir, "]", sep=''))
  message(paste("regexp: [", regexp, "]", sep=''))
  salmon_files <- fs::dir_ls(analysis_dir, regexp=regexp, recurse=T)
  salmon_tpm_df <- extract_salmon_quant_tpms(salmon_files, analysis_dir);
  salmon_numreads_df <- extract_salmon_quant_numreads(salmon_files, analysis_dir);
  write_out_tables(salmon_tpm_df, outdir, "STAR_Salmon", 'TPM_matrix', organism1, organism2);
  write_out_tables(salmon_numreads_df, outdir, "STAR_Salmon", 'NumReads_matrix', organism1, organism2);
}

## Function to extract TPM from each file
rsem_quant_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble::tibble(data$gene_id, data$TPM))
}

## Extract Bowtie2/RSEM quants
extract_bowtie2rsem_quants <- function(organism1, organism2, analysis_dir, outdir) {
  organisms = paste(organism1, organism2, sep='_')
  message("Extract Bowtie2/RSEM calculated quants")
  rsem_files <- fs::dir_ls(analysis_dir, glob="*.genes.results", recurse=T)

  # get the file list and pipe it into our extractor function
  rsem_df <- rsem_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(rsem_quant_extractor, .id="file.ID") %>%
    tidyr::pivot_wider(names_from=`file.ID`, values_from=`data$TPM`) %>%
    dplyr::rename_with(~ basename(.x)) %>%
    dplyr::rename_with(~ gsub(".genes.results", "", .x, fixed=TRUE)) %>%
    dplyr::rename(gene_id=`data$gene_id`)

  # filter for Acerv
  rsem_df_org1 <- rsem_df %>%
    dplyr::filter(grepl(organism1, gene_id))

  # filter for Smic
  rsem_df_org2 <- rsem_df %>%
    dplyr::filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  readr::write_csv(rsem_df, file=merged_file)
  readr::write_csv(rsem_df_org1, file=org1_file)
  readr::write_csv(rsem_df_org2, file=org2_file)
}

## Function to extract TPM from each file
kallisto_quant_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble::tibble(data$target_id, data$tpm))
}

## Extract Kallisto results
extract_kallisto_quants <- function(organism1, organism2, analysis_dir, outdir) {
  message("Extract Kallisto calculated quants")
  organisms = paste(organism1, organism2, sep='_')
  kallisto_files <- fs::dir_ls(analysis_dir, glob="*/abundance.tsv", recurse=T)

  # get the file list and pipe it into our extractor function
  kallisto_df <- kallisto_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(kallisto_quant_extractor, .id="file.ID") %>%
    tidyr::pivot_wider(names_from = `file.ID`, values_from=`data$tpm`) %>%
    dplyr::rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    dplyr::rename_with(~ gsub(KALLISTO_ABUNDANCE_PATH, "", .x, fixed=TRUE)) %>%
    dplyr::rename(gene_id=`data$target_id`)

  kallisto_df_org1 <- kallisto_df %>% dplyr::filter(grepl(organism1, gene_id))
  kallisto_df_org2 <- kallisto_df %>% dplyr::filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  readr::write_csv(kallisto_df, file=merged_file)
  readr::write_csv(kallisto_df_org1, file=org1_file)
  readr::write_csv(kallisto_df_org2, file=org2_file)
}

## Function to extract TPM from each file
bwa_quant_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble::tibble(data$Name, data$TPM))
}

## Extract bwa/Salmon results
extract_bwasalmon_results <- function(organism1, organism2, analysis_dir, outdir) {
  message("Extract bwa/Salmon results")
  organisms = paste(organism1, organism2, sep='_')
  bwa_files <- fs::dir_ls(analysis_dir, regexp=BWA_SALMON_REGEXP, recurse=T)

  # get the file list and pipe it into our extractor function
  bwa_df <- bwa_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(bwa_quant_extractor, .id="file.ID") %>%
    tidyr::pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    dplyr::rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    dplyr::rename_with(~ gsub(BWA_SALMON_QUANT_SF_PATH, "", .x, fixed=TRUE)) %>%
    dplyr::rename(gene_id = `data$Name`)

  bwa_df_org1 <- bwa_df %>% dplyr::filter(grepl(organism1, gene_id))
  bwa_df_org2 <- bwa_df %>% dplyr::filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')

  readr::write_csv(bwa_df, file=BWA_DF_MERGED_FILE)
  readr::write_csv(bwa_df_org1, file=BWA_DF_ORG1_FILE)
  readr::write_csv(bwa_df_org2, file=BWA_DF_ORG2_FILE)
}


##
## This is where the command line tool starts
##
#args = commandArgs(trailingOnly=TRUE)
#if (length(args) < 2) {
#  message("usage: Rscript importTPMS.R <indir> <outdir>");
#  q(save="no", status=1);
#}

#analysis_dir = normalizePath(args[1])
# input dir path *must* end with the path slash
#if (!endsWith(analysis_dir, '/')) {
#    analysis_dir <- paste(analysis_dir, '/', sep='');
#}
#out_dir = args[2]

#if (!dir.exists(out_dir)) {
#  dir.create(out_dir)
#}

#extract_salmon_quants(organism1, organism2, analysis_dir, out_dir, STAR_SALMON_REGEXP)

#select <- order(rowMeans(kallisto_df_Acerv)[2:126], decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds)[,c("condition","type")])
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
