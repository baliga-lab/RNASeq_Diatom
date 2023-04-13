# Make tidyverse package imports explicit so we have
# a cleaner package
library(tidyr)
library(dplyr)
library(tibble)
library(purrr)
library(readr)

library(fs)

# Input parameters
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

write_out_tables <- function(merged_df, outdir, algo, typename, org_vec) {
  organisms = paste(org_vec, collapse='_')
  # write into  file
  merged_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_Merged.csv", sep='')
  readr::write_csv(merged_df, file=merged_file)

  for (org in org_vec) {
      df_org <- merged_df %>% dplyr::filter(grepl(org, gene_id))
      org_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_", org, ".csv", sep='')
      readr::write_csv(df_org, file=org_file)
  }
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
extract_salmon_quants <- function(organisms, analysis_dir, outdir,
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
  tpm_dir <- paste(outdir, 'TPMs', sep='/')
  count_dir <- paste(outdir, 'Counts', sep='/')
  dir.create(tpm_dir, showWarnings=FALSE)
  dir.create(count_dir, showWarnings=FALSE)

  message("\nExtract Salmon calculated quants")
  message(paste("dir: [", analysis_dir, "]", sep=''))
  message(paste("regexp: [", regexp, "]", sep=''))
  salmon_files <- fs::dir_ls(analysis_dir, regexp=regexp, recurse=T)
  salmon_tpm_df <- extract_salmon_quant_tpms(salmon_files, analysis_dir);
  salmon_numreads_df <- extract_salmon_quant_numreads(salmon_files, analysis_dir);
  write_out_tables(salmon_tpm_df, tpm_dir, "STAR_Salmon", 'TPM_matrix', organisms);
  write_out_tables(salmon_numreads_df, count_dir, "STAR_Salmon", 'NumReads_matrix', organisms);
}

## Function to extract TPM from each file
rsem_quant_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble::tibble(data$gene_id, data$TPM))
}

## Extract Bowtie2/RSEM quants
extract_bowtie2rsem_quants <- function(org_vec, analysis_dir, outdir) {
  organisms = paste(org_vec, collapse='_')
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

  # write into  file
  merged_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Merged.csv", sep='')
  readr::write_csv(rsem_df, file=merged_file)

  for (org in org_vec) {
      # filter for Acerv
      rsem_df_org <- rsem_df %>% dplyr::filter(grepl(org, gene_id))
      org_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", org, ".csv", sep='')
      readr::write_csv(rsem_df_org, file=org_file)
  }
}

## Function to extract TPM from each file
kallisto_quant_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble::tibble(data$target_id, data$tpm))
}

## Extract Kallisto results
extract_kallisto_quants <- function(org_vec, analysis_dir, outdir) {
  message("Extract Kallisto calculated quants")
  organisms = paste(org_vec, collapse='_')
  kallisto_files <- fs::dir_ls(analysis_dir, glob="*/abundance.tsv", recurse=T)

  # get the file list and pipe it into our extractor function
  kallisto_df <- kallisto_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(kallisto_quant_extractor, .id="file.ID") %>%
    tidyr::pivot_wider(names_from = `file.ID`, values_from=`data$tpm`) %>%
    dplyr::rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    dplyr::rename_with(~ gsub(KALLISTO_ABUNDANCE_PATH, "", .x, fixed=TRUE)) %>%
    dplyr::rename(gene_id=`data$target_id`)

  # write into  file
  merged_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_Merged.csv", sep='')
  readr::write_csv(kallisto_df, file=merged_file)

  for (org in org_vec) {
      kallisto_df_org <- kallisto_df %>% dplyr::filter(grepl(org, gene_id))
      org_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", org, ".csv", sep='')
      readr::write_csv(kallisto_df_org, file=org_file)
  }
}

## Function to extract TPM from each file
bwa_quant_extractor = function(fname) {
  data = readr::read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble::tibble(data$Name, data$TPM))
}

## Extract bwa/Salmon results
extract_bwasalmon_results <- function(org_vec, analysis_dir, outdir) {
  message("Extract bwa/Salmon results")
  organisms = paste(org_vec, collapse='_')
  bwa_files <- fs::dir_ls(analysis_dir, regexp=BWA_SALMON_REGEXP, recurse=T)

  # get the file list and pipe it into our extractor function
  bwa_df <- bwa_files %>%
    purrr::set_names(.) %>%
    purrr::map_dfr(bwa_quant_extractor, .id="file.ID") %>%
    tidyr::pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    dplyr::rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    dplyr::rename_with(~ gsub(BWA_SALMON_QUANT_SF_PATH, "", .x, fixed=TRUE)) %>%
    dplyr::rename(gene_id = `data$Name`)

  merged_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
  readr::write_csv(bwa_df, file=merged_file)

  for (org in org_vec) {
      bwa_df_org <- bwa_df %>% dplyr::filter(grepl(org, gene_id))
      org_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", org, ".csv", sep='')
      readr::write_csv(bwa_df_org, file=org_file)
  }
}

#select <- order(rowMeans(kallisto_df_Acerv)[2:126], decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds)[,c("condition","type")])
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
