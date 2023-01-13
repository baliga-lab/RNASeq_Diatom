library(tidyverse)
library(fs)

# Input parameters
organism1 = 'Past'
organism2 = 'Smic'
organisms = paste(organism1, organism2, sep='_')
cat("ORGANISMS: ", organisms, "\n", sep='')

# derived constants
KALLISTO_ABUNDANCE_PATH = paste("/results_Kallisto_", organisms, "-reefGenomics/abundance.tsv", sep="")
BWA_SALMON_QUANT_SF_PATH = paste("/results_bwa_Salmon_", organisms, "-reefGenomics/salmon_quant/quant.sf", sep="")
BWA_SALMON_REGEXP = paste('*', BWA_SALMON_QUANT_SF_PATH, sep='')

## Function to extract TPM from each file
quant_tpm_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)
  # create and return a data frame like this:
  return(tibble(data$Name, data$TPM))
}
## Function to extract NumReads from each file
quant_numreads_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)
  # create and return a data frame like this:
  return(tibble(data$Name, data$NumReads))
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
    set_names(.) %>%
    map_dfr(quant_tpm_extractor, .id = "file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    # strip the input directory part out of the column name
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    # remove the last path part from the name
    rename_with(~ remove_quant_path(.)) %>%
    rename(gene_id=`data$Name`)
}

extract_salmon_quant_numreads = function(salmon_files, analysis_dir) {
  # get the file list and pipe it into our extractor function
  salmon_files %>%
    set_names(.) %>%
    map_dfr(quant_numreads_extractor, .id = "file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$NumReads`) %>%
    # strip the input directory part out of the column name
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    # remove the last path part from the name
    rename_with(~ remove_quant_path(.)) %>%
    rename(gene_id=`data$Name`)
}

write_out_tables <- function(merged_df, outdir, algo, typename, organism1, organism2) {
  organisms = paste(organism1, organism2, sep='_')
  df_org1 <- merged_df %>% filter(grepl(organism1, gene_id))
  df_org2 <- merged_df %>% filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_Merged.csv", sep='')
  org1_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/", algo, "_", organisms, "_", typename, "_", organism2, ".csv", sep='')
  write_csv(merged_df, file=merged_file)
  write_csv(df_org1, file=org1_file)
  write_csv(df_org2, file=org2_file)
}

## Extract Salmon calculated quants
extract_salmon_quants <- function(analysis_dir, outdir, regexp) {
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
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$gene_id, data$TPM))
}

## Extract Bowtie2/RSEM quants
extract_bowtie2rsem_quants <- function(analysis_dir, outdir) {
  message("Extract Bowtie2/RSEM calculated quants")
  rsem_files <- fs::dir_ls(analysis_dir, glob="*.genes.results", recurse=T)

  # get the file list and pipe it into our extractor function
  rsem_df <- rsem_files %>%
    set_names(.) %>%
    map_dfr(rsem_quant_extractor, .id="file.ID") %>%
    pivot_wider(names_from=`file.ID`, values_from=`data$TPM`) %>%
    rename_with(~ basename(.x)) %>%
    rename_with(~ gsub(".genes.results", "", .x, fixed=TRUE)) %>%
    rename(gene_id=`data$gene_id`)

  # filter for Acerv
  rsem_df_org1 <- rsem_df %>%
    filter(grepl(organism1, gene_id))

  # filter for Smic
  rsem_df_org2 <- rsem_df %>%
    filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  write_csv(rsem_df, file=merged_file)
  write_csv(rsem_df_org1, file=org1_file)
  write_csv(rsem_df_org2, file=org2_file)
}

## Function to extract TPM from each file
kallisto_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$target_id, data$tpm))
}

## Extract Kallisto results
extract_kallisto_quants <- function(analysis_dir, outdir) {
  message("Extract Kallisto calculated quants")
  kallisto_files <- fs::dir_ls(analysis_dir, glob="*/abundance.tsv", recurse=T)

  # get the file list and pipe it into our extractor function
  kallisto_df <- kallisto_files %>%
    set_names(.) %>%
    map_dfr(kallisto_quant_extractor, .id="file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$tpm`) %>%
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    rename_with(~ gsub(KALLISTO_ABUNDANCE_PATH, "", .x, fixed=TRUE)) %>%
    rename(gene_id=`data$target_id`)

  kallisto_df_org1 <- kallisto_df %>% filter(grepl(organism1, gene_id))
  kallisto_df_org2 <- kallisto_df %>% filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  write_csv(kallisto_df, file=merged_file)
  write_csv(kallisto_df_org1, file=org1_file)
  write_csv(kallisto_df_org2, file=org2_file)
}

## Function to extract TPM from each file
bwa_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$Name, data$TPM))
}

## Extract bwa/Salmon results
extract_bwasalmon_results <- function(analysis_dir, outdir) {
  message("Extract bwa/Salmon results")
  bwa_files <- fs::dir_ls(analysis_dir, regexp=BWA_SALMON_REGEXP, recurse=T)

  # get the file list and pipe it into our extractor function
  bwa_df <- bwa_files %>%
    set_names(.) %>%
    map_dfr(bwa_quant_extractor, .id="file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    rename_with(~ gsub(BWA_SALMON_QUANT_SF_PATH, "", .x, fixed=TRUE)) %>%
    rename(gene_id = `data$Name`)

  bwa_df_org1 <- bwa_df %>% filter(grepl(organism1, gene_id))
  bwa_df_org2 <- bwa_df %>% filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')

  write_csv(bwa_df, file=BWA_DF_MERGED_FILE)
  write_csv(bwa_df_org1, file=BWA_DF_ORG1_FILE)
  write_csv(bwa_df_org2, file=BWA_DF_ORG2_FILE)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  message("usage: Rscript importTPMS.R <indir> <outdir>");
  q(save="no", status=1);
}

STAR_SALMON_REGEXP = "*/*salmon_quant/quant.sf"
analysis_dir = normalizePath(args[1])
# input dir path *must* end with the path slash
if (!endsWith(analysis_dir, '/')) {
    analysis_dir <- paste(analysis_dir, '/', sep='');
}
out_dir = args[2]

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

extract_salmon_quants(analysis_dir, out_dir, STAR_SALMON_REGEXP)

#select <- order(rowMeans(kallisto_df_Acerv)[2:126], decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds)[,c("condition","type")])
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
