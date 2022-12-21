library(tidyverse)
library(fs)
library(purrr)
library(ComplexHeatmap)
library("pheatmap")

# Input parameters
organism1 = 'Past'
#organism1 = 'Acerv'
organism2 = 'Smic'
organisms = paste(organism1, organism2, sep='_')
cat("ORGANISMS: ", organisms, "\n", sep='')
#RNA_SEQ_ANALYSIS_DIR = "Pilot_Fail/X204SC21081158-Z01-F002/RNASeq_Analysis/"
RNA_SEQ_ANALYSIS_DIR = "past_smic_out/"
#RNA_SEQ_ANALYSIS_DIR = "/proj/omics4tb2/Global_Search/Pilot_Pass/X204SC21081158-Z01-F003/RNASeq_Analysis/"
# should I create this ?
#OUTDIR = '/proj/omics4tb2/wwu/GlobalSearch/Pilot_Pass/RNASeq_TPMs'
OUTDIR = 'past_smic_out/RNASeq_TPMs'
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
}

# derived constants
#STAR_SALMON_QUANT_SF_PATH = paste("/results_STAR_Salmon_", organisms, "-reefGenomics/salmon_quant/quant.sf", sep="")
STAR_SALMON_QUANT_SF_PATH = paste("/results_STAR_Salmon/salmon_quant/quant.sf", sep="")
#STAR_SALMON_REGEXP = paste('*', STAR_SALMON_QUANT_SF_PATH, sep='')
STAR_SALMON_REGEXP = "*/salmon_quant/quant.sf"
KALLISTO_ABUNDANCE_PATH = paste("/results_Kallisto_", organisms, "-reefGenomics/abundance.tsv", sep="")
BWA_SALMON_QUANT_SF_PATH = paste("/results_bwa_Salmon_", organisms, "-reefGenomics/salmon_quant/quant.sf", sep="")
BWA_SALMON_REGEXP = paste('*', BWA_SALMON_QUANT_SF_PATH, sep='')

BWA_DF_MERGED_FILE = paste(OUTDIR, "/bwa_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
BWA_DF_ORG1_FILE = paste(OUTDIR, "/bwa_Salmon_", organisms, "_TPM_matrix_Acerv.csv", sep='')
BWA_DF_ORG2_FILE = paste(OUTDIR, "/bwa_Salmon_", organisms, "_TPM_matrix_Smic.csv", sep='')

SALMON_DF_MERGED_FILE = paste(OUTDIR, "/STAR_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
SALMON_DF_ORG1_FILE = paste(OUTDIR, "/STAR_Salmon_", organisms, "_TPM_matrix_Acerv.csv", sep='')
SALMON_DF_ORG2_FILE = paste(OUTDIR, "/STAR_Salmon_", organisms, "_TPM_matrix_Smic.csv", sep='')

RSEM_DF_MERGED_FILE = paste(OUTDIR, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Merged.csv", sep='')
RSEM_DF_ORG1_FILE = paste(OUTDIR, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Acerv.csv", sep='')
RSEM_DF_ORG2_FILE = paste(OUTDIR, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Smic.csv", sep='')

KALLISTO_DF_MERGED_FILE = paste(OUTDIR, "/Kallisto_", organisms, "_TPM_matrix_Merged.csv", sep='')
KALLISTO_DF_ORG1_FILE = paste(OUTDIR, "/Kallisto_", organisms, "_TPM_matrix_Acerv.csv", sep='')
KALLISTO_DF_ORG2_FILE = paste(OUTDIR, "/Kallisto_", organisms, "_TPM_matrix_Smic.csv", sep='')

## Extract Salmon calculated quants
message("\nExtract Salmon calculated quants")
message(paste("dir: [", RNA_SEQ_ANALYSIS_DIR, "]", sep=''))
message(paste("regexp: [", STAR_SALMON_REGEXP, "]", sep=''))
salmon_files <- fs::dir_ls(RNA_SEQ_ANALYSIS_DIR, regexp=STAR_SALMON_REGEXP, recurse=T)
message("SALMON FILES FOUND: ")
print(salmon_files)
message("\n**** END SALMON FILES FOUND ****\n")

## Function to extract TPM from each file
quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t")

  # create and return a data frame like this:
  return(data_frame(data$Name, data$TPM))
}

# get the file list and pipe it into our extractor function
salmon_df <- salmon_files %>%
  set_names(.) %>%
  map_dfr(quant_extractor, .id = "file.ID") %>%
  pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
  rename_with(~ gsub(RNA_SEQ_ANALYSIS_DIR, "", .x, fixed=TRUE)) %>%
  rename_with(~ gsub(STAR_SALMON_QUANT_SF_PATH, "", .x, fixed=TRUE)) %>%
  rename(gene_id = `data$Name`)

salmon_df_org1 <- salmon_df %>% filter(grepl(organism1, gene_id))
salmon_df_org2 <- salmon_df %>% filter(grepl(organism2, gene_id))

# write into  file
write_csv(salmon_df, file=SALMON_DF_MERGED_FILE)
write_csv(salmon_df_org1, file=SALMON_DF_ORG1_FILE)
write_csv(salmon_df_org2, file=SALMON_DF_ORG2_FILE)

## Extract Bowtie2/RSEM quants
message("Extract Bowtie2/RSEM calculated quants")
rsem_files <- fs::dir_ls(RNA_SEQ_ANALYSIS_DIR, glob="*.genes.results", recurse=T)

## Function to extract TPM from each file
rsem_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t")

  # create and return a data frame like this:
  return(data_frame(data$gene_id, data$TPM))
}

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
write_csv(rsem_df, file=RSEM_DF_MERGED_FILE)
write_csv(rsem_df_org1, file=RSEM_DF_ORG1_FILE)
write_csv(rsem_df_org2, file=RSEM_DF_ORG2_FILE)

## Extract Kallisto results
message("Extract Kallisto calculated quants")
kallisto_files <- fs::dir_ls(RNA_SEQ_ANALYSIS_DIR, glob="*/abundance.tsv", recurse=T)

## Function to extract TPM from each file
kallisto_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t")

  # create and return a data frame like this:
  return(data_frame(data$target_id, data$tpm))
}

# get the file list and pipe it into our extractor function
kallisto_df <- kallisto_files %>%
  set_names(.) %>%
  map_dfr(kallisto_quant_extractor, .id="file.ID") %>%
  pivot_wider(names_from = `file.ID`, values_from=`data$tpm`) %>%
  rename_with(~ gsub(RNA_SEQ_ANALYSIS_DIR, "", .x, fixed=TRUE)) %>%
  rename_with(~ gsub(KALLISTO_ABUNDANCE_PATH, "", .x, fixed=TRUE)) %>%
  rename(gene_id=`data$target_id`)

kallisto_df_org1 <- kallisto_df %>% filter(grepl(organism1, gene_id))
kallisto_df_org2 <- kallisto_df %>% filter(grepl(organism2, gene_id))

# write into  file
write_csv(kallisto_df, file=KALLISTO_DF_MERGED_FILE)
write_csv(kallisto_df_org1, file=KALLISTO_DF_ORG1_FILE)
write_csv(kallisto_df_org2, file=KALLISTO_DF_ORG2_FILE)

## Extract bwa/Salmon results
message("Extract bwa/Salmon results")
bwa_files <- fs::dir_ls(RNA_SEQ_ANALYSIS_DIR, regexp=BWA_SALMON_REGEXP, recurse=T)

## Function to extract TPM from each file
bwa_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t")

  # create and return a data frame like this:
  return(data_frame(data$Name, data$TPM))
}

# get the file list and pipe it into our extractor function
bwa_df <- bwa_files %>%
  set_names(.) %>%
  map_dfr(bwa_quant_extractor, .id="file.ID") %>%
  pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
  rename_with(~ gsub(RNA_SEQ_ANALYSIS_DIR, "", .x, fixed=TRUE)) %>%
  rename_with(~ gsub(BWA_SALMON_QUANT_SF_PATH, "", .x, fixed=TRUE)) %>%
  rename(gene_id = `data$Name`)

bwa_df_org1 <- bwa_df %>% filter(grepl(organism1, gene_id))
bwa_df_org2 <- bwa_df %>% filter(grepl(organism2, gene_id))

# write into  file
write_csv(bwa_df, file=BWA_DF_MERGED_FILE)
write_csv(bwa_df_org1, file=BWA_DF_ORG1_FILE)
write_csv(bwa_df_org2, file=BWA_DF_ORG2_FILE)

select <- order(rowMeans(kallisto_df_Acerv)[2:126], decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
