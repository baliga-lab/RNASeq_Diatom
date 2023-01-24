## Imports kallisto results into DESEQ2 with tximport
import_kallisto2DESeq <- function(){
  library('tximport')
  library('rhdf5')
  library("DESeq2")
  library('DT')
  library('tidyverse')
  library('stringr')

  cat(" - Running kallisto import \n\n")

  ## Get list of abundance files from kallisto results
  cat("\t - Getting list of abundance files from kallisto \n\n")
  my.dir <- "/Volumes/omics4tb2/Global_Search/Pilot_Pass/X204SC21081158-Z01-F003/RNASeq_Analysis/"
  sample_id <- dir(file.path(my.dir))

  # remove folders for samples with issues and other script folders also remove ES1-36A since it seems to have very high number of reads across samples.
 sample_id <- sample_id[!(sample_id %in% c("Rawdata_Readme.pdf"))]


  # list of abundance files.
  files <- file.path(my.dir, sample_id, "results_Kallisto_Acerv_Smic-reefGenomics", "abundance.h5")
  names(files) <- sample_id
  cat("\t\t - Collected abundance files for ", length(files), "files\n\n")


  ## Collect sample metadata and prepare sample table
  cat("\t - Collecting meta-data informatioon for samples and creating sample table\n\n")
  # get metadata info from metadata worksheet and select relevant columns and delete NA rows
  meta_data <- read.delim("/Volumes/omics4tb2/Global_Search/Metadata/Mote_master_sample_name_sheet - Sheet1.csv", header=T, sep=",", stringsAsFactors = F)
  meta_data <- meta_data %>%
    dplyr::filter(Novogene_R. %in% sample_id)

  # remove samples from metadata that were earlier filtered and order
  #meta_data <- meta_data[grep("(S6_P_ST_36B|KS2_33A|ES3_33A|ES6_30A|ES1_36A)", meta_data$sample, invert=T),]
  # meta_data <- meta_data[order(meta_data$sample),]
  # meta_data$Temp <- as.factor(meta_data$Temp)
  # meta_data <- meta_data %>%
  #   mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "Exposed", "ExT")) %>%
  #   mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "Protected", "PrT"))


  # create a metadata column by combining temperature and Reef.site.name
 # meta_data$condition <- as.factor(paste(meta_data$Temp, meta_data$Reef.Site.Name, sep="_"))
  # meta_data <- meta_data %>%
    # mutate(condition = str_replace(condition, "KAUST", "AF")) %>%
    # mutate(condition = str_replace(condition, "Eilat", "ICN")) %>%
    # mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "KAUST", "AF")) %>%
    # mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "Eilat", "ICN")) %>%
    # mutate(Site = str_replace(Site, "KAUST", "AF")) %>%
    # mutate(Site = str_replace(Site, "Eilat", "ICN"))

  meta_data <- dplyr::mutate(meta_data, path = files)


  #rebuild_sample table
  sampleTable <- cbind(sampleName=sample_id, fileName=files, meta_data)
  # sampleTable <- sampleTable %>%
  #   mutate(condition = str_replace(condition, "KAUST", "AF")) %>%
  #   mutate(condition = str_replace(condition, "Eilat", "ICN")) %>%
  #   mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "KAUST", "AF")) %>%
  #   mutate(Reef.Site.Name = str_replace(Reef.Site.Name, "Eilat", "ICN")) %>%
  #   mutate(Site = str_replace(Site, "KAUST", "AF")) %>%
  #   mutate(Site = str_replace(Site, "Eilat", "ICN"))

  print(htmltools::tagList(datatable(sampleTable,rownames = F,escape = F)))
  datatable(sampleTable, caption = "Metadata sample table")


  ## Import kallisto abundances
  cat("\t - Importing count values from collected files\n\n")
  txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
  txi.kallisto.tpm <- tximport(files, type = "kallisto", txOut = TRUE,countsFromAbundance = "scaledTPM")

  print(htmltools::tagList(datatable(txi.kallisto$counts[1:6,1:6],rownames = F,escape = F)))
  datatable(txi.kallisto.tpm$counts[1:6,1:6], caption = "Kallisto count matrix")

  ## load data into DESEQ2
  cat("\t - Loading data into DESEQ2 and filtering for Spis vs Smic\n\n")
  dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
  # Seperate Spis vs Spis
  cat("\t\t - filtering for Spis\n\n")
  dds.Spis <- dds[grep("Spis", rownames(dds)),]
  #dds.Spis$condition <- relevel(dds.Spis$condition, ref = "30_AF")
  dds.Spis <- DESeq(dds.Spis)
  vsd.Spis <- vst(dds.Spis, blind = F)

  cat("\t\t - filtering for Smic\n\n")
  # Seperate Spis vs Spis
  dds.Smic <- dds[grep("Smic", rownames(dds)),]
  #dds.Smic$condition <- relevel(dds.Smic$condition, ref = "30_AF")
  dds.Smic <- DESeq(dds.Smic)
  vsd.Smic <- vst(dds.Smic, blind = F)

  return(list(meta_data=meta_data, dds.Spis = dds.Spis, vsd.Spis= vsd.Spis, dds.Smic = dds.Smic, vsd.Smic= vsd.Smic, sample_table= sampleTable))

}
