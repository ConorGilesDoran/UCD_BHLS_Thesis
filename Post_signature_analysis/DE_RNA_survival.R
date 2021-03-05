# Final Year Project 
# Bioinformatic Analysis of Copy Number Alteration Signature Signature Mechanisms in Ovarian Carcinoma
# Conor Giles Doran - UCD - Sept. 20 - Feb. 2021

# POST SIGNATURE ANALYSIS 
# Pairwise Differential Expression Analysis
# Gene set enrichment analysis
# RNA signature extraction
# Survival analysis 

#### PACKAGES ####

libs <- c("dplyr", "data.table", "edgeR", "limma", "DESeq2", "stringr", "ggplot2", 
          "VennDiagram", "RColorBrewer", "tidyr", "RNAseqR","sleuth", "fgsea", 
          "msigdbr", "magrittr", "nortest", "GSVA", "BBmisc", "survival", "survminer", 
          "reshape", "PCAtools")

lapply(libs, library, character.only = TRUE)


#### SETUP #### 

# script directory
SCRIPTDIR <- "G:/Shared drives/BHLS/2020/Conor Giles Doran/R Scripts/Thesis_scripts/Post_signature_analysis"

# data directories - TCGA data and Signature data from sig_extrapolation
DATADIR <- "G:/Shared drives/BHLS/data/TCGA_OV"
SIG_DATADIR <- "G:/Shared drives/BHLS/2020/Conor Giles Doran/Final_output/Signature_output"

# output directory
OUTDIR <- "G:/Shared drives/BHLS/2020/Conor Giles Doran/Final_output/Post_signature_analysis"


#### DATA ####

source(paste0(SCRIPTDIR, "/Post_sig_functions.R"))

# TCGA_OV clinical data
clinical_data.tsv <- paste0(DATADIR,'/clinical/ov_tcga_pancan_clinical.tsv')

# TCGA_OV raw count data
raw_count_data.RData <- paste0(DATADIR, '/RNAseq/ov_raw_counts.RData')

# CNA Signature results generated from sig_extrapolation script
DE_sig_data.rds <- paste0(SIG_DATADIR, '/DE_sig_data.rds')

# Clinical metadata
metadata.csv <- paste0(DATADIR, '/clinical/metadata.csv')

# Data provided from a previous kallisto alignment and 
# sleuth gene expression data (tpm)
counts_tpm_anno.RData <- paste0(DATADIR, "/RNAseq/count_tpm_anno.RData")

# outline design to be used in RNAseR DE functions
metadata_design <- "clinical_stage + histological_grade + max_signature"


#### PREPROCESSING ####

print("Preprocessing Input Data...")

# Read in data files 

if(!exists("clinical_data.tsv")){
  stop("Please specifiy clinical_data.tsv file directory")
} else {
  clinical_data <- read.table(clinical_data.tsv, sep = '\t', header = TRUE)
}

if(!exists("DE_sig_data.rds")){
  stop("Please specifiy signature data results (DE_sig_data.rds) file directory")
} else {
  DE_sig_data <- readRDS(DE_sig_data.rds)
  long_sig_results <- DE_sig_data$long_data
  clinical_sigs <- DE_sig_data$clinical_sigs  # sig results in format for joining with clinical 
}

if(!exists("counts_tpm_anno.RData")){
  stop("Please specifiy RData file directory with raw count data, 
  gene expression (tpm) data and gene annotation table")
} else {
  load(counts_tpm_anno.RData)
  raw_counts <- count_data
  
  # make log2tpm
  log2tpm <- dplyr::mutate(.data = tpm_tb, dplyr::across(where(is.numeric), log2))
  # aggregate ens IDs to get single ext ID for TPM
  # NB that ens IDs are unique, but map to multiple ext IDs
  agg_log2tpm_tb <- RNAseqR::group_agg_two(log2tpm, pattern = "_gene")
  
  # format for SSGSEA later
  agg_log2tpm_df <- as.data.frame(agg_log2tpm_tb) # as dataframe for ssgsea later
  # assign rownames as gene names
  rownames(agg_log2tpm_df) <- agg_log2tpm_df$external_gene_name
  # remove the other two columns, keep all tpm columns
  agg_log2tpm_df <- agg_log2tpm_df[, ! colnames(agg_log2tpm_df) %in% c("external_gene_name", "ensembl_gene_id")]
}


# get matching data across signature, clinical and raw count data
# place into a list of originals and matching

print("Identifying Matching Data...")

all_DE_data <- matching_data(clinical_sigs, clinical_data, raw_counts)

# get max signature tally for 311 matching samples

sig_tally_311 <- all_DE_data$matching_sig_data %>%
  group_by(max_signature) %>%
  tally()

saveRDS(all_DE_data, file = paste0(OUTDIR, "/all_DE_data.rds"))
saveRDS(sig_tally_311, file = paste0(OUTDIR, "/sig_tally_311.rds"))

# heatmap for matching 311 normalized signature by sample matrix

pdf(paste0(OUTDIR, "/311Norm_sig_x_sample.pdf"))
heatmap(as.matrix(all_DE_data$matching_sig_data[,-8]))
dev.off()

# write the joined clinical/signature data to file as 'metadata_csv' for RNAseqR functions 
# RNAseR functions will take the file path as input to locate the metadata

write.csv(all_DE_data$sig_clinical_join, file = metadata.csv, row.names = FALSE)


#### DE PROFILING ####

print("Performing DE analysis...")

dir.create(paste0(OUTDIR, "/combined"))

## DEseq2  

RNAseqR::DESeq2_module(raw_counts, anno_tb = anno_tb, tpm_tb = agg_log2tpm_tb, tag = 'final', 
                       metadata_csv = metadata.csv, metadata_design = metadata_design, output_dir = OUTDIR)

# copy main results list to a combined folder 
file.copy(file.path(paste0(OUTDIR, "/DEseq2/final.DESeq2_res_list.rds")), paste0(OUTDIR, "/combined"))

## edgeR 

RNAseqR::edgeR_module(raw_counts, anno_tb = anno_tb, tpm_tb = agg_log2tpm_tb, tag = 'final', 
                      metadata_csv = metadata.csv, metadata_design = metadata_design, 
                      output_dir = OUTDIR)

file.copy(file.path(paste0(OUTDIR, "/edgeR/final.edgeR_res_list.rds")), paste0(OUTDIR, "/combined"))

## limma 

# VOOM can be set to true or false
RNAseqR::limma_module(raw_counts, anno_tb = anno_tb, tpm_tb = agg_log2tpm_tb, tag = 'final', 
                      metadata_csv = metadata.csv, metadata_design = metadata_design, 
                      output_dir = OUTDIR, run_voom = TRUE)

file.copy(file.path(paste0(OUTDIR, "/limma/final.limma_res_list.rds")), paste0(OUTDIR, "/combined"))


#### JOINING ####

master_list <- RNAseqR::master_parse_join(input_dir = paste0(OUTDIR, '/combined'))

saveRDS(master_list, file = paste0(OUTDIR, '/combined/master_list.rds'))


#### OVERLAP ####

# getting overlap between the 2 and 3 packages used

# found in two
fitwo <- RNAseqR::found_in_two(master_list)

saveRDS(fitwo, file = paste0(OUTDIR, '/combined/fitwo.rds'))

# found in three
fithree <- RNAseqR::found_in_three(master_list)

saveRDS(fithree, file = paste0(OUTDIR, '/combined/fithree.rds'))


#### VENN DIAGRAMS ####

# across 3 packages 

RNAseqR::venn_3(master_list, tag = 'final', output_dir = paste0(OUTDIR, '/combined'))


#### GENE SET ENRICHMENT ####

print("Performing Gene Set Enrichment...")

# perform fgsea using master list and overlapped list objects
# res = DESeq2 master results
# sig res = overlap of edgeR and DESeq2 from master_list

fgsea_list <- lapply(names(master_list[["DESeq2"]]), function(f){
  fgsea_plot(res = master_list[["DESeq2"]][[f]], sig_res = fitwo[["DESeq2-edgeR"]][[f]], msigdb_species = "Homo sapiens", msigdb_cat = "H",
             gene_col = NULL, padj = 0.01, rank_col = NULL, output_dir = OUTDIR, tag = f, contrast = f)})

names(fgsea_list) <- names(master_list[["DESeq2"]])

saveRDS(fgsea_list, file = paste0(OUTDIR, '/fgsea/fgsea_list.rds'))

# create master fgsea result object with all fgsea_list objects combined into one

master_fgsea <- do.call(rbind, fgsea_list)
rownames(master_fgsea) <- NULL

saveRDS(master_fgsea, file = paste0(OUTDIR, '/fgsea/master_fgsea.rds'))


#### FGSEA BREAKDOWN ####

sig_pathway_data <- fgsea_breakdown(master_fgsea, num_sigs = 6)

saveRDS(sig_pathway_data, file = paste0(OUTDIR, '/fgsea/sig_pathway_data.rds'))


#### HALLMARK PLOTS ####

# colour blind friendly option can be set to TRUE
plot_pathway_enrichment(master_fgsea, sig_pathway_data, cb_friendly = FALSE)


#### SSGSEA GENE SETS ####

print("Extracting Unique Gene Sets...")

dir.create(paste0(OUTDIR, "/ssgsea"))

sig_genes <- get_sig_gene_sets(sig_pathway_data, master_fgsea)

# Master gene set with total number of genes from unique gene sets

all_genes <- unique(do.call(c, sig_genes))
agg_log2tpm_mat <- as.matrix(agg_log2tpm_df[all_genes,]) 

save(sig_genes, agg_log2tpm_mat, file = paste0(OUTDIR, '/ssgsea/ssgsea_input.Rdata'))

#### RNA SIGNATURES ####

print("Performing SSGSEA...")

RNA_sigs <- gsva(agg_log2tpm_mat, # master gene set
                 sig_genes,       # unique signature gene sets
                 method='ssgsea',
                 ssgsea.norm=T,
                 verbose = TRUE)

colnames(RNA_sigs) <- str_replace_all(colnames(RNA_sigs), "[[:punct:]]", "-")
rownames(RNA_sigs) <- paste0('r', c(1:length(rownames(RNA_sigs))))


# write to file
write.table(RNA_sigs,
            file = paste0(OUTDIR, "/ssgsea/RNA_sigs.tsv"),
            sep='\t',
            row.names=T,
            col.names=T)

# RNA signature heatmap
pdf(paste0(OUTDIR, "/ssgsea/RNA_sig_heatmap.pdf"))
heatmap(as.matrix(t(RNA_sigs)))
dev.off()

#### PCA/SURVIVAL ####

# join RNA sigs with clinical data

matching_clinical <- all_DE_data$matching_clinical

new_clinical <- join_clinical(matching_clinical, hall_sig_NES = RNA_sigs, geneset_names = paste0("r", c(1:4)))

saveRDS(new_clinical, file = paste0(OUTDIR, "/ssgsea/new_clinical.rds"))

## PCA ##

print("Performing PCA...")

RNA_PCA_list <- run_PCA(new_clinical, hall_sig_NES = RNA_sigs)
saveRDS(RNA_PCA_list, file = paste0(OUTDIR, "/PCA/RNA_PCA_list.rds"))

pdf(paste0(OUTDIR, "/PCA/RNA_sig_pca.pdf"), height = 10, width = 12)
RNA_PCA_list$sig_PCA
RNA_PCA_list$sig_clinical_PCA
RNA_PCA_list$scree_clinical
dev.off()

## SURVIVAL ##

print("Running Survival Analysis...")

surv_object <- Surv(time = new_clinical$OS.time, event = new_clinical$OS)

survival_res <- run_survival(new_clinical, surv_object, cb_friendly = FALSE) 

saveRDS(survival_res, file = paste0(OUTDIR, "/survival/surv_results.rds"))

#### session info ####
writeLines(capture.output(sessionInfo()), paste0(OUTDIR, "/Post_sig_sessionInfo.txt"))


