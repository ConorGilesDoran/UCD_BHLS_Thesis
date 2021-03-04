# Final Year Project - Helper functions for post signature analysis
# Bioinformatic Analysis of Copy Number Alteration Signature Signature Mechanisms in Ovarian Carcinoma
# Conor Giles Doran - UCD - Sept. 2020 - Feb. 2021

# Extra functions to be used within the 'DE_RNA_survival.R' script

#### MATCHING DATA ####
# combining normalized CNA signature results data with clinical data for DE analysis
# DE_sig_data - signature results in suitable format for merging with clinical
# clinical_data - clinical data for all samples
# count matrix - raw RNA_seq count matrix

matching_data <- function(DE_sig_data, clinical_data, count_matrix){
  # remove end of ID as clinical data does not have this
  rownames(DE_sig_data) <- substr(rownames(DE_sig_data), start = 1, stop = 12) 
  
  # clinical_data has 587 rows of data, signature table only has 415
  sig_clinical_data <- filter(clinical_data, submitter_id %in% rownames(DE_sig_data))
  
  DE_sig_clinical <- filter(DE_sig_data, rownames(DE_sig_data) %in% sig_clinical_data$submitter_id)
  
  # can join together after creating column for submitter id in the signature table
  DE_sig_clinical <- setDT(DE_sig_clinical, keep.rownames = 'submitter_id')[]
  
  full_clinical_data <- left_join(sig_clinical_data, DE_sig_clinical, by = c("submitter_id" = "submitter_id"))
  
  # have to replace fullstop with '-' in sample names, so filtering can be carried out
  colnames(count_matrix) <- str_replace_all(colnames(count_matrix), "[[:punct:]]", "-")
  
  # have to transpose to match up sample names
  count_matrix <- as.data.frame(t(count_matrix))
  
  # filter to match the samples selected from count matrix and clinical data
  match_count_matrix <- t(filter(count_matrix, rownames(count_matrix) %in% full_clinical_data$submitter_id))
  
  match_clinical <- filter(full_clinical_data, full_clinical_data$submitter_id %in% colnames(match_count_matrix))
  
  match_sig_data <- filter(DE_sig_data, rownames(DE_sig_data) %in%  match_clinical$submitter_id)
  
  full_clinical_data <- dplyr::rename(full_clinical_data, sample = submitter_id) # for RNAseqR package functions
  
  all_DE_data <- list(
    original_clinical = clinical_data,
    original_counts = count_matrix,
    original_sig_data = DE_sig_data,
    sig_clinical_join = full_clinical_data,
    matching_counts = match_count_matrix,
    matching_clinical = match_clinical,
    matching_sig_data = match_sig_data
  )
}


#### FGSEA BREAKDOWN ####
# processing fgsea results to calculate signature specific pathways and pathway tallies per signature
# master_fgsea - rds object containing all fgsea results merged in one large dataframe
# num_sigs - number of CNA signatures being analysed

fgsea_breakdown <- function(master_fgsea, num_sigs){
  
  genes_tally_tb <- as_tibble(master_fgsea$pathway) %>% 
    distinct() %>% 
    dplyr::rename(pathway = value)
  
  # signature tally per hallmark pathway 
  genes_tally <- list()
  for(i in 1:num_sigs){
    genes_tally[["pathway_enrichment"]][[paste0("sig_", i, "_tally")]] <- master_fgsea %>% 
      filter(grepl(paste0("max_signatures", i), contrast)) %>% 
      dplyr::select(-contrast) %>% 
      group_by(pathway) %>%    
      tally()
    
    genes_tally[["total_gene_sets"]][[paste0("signature", i)]] <- master_fgsea %>% 
      filter(grepl(paste0("max_signatures", i), contrast))
    
    gts <- genes_tally[["pathway_enrichment"]][[paste0("sig_", i, "_tally")]]
    sig <- paste0("max_signatures", i)
    genes_tally_tb <- left_join(genes_tally_tb, gts) %>% 
      dplyr::rename(!!quo_name(sig) := n)
    
    enrich_object <- genes_tally[["pathway_enrichment"]][[i]]
    
    genes_tally[["pathway_enrichment"]][[i]] <- enrich_object[rev(order(enrich_object[,2])),]
  }
  
  # assign na values to 0
  genes_tally_tb[is.na(genes_tally_tb)] <- 0
  
  # new tally object with pathway column and associated signature tallies
  sig_tally_mat <- as.matrix(genes_tally_tb[,c(2:length(colnames(genes_tally_tb)))])
  
  rownames(sig_tally_mat) <- unlist(genes_tally_tb[,1])
  
  # identify maximum CNA signature per hallmark pathway
  max_sigs <- colnames(sig_tally_mat)[max.col(sig_tally_mat, ties.method = "first")]
  
  genes_tally_tb$max_sig <- max_sigs
  
  gene_results <- list(
    fgsea_breakdown = genes_tally,
    sig_paths_tally = genes_tally_tb                         
  )
}


#### PLOT PATHWAY ENRICHMENT ####
# generates 3 individual plots for pathway enrichment
# - bar plot with top 30 pathways enriched across all signature comparisons
# - horizontal bar plot with signature enrichment % 
# - dot plot representing maximum signature involvement per pathway
# master_fgsea - rds object containing all fgsea results merged in one large dataframe
# sig_pathway_data - rds object generated from fgsea_breakdown function
# cb_friendly - colour blind friendly option, can be set to TRUE/FALSE

plot_pathway_enrichment <- function(master_fgsea, sig_pathway_data, cb_friendly = FALSE){

  fgsea_tally <- master_fgsea %>%                   
    group_by(pathway, contrast) %>%    
    tally()
  
  fgsea_tally <- fgsea_tally[rev(order(fgsea_tally[,3])),]
  
  # for plotting top 30 pathway - paired color brewer is cbf
  mycol <- c(brewer.pal(12, 'Paired'), "#66CDAA", "#7A7A7A", "#00EEEE")
  names(mycol) <- levels(factor(c(fgsea_tally$contrast)))

  dir.create(paste0(OUTDIR, "/pathway_counts"))
  
  total_tally <- ggplot(fgsea_tally[1:30,], aes(x = pathway, y = n, fill = contrast)) +
      geom_histogram(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = mycol,
                        name = "Contrasts") +                      
      labs(title = "Top 30 Results - Pathway Enrichment Across All Signature Comparisons", 
           x = "\n HMsigDB Pathways", y = "Gene Count \n") + 
      theme_bw() +
      theme(panel.grid = element_blank(), 
            axis.text.y = element_text(size = 10), 
            axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
            axis.title = element_text(size = 12), 
            plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
            plot.margin = unit(c(1,1,1,1), units = , "cm"), 
            legend.title = element_text(face = "bold", size = 8),
            legend.text = element_text(size = 12),
            legend.position = "bottom",
            legend.box.background = element_rect(color = "grey", size = 0.3)) +
      ggsave(paste0(OUTDIR, "/pathway_counts/top30_pathway_enriched.pdf"), height = 7, width = 15)

  # for max sig per hallmark plots
  signature_tally <- sig_pathway_data$sig_paths_tally
  max_sig_tally <- sig_pathway_data$sig_paths_tally[,-c(2:7)]
  max_sig_tally$max_sig <- substr(max_sig_tally$max_sig, start = 15, stop = 15) 
  factor(max_sig_tally$max_sig)
  
  # normalize to get 0:1 scale
  norm_max_tally <- YAPSA:::normalize_df_per_dim(signature_tally[,c(2:7)],1)
  norm_max_tally$pathway <- max_sig_tally$pathway
  norm_max_tally[,c(1:6)] <- norm_max_tally[,c(1:6)] * 100 
  sig_tally_df <- melt(norm_max_tally, id.vars = 'pathway', variable.name = 'signature')
  
  # need to make new dot plot same colour as enrichment plot
  if(cb_friendly == TRUE){
    enriched_col <- c(brewer.pal(6, 'Paired'))
    names(enriched_col) <- unique(sig_tally_df$signature)
    enriched_col2 <- c(brewer.pal(4, 'Paired'))
    names(enriched_col2) <- levels(max_sig_tally$max_sig)
  } else {
    enriched_col <- c(brewer.pal(6, 'Set1'))
    names(enriched_col) <- unique(sig_tally_df$signature)
    enriched_col2 <- c(brewer.pal(4, 'Set1'))
    names(enriched_col2) <- levels(max_sig_tally$max_sig)
  }
 
  # horizontal bar plot
  ggplot(sig_tally_df, aes(fill = signature, x = pathway, y = value)) +
    geom_col(position = "stack", width = .8) + 
    scale_fill_manual(values = enriched_col) +
    coord_flip() +
    labs(title = paste0("Hallmark Pathways - Signature Enrichment"), 
          x = "\n HMsigDB Pathways", y = "Signature Enrichment %\n") + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
          axis.title = element_text(size = 12), 
          plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
          plot.margin = unit(c(1,1,1,1), units = , "cm"), 
          legend.title = element_text(face = "bold", size = 8),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.box.background = element_rect(color = "grey", size = 0.3)) +
    ggsave(paste0(OUTDIR, "/pathway_counts/sig_pathway_enrichment.pdf"), height = 9, width = 10)
  
  # dot plot
  ggplot(max_sig_tally, aes(x = max_sig, y = pathway, colour = max_sig)) +
    geom_point(size = 3) +
    scale_colour_manual(values = enriched_col2) +
    labs(title = 'Max Signature Involvement per Hallmark pathway', 
          x = 'CNA Signature', y = 'HMSigDB Pathways \n') + 
    theme_bw() +
    theme(panel.grid = element_line(), 
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 12), 
          plot.title = element_text(size = 14, hjust = 0.5, face = 'bold'), 
          plot.margin = unit(c(0.5,0.5,0.5,0.5), units = , 'cm')) +
    ggsave(paste0(OUTDIR, "/pathway_counts/max_involvement_dot.pdf"), height = 10, width = 10)
}



#### GET GENE SETS ####

# extract unique-to-signature gene sets 
# master_fgsea - rds object containing all fgsea results merged in one large dataframe
# sig_pathway_data - rds object generated from fgsea_breakdown function
# uses the maximum signature identified per pathway to extract genes from master_fgsea 
# that are associated with that pathway

get_sig_gene_sets <- function(sig_pathway_data, master_fgsea){
  
  signature_tally <- sig_pathway_data$sig_paths_tally
  table(signature_tally$max_sig)
  max_num <- length(unique(signature_tally$max_sig))

  # get number of genes involved per signature
  genes_per_max_sig <- list()
  for(i in 1:max_num){
    max_genes <- filter(signature_tally, grepl(paste0("s", i), max_sig))
    gene_count <- sum(max_genes[,paste0("max_signatures",i)])
    genes_per_max_sig[[paste0("max_s", i, "_genes")]] <- gene_count
  }

  # get the signature gene sets
  sig_pathways <- list()
  for(i in 1:max_num){
    sig_paths <- filter(signature_tally, grepl(paste0("s", i), max_sig))
    sig_pathways[[paste0("s", i, "_pathways")]] <- filter(master_fgsea, master_fgsea$pathway %in% sig_paths$pathway & grepl(paste0("s", i), contrast))
  }

  # all in one step for gsva format
  # unique gene sets for the most abundant pathways of each signature
  sig_genes <- lapply(names(sig_pathways), function(x){
    t(as.matrix(unique(unlist(select(sig_pathways[[x]], external_gene_name)))))
  })

  names(sig_genes) <- paste0("s", c(1:max_num))
  
  return(sig_genes)
}


#### JOIN CLINICAL ####
# used for joining the ssgsea output and clinical data
# matching_clinical - the clinical data that was used for DE analysis
# hall_sig_NES - output of ssgsea (hallmark signature NES)
# geneset_names - must be specified for formatting

join_clinical <- function(matching_clinical, hall_sig_NES, geneset_names){
  
  clinical_NES <- as.data.frame(t(hall_sig_NES))
  
  colnames(clinical_NES) <- geneset_names
  
  max_NES <- colnames(clinical_NES)[max.col(clinical_NES, ties.method = "first")]
  
  clinical_NES$max_NES <- max_NES
  
  rownames(clinical_NES) <- str_replace_all(rownames(clinical_NES), "[[:punct:]]", "-")
  
  clinical_NES <- setDT(clinical_NES, keep.rownames = 'submitter_id')[]
  
  NES_clinical_data <- left_join(matching_clinical, clinical_NES, by = c("submitter_id" = "submitter_id"))
}


#### PCA ####
# used to run principal component analysis with the PCAtools package
# Generates data for RNA signature only PCA plot, clinical & RNA sig PCA plot and 
# scree plot for clinical & RNA sigs
# hall_sig_NES - output of ssgsea (hallmark signature NES)
# new_clinical - new clinical data generated from join_clinical function (RNA sigs and clinical joined together)

run_PCA <- function(new_clinical, hall_sig_NES){
  
  dir.create(paste0(OUTDIR, "/PCA"))
  
  rownames(new_clinical) <- new_clinical[,1]
  
  # multiple formatting adjustments must be made to allow proper joining 
  new_clinical <- as.matrix(select(new_clinical, -submitter_id))
  rownames(new_clinical) <- str_replace_all(rownames(new_clinical), "[[:punct:]]", "-")
  
  NES_mat <- as.matrix(hall_sig_NES)
  colnames(NES_mat) <- str_replace_all(colnames(NES_mat), "[[:punct:]]", "-")
  
  NES_mat <- hall_sig_NES[,which(colnames(NES_mat) %in% rownames(new_clinical))]
  colnames(NES_mat) <- str_replace_all(colnames(NES_mat), "[[:punct:]]", "-")
  
  colnames(NES_mat) <- sort(colnames(NES_mat))
  rownames(new_clinical) <- sort(rownames(new_clinical))
  
  p1 <- pca(t(NES_mat), removeVar = 0.1)
  p2 <- pca(NES_mat, metadata = new_clinical, removeVar = 0.1)
  
  # just signature points - transposed matrix
  biplot_1 <- biplot(p1, lab = c(colnames(t(NES_mat))), shape = c(colnames(t(NES_mat))))
  
  # PCA of matrix and clinical data PCA
  biplot_2 <- biplot(p2, drawConnectors = TRUE)
  
  # scree plot of matrix and clinical data PCA
  scree <- screeplot(p2, axisLabSize = 18, titleLabSize = 22)
  
  PCA_plot_list <- list(sig_PCA = biplot_1,
                        sig_clinical_PCA = biplot_2,
                        scree_clinical = scree)
}
  

#### SURVIVAL ####
# performs survival analysis on the RNA and CNA signature classifications
# univariate and multivariate data generated for all possible signature combinations
# Kaplan-Meier curves generated for CNA and RNA signatures
# new_clinical - new clinical data generated from join_clinical function (RNA sigs and clinical joined together)
# surv_object - survival object specifying Surv(time = , event =)
# cb_friendly - colour blind friendly option, can be set to TRUE/FALSE

run_survival <- function(new_clinical, surv_object, cb_friendly = FALSE){
  
  dir.create(paste0(OUTDIR, "/survival"))
  
  new_clinical$max_signature <- factor(new_clinical$max_signature)
  new_clinical$max_NES <- factor(new_clinical$max_NES)
  
  # univariate and multivariate for all CNA and RNA signature combinations
  
  CNA_sigs <- unique(new_clinical$max_signature)
  RNA_sigs <- unique(new_clinical$max_NES)
  
  CNA_cox <- list()
  for(i in 1:length(CNA_sigs)){
    baseline <- paste0(CNA_sigs[i])
    new_clinical$max_signature <- relevel(new_clinical$max_signature, ref = baseline)
    res.cox <- coxph(Surv(OS.time, OS) ~ max_signature, data = new_clinical)
    res.multicox <- coxph(Surv(OS.time, OS) ~ max_signature + histological_grade + clinical_stage + 
                       age_at_initial_pathologic_diagnosis, data = new_clinical)
    CNA_cox[["univariate"]][[paste0(CNA_sigs[i], "_baseline")]] <- res.cox
    CNA_cox[["multivariate"]][[paste0(CNA_sigs[i], "_baseline")]] <- res.multicox
  }
  
  RNA_cox <- list()
  for(i in 1:length(RNA_sigs)){
    baseline <- paste0(RNA_sigs[i])
    new_clinical$max_NES <- relevel(new_clinical$max_NES, ref = baseline)
    res.cox <- coxph(Surv(OS.time, OS) ~ max_NES, data = new_clinical)
    res.multicox <- coxph(Surv(OS.time, OS) ~ max_NES + histological_grade + clinical_stage + 
                            age_at_initial_pathologic_diagnosis, data = new_clinical)
    RNA_cox[["univariate"]][[paste0(RNA_sigs[i], "_baseline")]] <- res.cox
    RNA_cox[["multivariate"]][[paste0(RNA_sigs[i], "_baseline")]] <- res.multicox
  }
  
  new_clinical$max_signature <- as.character(new_clinical$max_signature)
  
  sig_fit <- survfit(surv_object ~ max_signature, data = new_clinical)
 
  exp_fit <- survfit(surv_object ~ max_NES, data = new_clinical)
  
  if(cb_friendly == TRUE){
    sig_col <- c(brewer.pal(6, 'Paired'))
    names(sig_col) <- levels(factor(new_clinical$max_signature))
    rna_sig_col <- c(brewer.pal(6, 'Paired'))
    names(rna_sig_col) <- paste0("r", c(1:length(CNA_sigs)))
  } else {
    sig_col <- c(brewer.pal(6, 'Set1'))
    names(sig_col) <- levels(factor(new_clinical$max_signature))
    rna_sig_col <- c(brewer.pal(6, 'Set1'))
    names(rna_sig_col) <- paste0("r", c(1:length(CNA_sigs)))   # RNA assigned to same colours as CNA
  }
  
  ggsurvplot(sig_fit, data = new_clinical, pval = TRUE,
             title = "Survival Curve (Max Sig)",
             legend = "bottom",
             xlab = "Time(Days)",
             ylab = "Survival Probability",
             legend.title = "Max sig",
             legend.labs = unique(new_clinical$max_signature),
             palette = sig_col,
             font.x = c(15, "plain", "black"),
             font.y = c(15, "plain", "black"),
             pval.size = 5,
             font.legend = c(20, "plain", "black")) 
  ggsave(filename = paste0(OUTDIR, "/survival/sig_survival.pdf"), height = 10, width = 12)
  
  
    
  ggsurvplot(exp_fit, data = new_clinical, pval = TRUE,
             title = "Survival Curve (Sig_NES)",
             legend = "bottom",
             xlab = "Time(Days)",
             ylab = "Survival Probability",
             legend.title = "Max NES",
             legend.labs = unique(new_clinical$max_NES),
             palette = rna_sig_col,
             font.x = c(15, "plain", "black"),
             font.y = c(15, "plain", "black"),
             pval.size = 5,
             font.legend = c(20, "plain", "black"))
  ggsave(filename = paste0(OUTDIR, "/survival/NES_survival.pdf"), height = 10, width = 12)
  
  survial_results <- list(CNA_signatures = CNA_cox,
                          RNA_signatures = RNA_cox,
                          sig_fit = sig_fit,
                          exp_fit = exp_fit
                          )
}
















