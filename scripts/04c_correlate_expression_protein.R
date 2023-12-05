library(qvalue)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
library(tidyr)
"%&%" = function(a,b) paste(a,b,sep="")

#Pull all genes that have a significant cis same association in any tissue, predicted or observed
all_sig_cis_same_results <- fread("output/3-9-2023/all_signficant_FDR0.05_cis_same_results_TOPMed.txt")
sig_pairs <- unique(all_sig_cis_same_results$transcript_ENSG)

#Prep protein and expression matrices
transposed_prot_expr <- fread("data/TOPMed_ALL_data/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt", header=T) %>% 
  separate(joint_id, sep="_", into = c("SLID","ENSG_ID")) %>%
  separate(ENSG_ID, sep="\\.", into = c("ENSG_ID","gene_version")) %>%
  select(!c("SLID","gene_version")) %>%
  distinct(ENSG_ID, .keep_all = T)

proteins <- unique(transposed_prot_expr$ENSG_ID)
IID <- colnames(transposed_prot_expr)
IID <- IID[2:length(IID)]
transposed_prot_expr <- transposed_prot_expr %>% select(!ENSG_ID)
prot_expr <- as.data.frame(t(transposed_prot_expr))
prot_expr <- cbind(IID, prot_expr)
colnames(prot_expr) <- c("IIDs",proteins)

prot_expr <- prot_expr %>% mutate(IID = as.numeric(IIDs), .before=IIDs, .keep="unused")

#Calculate correlation values
pred_tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta",
                  "Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24",
                  "Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum",
                  "Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus",
                  "Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
                  "Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue",
                  "Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid",
                  "Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa",
                  "Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle",
                  "Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal",
                  "Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate",
                  "Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg",
                  "Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus",
                  "Vagina","Whole_Blood")
obs_tissues <- c("Mono","PBMC","Tcell")
#Loop through every tissue and calculate correlation of expression with protein levels for every gene
corr_matrix <- data.frame(gene = sig_pairs)

for(tissue in pred_tissues) {
  print(tissue)
  corr_vec <- c() #Vector for storing correlation values as I loop through the genes
  
  pred_expr <- fread("output/4-7-2021/mashr_"%&%tissue%&%"/mashr_"%&%tissue%&%"_predict.txt")
  #Change column names
  gene_names <- colnames(pred_expr)
  gene_names <- gene_names[3:length(gene_names)]
  #Remove version numbers off of ENSG IDs
  shortened_gene_names <- unlist(lapply(gene_names, substr, start=1, stop=15))
  new_colnames <- c("FID","IID", shortened_gene_names)
  colnames(pred_expr) <- new_colnames
  
  pred_expr <- pred_expr %>% select(!FID)
  
  #Loop through genes
  for(gene in sig_pairs) {
    #print(gene)
    #Check if this tissue has a model for that gene.
    if(gene %in% colnames(pred_expr)) {
      #print("Gene in model")
      #Pull expression levels and protein levels
      expression_levels <- pred_expr %>% select(IID, gene)
      colnames(expression_levels) <- c("IID","predicted_expression")
      protein_levels <- prot_expr %>% select(IID, gene)
      colnames(protein_levels) <- c("IID","protein_abundance")
      
      combined_levels <- inner_join(expression_levels, protein_levels, by=c("IID"="IID"))
      corr_val <- cor(combined_levels$predicted_expression, combined_levels$protein_abundance, method="pearson")
      
      corr_vec <- append(corr_vec, corr_val)
    } else {
      #print("Gene not in model")
      #Store NA as correlation value
      corr_vec <- append(corr_vec, NA)
    }
  }
  corr_matrix[,tissue] <- corr_vec
}

#Repeat for observed expression tissues
for(tissue in obs_tissues) {
  print(tissue)
  corr_vec <- c() #Vector for storing correlation values as I loop through the genes
  
  obs_expr <- fread("/home/chris/topmed_expression_whole_genome/expression/"%&%tissue%&%"_expression_ALLage_sex_adj_rinv_PC10_gene_id_tpm_0.1_expression_PCs10.txt", header = T) %>%
    mutate(IID = sidno, .before=sidno, .keep="unused")
  
  #Loop through genes
  for(gene in sig_pairs) {
    #Check if this tissue has a model for that gene.
    if(gene %in% colnames(obs_expr)) {
      #print("Gene in model")
      #Pull expression levels and protein levels
      expression_levels <- obs_expr %>% select(IID, gene)
      colnames(expression_levels) <- c("IID","observed_expression")
      protein_levels <- prot_expr %>% select(IID, gene)
      colnames(protein_levels) <- c("IID","protein_abundance")
      
      combined_levels <- inner_join(expression_levels, protein_levels, by=c("IID"="IID"))
      corr_val <- cor(combined_levels$observed_expression, combined_levels$protein_abundance, method="pearson")
      
      corr_vec <- append(corr_vec, corr_val)
    } else {
      #print("Gene not in model")
      #Store NA as correlation value
      corr_vec <- append(corr_vec, NA)
    }
  }
  corr_matrix[,tissue] <- corr_vec
}

fwrite(corr_matrix, "output/4-24-2023/significant_cis_same_genes_expression_protein_levels_correlation_matrix.txt")

#Pull out max correlation values
obs_sig_correlations <- corr_matrix %>% select(all_of(obs_tissues))
pred_sig_correlations <- corr_matrix %>% select(all_of(pred_tissues))

max_corr_vals <- data.frame(gene = sig_pairs)
max_corr_vals$observed_expression <- apply(obs_sig_correlations, 1, max, na.rm=T)
max_corr_vals$observed_expression[max_corr_vals$observed_expression == -Inf] <- NA
max_corr_vals$predicted_expression <- apply(pred_sig_correlations, 1, max, na.rm=T)
max_corr_vals$predicted_expression[max_corr_vals$predicted_expression == -Inf] <- NA

fwrite(max_corr_vals, "output/4-24-2023/significant_cis_same_genes_expression_protein_levels_max_correlations_matrix.txt")

