library(qvalue)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
library(tidyr)
"%&%" = function(a,b) paste(a,b,sep="")

#Loop through each tissue and pull out FDR significant hits
#Also keep track of how many tissues each significant pair is found in
tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood")

#Do it for one tissue before looping through rest
tissue <- tissues[1]
all_results <- fread("output/9-16-2021/mashr_" %&% tissue %&% "/mashr_" %&% tissue %&%"_every_protein_association.txt") %>%
  separate(gene,sep="\\.",into=c("gene","gene_version")) %>%
  mutate(pair = protein %&% "_" %&% gene, tissue = tissue)

## Perform this if calculating within tissue
#Calculate q values
pvalues <- all_results$pvalue
qobj <- qvalue(p=pvalues)

#all_results$qvalue <- qobj$qvalues

#Filter by FDR 0.05
#all_significant_results <- all_results %>%
#filter(qvalue < 0.05)

for (j in 2:length(tissues)) {
  tissue <- tissues[j]
  print(tissue)
  results <- fread("output/9-16-2021/mashr_" %&% tissue %&% "/mashr_" %&% tissue %&%"_every_protein_association.txt") %>%
    separate(gene,sep="\\.",into=c("gene","gene_version")) %>%
    mutate(pair = protein %&% "_" %&% gene, tissue = tissue)
  ## Perform this is calculating within tissue
  #Calculate q values
  #pvalues <- results$pvalue
  #qobj <- qvalue(p=pvalues)
  
  #results$qvalue <- qobj$qvalues
  
  #Filter by FDR 0.05
  #significant_results <- results %>%
  #filter(qvalue < 0.05)
  
  #all_significant_results <- rbind(all_significant_results, significant_results)
  all_results <- rbind(all_results,results)
}
#fwrite(all_significant_results, "output/9-23-2022/INTERVAL_all_tissues_FDR0.05_within_tissue_pairs.txt")

##Run this is calculating qvalue across tissue
#Calculate q values
pvalues <- all_results$pvalue
qobj <- qvalue(p=pvalues)

all_results$qvalue <- qobj$qvalues

fwrite(all_results, "output/9-23-2022/INTERVAL_all_tissues_across_tissue_qvals_pairs.txt")

#Filter by FDR 0.05
all_significant_results <- all_results %>%
  filter(qvalue < 0.05)

fwrite(all_significant_results, "output/9-23-2022/INTERVAL_all_tissues_FDR0.05_across_tissue_qvals_pairs.txt")
