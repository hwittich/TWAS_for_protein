#We have a list of INTERVAL FDR0.05 hits (qvalues calculated within tissue) that were tested in TOPMed
#Loop through those, calculate qvalues and filter out hits that aren't FDR 0.05 significant
#Compile these lists across tissues for plotting

library(qvalue)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
"%&%" = function(a,b) paste(a,b,sep="")
tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood")

#Do it for one tissue before looping through rest
tissue <- tissues[1]
all_significant_results <- fread("output/9-16-2022/"%&%tissue%&%"_INTERVAL_TOPMed_qvalues.txt") %>%
  filter(TOPMed_qvalue<0.05) %>%
  mutate(tissue=tissue)

for (j in 2:length(tissues)) {
  tissue <- tissues[j]
  print(tissue)
  
  significant_results <- fread("output/9-16-2022/"%&%tissue%&%"_INTERVAL_TOPMed_qvalues.txt") %>%
    filter(TOPMed_qvalue<0.05) %>%
    mutate(tissue=tissue)
  
  all_significant_results <- rbind(all_significant_results, significant_results)
}

fwrite(all_significant_results, "output/9-16-2022/all_INTERVAL_FDR0.05_within_tissue_TOPMed_replicated_FDR0.05_within_tissue_pairs.txt")
