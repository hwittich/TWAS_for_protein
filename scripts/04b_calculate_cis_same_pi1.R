library(qvalue)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
library(tidyr)
"%&%" = function(a,b) paste(a,b,sep="")

#First, TOPMed
tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta",
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
pi1s <- data.frame(tissue = character(), trans_acting = numeric(), cis_acting = numeric(),
                   cis_same = numeric(), cis_different = numeric()) # dataframe to store pi1 results
#all_results <- fread("output/3-9-2023/TOPMed_all_tissues_within_tissue_qvals_pairs.txt")

end_time <- Sys.time() #to feed as start time into the loop
for (tissue in tissues) {
  start_time <- end_time
  print(tissue)
  pi1_vec <- c(tissue)
  
  #Check to see if joined dataframe already created
  if(!file.exists("output/3-9-2023/"%&%tissue%&%"pred_expr_annotated_TOPMed_qvalues.txt")) {
    print("Annotating for cis_trans")
    #If not, create it
    tissue_results <- fread("output/4-7-2021/mashr_" %&% tissue %&% "/mashr_" %&% tissue %&%"_every_protein_association.txt") %>%
      separate(gene,sep="\\.",into=c("gene","gene_version")) %>%
      mutate(transcript_ENSG = gene) %>%
      select(protein, transcript_ENSG, pvalue)
    
    #Annotate TOPMed proteins
    TOPMed_protein_annotation <- fread("data/TOPMed_ALL_data/somascan1.3k_gencode.v32_annotation.txt", header=T)
    TOPMed_annotated_proteins <- tissue_results %>%
      left_join(TOPMed_protein_annotation, by=c("protein"="SomaId")) %>%
      separate(ENSG_id, sep="\\.",into=c("ENSG_id","gene_version")) %>%
      mutate(protein_ENSG=ENSG_id, protein_chr=chr, protein_start=start) %>%
      select(protein, protein_ENSG, protein_chr, protein_start, transcript_ENSG, pvalue)
    rm(tissue_results)
    
    #Now annotate the location of the transcript
    predicted_transcript_annotation <- fread("data/TOPMed_ALL_data/gene_annotation_gencode.v37_hg38_sanitized_BP_CHROM.txt") %>%
      separate(gene_id,sep="\\.",into=c("ENSG_id","gene_version"))
    #Remove version number off of transcript ENSG IDs
    TOPMed_annotated_pairs <- TOPMed_annotated_proteins %>%
      left_join(predicted_transcript_annotation, by=c("transcript_ENSG"="ENSG_id")) %>%
      mutate(transcript_chr = chrchr, transcript_start = start) %>%
      select(protein, protein_ENSG, protein_chr, protein_start, transcript_ENSG, transcript_chr, transcript_start, pvalue)
    rm(TOPMed_annotated_proteins)
    
    #Finally annotate cis_trans
    #If different chromosome, trans
    #If same chromosome but greater than 1000000 bp away, trans
    #If same chromosome but less than 1000000 bp away, cis
    #Annotate for cis_same or cis-different
    results_cis_trans_annotated <- TOPMed_annotated_pairs %>%
      mutate(cis_trans = ifelse(!protein_chr==transcript_chr, "trans_acting",
                                ifelse((protein_start < (transcript_start-1000000)) | (protein_start > (transcript_start+1000000)), "trans_acting","cis_acting"))) %>%
      mutate(cis_same = ifelse(transcript_ENSG == protein_ENSG, "cis_same",
                               ifelse(cis_trans == "trans_acting", "trans_acting", "cis_different")))
    rm(TOPMed_annotated_pairs)
    print(table(results_cis_trans_annotated$cis_trans))
    print(table(results_cis_trans_annotated$cis_same))
    
    fwrite(results_cis_trans_annotated, "output/3-9-2023/"%&%tissue%&%"pred_expr_annotated_TOPMed_qvalues.txt",sep="\t")
  } else {
    #Read it in
    print("Loading annotated pairs")
    results_cis_trans_annotated <- fread("output/3-9-2023/"%&%tissue%&%"pred_expr_annotated_TOPMed_qvalues.txt")
  }
  
  #Now pull out cis_same results, calculate qvalues, and save significant pairs
  filtered_results <- results_cis_trans_annotated %>%
    filter(cis_same == "cis_same")
  
  pvalues <- filtered_results$pvalue
  qobj <- try(qvalue(p=pvalues, lambda = seq(0.05, 0.75, 0.05)))
  
  ### Save significant results
  filtered_results$cis_same_qvalue <- qobj$qvalues
  fwrite(filtered_results, "output/3-9-2023/"%&%tissue%&%"TOPMed_cis_same_results.txt")
  significant_cis_same_results <- filtered_results %>% filter(cis_same_qvalue<0.05)
  fwrite(significant_cis_same_results, "output/3-9-2023/"%&%tissue%&%"TOPMed_significant_cis_same_results.txt")
  
  rm(filtered_results)
  rm(results_cis_trans_annotated)
  
  end_time <- Sys.time()
  print(end_time-start_time)
}
rm(all_results)

# Now observed expression
tissues <- c("Mono","PBMC","Tcell")
pi1s <- data.frame(tissue = character(), trans_acting = numeric(), cis_acting = numeric(),
                   cis_same = numeric(), cis_different = numeric()) # dataframe to store pi1 results
#all_results <- fread("output/3-7-2023/every_tissue_obs_expr_results")

end_time <- Sys.time() #to feed as start time into the loop
for (tissue in tissues) {
  start_time <- end_time
  print(tissue)
  pi1_vec <- c(tissue)
  
  #Check to see if joined dataframe already created
  if(!file.exists("output/3-9-2023/"%&%tissue%&%"obs_expr_annotated_TOPMed_qvalues.txt")) {
    print("Annotating for cis_trans")
    #If not, create it
    tissue_results <- fread("output/3-7-2023/"%&%tissue%&%"/obs_expr_"%&%tissue%&%"_every_protein_association.txt") %>%
      mutate(transcript_ENSG = gene) %>%
      select(protein, transcript_ENSG, pvalue)
    
    #Annotate TOPMed proteins
    TOPMed_protein_annotation <- fread("data/TOPMed_ALL_data/somascan1.3k_gencode.v32_annotation.txt", header=T)
    TOPMed_annotated_proteins <- tissue_results %>%
      left_join(TOPMed_protein_annotation, by=c("protein"="SomaId")) %>%
      separate(ENSG_id, sep="\\.",into=c("ENSG_id","gene_version")) %>%
      mutate(protein_ENSG=ENSG_id, protein_chr=chr, protein_start=start) %>%
      select(protein, protein_ENSG, protein_chr, protein_start, transcript_ENSG, pvalue)
    rm(tissue_results)
    
    #Now annotate the location of the transcript
    predicted_transcript_annotation <- fread("data/TOPMed_ALL_data/gene_annotation_gencode.v37_hg38_sanitized_BP_CHROM.txt") %>%
      separate(gene_id,sep="\\.",into=c("ENSG_id","gene_version"))
    #Remove version number off of transcript ENSG IDs
    TOPMed_annotated_pairs <- TOPMed_annotated_proteins %>%
      left_join(predicted_transcript_annotation, by=c("transcript_ENSG"="ENSG_id")) %>%
      mutate(transcript_chr = chrchr, transcript_start = start) %>%
      select(protein, protein_ENSG, protein_chr, protein_start, transcript_ENSG, transcript_chr, transcript_start, pvalue)
    rm(TOPMed_annotated_proteins)
    
    #Finally annotate cis_trans
    #If different chromosome, trans
    #If same chromosome but greater than 1000000 bp away, trans
    #If same chromosome but less than 1000000 bp away, cis
    #Annotate for cis_same or cis-different
    results_cis_trans_annotated <- TOPMed_annotated_pairs %>%
      mutate(cis_trans = ifelse(!protein_chr==transcript_chr, "trans_acting",
                                ifelse((protein_start < (transcript_start-1000000)) | (protein_start > (transcript_start+1000000)), "trans_acting","cis_acting"))) %>%
      mutate(cis_same = ifelse(transcript_ENSG == protein_ENSG, "cis_same",
                               ifelse(cis_trans == "trans_acting", "trans_acting", "cis_different")))
    rm(TOPMed_annotated_pairs)
    print(table(results_cis_trans_annotated$cis_trans))
    print(table(results_cis_trans_annotated$cis_same))
    
    fwrite(results_cis_trans_annotated, "output/3-9-2023/"%&%tissue%&%"obs_expr_annotated_TOPMed_qvalues.txt",sep="\t")
  } else {
    #Read it in
    print("Loading annotated pairs")
    results_cis_trans_annotated <- fread("output/3-9-2023/"%&%tissue%&%"obs_expr_annotated_TOPMed_qvalues.txt")
  }
  
  #Now pull out cis_same results, calculate qvalues, and save significant pairs
  filtered_results <- results_cis_trans_annotated %>%
    filter(cis_same == "cis_same")
  
  pvalues <- filtered_results$pvalue
  qobj <- try(qvalue(p=pvalues, lambda = seq(0.05, 0.75, 0.05)))
  
  ### Save significant results
  filtered_results$cis_same_qvalue <- qobj$qvalues
  fwrite(filtered_results, "output/3-9-2023/"%&%tissue%&%"TOPMed_cis_same_results.txt")
  significant_cis_same_results <- filtered_results %>% filter(cis_same_qvalue<0.05)
  fwrite(significant_cis_same_results, "output/3-9-2023/"%&%tissue%&%"TOPMed_significant_cis_same_results.txt")
  
  rm(filtered_results)
  rm(results_cis_trans_annotated)
  
  end_time <- Sys.time()
  print(end_time-start_time)
  
}
rm(all_results)