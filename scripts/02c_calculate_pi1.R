library(qvalue)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
"%&%" = function(a,b) paste(a,b,sep="")
tissues <- c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood")
pi1s <- data.frame(tissue = character(), trans_acting = numeric(), cis_acting = numeric(),
                   cis_same = numeric(), cis_different = numeric()) # dataframe to store pi1 results
for (tissue in tissues) {
  print(tissue)
  pi1_vec <- c(tissue)
  results <- fread("output/9-16-2021/mashr_"%&%tissue%&%"/mashr_"%&%tissue%&%"_every_protein_association.txt") %>%
    select(protein, gene, pvalue)
  
  #Divide by cis and trans
  #First, annotate ENSG ID, chr, and bp for the target protein
  INTERVAL_protein_annotation <- fread("data/INTERVAL_data/INTERVAL_annotation_file.txt")
  annotated_proteins_results <- results %>%
    left_join(INTERVAL_protein_annotation, by=c("protein"="SomaId")) %>%
    mutate(protein_ENSG=ENSG_id, protein_chr=chr, protein_start=start) %>%
    select(protein, protein_ENSG, protein_chr, protein_start, gene, pvalue)
  rm(results)
  
  #Next, annotate chr and bp for the predicted transcript
  predicted_transcript_annotation <- fread("data/TOPMed_ALL_data/gene_annotation_gencode.v37_hg38_sanitized_BP_CHROM.txt") %>%
    separate(gene_id,sep="\\.",into=c("ENSG_id","gene_version"))
  #Remove version number off of transcript ENSG IDs
  annotated_proteins_transcripts_results <- annotated_proteins_results %>%
    separate(gene,sep="\\.",into=c("transcript_ENSG","gene_version")) %>%
    left_join(predicted_transcript_annotation, by=c("transcript_ENSG"="ENSG_id")) %>%
    mutate(transcript_chr = chrchr, transcript_start = start) %>%
    select(protein, protein_ENSG, protein_chr, protein_start, transcript_ENSG, transcript_chr, transcript_start, pvalue)
  rm(annotated_proteins_results)
  
  #Finally, annotate for cis and trans
  #If different chromosome, trans
  #If same chromosome but greater than 1000000 bp away, trans
  #If same chromosome but less than 1000000 bp away, cis
  #Annotate for cis_same or cis-different
  results_cis_trans_annotated <- annotated_proteins_transcripts_results %>%
    mutate(cis_trans = ifelse(!protein_chr==transcript_chr, "trans_acting",
                              ifelse((protein_start < (transcript_start-1000000)) | (protein_start > (transcript_start+1000000)), "trans_acting","cis_acting"))) %>%
    mutate(cis_same = ifelse(transcript_ENSG == protein_ENSG, "cis_same",
                             ifelse(cis_trans == "trans_acting", "trans_acting", "cis_different")))
  rm(annotated_proteins_transcripts_results)
  print(table(results_cis_trans_annotated$cis_trans))
  print(table(results_cis_trans_annotated$cis_same))
  
  #Now loop through and calculate pi1 and make histogram
  groups <- c("trans_acting","cis_acting", "cis_same","cis_different")
  for(i in 1:length(groups)) {
    print(groups[i])
    if (i == 2) {
      filtered_results <- results_cis_trans_annotated %>%
        filter(cis_trans == groups[i])
    } else {
      filtered_results <- results_cis_trans_annotated %>%
        filter(cis_same == groups[i])
    }
    
    pvalues <- filtered_results$pvalue
    qobj <- qvalue(p=pvalues)
    
    pi0 <- qobj$pi0 #proportion of true null hypotheses
    pi1 <- 1-pi0
    pi1_vec <- append(pi1_vec,pi1)
    
    png(file="output/9-9-2022/"%&%tissue%&%"_"%&%groups[i]%&%"pvalue_histogram.png")
    print(hist(qobj))
    dev.off()
  }
  pi1s <- rbind(pi1s, pi1_vec)
  rm(filtered_results)
  rm(results_cis_trans_annotated)
}
colnames(pi1s) <- c("tissue","trans_acting","cis_acting","cis_same","cis_different")
fwrite(pi1s, "output/9-9-2022/pi1s.txt",sep="\t",col.names=T)
