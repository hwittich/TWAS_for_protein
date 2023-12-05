suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
"%&%" = function(a,b) paste(a,b,sep="")

INTERVAL_sig_hits <- fread("output/9-23-2022/INTERVAL_all_tissues_FDR0.05_pairs.txt")

unique_significant_INTERVAL_results <- INTERVAL_sig_hits %>% distinct(pair, .keep_all = T)

num_associations <- as.data.frame(table(unique_significant_INTERVAL_results$gene))

#num_associations <- as.data.frame(table(INTERVAL_sig_hits$gene))
colnames(num_associations) <- c("transcript_ENSG","num_targets")
master_regulators <- num_associations[order(-num_associations$num_targets),]

predicted_transcript_annotation <- fread("data/TOPMed_ALL_data/gene_annotation_gencode.v37_hg38_sanitized_BP_CHROM.txt") %>%
  separate(gene_id,sep="\\.",into=c("ENSG_id","gene_version"))
#Remove version number off of transcript ENSG IDs

annotated_INTERVAL_master_regulators <- master_regulators %>%
  left_join(predicted_transcript_annotation, by=c("transcript_ENSG"="ENSG_id"))

#Now I want to cluster the master regulators into loci
#If two master regulators are within 200,000 bp of each other, combine into one locus

INTERVAL_master_regulators <- annotated_INTERVAL_master_regulators %>% filter(num_targets >= 50) %>%
  arrange(desc(start))
master_regulatory_genes <- c(INTERVAL_master_regulators$transcript_ENSG)

#Vector for storing which locus each gene gets assigned to
locus_annotations <- c()

#Initialize a locus based on the first gene
gene <- INTERVAL_master_regulators %>% filter(transcript_ENSG == master_regulatory_genes[1])
print(gene$gene_name)
print(gene$chrchr)
print(gene$start)

master_regulatory_loci <- data.frame(locus = "locus1", chr = gene$chrchr, start = gene$start, end = gene$start)
locus_annotations <- append(locus_annotations, "locus1")
print("locus1")
master_regulatory_genes <- master_regulatory_genes[2:length(master_regulatory_genes)]
#Loop through the master regulators and figure out the bounds of the master regulatory loci
while(length(master_regulatory_genes != 0)) {
  #Loop through remaining genes, if they fit in an existing locus, edit the bounds of the locus if necessary
  #Else, make a new locus centered around the gene
  gene <- INTERVAL_master_regulators %>% filter(transcript_ENSG == master_regulatory_genes[1])
  print(gene$gene_name)
  print(gene$chrchr)
  print(gene$start)
  
  locus_match = F
  #Loop through loci we have identified
  for (row in 1:nrow(master_regulatory_loci)) {
    locus <- master_regulatory_loci[row,]
    
    #If same chromosome and +- 200,000bp away from ends of locus, add
    if (gene$chrchr == locus$chr & (gene$start > locus$start-200000 & gene$start < locus$end+200000)) {
      locus_match = T
      locus_annotations <- append(locus_annotations, locus$locus)
      print(locus$locus)
      
      #Update bounds of locus if necessary
      if (gene$start < locus$start) {
        locus$start <- gene$start
      }
      if (gene$start > locus$end) {
        locus$end <- gene$start
      }
      
      master_regulatory_loci[row,] <- locus
    }
  }
  
  print(locus_match)
  #If not added to any loci, create new locus
  if (locus_match == F) {
    new_locus <- data.frame(locus = "locus"%&%(nrow(master_regulatory_loci)+1), chr = gene$chrchr, start = gene$start, end = gene$start)
    locus_annotations <- append(locus_annotations, new_locus$locus)
    print(locus$locus)
    master_regulatory_loci <- rbind(master_regulatory_loci, new_locus)
  }
  
  #Remove master regulatory gene from loop
  if(length(master_regulatory_genes) > 1) {
    master_regulatory_genes <- master_regulatory_genes[2:length(master_regulatory_genes)]
  } else {
    master_regulatory_genes <- character()
  }
}