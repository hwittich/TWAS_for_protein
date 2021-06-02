suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
"%&%" = function(a,b) paste(a,b,sep="")

args <- commandArgs(trailingOnly=TRUE)
results <- args[1] #Protein association results
model <- args[2] #GTEx tissue type
gene_anno <- args[3] #Gene annotation file from TOPMed
protein_anno <- args[4] #Protein annotation file from GENCODE
multigene <- args[5] #List of multigene aptamers
outdir <- args[6] #Directory to put all plots

protein_associations <- read.table(results,header=TRUE,sep="\t")

#Calculate the bonferroni threshold
#=0.5/((#proteins tested * #transcripts predicted)/total # of tests)
###These values are specific to my input data and may need to be adjusted
n_proteins<-1279 #length of protein matrix (sort of, duplicates were removed)
n_transcripts<-length(protein_associations$gene)/n_proteins #length of association output
alpha<-0.05
bonferroni_threshold=alpha/(n_proteins*n_transcripts)
#Filter out non-significant pairs
sig_hits <- filter(protein_associations,pvalue<bonferroni_threshold) %>% separate(gene,sep="\\.",into=c("gene_id","gene_version"))

#Now annotate the positions of the predicted transcript and the protein
gene_annotation <- read.table(gene_anno,header=T,sep=" ") %>% select(gene_id,chrchr,start)
sig_hits <- left_join(sig_hits, gene_annotation, by = c("gene_id"="gene_id"), copy=TRUE) %>% rename("gene_chr"="chrchr","gene_start"="start")

protein_annotation <- read.table(protein_anno,header=T,sep="\t") %>% select(SomaId,ENSG_id,chr,start) %>% rename("protein_ID"="ENSG_id","protein_chr"="chr","protein_start"="start")
sig_hits <- left_join(sig_hits,protein_annotation, by = c("protein"="SomaId"), copy=TRUE)

#Label the multigene aptamers
multigene_aptamers <- read.table(multigene,header=T) %>% separate(joint_ID,sep="_",into=c("SL_ID","ENSG_ID"))
sig_hits <- mutate(sig_hits, multigene = ifelse(protein%in%multigene_aptamers$SL_ID,TRUE,FALSE))

#Now adjust for relative chromosomal positions
chromosomes<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
lengths<-c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415)
#Calculate chr start positions by summing length of prior chromosomes
starts <- rep(NA,24)
sum=0
for(i in 1:24){
  starts[i] = sum
  sum = sum + lengths[i]
}
chr_len<-data.frame(chromosomes,starts)

#Now join the lengths to the sig_hits dataframe
sig_hits<-left_join(sig_hits,chr_len,by=c("gene_chr"="chromosomes"),copy=TRUE) %>% rename("gene_chr_start"="starts")
sig_hits<-left_join(sig_hits,chr_len,by=c("protein_chr"="chromosomes"),copy=TRUE) %>% rename("protein_chr_start"="starts")

#Finally, add the chromosome starting position to the gene starting positions
sig_hits<-mutate(sig_hits,new_protein_start = protein_chr_start + protein_start) %>% mutate(new_gene_start = gene_chr_start+gene_start)

#Then, let's compare the positions to determine if a gene is trans-acting or cis-acting
sig_hits<-mutate(sig_hits,cis_trans = ifelse((new_gene_start-1000000)<=new_protein_start & new_protein_start<=(new_gene_start+1000000),"cis-acting","trans-acting"))

#Finally, plot!
fig1 <- ggplot(sig_hits,aes(x=new_gene_start,y=new_protein_start,size=pvalue,color=cis_trans,group=multigene)) +
  geom_point(aes(shape=multigene)) +
  scale_shape_manual(values=c(21,22)) +
  scale_size_continuous(guide=FALSE) +
  scale_color_manual(values=c("dark gray","dark red"),na.value="dark gray") +
  geom_hline(yintercept=starts[1:22],size=0.1) +
  geom_vline(xintercept=starts[1:22],size=0.1) +
  coord_cartesian(xlim=c(starts[1],starts[22]),ylim=c(starts[1],starts[22])) +
  theme_classic(10) +
  labs(x="Gene",y="Target Protein",title=model) +
  scale_x_continuous(breaks=starts[1:22],labels=c(1:14,"",16,"",18,"",20,"",22)) +
  scale_y_continuous(breaks=starts[1:22],labels=c(1:16,"",18,"",20,"",22))

png(filename=outdir %&% model %&% "_chromosomal_position_of_significant_pairs.png",width=600,height=600)
fig1
dev.off()
