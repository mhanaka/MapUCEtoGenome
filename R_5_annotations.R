# Script info ----
#Script to data wrangle the outputs from mapping UCEs to genome
#Requirements: datawrangle.R outputs and "Synthetic chromosome" section in chromomap.R
#Written with R version 4.3.3 (2024-02-29)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: April 2024

## Load prerequisite ----
source(file="R_2_datawrangle.R")

## Read files ----
gff_Amil <- read.delim("data/genomesToMatch/Acropora_millepora_GCF_013753865_1/Acropora_millepora_GCF_013753865_1_with_introns.gff",skip=8,header=FALSE)
gaf_Amil <- read.delim("data/genomesToMatch/Acropora_millepora_GCF_013753865_1/GCF_013753865.1_Amil_v2.1_gene_ontology.gaf",skip=8,header=TRUE)
gff_Nvect <- read.delim("data/genomesToMatch/Nematostella_vectensis_GCF_932526225_1/Nematostella_vectensis_GCF_932526225_1_with_introns.gff",skip=8,header=FALSE)
gaf_Nvect <- read.delim("data/genomesToMatch/Nematostella_vectensis_GCF_932526225_1/GCF_932526225.1_jaNemVect1.1_gene_ontology.gaf",skip=8,header=TRUE)
  #Ahya and Cjar are only annotated to type level, not gene description level. Removed from analyses

gaf_Amil = gaf_Amil %>% mutate(genes=paste0('gene-',Symbol)) %>% filter(Aspect=='P') %>% select(genes,GO_ID) %>% distinct() %>% 
  group_by(genes) %>% summarise(GO_ID=paste0(GO_ID,collapse=', ')) %>% ungroup()    ##P=Biological processes
gaf_Nvect = gaf_Nvect %>% mutate(genes=paste0('gene-',Symbol)) %>% filter(Aspect=='P') %>% select(genes,GO_ID) %>% distinct() %>% 
  group_by(genes) %>% summarise(GO_ID=paste0(GO_ID,collapse=', ')) %>% ungroup() 


# All UCE identified with annotations ----
tmp_Amil <- Acropora_millepora_G_chromo_genome_relaxed %>% select(UCE,`GFF type list`) %>% 
  mutate(genes=case_when(grepl("ID=",`GFF type list`) ~ 
                           case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ #If it matches this then there are multiple gene hits
                                       paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                     TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                         TRUE ~ ""), genome="Amil") %>% select(-`GFF type list`) %>% 
  mutate(genes=strsplit(genes,";")) %>% unnest(genes)
uce_relaxed_gff_Amil <- left_join(tmp_Amil, 
                                gff_Amil %>% select(V9) %>% mutate(genes=case_when(grepl("Parent=",V9) ~ str_extract(V9,"Parent=(.*)\\;Dbxref.*",group=1), TRUE ~ ''),
                                                                   annot=case_when(grepl("Parent=",V9) ~ str_extract(V9,"product=(.*)\\;transcript.*",group=1), TRUE ~ '')) %>% 
                                  select(-V9) %>% filter(grepl("gene-LOC",genes)==T) %>% group_by(genes) %>% mutate(annot=paste0(annot, collapse=" OR ")) %>% ungroup() %>% distinct(),
                                by='genes',unmatched='drop')
uce_relaxed_gff_gaf_Amil <- left_join(uce_relaxed_gff_Amil,gaf_Amil,by='genes') %>%
  group_by(UCE) %>% mutate(genes=paste0(genes, collapse=" AND "),
                                 annot=paste0(annot, collapse=" AND "),
                                 GO_ID=paste0(GO_ID, collapse=" AND ")) %>% ungroup() %>% distinct()

tmp_Nvect <- Nematostella_vectensis_G_chromo_genome_relaxed %>% select(UCE,`GFF type list`) %>% 
  mutate(genes=case_when(grepl("ID=",`GFF type list`) ~ 
                           case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ #If it matches this then there are multiple gene hits
                                       paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                     TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                         TRUE ~ ""), genome="Nvect") %>% select(-`GFF type list`) %>% 
  mutate(genes=strsplit(genes,";")) %>% unnest(genes)
uce_relaxed_gff_Nvect <- left_join(tmp_Nvect,
                                 gff_Nvect %>% select(V9) %>% mutate(genes=case_when(grepl("Parent=",V9) ~ str_extract(V9,"Parent=(.*)\\;Dbxref.*",group=1), TRUE ~ ''),
                                                                     annot=case_when(grepl("Parent=",V9) ~ str_extract(V9,"product=(.*)\\;transcript.*",group=1), TRUE ~ '')) %>% 
                                   select(-V9) %>% filter(grepl("gene-LOC",genes)==T) %>% group_by(genes) %>% mutate(annot=paste0(annot, collapse=" OR ")) %>% ungroup() %>% distinct(),
                                 by='genes',unmatched='drop')
uce_relaxed_gff_gaf_Nvect <- left_join(uce_relaxed_gff_Nvect,gaf_Nvect,by='genes') %>%
  group_by(UCE) %>% mutate(genes=paste0(genes, collapse=" AND "),
                                 annot=paste0(annot, collapse=" AND "),
                                 GO_ID=paste0(GO_ID, collapse=" AND ")) %>% ungroup() %>% distinct()

uce_relaxed_gff <- rbind(uce_relaxed_gff_gaf_Amil, uce_relaxed_gff_gaf_Nvect) %>% filter(!grepl("NA|NA AND NA",GO_ID)) %>% mutate(UCE=paste0('uce-',UCE))
uce_relaxed_gff_summary = uce_relaxed_gff %>% pivot_wider(names_from=genome, values_from=c(genes,annot,GO_ID)) 

uce_relaxed_gff_summary_n = rbind(uce_relaxed_gff_gaf_Amil %>% select(UCE,genome,GO_ID) %>% mutate(GO_ID=strsplit(GO_ID,', ')) %>% unnest(GO_ID),
                                uce_relaxed_gff_gaf_Nvect %>% select(UCE,genome,GO_ID) %>% mutate(GO_ID=strsplit(GO_ID,', ')) %>% unnest(GO_ID)) %>% 
  group_by(GO_ID,genome) %>% summarise(n=n()) %>% ungroup() %>% pivot_wider(names_from=genome,values_from=n) %>% arrange(desc(Amil))


# Output ----
write_tsv(uce_relaxed_gff_summary, file="export/mappedUCE_annotations.tsv")
write_tsv(uce_relaxed_gff_summary_n, file="export/mappedUCE_annotations_summary.tsv")
