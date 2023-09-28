# Script info ----
#Script to data wrangle the outputs from mapping UCEs to genome
#Requirements: datawrangle.R outputs and "Synthetic chromosome" section in chromomap.R
#Written with R-4.2.2
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: August 2023

## Load prerequisite ----
source(file="datawrangle.R")
  #Also run # Synthetic chromosome ---- section in chromomap.R

## Read files ----
#gff_Ahya <- read.delim("data/genomesToMatch/Acropora_hyacinthus_GCA_020536085_1/Acropora_hyacinthus_GCA_020536085_1_with_introns_edited.gff",header=FALSE)
gff_Amil <- read.delim("data/genomesToMatch/Acropora_millepora_GCF_013753865_1/Acropora_millepora_GCF_013753865_1_with_introns.gff",skip=8,header=FALSE)
#gff_Cjar <- read.delim("data/genomesToMatch/Catalaphyllia_jardinei_GCA_022496165_2/Catalaphyllia_jardinei_GCA_022496165_2_with_introns_edited.gff",header=FALSE)
gff_Nvec <- read.delim("data/genomesToMatch/Nematostella_vectensis_GCF_932526225_1/Nematostella_vectensis_GCF_932526225_1_with_introns.gff",skip=8,header=FALSE)
  #Ahya and Cjar aren't really annotated to the gene description level

# v1. All 215 loci used in the phylogeny ----
uce_inphylogeny_loose_gff <- 
  left_join(uce_inphylogeny_loose %>% select(UCE),
            Acropora_millepora_G_chromo_genome_loose %>% select(UCE,`GFF type list`,Chr),by="UCE") %>% 
  mutate(genes_Amil=case_when(grepl("ID=",`GFF type list`) ~ 
                      case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ 
                                  paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                      TRUE ~ `GFF type list`)) %>% 
  rename(Chr_Amil=Chr.x) %>% select(-`GFF type list`)
uce_inphylogeny_loose_gff <- 
  left_join(uce_inphylogeny_loose_gff,
          Nematostella_vectensis_G_chromo_genome_loose %>% select(UCE,`GFF type list`,Chr),by="UCE") %>% 
  mutate(genes_Nvec=case_when(grepl("ID=",`GFF type list`) ~ 
                      case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ 
                                  paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                      TRUE ~ `GFF type list`)) %>% 
  rename(Chr_Nvec=Chr.x) %>% select(-`GFF type list`)
#This takes awhile: 
tmp_Amil <- filter(gff_Amil,grepl(paste(uce_inphylogeny_loose_gff$genes_Amil,collapse='|'),V9))
uce_inphylogeny_loose_gff_Amil <- tmp_Amil %>% select(V3:V5,V9) %>% 
  mutate(genes_Amil=case_when(grepl("Parent=",V9) ~ str_extract(V9,"Parent=(.*)\\;Dbxref.*",group=1),
                         TRUE ~ "remove"),
         annot_Amil=case_when(grepl("Parent=",V9) ~ str_extract(V9,"product=(.*)\\;transcript.*",group=1),
                         TRUE ~ "remove")) %>% 
  filter(genes_Amil!="remove") %>% select(genes_Amil,annot_Amil) %>% group_by(genes_Amil) %>% slice(1) %>% 
  mutate(annot_Amil=sub("%2C",",",annot_Amil)) 
uce_inphylogeny_loose_gff_annotated <- 
  left_join(uce_inphylogeny_loose_gff,uce_inphylogeny_loose_gff_Amil,by="genes_Amil")
#This takes awhile again: 
tmp_Nvec <- filter(gff_Nvec,grepl(paste(uce_inphylogeny_loose_gff$genes_Nvec,collapse='|'),V9))
uce_inphylogeny_loose_gff_Nvec <- tmp_Nvec %>% select(V3:V5,V9) %>% 
  mutate(genes_Nvec=case_when(grepl("Parent=",V9) ~ str_extract(V9,"Parent=(.*)\\;Dbxref.*",group=1),
                              TRUE ~ "remove"),
         annot_Nvec=case_when(grepl("Parent=",V9) ~ str_extract(V9,"product=(.*)\\;transcript.*",group=1),
                              TRUE ~ "remove")) %>% 
  filter(genes_Nvec!="remove") %>% select(genes_Nvec,annot_Nvec) %>% group_by(genes_Nvec) %>% slice(1) %>% 
  mutate(annot_Nvec=sub("%2C",",",annot_Nvec)) 
uce_inphylogeny_loose_gff_annotated <- 
  left_join(uce_inphylogeny_loose_gff_annotated,uce_inphylogeny_loose_gff_Nvec,by="genes_Nvec") %>% 
  mutate(UCE=paste0("uce-",UCE)) 

  
# v2. Look specifically in the selected area of AhyaChr1 ----
Ahya_Chr1cluster <- 
  left_join(synth_chromo1_loose_toplot %>% filter(V3>25000000 & V2=="Chr1_Ahya") %>% select(-V2) %>% rename(UCE=V1) %>% mutate(UCE=sub("_Ahya","",UCE)),
            Acropora_millepora_G_chromo_genome_loose %>% select(UCE,`GFF type list`), by="UCE")
Ahya_Chr1cluster = Ahya_Chr1cluster %>% 
  mutate(genes=case_when(grepl("ID=",`GFF type list`) ~ 
                           case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                     TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                         TRUE ~ `GFF type list`),
         genes=as.factor(genes)) %>% 
  select(-`GFF type list`)
#tmp <- filter(gff_Amil,grepl(paste(Ahya_Chr1cluster$genes,collapse='|'),V9)) only run once
Ahya_Chr1cluster_gff = tmp %>% select(V3:V5,V9) %>% 
  mutate(genes=case_when(grepl("Parent=",V9) ~ str_extract(V9,"Parent=(.*)\\;Dbxref.*",group=1),
                         TRUE ~ "remove"),
         annot=case_when(grepl("Parent=",V9) ~ str_extract(V9,"product=(.*)\\;transcript.*",group=1),
                         TRUE ~ "remove")) %>% 
  filter(genes!="remove") %>% select(genes,annot) %>%
  group_by(genes) %>% slice(1) %>% 
  mutate(annot=sub("%2C",",",annot)) 
left_join(Ahya_Chr1cluster,Ahya_Chr1cluster_gff,by="genes") %>% 
  mutate(UCE=paste0("uce-",UCE)) %>% rename(UCEstart=V3) %>% rename(UCEend=V4)

# Output ----
write_tsv(uce_inphylogeny_loose_gff_annotated,file="export/mappedUCE_inphylogeny_annotations.tsv")
write_tsv(left_join(Ahya_Chr1cluster,Ahya_Chr1cluster_gff,by="genes") %>% mutate(UCE=paste0("uce-",UCE)) %>% rename(UCEstart=V3) %>% rename(UCEend=V4), 
          file="export/AhyaChr1cluster_annotations.tsv")
