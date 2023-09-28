# Script info ----
#Script to data wrangle the outputs from mapping UCEs to genome
#Requirements: *_uce_type_summary.txt, *chromoinfo.txt, taxon.tsv (Order Taxon identifiers(like GCF_####))
 #Optional: hexatrans.tsv (list of UCE loci that are from transcriptomics), infosites.csv (informative sites output from phylogeny.sh)
#Written with R-4.2.2
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: August 2023

## Load necessary libraries ----
library(tidyverse)
#Run below once
#system('mkdir export')

## Read all files ----
  #Excluding Octocorallias that didn't return any results
  #Unhash if you also have those files
#hexatrans <- read.delim("../Cowman etal 2020_phylogenies/Cowman_etal_Hexa_v2_PROBE_SETS/hexatrans.tsv",header=FALSE)
taxon <- read.delim("data/taxon.tsv",header=TRUE) 
#infosites <- read.delim("data/infosites-i50.csv",header=TRUE,sep=',') 
uce_type_summaries <- list.files(path="data/genomesToMatch/",pattern="*uce_type_summary.txt",recursive=TRUE)
chromo_infos <- list.files(path="data/genomesToMatch/",pattern="*chromoinfo.txt",recursive=TRUE)

all <- lapply(uce_type_summaries, function(i){
  i <- paste0("data/genomesToMatch/",i)
  read_delim(i,col_names=FALSE,show_col_types=FALSE,delim="\t")
})
names(all) <- sub("([A-Za-z]+_[a-z]+_[A-Z]).*","\\1",uce_type_summaries)

for (i in chromo_infos){
  genome <- sub("([A-Za-z]+_[a-z]+_[A-Z]).*","\\1",i)
  genome <- paste0(genome,"_chromo")
  i <- paste0("data/genomesToMatch/",i)
  assign(genome,read_delim(i,col_names=FALSE,show_col_types=FALSE))
}

# Data wrangle ----
colnames(hexatrans) <- "UCE"
hexatrans = hexatrans %>% mutate(origin='transcriptome')
infosites = infosites %>% select(locus,percent_variation_infositesDIVlength) %>% 
  rename(UCE=locus) %>% mutate(UCE=sub("uce-","",sub(".nexus","",UCE)))

taxon = taxon %>% 
  mutate(NCBI.RefSeq...other.identifiers=substr(NCBI.RefSeq...other.identifiers,1,1),
         genome=paste0(sub(" ","_",Taxon),"_",NCBI.RefSeq...other.identifiers)) %>% 
  select(Order,genome)

df_all <- do.call(rbind.data.frame,all) %>% 
  mutate(genome=sub("\\..*","",row.names(.)))
colnames(df_all) <- c("UCE","Scaf","UCEstart","Type","Distance","GFF type list","genome")
df_all = left_join(df_all,hexatrans,by="UCE") %>% 
  mutate(origin=case_when(is.na(origin) ~ "genome",  TRUE ~ origin),
         UCE=sub("uce-","",UCE))
  #"Distance is from previous UCE on scaffold" - but when you group_by scaffold, this is not always the case. So add Corrected_Distance column
df_all = left_join(df_all,taxon, by="genome") %>% 
  arrange(UCEstart) %>% group_by(Scaf,genome,Order) %>% 
  mutate(corrected_Distance=UCEstart-lag(UCEstart,default=first(UCEstart))) %>% ungroup() %>% 
  mutate(corrected_Distance=case_when(corrected_Distance==0 ~ as.numeric("NA"),
                                      TRUE ~ corrected_Distance)) %>% 
  arrange(Order,origin,genome,UCE)
#Quality control (remove row with double UCE match)
rowsToRemove <- df_all %>% group_by(UCE,Scaf,genome,origin,Order) %>% add_tally() %>% filter(n>1) %>% filter(corrected_Distance!="NA") %>% select(-n) %>% ungroup()
df_all <- anti_join(df_all,rowsToRemove)

#Group by scaffolds
all_scaf <- df_all %>% 
  group_by(Scaf,genome,Order) %>% add_tally() %>% ungroup() %>% rename(n.UCEonScaf=n)

#Group by genes
all_genes <- df_all %>% 
  mutate(genes=case_when(grepl("ID=",`GFF type list`) ~ 
                           case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                     TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                         TRUE ~ `GFF type list`),
         genes=as.factor(genes)) %>% 
  group_by(Scaf,genome,Order,genes) %>% add_tally() %>% ungroup() %>% rename(n.UCEonGene=n)

#Group by identified loci types
all_type <- df_all %>% 
  group_by(Type,genome,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup()

#Group by UCE (to find out commonly found UCE loci - for genomes with chromoinfo only)
common_uce_chromo <- df_all %>% 
  filter(genome=="Acropora_hyacinthus_G"| genome=="Acropora_millepora_G" | 
           genome=="Catalaphyllia_jardinei_G" | genome=="Nematostella_vectensis_G") %>% droplevels() %>%  #genome=="Montipora_capitata_C" | 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% droplevels()

#Group by UCE (to find out commonly found UCE loci - all 29 genomes. tried a few, n>=15 is a good cutoff for exploratory purposes
common_uce_all <- df_all %>% 
  group_by(UCE) %>% summarise(n=n())  %>% ungroup() %>% filter(n >= 15) %>% droplevels() 

#Take transcriptome based baits into consideration
all_type_trans <- df_all %>% 
  group_by(Type,genome,origin,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n,na.rm=TRUE),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup()

## Clean up chromo df for chromoMap ----
#chromosome name,start,end (but we also need the scaffold name for later so keeping that as the 4th column)
#chromosome name: ‘chr1’ or ‘1’ or ‘ch1’
#chromosome start: a numeric value to specify chromosome start position. Typically '1'.
#chromosome end: a numeric value specifying chromosome/contig/region end position. Typically length of chromosome.
Acropora_hyacinthus_G_chromo = Acropora_hyacinthus_G_chromo %>% 
  mutate(X1=paste0("Chr",X1),  X4=1, X2=sub("\\.1","",X2)) %>% 
  relocate(X4,.before=X3) %>% relocate(X2,.after=X3)
colnames(Acropora_hyacinthus_G_chromo) <- c("Chr","Chr_start","Chr_end","Scaf")
Catalaphyllia_jardinei_G_chromo = Catalaphyllia_jardinei_G_chromo %>% 
mutate(X1=paste0("Chr",X1),  X4=1,  X2=sub("\\.1","",X2)) %>% 
  relocate(X4,.before=X3) %>% relocate(X2,.after=X3)
colnames(Catalaphyllia_jardinei_G_chromo) <- c("Chr","Chr_start","Chr_end","Scaf")
Montipora_capitata_C_chromo = Montipora_capitata_C_chromo %>% 
  mutate(X4=1) %>% 
  relocate(X4,.before=X2)
colnames(Montipora_capitata_C_chromo) <- c("Chr","Chr_start","Chr_end","Scaf")
Acropora_millepora_G_chromo = Acropora_millepora_G_chromo %>% 
  mutate(X3=sub("chromosome","Chr",X3), X4=1) %>% 
  relocate(X3,.before=X1) %>% relocate(X1,.after=X4) %>% relocate(X2,.after=X4)
colnames(Acropora_millepora_G_chromo) <- c("Chr","Chr_start","Chr_end","Scaf")
Nematostella_vectensis_G_chromo = Nematostella_vectensis_G_chromo %>% 
mutate(X3=sub("chromosome","Chr",X3), X4=1) %>% 
  relocate(X3,.before=X1) %>% relocate(X1,.after=X4) %>% relocate(X2,.after=X4)
colnames(Nematostella_vectensis_G_chromo) <- c("Chr","Chr_start","Chr_end","Scaf")


# Genome with chromosome info ----
## Acropora_hyacinthus_G ----
Acropora_hyacinthus_G_chromo_genome <- 
  left_join(df_all %>% filter(genome=="Acropora_hyacinthus_G"),
            Acropora_hyacinthus_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_hyacinthus_G_chromo_genome = Acropora_hyacinthus_G_chromo_genome %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_hyacinthus_G_chromo[match(Acropora_hyacinthus_G_chromo_genome$Chr, Acropora_hyacinthus_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo <- left_join(common_uce_chromo,Acropora_hyacinthus_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)
uce_inphylogeny <- left_join(infosites,Acropora_hyacinthus_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)

## Acropora_millepora_G ----
Acropora_millepora_G_chromo_genome <- 
  left_join(df_all %>% filter(genome=="Acropora_millepora_G"),
            Acropora_millepora_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_millepora_G_chromo_genome = Acropora_millepora_G_chromo_genome %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_millepora_G_chromo[match(Acropora_millepora_G_chromo_genome$Chr, Acropora_millepora_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo <- 
  left_join(common_uce_chromo,Acropora_millepora_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)
uce_inphylogeny <- 
  left_join(uce_inphylogeny,Acropora_millepora_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)

## Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo_genome <- 
  left_join(df_all %>% filter(genome=="Catalaphyllia_jardinei_G"),
            Catalaphyllia_jardinei_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Catalaphyllia_jardinei_G_chromo_genome = Catalaphyllia_jardinei_G_chromo_genome %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Catalaphyllia_jardinei_G_chromo[match(Catalaphyllia_jardinei_G_chromo_genome$Chr, Catalaphyllia_jardinei_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo <- 
  left_join(common_uce_chromo,Catalaphyllia_jardinei_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)
uce_inphylogeny <- 
  left_join(uce_inphylogeny,Catalaphyllia_jardinei_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)

## Montipora_capitata_C ----
Montipora_capitata_C_chromo_genome <- 
  left_join(df_all %>% filter(genome=="Montipora_capitata_C"),
            Montipora_capitata_C_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Montipora_capitata_C_chromo_genome = Montipora_capitata_C_chromo_genome %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Montipora_capitata_C_chromo[match(Montipora_capitata_C_chromo_genome$Chr, Montipora_capitata_C_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo <- 
  left_join(common_uce_chromo,Montipora_capitata_C_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)
uce_inphylogeny <- 
  left_join(uce_inphylogeny,Montipora_capitata_C_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)

## Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo_genome <- 
  left_join(df_all %>% filter(genome=="Nematostella_vectensis_G"),
            Nematostella_vectensis_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Nematostella_vectensis_G_chromo_genome = Nematostella_vectensis_G_chromo_genome %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Nematostella_vectensis_G_chromo[match(Nematostella_vectensis_G_chromo_genome$Chr, Nematostella_vectensis_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo <- 
  left_join(common_uce_chromo,Nematostella_vectensis_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)
uce_inphylogeny <- 
  left_join(uce_inphylogeny,Nematostella_vectensis_G_chromo_genome %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)

#~~~~~~~~~~ ----

# Repeat for loose ----

## Read all files ----
uce_type_summaries_loose <- list.files(path="data/genomesToMatch_loose/",pattern="*uce_type_summary.txt",recursive=TRUE)

all_loose <- lapply(uce_type_summaries_loose, function(i){
  i <- paste0("data/genomesToMatch_loose/",i)
  read_delim(i,col_names=FALSE,show_col_types=FALSE,delim="\t")
})
names(all_loose) <- sub("([A-Za-z]+_[a-z]+_[A-Z]).*","\\1",uce_type_summaries_loose)

## Data wrangle ----
df_all_loose <- do.call(rbind.data.frame, all_loose) %>% 
  mutate(genome=sub("\\..*","",row.names(.)))
colnames(df_all_loose) <- c("UCE","Scaf","UCEstart","Type","Distance","GFF type list","genome")
df_all_loose = left_join(df_all_loose,hexatrans,by="UCE") %>% 
  mutate(origin=case_when(is.na(origin) ~ "genome",  TRUE ~ origin),
         UCE=sub("uce-","",UCE))
  #"Distance is from previous UCE on scaffold" - but when you group_by scaffold, this is not always the case. So add Corrected_Distance column
df_all_loose = left_join(df_all_loose,taxon, by="genome") %>% 
  arrange(UCEstart) %>% group_by(Scaf,genome,Order) %>% 
  mutate(corrected_Distance=UCEstart-lag(UCEstart,default=first(UCEstart))) %>% ungroup() %>% 
  mutate(corrected_Distance=case_when(corrected_Distance==0 ~ as.numeric("NA"),
                                      TRUE ~ corrected_Distance)) %>% 
  arrange(Order,origin,genome,UCE)
#Quality control (remove row with multiple match)
rowsToRemove <- df_all_loose %>% group_by(UCE,Scaf,genome) %>% add_tally() %>% filter(n>1) %>% select(-n) %>% ungroup() %>% arrange(Scaf)
df_all_loose <- anti_join(df_all_loose,rowsToRemove[seq(2,nrow(rowsToRemove),2),])

#Group by scaffolds
all_loose_scaf <- df_all_loose %>% 
  group_by(Scaf,genome,Order) %>% add_tally() %>% ungroup() %>% rename(n.UCEonScaf=n)

#Group by genes
all_loose_genes <- df_all_loose %>% 
  mutate(genes=case_when(grepl("ID=",`GFF type list`) ~ 
                           case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                     TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                         TRUE ~ `GFF type list`),
         genes=as.factor(genes)) %>% 
  group_by(Scaf,genome,Order,genes) %>% add_tally() %>% ungroup() %>% rename(n.UCEonGene=n)

#Group by identified loci types
all_loose_type <- df_all_loose %>% 
  group_by(Type,genome,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup()
#Group by UCE (to find out commonly found UCE loci - for genomes with chromoinfo only)
common_uce_chromo_loose <- df_all_loose %>% 
  filter(genome=="Acropora_hyacinthus_G"| genome=="Acropora_millepora_G" | genome=="Catalaphyllia_jardinei_G" |
           genome=="Montipora_capitata_C" | genome=="Nematostella_vectensis_G") %>% droplevels() %>% 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% droplevels() #filter(n >= 4) %>%
#Group by UCE (to find out commonly found UCE loci - all 29 genomes. tried a few, n>=15 is a good cutoff for exploratory purposes
common_uce_all_loose <- df_all_loose %>% 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% filter(n >= 15) %>% droplevels() 

#Take transcriptome based baits into consideration
all_loose_type_trans <- df_all_loose %>% 
  group_by(Type,genome,origin,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n,na.rm=TRUE),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup() 

## Genome with chromosome info ----
### Acropora_hyacinthus_G ----
Acropora_hyacinthus_G_chromo_genome_loose <- 
  left_join(df_all_loose %>% filter(genome=="Acropora_hyacinthus_G"),
            Acropora_hyacinthus_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_hyacinthus_G_chromo_genome_loose = Acropora_hyacinthus_G_chromo_genome_loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_hyacinthus_G_chromo[match(Acropora_hyacinthus_G_chromo_genome_loose$Chr, Acropora_hyacinthus_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_loose <- 
  left_join(common_uce_chromo_loose,Acropora_hyacinthus_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)
uce_inphylogeny_loose <- left_join(infosites,Acropora_hyacinthus_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)

### Acropora_millepora_G ----
Acropora_millepora_G_chromo_genome_loose <- 
  left_join(df_all_loose %>% filter(genome=="Acropora_millepora_G"),
            Acropora_millepora_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_millepora_G_chromo_genome_loose = Acropora_millepora_G_chromo_genome_loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_millepora_G_chromo[match(Acropora_millepora_G_chromo_genome_loose$Chr, Acropora_millepora_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_loose <- 
  left_join(common_uce_chromo_loose,Acropora_millepora_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)
uce_inphylogeny_loose <- 
  left_join(uce_inphylogeny_loose,Acropora_millepora_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)

### Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo_genome_loose <- 
  left_join(df_all_loose %>% filter(genome=="Catalaphyllia_jardinei_G"),
            Catalaphyllia_jardinei_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Catalaphyllia_jardinei_G_chromo_genome_loose = Catalaphyllia_jardinei_G_chromo_genome_loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Catalaphyllia_jardinei_G_chromo[match(Catalaphyllia_jardinei_G_chromo_genome_loose$Chr, Catalaphyllia_jardinei_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_loose <- 
  left_join(common_uce_chromo_loose,Catalaphyllia_jardinei_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)
uce_inphylogeny_loose <- 
  left_join(uce_inphylogeny_loose,Catalaphyllia_jardinei_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)

### Montipora_capitata_C ----
Montipora_capitata_C_chromo_genome_loose <- 
  left_join(df_all_loose %>% filter(genome=="Montipora_capitata_C"),
            Montipora_capitata_C_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Montipora_capitata_C_chromo_genome_loose = Montipora_capitata_C_chromo_genome_loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Montipora_capitata_C_chromo[match(Montipora_capitata_C_chromo_genome_loose$Chr, Montipora_capitata_C_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_loose <- 
  left_join(common_uce_chromo_loose,Montipora_capitata_C_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)
uce_inphylogeny_loose <- 
  left_join(uce_inphylogeny_loose,Montipora_capitata_C_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)

### Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo_genome_loose <- 
  left_join(df_all_loose %>% filter(genome=="Nematostella_vectensis_G"),
            Nematostella_vectensis_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Nematostella_vectensis_G_chromo_genome_loose = Nematostella_vectensis_G_chromo_genome_loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Nematostella_vectensis_G_chromo[match(Nematostella_vectensis_G_chromo_genome_loose$Chr, Nematostella_vectensis_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_loose <- 
  left_join(common_uce_chromo_loose,Nematostella_vectensis_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)
uce_inphylogeny_loose <- 
  left_join(uce_inphylogeny_loose,Nematostella_vectensis_G_chromo_genome_loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)

###---###---###---###---###---###---###---###---###---###---###---###---###

# Repeat for C99 loose ----

## Read all files ----
uce_type_summaries_C99loose <- list.files(path="data/genomesToMatch_C99loose/",pattern="*uce_type_summary.txt",recursive=TRUE)

all_C99loose <- lapply(uce_type_summaries_C99loose, function(i){
  i <- paste0("data/genomesToMatch_C99loose/",i)
  read_delim(i,col_names=FALSE,show_col_types=FALSE,delim="\t")
})
names(all_C99loose) <- sub("([A-Za-z]+_[a-z]+_[A-Z]).*","\\1",uce_type_summaries_C99loose)

## Data wrangle ----
df_all_C99loose <- do.call(rbind.data.frame, all_C99loose) %>% 
  mutate(genome=sub("\\..*","",row.names(.)))
colnames(df_all_C99loose) <- c("UCE","Scaf","UCEstart","Type","Distance","GFF type list","genome")
df_all_C99loose = left_join(df_all_C99loose,hexatrans,by="UCE") %>% 
  mutate(origin=case_when(is.na(origin) ~ "genome",  TRUE ~ origin),
         UCE=sub("_C99","",sub("uce-","",UCE)))
#"Distance" is from "previous UCE on scaffold".
df_all_C99loose = left_join(df_all_C99loose,taxon, by="genome")

#Group by scaffolds
all_C99loose_scaf <- df_all_C99loose %>% 
  group_by(Scaf,genome,Order) %>% add_tally() %>% ungroup() %>% rename(n.UCEonScaf=n)
#Group by identified loci types
all_C99loose_type <- df_all_C99loose %>% 
  group_by(Type,genome,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup()
#Group by UCE (to find out commonly found UCE loci - for genomes with chromoinfo only)
common_uce_chromo_C99loose <- df_all_C99loose %>% 
  filter(genome=="Acropora_hyacinthus_G"| genome=="Acropora_millepora_G" | genome=="Catalaphyllia_jardinei_G" |
           genome=="Montipora_capitata_C" | genome=="Nematostella_vectensis_G") %>% droplevels() %>% 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% droplevels() #filter(n >= 4) %>%
#Group by UCE (to find out commonly found UCE loci - all 29 genomes. tried a few, n>=15 is a good cutoff for exploratory purposes
common_uce_all_C99loose <- df_all_C99loose %>% 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% filter(n >= 15) %>% droplevels() 

#Take transcriptome based baits into consideration
all_C99loose_type_trans <- df_all_C99loose %>% 
  group_by(Type,genome,origin,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup() 
  #Nothing mapped to transcriptome 

## Genome with chromosome info ----
### Acropora_hyacinthus_G ----
Acropora_hyacinthus_G_chromo_genome_C99loose <- 
  left_join(df_all_C99loose %>% filter(genome=="Acropora_hyacinthus_G"),
            Acropora_hyacinthus_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_hyacinthus_G_chromo_genome_C99loose = Acropora_hyacinthus_G_chromo_genome_C99loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_hyacinthus_G_chromo[match(Acropora_hyacinthus_G_chromo_genome_C99loose$Chr, Acropora_hyacinthus_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_C99loose <- 
  left_join(common_uce_chromo_C99loose,Acropora_hyacinthus_G_chromo_genome_C99loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)

### Acropora_millepora_G ----
Acropora_millepora_G_chromo_genome_C99loose <- 
  left_join(df_all_C99loose %>% filter(genome=="Acropora_millepora_G"),
            Acropora_millepora_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_millepora_G_chromo_genome_C99loose = Acropora_millepora_G_chromo_genome_C99loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_millepora_G_chromo[match(Acropora_millepora_G_chromo_genome_C99loose$Chr, Acropora_millepora_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_C99loose <- 
  left_join(common_uce_chromo_C99loose,Acropora_millepora_G_chromo_genome_C99loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)

### Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo_genome_C99loose <- 
  left_join(df_all_C99loose %>% filter(genome=="Catalaphyllia_jardinei_G"),
            Catalaphyllia_jardinei_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Catalaphyllia_jardinei_G_chromo_genome_C99loose = Catalaphyllia_jardinei_G_chromo_genome_C99loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Catalaphyllia_jardinei_G_chromo[match(Catalaphyllia_jardinei_G_chromo_genome_C99loose$Chr, Catalaphyllia_jardinei_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_C99loose <- 
  left_join(common_uce_chromo_C99loose,Catalaphyllia_jardinei_G_chromo_genome_C99loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)

### Montipora_capitata_C ----
Montipora_capitata_C_chromo_genome_C99loose <- 
  left_join(df_all_C99loose %>% filter(genome=="Montipora_capitata_C"),
            Montipora_capitata_C_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Montipora_capitata_C_chromo_genome_C99loose = Montipora_capitata_C_chromo_genome_C99loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Montipora_capitata_C_chromo[match(Montipora_capitata_C_chromo_genome_C99loose$Chr, Montipora_capitata_C_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_C99loose <- 
  left_join(common_uce_chromo_C99loose,Montipora_capitata_C_chromo_genome_C99loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)

### Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo_genome_C99loose <- 
  left_join(df_all_C99loose %>% filter(genome=="Nematostella_vectensis_G"),
            Nematostella_vectensis_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Nematostella_vectensis_G_chromo_genome_C99loose = Nematostella_vectensis_G_chromo_genome_C99loose %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Nematostella_vectensis_G_chromo[match(Nematostella_vectensis_G_chromo_genome_C99loose$Chr, Nematostella_vectensis_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_C99loose <- 
  left_join(common_uce_chromo_C99loose,Nematostella_vectensis_G_chromo_genome_C99loose %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)

###---###---###---###---###---###---###---###---###---###---###---###---###

# Output ----
#csv file with each type and percentage 
write.csv(all_type %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome.csv",row.names=FALSE)
write.csv(all_type_trans %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome_transcriptomeSeparated.csv",row.names=FALSE)
write.csv(all_loose_type %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome_loose.csv",row.names=FALSE)
write.csv(all_loose_type_trans %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome_loose_transcriptomeSeparated.csv",row.names=FALSE)
#csv file with UCE locus closer than 1000bp
write.csv(all_scaf %>% filter(corrected_Distance <= 1000),
          file="export/uce_distances_less1000.csv",row.names=FALSE)
write.csv(all_loose_scaf %>% filter(corrected_Distance <= 1000),
          file="export/uce_distances_loose_less1000.csv",row.names=FALSE)

#csv file with number of loci mapped to each chromosome
chromos_count <- rbind(
  Acropora_hyacinthus_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_hyacinthus_G") %>% arrange(Chr),
  Acropora_millepora_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_millepora_G") %>% arrange(Chr),
  Catalaphyllia_jardinei_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Catalaphyllia_jardinei_G") %>% arrange(Chr),
  Montipora_capitata_C_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Montipora_capitata_C") %>% arrange(Chr),
  Nematostella_vectensis_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:15),"others")),genome="Nematostella_vectensis_G") %>% arrange(Chr))
write.csv(chromos_count,file="export/chromosome_uce_loci_count.csv",row.names=FALSE)
chromos_count_loose <- rbind(
  Acropora_hyacinthus_G_chromo_genome_loose %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_hyacinthus_G") %>% arrange(Chr),
  Acropora_millepora_G_chromo_genome_loose %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_millepora_G") %>% arrange(Chr),
  Catalaphyllia_jardinei_G_chromo_genome_loose %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Catalaphyllia_jardinei_G") %>% arrange(Chr),
  Montipora_capitata_C_chromo_genome_loose %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Montipora_capitata_C") %>% arrange(Chr),
  Nematostella_vectensis_G_chromo_genome_loose %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:15),"others")),genome="Nematostella_vectensis_G") %>% arrange(Chr))
write.csv(chromos_count_loose,file="export/chromosome_uce_loci_count_loose.csv",row.names=FALSE)
