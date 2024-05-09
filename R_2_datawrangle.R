# Script info ----
#Script to data wrangle the outputs from mapping UCEs to genome
#Requirements: *_uce_type_summary.txt, *chromoinfo.txt, taxon.tsv (Order Taxon identifiers(like GCF_####))
 #Optional: hexatrans.tsv (list of UCE loci that are from transcriptomics), infosites.csv (informative sites output from phylogeny.sh)
#Written with R version 4.3.3 (2024-02-29)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: April 2024

## Load necessary libraries ----
library(tidyverse)
#Run below once
#system('mkdir export')

## Read all files ----
  #Note: Excluding Octocorallias that didn't return any results
  #Hash 'hexatrans' and 'infosites' if you don't have those files. Note you will have to also comment out further down the script
hexatrans <- read.delim("../Cowman etal 2020_phylogenies/Cowman_etal_Hexa_v2_PROBE_SETS/hexatrans.tsv",header=FALSE)
taxon <- read.delim("data/taxon.tsv",header=TRUE) 
infosites <- read.delim("data/infosites-i50_id70.csv",header=TRUE,sep=',') 
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
infosites = infosites %>% mutate(percent_variation_infositesDIVlength=informative_sites/length) %>% 
  select(locus,percent_variation_infositesDIVlength) %>% 
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

#~~~~~~~~~~~~~~~~ ----

# Repeat for relaxed filtering----

## Read all files ----
uce_type_summaries_relaxed <- list.files(path="data/genomesToMatch_extreme/",pattern="*uce_type_summary.txt",recursive=TRUE)

all_relaxed <- lapply(uce_type_summaries_relaxed, function(i){
  i <- paste0("data/genomesToMatch_extreme/",i)
  read_delim(i,col_names=FALSE,show_col_types=FALSE,delim="\t")
})
names(all_relaxed) <- sub("([A-Za-z]+_[a-z]+_[A-Z]).*","\\1",uce_type_summaries_relaxed)

## Data wrangle ----
df_all_relaxed <- do.call(rbind.data.frame, all_relaxed) %>% mutate(genome=sub("\\..*","",row.names(.)))
colnames(df_all_relaxed) <- c("UCE","Scaf","UCEstart","Type","Distance","GFF type list","genome")
df_all_relaxed = left_join(df_all_relaxed,hexatrans,by="UCE") %>% 
  mutate(origin=case_when(is.na(origin) ~ "genome",  TRUE ~ origin), UCE=sub("uce-","",UCE))
df_all_relaxed = left_join(df_all_relaxed,taxon, by="genome")
#Quality control (remove row with multiple match)
rowsToRemove <- df_all_relaxed %>% group_by(UCE,Scaf,genome) %>% add_tally() %>% filter(n>1) %>% slice(-1) %>% select(-n) %>% arrange(UCE,Scaf,genome) %>% ungroup()
df_all_relaxed <- anti_join(df_all_relaxed,rowsToRemove)
#"Distance is from previous UCE on scaffold" - but when you group_by scaffold, this is not always the case. So add Corrected_Distance column
df_all_relaxed = df_all_relaxed %>% 
  arrange(UCEstart) %>% group_by(Scaf,genome,Order) %>% 
  mutate(corrected_Distance=UCEstart-lag(UCEstart,default=first(UCEstart))) %>% ungroup() %>% 
  mutate(corrected_Distance=case_when(corrected_Distance==0 ~ as.numeric("NA"), TRUE ~ corrected_Distance)) %>% 
  arrange(Order,origin,genome,UCE)

#Group by scaffolds
all_relaxed_scaf <- df_all_relaxed %>% 
  group_by(Scaf,genome,Order) %>% add_tally() %>% ungroup() %>% rename(n.UCEonScaf=n)
#Extract close loci and the neighbours
  #The idea is: if one locus has Distance<1100, the next neighbour is also 'close' to another locus. So for each locus with distance<1100, at least a second locus needs to be extracted.
close_loci_relaxed_groups <- rbind(  #loci groups that COULD be close (but may not be for some taxa)
  df_all_relaxed %>% arrange(genome,Scaf,UCEstart) |> filter(lead(corrected_Distance<=1100)),
  df_all_relaxed %>% arrange(genome,Scaf,UCEstart) %>% filter(corrected_Distance <= 1100),
  df_all_relaxed %>% arrange(genome,Scaf,UCEstart) |> filter(lag(corrected_Distance<=1100))) %>% 
  arrange(Order,genome,Scaf,UCEstart) %>% select(-`GFF type list`) %>% unique() %>% 
  mutate(neighbourGroups=case_when(corrected_Distance<=1100 ~ paste0(lag(UCE),",",UCE), TRUE ~ "")) %>% 
  select(neighbourGroups) %>% unique() %>% separate(neighbourGroups,into=c("UCE1","UCE2"),sep=",") %>% filter(!is.na(UCE2)) %>% 
  mutate(neighbourGroups_new=case_when(as.numeric(UCE1) < as.numeric(UCE2) ~ paste0(UCE1,", ",UCE2), TRUE ~ paste0(UCE2,", ",UCE1))) %>% 
  select(neighbourGroups_new) %>% unique() %>%  #---up to here you have a list of loci pairs, but some loci may appear multiple times
  mutate(UCE1=str_extract(neighbourGroups_new,"(.*), (.*)",group=1),
         UCE2=str_extract(neighbourGroups_new,"(.*), (.*)",group=2)) %>% pivot_longer(cols=UCE1:UCE2,values_to='UCE') %>% select(-name) %>% 
  group_by(UCE) %>% add_tally() %>% ungroup()
#There is one group with 4loci grouping: 116142, 3911419, 9901734, 9901736
dput(which(grepl("116142|3911419|9901734|9901736", close_loci_relaxed_groups$neighbourGroups_new)))
close_loci_relaxed_groups_n4=close_loci_relaxed_groups[c(79,80,81,82,245,246),] %>% select(UCE) %>% 
  unique() %>% arrange(UCE) %>% mutate(neighbourGroups=paste0(UCE,collapse=", ")) %>% relocate(neighbourGroups,.before=UCE)
#Now sort out 3loci groups
close_loci_relaxed_groups_n3 = close_loci_relaxed_groups[-c(79,80,81,82,245,246),] %>% filter(n>1) %>% arrange(as.numeric(UCE)) %>%
  mutate(UCE1=str_extract(neighbourGroups_new,"(.*), (.*)",group=1),UCE2=str_extract(neighbourGroups_new,"(.*), (.*)",group=2)) %>% 
  mutate(neighbourGroups=case_when(UCE1==lead(UCE1) ~ paste0(UCE1,", ",UCE2,", ",lead(UCE2)),
                                   UCE2==lead(UCE2) ~ paste0(UCE1,", ",UCE2,", ",lead(UCE1)),
                                   UCE1==lead(UCE2) ~ paste0(lead(UCE1),", ",UCE1,", ",UCE2))) %>% 
  select(neighbourGroups) %>% filter(!is.na(neighbourGroups)) %>% unique() %>% 
  mutate(UCE1=str_extract(neighbourGroups,"(.*), (.*), (.*)",group=1),
         UCE2=str_extract(neighbourGroups,"(.*), (.*), (.*)",group=2),
         UCE3=str_extract(neighbourGroups,"(.*), (.*), (.*)",group=3)) %>% pivot_longer(cols=UCE1:UCE3,values_to='UCE') %>% select(-name)
#Now sort out pairs
close_loci_relaxed_groups_n2 = close_loci_relaxed_groups %>% filter(!(UCE %in% close_loci_relaxed_groups_n4$UCE)) %>% filter(!(UCE %in% close_loci_relaxed_groups_n3$UCE)) %>% 
  rename(neighbourGroups=neighbourGroups_new) %>% select(-n)
#Combine
close_loci_relaxed_groups = rbind(close_loci_relaxed_groups_n2,close_loci_relaxed_groups_n3,close_loci_relaxed_groups_n4)
#Get all distance info for these loci to see how many taxa are truly close
close_loci_relaxed = left_join(close_loci_relaxed_groups,df_all_relaxed,by="UCE") %>% select(-`GFF type list`) %>% 
  mutate(neighbourGroups=paste0("uce-",gsub(", ",", uce-",neighbourGroups))) %>% arrange(Order,genome,Scaf,UCEstart) %>% 
  group_by(Order,genome,Scaf,neighbourGroups) %>% add_tally() %>% filter(n!=1) %>% ungroup() %>% select(-n)
  #Then manually go through

#Group by genes
all_relaxed_genes <- df_all_relaxed %>% 
  mutate(genes=case_when(grepl("ID=",`GFF type list`) ~ 
                           case_when(grepl("ID=",str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)) ~ paste0(str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=1),";",str_extract(`GFF type list`,"ID=(.*)\\)(.*)ID=(.*)\\)(.*)",group=3)),
                                     TRUE ~ str_extract(`GFF type list`,"ID=(.*)\\)(.*)",group=1)),
                         TRUE ~ `GFF type list`),
         genes=as.factor(genes)) %>% 
  group_by(Scaf,genome,Order,genes) %>% add_tally() %>% ungroup() %>% rename(n.UCEonGene=n)

#Group by identified loci types
all_relaxed_type <- df_all_relaxed %>% 
  group_by(Type,genome,Order) %>% summarise(n=n()) %>% ungroup() %>% 
  group_by(genome) %>% 
  mutate(percent=n/sum(n),
         percent=round(percent,digits=3),
         Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified"))) %>% 
  ungroup()
#Group by UCE (to find out commonly found UCE loci - for genomes with chromoinfo only)
common_uce_chromo_relaxed <- df_all_relaxed %>% 
  filter(genome=="Acropora_hyacinthus_G"| genome=="Acropora_millepora_G" | genome=="Catalaphyllia_jardinei_G" |
           genome=="Montipora_capitata_C" | genome=="Nematostella_vectensis_G") %>% droplevels() %>% 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% droplevels() #filter(n >= 4) %>%
#Group by UCE (to find out commonly found UCE loci - all 29 genomes. tried a few, n>=15 is a good cutoff for exploratory purposes
common_uce_all_relaxed <- df_all_relaxed %>% 
  group_by(UCE) %>% summarise(n=n()) %>% ungroup() %>% filter(n >= 15) %>% droplevels() 

#Take transcriptome based baits into consideration
all_relaxed_type_trans <- df_all_relaxed %>% 
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
Acropora_hyacinthus_G_chromo_genome_relaxed <- 
  left_join(df_all_relaxed %>% filter(genome=="Acropora_hyacinthus_G"),
            Acropora_hyacinthus_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_hyacinthus_G_chromo_genome_relaxed = Acropora_hyacinthus_G_chromo_genome_relaxed %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_hyacinthus_G_chromo[match(Acropora_hyacinthus_G_chromo_genome_relaxed$Chr, Acropora_hyacinthus_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_relaxed <- 
  left_join(common_uce_chromo_relaxed,Acropora_hyacinthus_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)
uce_inphylogeny_relaxed <- left_join(infosites,Acropora_hyacinthus_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_hyacinthus_G=Chr)

### Acropora_millepora_G ----
Acropora_millepora_G_chromo_genome_relaxed <- 
  left_join(df_all_relaxed %>% filter(genome=="Acropora_millepora_G"),
            Acropora_millepora_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Acropora_millepora_G_chromo_genome_relaxed = Acropora_millepora_G_chromo_genome_relaxed %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Acropora_millepora_G_chromo[match(Acropora_millepora_G_chromo_genome_relaxed$Chr, Acropora_millepora_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_relaxed <- 
  left_join(common_uce_chromo_relaxed,Acropora_millepora_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)
uce_inphylogeny_relaxed <- 
  left_join(uce_inphylogeny_relaxed,Acropora_millepora_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Acropora_millepora_G=Chr)

### Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo_genome_relaxed <- 
  left_join(df_all_relaxed %>% filter(genome=="Catalaphyllia_jardinei_G"),
            Catalaphyllia_jardinei_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Catalaphyllia_jardinei_G_chromo_genome_relaxed = Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Catalaphyllia_jardinei_G_chromo[match(Catalaphyllia_jardinei_G_chromo_genome_relaxed$Chr, Catalaphyllia_jardinei_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_relaxed <- 
  left_join(common_uce_chromo_relaxed,Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)
uce_inphylogeny_relaxed <- 
  left_join(uce_inphylogeny_relaxed,Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Catalaphyllia_jardinei_G=Chr)

### Montipora_capitata_C ----
Montipora_capitata_C_chromo_genome_relaxed <- 
  left_join(df_all_relaxed %>% filter(genome=="Montipora_capitata_C"),
            Montipora_capitata_C_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Montipora_capitata_C_chromo_genome_relaxed = Montipora_capitata_C_chromo_genome_relaxed %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Montipora_capitata_C_chromo[match(Montipora_capitata_C_chromo_genome_relaxed$Chr, Montipora_capitata_C_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_relaxed <- 
  left_join(common_uce_chromo_relaxed,Montipora_capitata_C_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)
uce_inphylogeny_relaxed <- 
  left_join(uce_inphylogeny_relaxed,Montipora_capitata_C_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Montipora_capitata_C=Chr)

### Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo_genome_relaxed <- 
  left_join(df_all_relaxed %>% filter(genome=="Nematostella_vectensis_G"),
            Nematostella_vectensis_G_chromo %>% select(c(Chr,Scaf)), by="Scaf")
Nematostella_vectensis_G_chromo_genome_relaxed = Nematostella_vectensis_G_chromo_genome_relaxed %>% 
  mutate(UCEend=UCEstart+120,
         Chr = case_when(is.na(Chr)==FALSE ~ Chr,  TRUE ~ 'others'),
         max = unlist(Nematostella_vectensis_G_chromo[match(Nematostella_vectensis_G_chromo_genome_relaxed$Chr, Nematostella_vectensis_G_chromo$Chr),3]),
         remove=case_when(max == '0' ~ 'yes',  UCEend > max ~ 'yes',  TRUE ~ 'no')) %>% 
  relocate(UCEend, .after=UCEstart)
common_uce_chromo_relaxed <- 
  left_join(common_uce_chromo_relaxed,Nematostella_vectensis_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)
uce_inphylogeny_relaxed <- 
  left_join(uce_inphylogeny_relaxed,Nematostella_vectensis_G_chromo_genome_relaxed %>% select(UCE,Chr),by="UCE") %>% 
  rename(Nematostella_vectensis_G=Chr)

###---###---###---###---###---###---###---###---###---###---###---###---###

# Output ----
#csv file with each type and percentage 
write.csv(all_type %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome.csv",row.names=FALSE)
write.csv(all_type_trans %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome_transcriptomeSeparated.csv",row.names=FALSE)
write.csv(all_relaxed_type %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome_relaxed.csv",row.names=FALSE)
write.csv(all_relaxed_type_trans %>% pivot_wider(names_from="Type",values_from=c("n","percent")), 
          file="export/uce_type_on_genome_relaxed_transcriptomeSeparated.csv",row.names=FALSE)

#csv file with UCE locus less than 1100bp apart consistently across genomes
write.csv(close_loci_relaxed %>% mutate(Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic")),
          file="export/uce_distances_relaxed_consistentlyclose.csv",row.names=FALSE)

#csv file of nearly-raw input
write.csv(df_all %>% mutate(Type=as.factor(replace_na(Type,"Unclassified")),
                            Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic")),
          file="export/uce_mapresults_all.csv",row.names=FALSE)
write.csv(df_all_relaxed %>% mutate(Type=as.factor(replace_na(Type,"Unclassified")),
                                  Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic")),
          file="export/uce_mapresults_all_relaxed.csv",row.names=FALSE)

#csv file with number of loci mapped to each chromosome
chromos_count <- rbind(
  Acropora_hyacinthus_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_hyacinthus_G") %>% arrange(Chr),
  Acropora_millepora_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_millepora_G") %>% arrange(Chr),
  Catalaphyllia_jardinei_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Catalaphyllia_jardinei_G") %>% arrange(Chr),
  Montipora_capitata_C_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Montipora_capitata_C") %>% arrange(Chr),
  Nematostella_vectensis_G_chromo_genome %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:15),"others")),genome="Nematostella_vectensis_G") %>% arrange(Chr))
write.csv(chromos_count,file="export/chromosome_uce_loci_count.csv",row.names=FALSE)
chromos_count_relaxed <- rbind(
  Acropora_hyacinthus_G_chromo_genome_relaxed %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_hyacinthus_G") %>% arrange(Chr),
  Acropora_millepora_G_chromo_genome_relaxed %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Acropora_millepora_G") %>% arrange(Chr),
  Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Catalaphyllia_jardinei_G") %>% arrange(Chr),
  Montipora_capitata_C_chromo_genome_relaxed %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:14),"others")),genome="Montipora_capitata_C") %>% arrange(Chr),
  Nematostella_vectensis_G_chromo_genome_relaxed %>% group_by(Chr) %>% summarise(n=n()) %>% ungroup() %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:15),"others")),genome="Nematostella_vectensis_G") %>% arrange(Chr))
write.csv(chromos_count_relaxed,file="export/chromosome_uce_loci_count_relaxed.csv",row.names=FALSE)

uce_inphylogeny_count <- uce_inphylogeny %>% select(-percent_variation_infositesDIVlength) %>% 
  pivot_longer(cols=Acropora_hyacinthus_G:Nematostella_vectensis_G,names_to='genome',values_to='Chr') %>% 
  group_by(genome,Chr) %>% summarise(n=n()) %>% ungroup() %>% filter(Chr!="NA") %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:15),"others"))) %>% arrange(genome,Chr)
write.csv(uce_inphylogeny_count,file="export/chromosome_uce_loci_inphylogeny_count.csv",row.names=FALSE)
uce_inphylogeny_relaxed_count <- uce_inphylogeny_relaxed %>% select(-percent_variation_infositesDIVlength) %>% 
  pivot_longer(cols=Acropora_hyacinthus_G:Nematostella_vectensis_G,names_to='genome',values_to='Chr') %>% 
  group_by(genome,Chr) %>% summarise(n=n()) %>% ungroup() %>% filter(Chr!="NA") %>% mutate(Chr=factor(Chr,levels=c(paste0("Chr",1:15),"others"))) %>% arrange(genome,Chr)
write.csv(uce_inphylogeny_relaxed_count,file="export/chromosome_uce_loci_inphylogeny_count_relaxed.csv",row.names=FALSE)

