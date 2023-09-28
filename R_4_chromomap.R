# Script info ----
#Script to plot data wrangled outputs from mapping UCEs to genome
#Requirements: datawrangle.R, plot.R outputs
#Written with R-4.2.2
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: May 2023

## Load prerequisite ----
source(file="datawrangle.R") 
#source(file="plot.R") but it'll spit errors so better to run line by line
library(grid)
library(gridExtra)
library(showtext)
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE) #Change to FALSE to turn it off

#install.packages("chromoMap")
library(chromoMap)
packageDescription("chromoMap",fields="Version")
#Note: Legend aesthetics needs to be edited in source:
#.libPaths() [1]/chromoMap/htmlwidgets/lib/chromoMap-3.0/chromoMap-3.0.js, 
#go to 'glgnd' and edit "font-size"
#go to 'w=rng.length' and edit *20 (height of legend key)


# Map to chromosome ----

# default ----
## Acropora_hyacinthus_G ----
Acropora_hyacinthus_G_chromo = Acropora_hyacinthus_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start),
         Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Acropora_hyacinthus_G_chromo_genome_toplot = Acropora_hyacinthus_G_chromo_genome %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Acropora_hyacinthus_G_chromo,file="data/Acropora_hyacinthus_G_chromo.csv",row.names=FALSE)
write.csv(Acropora_hyacinthus_G_chromo_genome_toplot,file="data/Acropora_hyacinthus_G_chromo_genome.csv",row.names=FALSE)
Acropora_hyacinthus_G_chromo <- read.csv("data/Acropora_hyacinthus_G_chromo.csv",header=FALSE,skip=1)
Acropora_hyacinthus_G_chromo_genome_toplot <- read.csv("data/Acropora_hyacinthus_G_chromo_genome.csv",header=FALSE,skip=1)

chromoMap(list(Acropora_hyacinthus_G_chromo),list(Acropora_hyacinthus_G_chromo_genome_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical",data_colors = list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color=c("lightgrey"), text_font_size=11,
          title="UCEs on Acropora hyacinthus genome", title_font_size=12,
          legend=T,lg_x=165,lg_y=270,interactivity=F,export.options=T) #n_win.factor=2,labels=T,label_font=3,label_angle=-90,
#Save in viewer pane. 

## Acropora_millepora_G ----
Acropora_millepora_G_chromo = Acropora_millepora_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start),
         Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Acropora_millepora_G_chromo_genome_toplot = Acropora_millepora_G_chromo_genome %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Acropora_millepora_G_chromo,file="data/Acropora_millepora_G_chromo.csv",row.names=FALSE)
write.csv(Acropora_millepora_G_chromo_genome_toplot,file="data/Acropora_millepora_G_chromo_genome.csv",row.names=FALSE)
Acropora_millepora_G_chromo <- read.csv("data/Acropora_millepora_G_chromo.csv",header=FALSE,skip=1)
Acropora_millepora_G_chromo_genome_toplot <- read.csv("data/Acropora_millepora_G_chromo_genome.csv",header=FALSE,skip=1)

chromoMap(list(Acropora_millepora_G_chromo),list(Acropora_millepora_G_chromo_genome_toplot),
          data_based_color_map = T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type = "categorical", data_colors = list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color = c("lightgrey"),text_font_size=11,
          title="UCEs on Acropora millepora genome",title_font_size=12,
          legend=T,lg_x=200,lg_y=270,interactivity=F,export.options=T)

## Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo = Catalaphyllia_jardinei_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start),
         Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Catalaphyllia_jardinei_G_chromo_genome_toplot = Catalaphyllia_jardinei_G_chromo_genome %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Catalaphyllia_jardinei_G_chromo,file="data/Catalaphyllia_jardinei_G_chromo.csv",row.names=FALSE)
write.csv(Catalaphyllia_jardinei_G_chromo_genome_toplot,file="data/Catalaphyllia_jardinei_G_chromo_genome.csv",row.names=FALSE)
Catalaphyllia_jardinei_G_chromo <- read.csv("data/Catalaphyllia_jardinei_G_chromo.csv",header=FALSE,skip=1)
Catalaphyllia_jardinei_G_chromo_genome_toplot <- read.csv("data/Catalaphyllia_jardinei_G_chromo_genome.csv",header=FALSE,skip=1)

chromoMap(list(Catalaphyllia_jardinei_G_chromo),list(Catalaphyllia_jardinei_G_chromo_genome_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors = list(hex_values[c(5,4,2)]),discrete.domain=list("Exon","Intron","Intergenic"),
          chr_color=c("lightgrey"), text_font_size=11,
          title="UCEs on Catalaphyllia jardinei genome",title_font_size=12,
          legend=T,lg_x=200,lg_y=270,interactivity=F,export.options=T)

## Montipora_capitata_C ----
Montipora_capitata_C_chromo = Montipora_capitata_C_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start),
         Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Montipora_capitata_C_chromo_genome_toplot = Montipora_capitata_C_chromo_genome %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"))
write.csv(Montipora_capitata_C_chromo,file="data/Montipora_capitata_C_chromo.csv",row.names=FALSE)
write.csv(Montipora_capitata_C_chromo_genome_toplot,file="data/Montipora_capitata_C_chromo_genome.csv",row.names=FALSE)
Montipora_capitata_C_chromo <- read.csv("data/Montipora_capitata_C_chromo.csv",header=FALSE,skip=1)
Montipora_capitata_C_chromo_genome_toplot <- read.csv("data/Montipora_capitata_C_chromo_genome.csv",header=FALSE,skip=1) %>% 
  mutate(V5=fct_relevel(V5,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))

chromoMap(list(Montipora_capitata_C_chromo),list(Montipora_capitata_C_chromo_genome_toplot),
          data_based_color_map = T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type = "categorical", data_colors = list(rev(hex_values[c(1,2,5)])),discrete.domain=list("Exon","Intergenic","Unclassified"),
          chr_color = c("lightgrey"),text_font_size=11,
          title="UCEs on Montipora capitata genome",title_font_size=12,
          legend=T,lg_x = 190,lg_y =260,interactivity=F,export.options=T)
#Save in viewer pane

## Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo = Nematostella_vectensis_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start),
         Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Nematostella_vectensis_G_chromo_genome_toplot = Nematostella_vectensis_G_chromo_genome %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Nematostella_vectensis_G_chromo,file="data/Nematostella_vectensis_G_chromo.csv",row.names=FALSE)
write.csv(Nematostella_vectensis_G_chromo_genome_toplot,file="data/Nematostella_vectensis_G_chromo_genome.csv",row.names=FALSE)
Nematostella_vectensis_G_chromo <- read.csv("data/Nematostella_vectensis_G_chromo.csv",header=FALSE,skip=1)
Nematostella_vectensis_G_chromo_genome_toplot <- read.csv("data/Nematostella_vectensis_G_chromo_genome.csv",header=FALSE,skip=1)

chromoMap(list(Nematostella_vectensis_G_chromo),list(Nematostella_vectensis_G_chromo_genome_toplot),
          data_based_color_map=T,chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color=c("lightgrey"),text_font_size=11, 
          title="UCEs on Nematostella vectensis genome",title_font_size=12,
          legend=T,lg_x=155,lg_y=270,interactivity=F,export.options=T)

###---###---###---###---###---###---###---###---###---###---###---###---###

# loose ----
## Acropora_hyacinthus_G ----
Acropora_hyacinthus_G_chromo_genome_loose_toplot = Acropora_hyacinthus_G_chromo_genome_loose %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Acropora_hyacinthus_G_chromo_genome_loose_toplot,file="data/Acropora_hyacinthus_G_chromo_genome_loose.csv",row.names=FALSE)
Acropora_hyacinthus_G_chromo_genome_loose_toplot <- read.csv("data/Acropora_hyacinthus_G_chromo_genome_loose.csv",header=FALSE,skip=1)

chromoMap(list(Acropora_hyacinthus_G_chromo),list(Acropora_hyacinthus_G_chromo_genome_loose_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical",data_colors = list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color=c("lightgrey"), text_font_size=11,
          title="UCEs on Acropora hyacinthus genome", title_font_size=12,
          legend=T,lg_x=165,lg_y=270,interactivity=F,export.options=T)

## Acropora_millepora_G ----
Acropora_millepora_G_chromo_genome_loose_toplot = Acropora_millepora_G_chromo_genome_loose %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Acropora_millepora_G_chromo_genome_loose_toplot,file="data/Acropora_millepora_G_chromo_genome_loose.csv",row.names=FALSE)
Acropora_millepora_G_chromo_genome_loose_toplot <- read.csv("data/Acropora_millepora_G_chromo_genome_loose.csv",header=FALSE,skip=1)

chromoMap(list(Acropora_millepora_G_chromo),list(Acropora_millepora_G_chromo_genome_loose_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color = c("lightgrey"),text_font_size=11,#chr_curve=5,
          title="UCEs on Acropora millepora genome",title_font_size=12,
          legend=T,lg_x=200,lg_y=270,interactivity=F,export.options=T)

## Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo_genome_loose_toplot = Catalaphyllia_jardinei_G_chromo_genome_loose %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Catalaphyllia_jardinei_G_chromo_genome_loose_toplot,file="data/Catalaphyllia_jardinei_G_chromo_genome_loose.csv",row.names=FALSE)
Catalaphyllia_jardinei_G_chromo_genome_loose_toplot <- read.csv("data/Catalaphyllia_jardinei_G_chromo_genome_loose.csv",header=FALSE,skip=1)

chromoMap(list(Catalaphyllia_jardinei_G_chromo),list(Catalaphyllia_jardinei_G_chromo_genome_loose_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(hex_values[c(5,4,2)]),discrete.domain=list("Exon","Intron","Intergenic"),
          chr_color=c("lightgrey"), text_font_size=11,
          title="UCEs on Catalaphyllia jardinei genome",title_font_size=12,
          legend=T,lg_x=200,lg_y=270,interactivity=F,export.options=T)

## Montipora_capitata_C ----
Montipora_capitata_C_chromo_genome_loose_toplot = Montipora_capitata_C_chromo_genome_loose %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"))
write.csv(Montipora_capitata_C_chromo_genome_loose_toplot,file="data/Montipora_capitata_C_chromo_genome_loose.csv",row.names=FALSE)
Montipora_capitata_C_chromo_genome_loose_toplot <- read.csv("data/Montipora_capitata_C_chromo_genome_loose.csv",header=FALSE,skip=1) %>% 
  mutate(V5=fct_relevel(V5,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))

chromoMap(list(Montipora_capitata_C_chromo),list(Montipora_capitata_C_chromo_genome_loose_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[c(1,2,5)])),discrete.domain=list("Exon","Intergenic","Unclassified"),
          chr_color=c("lightgrey"),text_font_size=11,
          title="UCEs on Montipora capitata genome",title_font_size=12,
          legend=T,lg_x = 190,lg_y =260,interactivity=F,export.options=T)

## Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo_genome_loose_toplot = Nematostella_vectensis_G_chromo_genome_loose %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Nematostella_vectensis_G_chromo_genome_loose_toplot,file="data/Nematostella_vectensis_G_chromo_genome_loose.csv",row.names=FALSE)
Nematostella_vectensis_G_chromo_genome_loose_toplot <- read.csv("data/Nematostella_vectensis_G_chromo_genome_loose.csv",header=FALSE,skip=1)

chromoMap(list(Nematostella_vectensis_G_chromo),list(Nematostella_vectensis_G_chromo_genome_loose_toplot),
          data_based_color_map=T,chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color=c("lightgrey"),text_font_size=11, 
          title="UCEs on Nematostella vectensis genome",title_font_size=12,
          legend=T,lg_x=155,lg_y=270,interactivity=F,export.options=T)

#after saving all pngs in the viewer pane:
#Remove all white spaces
system("mogrify -trim export/chromoMap_*.png")
system("montage export/chromoMap_Ahyacinthus.png export/chromoMap_Amillepora.png export/chromoMap_Cjardinei.png export/chromomap_Nvectensis.png -tile 2x2 -geometry +20+20 export/mappedUCE_chromomap_All.png")  #export/chromoMap_Mcapitata.png 
system("montage export/chromoMap_Ahyacinthus_loose.png export/chromoMap_Amillepora_loose.png export/chromoMap_Cjardinei_loose.png export/chromomap_Nvectensis_loose.png -tile 2x2 -geometry +20+20 export/mappedUCE_chromomap_All_loose.png") #export/chromoMap_Mcapitata_loose.png 

# ~~~~~~~~~~~~~~~~~ ----

# Synthetic chromosome ----
#common_uce_chromo was just so we can visualise a potential pattern. 
#Now we know there are some patterns, we want to place 4 taxa into one page. We can do this in chromomap by artificially placing all 4 taxa as '4 chromosomes'

## Ahya Chr1 ----
synth_chromo1_loose_tmp <- common_uce_chromo_loose %>% filter(n>=3) %>% select(-Montipora_capitata_C) %>% 
  filter(Acropora_hyacinthus_G=="Chr1") %>% select(UCE)
synth_chromo1_loose_toplot <-  #Amil and Nvec filtered because of multiple chromosomes
  rbind(left_join(synth_chromo1_loose_tmp,Acropora_hyacinthus_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromo1_loose_tmp,Acropora_millepora_G_chromo_genome_loose %>% filter(Chr=="Chr8") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromo1_loose_tmp,Catalaphyllia_jardinei_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromo1_loose_tmp,Nematostella_vectensis_G_chromo_genome_loose %>% filter(Chr=="Chr9") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% 
  mutate(UCE=paste0(UCE,"_",genome),
         Chr=paste0(Chr,"_",genome)) %>% select(-genome) %>% 
  na.exclude(Chr) %>% droplevels()
synth_chromo1 <- 
  rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr1") %>% mutate(genome="Ahya"),
        Acropora_millepora_G_chromo %>% filter(V1=="Chr8") %>% mutate(genome="Amil"),
        Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr10") %>% mutate(genome="Cjar"),
        Nematostella_vectensis_G_chromo %>% filter(V1=="Chr9") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% 
  mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)

write.csv(synth_chromo1,file="data/synthetic_chromosome1.csv",row.names=FALSE)
write.csv(synth_chromo1_loose_toplot,file="data/synthetic_chromosome1_loose_toplot.csv",row.names=FALSE)
synth_chromo1 <- read.csv("data/synthetic_chromosome1.csv",header=FALSE,skip=1)
synth_chromo1_loose_toplot <- read.csv("data/synthetic_chromosome1_loose_toplot.csv",header=FALSE,skip=1)
write.csv(cbind(synth_chromo1_loose_tmp,1,synth_chromo1_loose_tmp,1),file="data/synthetic_chromosome1_links.csv",row.names=FALSE)
synth_chromo1_links <- 
  rbind(read.csv("data/synthetic_chromosome1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Amil"),V5=1),
        read.csv("data/synthetic_chromosome1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Nvec"),V5=1))
        
chromoMap(list(synth_chromo1[c(3,1,2,4),]), 
          list(synth_chromo1_loose_toplot), # %>% filter(V2!="Chr14_Cjar")
          left_margin=80,n_win.factor=2,
          show.links=T, loci_links=synth_chromo1_links,
          chr_color=c("lightgrey"),ch_gap=40,#legend=T,labels=T,label_font=7,
          interactivity=T,export.options=T) #links.colors doesn't work...

## Ahya Chr14 ----
synth_chromo14_loose_tmp <- common_uce_chromo_loose %>% filter(n>=3) %>% select(-Montipora_capitata_C) %>% 
  filter(Acropora_hyacinthus_G=="Chr14") %>% select(UCE)
synth_chromo14_loose_toplot <- #Amil and Cjar filtered because of multiple chromosomes
  rbind(left_join(synth_chromo14_loose_tmp,Acropora_hyacinthus_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromo14_loose_tmp,Acropora_millepora_G_chromo_genome_loose %>% filter(Chr=="Chr11") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromo14_loose_tmp,Catalaphyllia_jardinei_G_chromo_genome_loose %>% filter(Chr!="others") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromo14_loose_tmp,Nematostella_vectensis_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% 
  mutate(UCE=paste0(UCE,"_",genome),
         Chr=paste0(Chr,"_",genome)) %>% select(-genome) %>% 
  na.exclude(Chr) %>% droplevels()
synth_chromo14 <- 
  rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Ahya"),
        Acropora_millepora_G_chromo %>% filter(V1=="Chr11") %>% mutate(genome="Amil"),
        Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Cjar"),
        Nematostella_vectensis_G_chromo %>% filter(V1=="Chr2") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% 
  mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)

write.csv(synth_chromo14,file="data/synthetic_chromosome14.csv",row.names=FALSE)
write.csv(synth_chromo14_loose_toplot,file="data/synthetic_chromosome14_loose_toplot.csv",row.names=FALSE)
synth_chromo14 <- read.csv("data/synthetic_chromosome14.csv",header=FALSE,skip=1)
synth_chromo14_loose_toplot <- read.csv("data/synthetic_chromosome14_loose_toplot.csv",header=FALSE,skip=1)

write.csv(cbind(synth_chromo14_loose_tmp,1,synth_chromo14_loose_tmp,1),file="data/synthetic_chromosome14_links.csv",row.names=FALSE)

synth_chromo14_links <- 
  rbind(read.csv("data/synthetic_chromosome14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Amil"),V5=1),
        read.csv("data/synthetic_chromosome14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Nvec"),V5=1))
chromoMap(list(synth_chromo14[c(3,1,2,4),]), 
          list(synth_chromo14_loose_toplot), # %>% filter(V2!="Chr14_Cjar")
          left_margin=80,n_win.factor=2,
          show.links=T, loci_links=synth_chromo14_links,
          chr_color=c("lightgrey"),ch_gap=50,#legend=T,labels=T,label_font=7,
          interactivity=F,export.options=T) #links.colors doesn't work...


## Use uces that were used in the i50 phylogeny instead of common ones (similar anyway) ----
### Ahya Chr1 ----
synth_chromo1_loose_tmp <- uce_inphylogeny_loose %>% select(-Montipora_capitata_C) %>% 
  filter(Acropora_hyacinthus_G=="Chr1") %>% select(UCE)
synth_chromo1_loose_toplot <-  #Amil and Nvec filtered because of multiple chromosomes
  rbind(left_join(synth_chromo1_loose_tmp,Acropora_hyacinthus_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromo1_loose_tmp,Acropora_millepora_G_chromo_genome_loose %>% filter(Chr=="Chr8") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromo1_loose_tmp,Catalaphyllia_jardinei_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromo1_loose_tmp,Nematostella_vectensis_G_chromo_genome_loose %>% filter(Chr=="Chr9") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% 
  mutate(UCE=paste0(UCE,"_",genome),
         Chr=paste0(Chr,"_",genome)) %>% select(-genome) %>% 
  na.exclude(Chr) %>% droplevels()
synth_chromo1 <- 
  rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr1") %>% mutate(genome="Ahya"),
        Acropora_millepora_G_chromo %>% filter(V1=="Chr8") %>% mutate(genome="Amil"),
        Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr10") %>% mutate(genome="Cjar"),
        Nematostella_vectensis_G_chromo %>% filter(V1=="Chr9") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% 
  mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)

write.csv(synth_chromo1,file="data/synthetic_chromosome1.csv",row.names=FALSE)
write.csv(synth_chromo1_loose_toplot,file="data/synthetic_chromosome1_loose_toplot.csv",row.names=FALSE)
synth_chromo1 <- read.csv("data/synthetic_chromosome1.csv",header=FALSE,skip=1)
synth_chromo1_loose_toplot <- read.csv("data/synthetic_chromosome1_loose_toplot.csv",header=FALSE,skip=1)
write.csv(cbind(synth_chromo1_loose_tmp,1,synth_chromo1_loose_tmp,1),file="data/synthetic_chromosome1_links.csv",row.names=FALSE)
synth_chromo1_links <- 
  rbind(read.csv("data/synthetic_chromosome1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Amil"),V5=1),
        read.csv("data/synthetic_chromosome1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Nvec"),V5=1))

chromoMap(list(synth_chromo1[c(3,1,2,4),]), 
          list(synth_chromo1_loose_toplot), # %>% filter(V2!="Chr14_Cjar")
          left_margin=80,n_win.factor=2,
          show.links=T, loci_links=synth_chromo1_links,
          chr_color=c("lightgrey"),ch_gap=40,#legend=T,labels=T,label_font=7,
          interactivity=T,export.options=T) #links.colors doesn't work...

### Ahya Chr14 ----
synth_chromo14_loose_tmp <- uce_inphylogeny_loose %>% select(-Montipora_capitata_C) %>% 
  filter(Acropora_hyacinthus_G=="Chr14") %>% select(UCE)
synth_chromo14_loose_toplot <- #Amil and Cjar filtered because of multiple chromosomes
  rbind(left_join(synth_chromo14_loose_tmp,Acropora_hyacinthus_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromo14_loose_tmp,Acropora_millepora_G_chromo_genome_loose %>% filter(Chr=="Chr11") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromo14_loose_tmp,Catalaphyllia_jardinei_G_chromo_genome_loose %>% filter(Chr!="others") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromo14_loose_tmp,Nematostella_vectensis_G_chromo_genome_loose %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% 
  mutate(UCE=paste0(UCE,"_",genome),
         Chr=paste0(Chr,"_",genome)) %>% select(-genome) %>% 
  na.exclude(Chr) %>% droplevels()
synth_chromo14 <- 
  rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Ahya"),
        Acropora_millepora_G_chromo %>% filter(V1=="Chr11") %>% mutate(genome="Amil"),
        Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Cjar"),
        Nematostella_vectensis_G_chromo %>% filter(V1=="Chr2") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% 
  mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)

write.csv(synth_chromo14,file="data/synthetic_chromosome14.csv",row.names=FALSE)
write.csv(synth_chromo14_loose_toplot,file="data/synthetic_chromosome14_loose_toplot.csv",row.names=FALSE)
synth_chromo14 <- read.csv("data/synthetic_chromosome14.csv",header=FALSE,skip=1)
synth_chromo14_loose_toplot <- read.csv("data/synthetic_chromosome14_loose_toplot.csv",header=FALSE,skip=1)

write.csv(cbind(synth_chromo14_loose_tmp,1,synth_chromo14_loose_tmp,1),file="data/synthetic_chromosome14_links.csv",row.names=FALSE)

synth_chromo14_links <- 
  rbind(read.csv("data/synthetic_chromosome14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Amil"),V5=1),
        read.csv("data/synthetic_chromosome14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Nvec"),V5=1))
chromoMap(list(synth_chromo14[c(3,1,2,4),]), 
          list(synth_chromo14_loose_toplot), # %>% filter(V2!="Chr14_Cjar")
          left_margin=80,n_win.factor=2,
          show.links=T, loci_links=synth_chromo14_links,
          chr_color=c("lightgrey"),ch_gap=40,#legend=T,labels=T,label_font=7,
          interactivity=F,export.options=T) #links.colors doesn't work...

# After saving all pngs in the viewer pane ----
#Remove all white spaces
system("mogrify -trim -border 10 -bordercolor none export/chromoMap_*.png")
system("montage export/chromoMap_Ahyacinthus.png export/chromoMap_Amillepora.png export/chromoMap_Cjardinei.png export/chromomap_Nvectensis.png -tile 2x2 -geometry +20+20 export/mappedUCE_chromomap_All.png")  #export/chromoMap_Mcapitata.png 
system("montage export/chromoMap_Ahyacinthus_loose.png export/chromoMap_Amillepora_loose.png export/chromoMap_Cjardinei_loose.png export/chromomap_Nvectensis_loose.png -tile 2x2 -geometry +20+20 export/mappedUCE_chromomap_All_loose.png") #export/chromoMap_Mcapitata_loose.png 
