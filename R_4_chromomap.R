# Script info ----
#Script to plot data wrangled outputs from mapping UCEs to genome
#Requirements: datawrangle.R, plot.R outputs
#Written with R version 4.3.3 (2024-02-29)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: April 2024

## Load prerequisite ----
source(file="R_2_datawrangle.R") 
  #also R_3_plot.R, better to run line by line than to source() it
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
  #go to 'glgnd' and edit "font-size","11px"
  #go to 'w=rng.length' and edit *20 (height of legend key)-> *15


# Map to chromosome ----

# relaxed ----
## Acropora_hyacinthus_G ----
Acropora_hyacinthus_G_chromo = Acropora_hyacinthus_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start), Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Acropora_hyacinthus_G_chromo_genome_relaxed_toplot = Acropora_hyacinthus_G_chromo_genome_relaxed %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Acropora_hyacinthus_G_chromo,file="data/Acropora_hyacinthus_G_chromo.csv",row.names=FALSE)
write.csv(Acropora_hyacinthus_G_chromo_genome_relaxed_toplot,file="data/Acropora_hyacinthus_G_chromo_genome_relaxed.csv",row.names=FALSE)
Acropora_hyacinthus_G_chromo <- read.csv("data/Acropora_hyacinthus_G_chromo.csv",header=FALSE,skip=1)
Acropora_hyacinthus_G_chromo_genome_relaxed_toplot <- read.csv("data/Acropora_hyacinthus_G_chromo_genome_relaxed.csv",header=FALSE,skip=1)

chromoMap(list(Acropora_hyacinthus_G_chromo),list(Acropora_hyacinthus_G_chromo_genome_relaxed_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical",data_colors = list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color=c("lightgrey"), text_font_size=11,
          title="Acropora hyacinthus", title_font_size=12,
          legend=T,lg_x=165,lg_y=270,interactivity=F,export.options=T)
  #Save with bottom button of viewer pane : chromoMap_Ahyacinthus_relaxed (.png if using mogrify, and .svg if using AdobeIllustrator)

## Acropora_millepora_G ----
Acropora_millepora_G_chromo = Acropora_millepora_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start), Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Acropora_millepora_G_chromo_genome_relaxed_toplot = Acropora_millepora_G_chromo_genome_relaxed %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Acropora_millepora_G_chromo,file="data/Acropora_millepora_G_chromo.csv",row.names=FALSE)
write.csv(Acropora_millepora_G_chromo_genome_relaxed_toplot,file="data/Acropora_millepora_G_chromo_genome_relaxed.csv",row.names=FALSE)
Acropora_millepora_G_chromo <- read.csv("data/Acropora_millepora_G_chromo.csv",header=FALSE,skip=1)
Acropora_millepora_G_chromo_genome_relaxed_toplot <- read.csv("data/Acropora_millepora_G_chromo_genome_relaxed.csv",header=FALSE,skip=1)

chromoMap(list(Acropora_millepora_G_chromo),list(Acropora_millepora_G_chromo_genome_relaxed_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color = c("lightgrey"),text_font_size=11,#chr_curve=5,
          title="Acropora millepora",title_font_size=12,
          legend=T,lg_x=200,lg_y=270,interactivity=F,export.options=T)
  #chromoMap_Amillepora_relaxed

## Catalaphyllia_jardinei_G ----
Catalaphyllia_jardinei_G_chromo = Catalaphyllia_jardinei_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start), Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Catalaphyllia_jardinei_G_chromo_genome_relaxed_toplot = Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% 
  droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Catalaphyllia_jardinei_G_chromo,file="data/Catalaphyllia_jardinei_G_chromo.csv",row.names=FALSE)
write.csv(Catalaphyllia_jardinei_G_chromo_genome_relaxed_toplot,file="data/Catalaphyllia_jardinei_G_chromo_genome_relaxed.csv",row.names=FALSE)
Catalaphyllia_jardinei_G_chromo <- read.csv("data/Catalaphyllia_jardinei_G_chromo.csv",header=FALSE,skip=1)
Catalaphyllia_jardinei_G_chromo_genome_relaxed_toplot <- read.csv("data/Catalaphyllia_jardinei_G_chromo_genome_relaxed.csv",header=FALSE,skip=1)

chromoMap(list(Catalaphyllia_jardinei_G_chromo),list(Catalaphyllia_jardinei_G_chromo_genome_relaxed_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(hex_values[c(5,4,2)]),discrete.domain=list("Exon","Intron","Intergenic"),
          chr_color=c("lightgrey"), text_font_size=11,
          title="Catalaphyllia jardinei",title_font_size=12,
          legend=T,lg_x=200,lg_y=270,interactivity=F,export.options=T)
  #chromoMap_Cjardinei_relaxed

## Montipora_capitata_C ----
Montipora_capitata_C_chromo = Montipora_capitata_C_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start), Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Montipora_capitata_C_chromo_genome_relaxed_toplot = Montipora_capitata_C_chromo_genome_relaxed %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"))
write.csv(Montipora_capitata_C_chromo,file="data/Montipora_capitata_C_chromo.csv",row.names=FALSE)
write.csv(Montipora_capitata_C_chromo_genome_relaxed_toplot,file="data/Montipora_capitata_C_chromo_genome_relaxed.csv",row.names=FALSE)
Montipora_capitata_C_chromo <- read.csv("data/Montipora_capitata_C_chromo.csv",header=FALSE,skip=1)
Montipora_capitata_C_chromo_genome_relaxed_toplot <- read.csv("data/Montipora_capitata_C_chromo_genome_relaxed.csv",header=FALSE,skip=1) %>% 
  mutate(V5=fct_relevel(V5,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))

chromoMap(list(Montipora_capitata_C_chromo),list(Montipora_capitata_C_chromo_genome_relaxed_toplot),
          data_based_color_map=T, chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[c(1,2,5)])),discrete.domain=list("Exon","Intergenic","Unclassified"),
          chr_color=c("lightgrey"),text_font_size=11,
          title="Montipora capitata",title_font_size=12,
          legend=T,lg_x = 190,lg_y =260,interactivity=F,export.options=T)
  #chromoMap_Mcapitata_relaxed

## Nematostella_vectensis_G ----
Nematostella_vectensis_G_chromo = Nematostella_vectensis_G_chromo %>% 
  mutate(Chr_start=as.integer(Chr_start), Chr_end=as.integer(Chr_end)) %>% select(-Scaf)
Nematostella_vectensis_G_chromo_genome_relaxed_toplot = Nematostella_vectensis_G_chromo_genome_relaxed %>% 
  filter(Chr != 'others' & remove == 'no') %>% select(c(UCE,Chr,UCEstart,UCEend,Type)) %>% droplevels() %>% 
  mutate(Type=as.factor(replace_na(Type,"Unclassified")),
         Type=recode(Type,"E"="Exon","I"="Intron","EI"="Exon and Intron","N"="Intergenic"),
         Type=fct_relevel(Type,c("Exon","Intron","Exon and Intron","Intergenic","Unclassified")))
write.csv(Nematostella_vectensis_G_chromo,file="data/Nematostella_vectensis_G_chromo.csv",row.names=FALSE)
write.csv(Nematostella_vectensis_G_chromo_genome_relaxed_toplot,file="data/Nematostella_vectensis_G_chromo_genome_relaxed.csv",row.names=FALSE)
Nematostella_vectensis_G_chromo <- read.csv("data/Nematostella_vectensis_G_chromo.csv",header=FALSE,skip=1)
Nematostella_vectensis_G_chromo_genome_relaxed_toplot <- read.csv("data/Nematostella_vectensis_G_chromo_genome_relaxed.csv",header=FALSE,skip=1)

chromoMap(list(Nematostella_vectensis_G_chromo),list(Nematostella_vectensis_G_chromo_genome_relaxed_toplot),
          data_based_color_map=T,chr.scale.ticks=6,left_margin=40,ch_gap=1,
          data_type="categorical", data_colors=list(rev(hex_values[2:5])),discrete.domain=list("Exon","Intron","Exon and Intron","Intergenic"),
          chr_color=c("lightgrey"),text_font_size=11, 
          title="Nematostella vectensis",title_font_size=12,
          legend=T,lg_x=155,lg_y=270,interactivity=F,export.options=T)
  #chromoMap_Nvectensis_relaxed

#After saving all pngs in the viewer pane: ----
#Remove all white spaces
system("mogrify -trim export/chromoMap_*.png")
system("montage export/chromoMap_Ahyacinthus.png export/chromoMap_Amillepora.png export/chromoMap_Cjardinei.png export/chromomap_Nvectensis.png -tile 2x2 -geometry +20+20 export/mappedUCE_chromomap_All.png")  #export/chromoMap_Mcapitata.png 
#system("montage export/chromoMap_Ahyacinthus_relaxed.png export/chromoMap_Amillepora_relaxed.png export/chromoMap_Cjardinei_relaxed.png export/chromomap_Nvectensis_relaxed.png -tile 2x2 -geometry +20+20 export/mappedUCE_chromomap_All_relaxed.png") #export/chromoMap_Mcapitata_relaxed.png 
system("montage export/chromoMap_Ahyacinthus_relaxed.png export/chromoMap_Amillepora_relaxed.png export/chromoMap_Cjardinei_relaxed.png export/chromoMap_Mcapitata_relaxed.png export/chromomap_Nvectensis_relaxed.png -tile 2x3 -geometry +20+20 export/mappedUCE_chromomap_All_relaxed.png")

# ~~~~~~~~~~~~~~~~~ ----

# Synthetic chromosome ----
#common_uce_chromo was just so we can visualise a potential pattern. 
#Now we know there are some patterns, we want to place all taxa into one page. We can do this in chromomap by artificially placing all taxa as '4/5 chromosomes'
#Use uces that were used in the i50 phylogeny instead of common ones (similar anyway)

### Ahya Chr1 ----
synth_chromoAhya1_relaxed_tmp <- uce_inphylogeny_relaxed %>% filter(Acropora_hyacinthus_G=="Chr1") %>% select(UCE)
synth_chromoAhya1_toplot <-  #filter for the majority chromosome
  rbind(left_join(synth_chromoAhya1_relaxed_tmp,Acropora_hyacinthus_G_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromoAhya1_relaxed_tmp,Acropora_millepora_G_chromo_genome_relaxed %>% filter(Chr=="Chr8") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromoAhya1_relaxed_tmp,Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% filter(Chr=="Chr10") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromoAhya1_relaxed_tmp,Montipora_capitata_C_chromo_genome_relaxed %>% filter(Chr=="Chr14") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Mcap"),
        left_join(synth_chromoAhya1_relaxed_tmp,Nematostella_vectensis_G_chromo_genome_relaxed %>% filter(Chr=="Chr9") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% mutate(UCE=paste0(UCE,"_",genome), Chr=paste0(Chr,"_",genome)) %>% 
  select(-genome) %>% na.exclude(Chr) %>% droplevels()
synth_chromoAhya1 <- 
  rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr1") %>% mutate(genome="Ahya"),
        Acropora_millepora_G_chromo %>% filter(V1=="Chr8") %>% mutate(genome="Amil"),
        Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr10") %>% mutate(genome="Cjar"),
        Montipora_capitata_C_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Mcap"),
        Nematostella_vectensis_G_chromo %>% filter(V1=="Chr9") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)
write.csv(synth_chromoAhya1,file="data/synthetic_chromosome_Ahya1.csv",row.names=FALSE)
write.csv(synth_chromoAhya1_toplot,file="data/synthetic_chromosome_Ahya1_toplot.csv",row.names=FALSE)
synth_chromoAhya1 <- read.csv("data/synthetic_chromosome_Ahya1.csv",header=FALSE,skip=1)
synth_chromoAhya1_toplot <- read.csv("data/synthetic_chromosome_Ahya1_toplot.csv",header=FALSE,skip=1)
write.csv(cbind(synth_chromoAhya1_relaxed_tmp,1,synth_chromoAhya1_relaxed_tmp,1),file="data/synthetic_chromosome_Ahya1_links.csv",row.names=FALSE)
synth_chromoAhya1_links <- 
  rbind(read.csv("data/synthetic_chromosome_Ahya1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome_Ahya1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Nvec"),V5=1),
        read.csv("data/synthetic_chromosome_Ahya1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome_Ahya1_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Mcap"),V5=1))

chromoMap(list(synth_chromoAhya1[c(2,1,5),-4]), list(synth_chromoAhya1_toplot),
          left_margin=80,n_win.factor=2, show.links=T, loci_links=synth_chromoAhya1_links,
          chr_color=c("lightgrey"),ch_gap=30,links.colors=c('#27445C','black'),#for some reason you have to specify two colours to work (first colour actually goes into plot)
          legend=T,interactivity=T,export.options=T) #set legend=T to make it easier to remove
  #save: chromoMap_link_inphylogeny_AhyaChr1_a
system('mogrify -crop 2940x900 export/chromoMap_link_inphylogeny_AhyaChr1_a.png')
chromoMap(list(synth_chromoAhya1[c(3,1,4),-4]), list(synth_chromoAhya1_toplot),
          left_margin=80,n_win.factor=2, show.links=T, loci_links=synth_chromoAhya1_links,
          chr_color=c("lightgrey"),ch_gap=30,links.colors=c('#27445C','black'),
          legend=T,interactivity=T,export.options=T) 
  #save: chromoMap_link_inphylogeny_AhyaChr1_b
system('mogrify -crop 2940x900 export/chromoMap_link_inphylogeny_AhyaChr1_b.png')
system('rm export/chromoMap_link_inphylogeny_AhyaChr1_*-1.png')


### Ahya Chr14 ----
synth_chromoAhya14_relaxed_tmp <- uce_inphylogeny_relaxed %>% filter(Acropora_hyacinthus_G=="Chr14") %>% select(UCE)
synth_chromoAhya14_toplot <-  #filter for the majority chromosome
  rbind(left_join(synth_chromoAhya14_relaxed_tmp,Acropora_hyacinthus_G_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromoAhya14_relaxed_tmp,Acropora_millepora_G_chromo_genome_relaxed %>% filter(Chr=="Chr11") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromoAhya14_relaxed_tmp,Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromoAhya14_relaxed_tmp,Montipora_capitata_C_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Mcap"),
        left_join(synth_chromoAhya14_relaxed_tmp,Nematostella_vectensis_G_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% mutate(UCE=paste0(UCE,"_",genome), Chr=paste0(Chr,"_",genome)) %>% 
  select(-genome) %>% na.exclude(Chr) %>% droplevels()
synth_chromoAhya14 <- 
  rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Ahya"),
        Acropora_millepora_G_chromo %>% filter(V1=="Chr11") %>% mutate(genome="Amil"),
        Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr14") %>% mutate(genome="Cjar"),
        Montipora_capitata_C_chromo %>% filter(V1=="Chr7") %>% mutate(genome="Mcap"),
        Nematostella_vectensis_G_chromo %>% filter(V1=="Chr2") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)
write.csv(synth_chromoAhya14,file="data/synthetic_chromosome_Ahya14.csv",row.names=FALSE)
write.csv(synth_chromoAhya14_toplot,file="data/synthetic_chromosome_Ahya14_toplot.csv",row.names=FALSE)
synth_chromoAhya14 <- read.csv("data/synthetic_chromosome_Ahya14.csv",header=FALSE,skip=1)
synth_chromoAhya14_toplot <- read.csv("data/synthetic_chromosome_Ahya14_toplot.csv",header=FALSE,skip=1)
write.csv(cbind(synth_chromoAhya14_relaxed_tmp,1,synth_chromoAhya14_relaxed_tmp,1),file="data/synthetic_chromosome_Ahya14_links.csv",row.names=FALSE)
synth_chromoAhya14_links <- 
  rbind(read.csv("data/synthetic_chromosome_Ahya14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome_Ahya14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Nvec"),V5=1),
        read.csv("data/synthetic_chromosome_Ahya14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome_Ahya14_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Ahya"),V3=paste0(V3,"_Mcap"),V5=1))

chromoMap(list(synth_chromoAhya14[c(2,1,5),-4]), list(synth_chromoAhya14_toplot),
          left_margin=80,n_win.factor=2, show.links=T, loci_links=synth_chromoAhya14_links,
          chr_color=c("lightgrey"),ch_gap=30,links.colors=c('#27445C','black'),#for some reason you have to specify two colours to work (first colour actually goes into plot)
          legend=T,interactivity=T,export.options=T) #set legend=T to make it easier to remove
  #save: chromoMap_link_inphylogeny_AhyaChr14_a
system('mogrify -crop 2940x900 export/chromoMap_link_inphylogeny_AhyaChr14_a.png')
chromoMap(list(synth_chromoAhya14[c(3,1,4),-4]), list(synth_chromoAhya14_toplot),
          left_margin=80,n_win.factor=2, show.links=T, loci_links=synth_chromoAhya14_links,
          chr_color=c("lightgrey"),ch_gap=30,links.colors=c('#27445C','black'),
          legend=T,interactivity=T,export.options=T) 
  #save: chromoMap_link_inphylogeny_AhyaChr14_b
system('mogrify -crop 2940x900 export/chromoMap_link_inphylogeny_AhyaChr14_b.png')
system('rm export/chromoMap_link_inphylogeny_AhyaChr14_*-1.png')



## Nvec Chr13 ----
synth_chromoNvec13_relaxed_tmp <- uce_inphylogeny_relaxed %>% filter(Nematostella_vectensis_G=="Chr13") %>% select(UCE)
synth_chromoNvec13_toplot <-  #filter for the majority chromosome
  rbind(left_join(synth_chromoNvec13_relaxed_tmp,Acropora_hyacinthus_G_chromo_genome_relaxed %>% filter(Chr=="Chr9") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Ahya"),
        left_join(synth_chromoNvec13_relaxed_tmp,Acropora_millepora_G_chromo_genome_relaxed %>% filter(Chr=="Chr4") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Amil"),
        left_join(synth_chromoNvec13_relaxed_tmp,Catalaphyllia_jardinei_G_chromo_genome_relaxed %>% filter(Chr=="Chr1") %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Cjar"),
        left_join(synth_chromoNvec13_relaxed_tmp,Montipora_capitata_C_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Mcap"),
        left_join(synth_chromoNvec13_relaxed_tmp,Nematostella_vectensis_G_chromo_genome_relaxed %>% select(UCE,UCEstart,UCEend,Chr),by="UCE") %>% mutate(genome="Nvec")) %>% 
  relocate(Chr,.before="UCEstart") %>% mutate(UCE=paste0(UCE,"_",genome), Chr=paste0(Chr,"_",genome)) %>% 
  select(-genome) %>% na.exclude(Chr) %>% droplevels()
synth_chromoNvec13 <- rbind(Acropora_hyacinthus_G_chromo %>% filter(V1=="Chr9") %>% mutate(genome="Ahya"),
                            Acropora_millepora_G_chromo %>% filter(V1=="Chr4") %>% mutate(genome="Amil"),
                            Catalaphyllia_jardinei_G_chromo %>% filter(V1=="Chr1") %>% mutate(genome="Cjar"),
                            Montipora_capitata_C_chromo %>% filter(V1=="Chr11") %>% mutate(genome="Mcap"),
                            Nematostella_vectensis_G_chromo %>% filter(V1=="Chr13") %>% mutate(genome="Nvec")) %>% 
  arrange(genome) %>% mutate(V1=paste0(V1,"_",genome)) %>% select(-genome)
write.csv(synth_chromoNvec13_toplot,file="data/synthetic_chromosome_Nvec13_toplot.csv",row.names=FALSE)
synth_chromoNvec13_toplot <- read.csv("data/synthetic_chromosome_Nvec13_toplot.csv",header=FALSE,skip=1)
write.csv(cbind(synth_chromoNvec13_relaxed_tmp,1,synth_chromoNvec13_relaxed_tmp,1),file="data/synthetic_chromosome_Nvec13_links.csv",row.names=FALSE)
synth_chromoNvec13_links <- 
  rbind(read.csv("data/synthetic_chromosome_Nvec13_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Amil"),V3=paste0(V3,"_Nvec"),V5=1),
        read.csv("data/synthetic_chromosome_Nvec13_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Nvec"),V3=paste0(V3,"_Ahya"),V5=1),
        read.csv("data/synthetic_chromosome_Nvec13_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Cjar"),V3=paste0(V3,"_Nvec"),V5=1),
        read.csv("data/synthetic_chromosome_Nvec13_links.csv",header=FALSE,skip=1) %>% mutate(V1=paste0(V1,"_Nvec"),V3=paste0(V3,"_Mcap"),V5=1))
chromoMap(list(synth_chromoNvec13[c(1,5,2),-4]), list(synth_chromoNvec13_toplot),
          left_margin=80,n_win.factor=2, show.links=T, loci_links=synth_chromoNvec13_links,
          chr_color=c("lightgrey"),ch_gap=30,links.colors=c('#27445C','black'),#for some reason you have to specify two colours to work (first colour actually goes into plot)
          legend=T,interactivity=T,export.options=T) #set legend=T to make it easier to remove
  #save: chromoMap_link_inphylogeny_NvecChr13_a
system('mogrify -crop 2940x900 export/chromoMap_link_inphylogeny_NvecChr13_a.png')
chromoMap(list(synth_chromoNvec13[c(3,5,4),-4]), list(synth_chromoNvec13_toplot),
          left_margin=80,n_win.factor=2, show.links=T, loci_links=synth_chromoNvec13_links,
          chr_color=c("lightgrey"),ch_gap=30,links.colors=c('#27445C','black'),
          legend=T,interactivity=T,export.options=T) 
  #save: chromoMap_link_inphylogeny_NvecChr13_b
system('mogrify -crop 2940x900 export/chromoMap_link_inphylogeny_NvecChr13_b.png')
system('rm export/chromoMap_link_inphylogeny_NvecChr13_*-1.png')


# After saving all pngs in the viewer pane ----
system("mogrify -gravity NorthWest -strokewidth 5 -pointsize 40 -annotate +10+10 'A' export/chromoMap_link_inphylogeny_AhyaChr1_a-0.png")
system("mogrify -gravity NorthWest -strokewidth 5 -pointsize 40 -annotate +10+10 'B' export/chromoMap_link_inphylogeny_AhyaChr14_a-0.png")
system("montage export/chromoMap_link_inphylogeny_AhyaChr1_a-0.png export/chromoMap_link_inphylogeny_AhyaChr14_a-0.png -tile 1x2 -geometry +20+20 export/mappedUCE_chromomap_link_inphylogeny_AhyaChr1-14_AmilAhyaNvec.png")

system("mogrify -gravity NorthWest -strokewidth 5 -pointsize 40 -annotate +10+10 'A' export/chromoMap_link_inphylogeny_AhyaChr1_b-0.png")
system("mogrify -gravity NorthWest -strokewidth 5 -pointsize 40 -annotate +10+10 'B' export/chromoMap_link_inphylogeny_AhyaChr14_b-0.png")
system("montage export/chromoMap_link_inphylogeny_AhyaChr1_b-0.png export/chromoMap_link_inphylogeny_AhyaChr14_b-0.png -tile 1x2 -geometry +20+20 export/mappedUCE_chromomap_link_inphylogeny_AhyaChr1-14_CjarAhyaMcap.png")
