# Script info ----
#Script to plot trees with genomes used for mapping UCEs
#file type: .tree or .treefile, metadata in .csv if you want to further rename tips
#Written with R version 4.2.2 (2022-10-31)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: August 2023

## Load necessary libraries ----
library(tidyverse)
library(ggtree)
library(treeio)
#package:ape needed for root()
packageDescription("ggtree", fields="Version")  #Check version

library(showtext)
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE) #Change to FALSE to turn it off

# Data prep ----
## Read files ----
tree <- read.tree("data/siteconcord.ML.i50.map.215loci.cf.tree")  #IQ-TREE output with site concordance factors
#info <- read.csv("../metadata_Acro10.csv")

#Run below once
#system('mkdir export')

# Root tree ----
tree_rooted <- tree
MRCA(tree_rooted,c("Actinia_tenebrosa_GCF_009602425_1","Nematostella_vectensis_GCF_932526225_1"))
tree_rooted <- ape::root(tree_rooted,node=41,resolve.root=TRUE)   #Actiniaria
#Add info to the plot
treedf <- as_tibble(tree_rooted)
treedf = treedf %>% #left_join(treedf, info, by='label') %>% 
  mutate(UfBS=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub("\\/.*","",label)),
         gCF=case_when(grepl("C",label) ~ "NA",
                       TRUE ~ gsub("\\/.*","",str_extract(label,"([^/]+)(?:/[^/]+){1}$"))),
         sCF=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub(".*\\/","",label)),
         UfBS=as.numeric(UfBS),
         gCf=as.numeric(gCF),
         sCF=as.numeric(sCF)) #IGNORE WARNING
tree_rooted <- treeio::as.treedata(treedf) #convert it back from tibble to tree object
#Plot with node suppport values and PSH labels
MRCA(tree_rooted,c("Acropora_millepora_GCF_013753865_1","Desmophyllum_pertusum_GCA_029204205_1"))
MRCA(tree_rooted,c("Discosoma_sp_ReefGenomics_v1","Amplexidiscus_fenestrafer_ReefGenomics_v1"))
MRCA(tree_rooted,c("Actinia_tenebrosa_GCF_009602425_1","Nematostella_vectensis_GCF_932526225_1"))
g <- ggtree(tree_rooted)+
  geom_tiplab(size=8,linetype=NULL)+
  xlim(0,2.7)+
  geom_treescale(x=0,y=20,fontsize=10)+
  geom_cladelab(node=39,label=" Scleractinia",offset=0.85,fontface=2,angle=90,vjust=1.5,hjust=0.5,
                barcolour='navy',barsize=2,fontsize=12)+
  geom_cladelab(node=41,label=" Corallimorpharia",offset=1.2,fontface=2,angle=90,vjust=1.5,hjust=0.5,
                barcolour='darkgreen',barsize=2,fontsize=12)+
  geom_cladelab(node=59,label=" Actiniaria",offset=1.3,fontface=2,angle=90,vjust=1.5,hjust=0.5,
                barcolour='salmon',barsize=2,fontsize=12)+
  geom_nodepoint(aes(subset=UfBS>=95),fill="white",size=4,shape=21)+
  geom_nodepoint(aes(subset=sCF >= 34),fill="orange",size=2,shape=21)+ #UfBS >= 95 & 
  geom_nodepoint(aes(subset=sCF >= 50),fill="white",size=2,shape=21)+ #UfBS == 100 & 
  #geom_nodepoint(aes(subset=sCF < 34),fill="red",size=2,shape=21)+
  theme(text=element_text(family="work",size=20))
#Find the long branch
#which(grepl(max(g$data$branch.length),g$data$branch.length))
#g$data[g$data$node %in% c(1),"x"] <- 0.01
#g = g %>% 
#  collapse(205,mode="none") %>% 
#  collapse(241,mode="none")
ggsave("export/mapgenome_tree_i50_215loci.png",plot=g,
       device="png",width=8,height=6,units="in",dpi=300)
ggtree::get_taxa_name(g)
