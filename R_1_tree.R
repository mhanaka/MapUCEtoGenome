# Script info ----
#Script to plot trees with genomes used for mapping UCEs
#file type: .tree or .treefile, (optional) metadata in .csv
#Written with R version 4.2.2 (2022-10-31)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: November 2023

## Load necessary libraries ----
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
packageDescription("ggtree", fields="Version")  #Check version

library(ggtext)
library(showtext)
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE) #Change to FALSE to turn it off

# Data prep ----
## Read files ----
tree <- read.tree("data/siteconcord.i50.gen.map.215loci.cf.tree")  #IQ-TREE output with site concordance factors
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
  mutate(SHalrt=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub("\\/.*","",label)),
         UfBS=case_when(grepl("C",label) ~ "NA",
                       TRUE ~ gsub("\\/.*","",str_extract(label,"([^/]+)(?:/[^/]+){1}$"))),
         sCF=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub(".*\\/","",label)),
         SHalrt=as.numeric(SHalrt),
         UfBS=as.numeric(UfBS),
         sCF=as.numeric(sCF),
         UfBS_cat=factor(case_when(UfBS>=95 ~ 'UfBS ≥ 95',  TRUE ~ 'UfBS < 95')),
         sCF_cat=factor(case_when(sCF>=34 & sCF<50 ~ 'sCF ≥ 34',  sCF>=50 ~ 'sCF ≥ 50',  TRUE ~ 'sCF < 34')),
         inProbe=factor(case_when(grepl("digitifera|Nematostella",label) ~ 'Both',
                           grepl("Discosoma|Exaiptasia|Montastraea",label) ~ 'Genome',
                           grepl("Orbicella|damicornis",label) ~ 'Transcriptome',
                           TRUE ~ 'Not used'),levels=c('Genome','Transcriptome','Both','Not used'))) #IGNORE WARNING
tree_rooted <- treeio::as.treedata(treedf) #convert it back from tibble to tree object
#Plot with node suppport values and PSH labels
MRCA(tree_rooted,c("Acropora_millepora_GCF_013753865_1","Desmophyllum_pertusum_GCA_029204205_1")) #all Scleractinia
MRCA(tree_rooted,c("Pocillopora_meandrina_Cyanophora_v1","Desmophyllum_pertusum_GCA_029204205_1")) #Robust
MRCA(tree_rooted,c("Acropora_millepora_GCF_013753865_1","Porites_evermanni_GCA_942486025_1")) #Complex
MRCA(tree_rooted,c("Discosoma_sp_ReefGenomics_v1","Amplexidiscus_fenestrafer_ReefGenomics_v1")) #Corallimorpharia
MRCA(tree_rooted,c("Actinia_tenebrosa_GCF_009602425_1","Nematostella_vectensis_GCF_932526225_1")) #Actiniaria

#Edit tip label again so that latin names are italicised
tree_rooted@phylo$tip.label <- 
  case_when(grepl('Pocillopora_meandrina_GCA',tree_rooted@phylo$tip.label) ~ 'italic("Pocillopora meandrina")~-~NCBI',
            grepl('sp_',tree_rooted@phylo$tip.label) ~ 
              paste0('italic("',str_extract(tree_rooted@phylo$tip.label,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
            TRUE ~ paste0('italic("',gsub("_"," ",str_extract(tree_rooted@phylo$tip.label,"([A-Z].*)_([A-Z].*)",group=1)),'")'))

#Individual tree plot
g <- ggtree(tree_rooted)+
  geom_tiplab(size=8,linetype=NULL,parse=TRUE)+ #parse=T needed for partial italic of tip labels
  xlim(-0.05,2.1)+
  geom_treescale(x=0,y=20,fontsize=10)+
  geom_rootedge(rootedge=0.05)+
  geom_cladelab(node=44,label=" Scleractinia - Robust",offset=0.63,fontface=2,angle=90,vjust=1.2,hjust=0.5,
                barcolour='navy',barsize=2,fontsize=11)+
  geom_cladelab(node=38,label=" Scleractinia - Complex",offset=0.415,fontface=2,angle=90,vjust=1.2,hjust=0.5,
                barcolour='navy',barsize=2,fontsize=11)+
  geom_cladelab(node=41,label=" Corallimorpharia",offset=0.5,fontface=2,angle=90,vjust=1.2,hjust=0.5,
                barcolour='darkgreen',barsize=2,fontsize=11)+
  geom_cladelab(node=59,label=" Actiniaria",offset=0.4,fontface=2,angle=90,vjust=1.2,hjust=0.5,
                barcolour='salmon',barsize=2,fontsize=11)+
  geom_nodepoint(aes(fill=UfBS_cat, subset=UfBS_cat == 'UfBS ≥ 95'),size=3,shape=21)+ 
  geom_nodepoint(aes(colour=sCF_cat, subset=sCF_cat == 'sCF ≥ 34'),size=1.5,shape=21,fill='white')+  
  geom_nodepoint(aes(colour=sCF_cat, subset=sCF_cat == 'sCF ≥ 50'),size=1.5,shape=21,fill='green')+
  scale_fill_manual(values=c('white','white','green'))+
  scale_colour_manual(values=c('black','black'))+
  guides(fill=guide_legend(title='Ultrafast Bootstrap',order=1,keywidth=0.5,keyheight=0.5),
         colour=guide_legend(title='Site Concordance Factor',order=2,keywidth=0.5,keyheight=0.5,
                             override.aes=list(fill=c('white','green'),colour=c('black','black'))))+
  theme(text=element_text(family="work",size=10), 
        legend.position=c(0.2,0.8),legend.margin=margin(0,0,0,0),
        legend.text=element_text(size=16,margin=margin(l=-5)), 
        legend.title=element_text(size=18,margin=margin(b=-5)),
        legend.box.margin=margin(3,3,3,3),legend.spacing=unit(3,'mm'),
        legend.box.background=element_rect(colour='black',linewidth=0.5))
g
#Find the long branch
#which(grepl(max(g$data$branch.length),g$data$branch.length))
#g$data[g$data$node %in% c(1),"x"] <- 0.01
#g = g %>% collapse(205,mode="none") %>% collapse(241,mode="none")
ggsave("export/mapgenome_tree_i50_215loci.png",plot=g,device="png",width=6,height=5,units="in",dpi=300)
tree_tiporder <- ggtree::get_taxa_name(g) 

#Tree plot for composite
g_tree <- ggtree(tree_rooted)+
  geom_tiplab(aes(colour=inProbe),size=5,linetype=NULL,parse=T,offset=0.02)+ #parse=T needed for partial italic of tip labels
  xlim(-0.03,2.1)+
  geom_treescale(x=0,y=15,fontsize=5)+
  geom_rootedge(rootedge=0.03)+
  geom_cladelab(node=44,label=" Scleractinia - Robust",offset=0.645,fontface=2,angle=90,vjust=1.4,hjust=0.5,
                barcolour='navy',barsize=2,fontsize=7)+
  geom_cladelab(node=38,label=" Scleractinia - Complex",offset=0.43,fontface=2,angle=90,vjust=1.4,hjust=0.5,
                barcolour='navy',barsize=2,fontsize=7)+
  geom_cladelab(node=41,label=" Corallimorpharia",offset=0.5,fontface=2,angle=90,vjust=1.4,hjust=0.5,
                barcolour='darkgreen',barsize=2,fontsize=7)+
  geom_cladelab(node=59,label=" Actiniaria",offset=0.4,fontface=2,angle=90,vjust=1.4,hjust=0.5,
                barcolour='salmon',barsize=2,fontsize=7)+
  geom_nodepoint(aes(shape=UfBS_cat, subset=UfBS_cat == 'UfBS ≥ 95'),size=3,fill='white')+ 
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 34'),size=1.5,shape=21)+  
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 50'),size=1.5,shape=21)+
  scale_shape_manual(values=21)+
  scale_fill_manual(values=c('white','green'))+
  scale_colour_manual(values=c('blue','red','purple','black'))+  #,labels=paste0("<span style='margin-top:0px; color:",c('blue','red','purple','black'),"'>",c('Genome','Transcriptome','Both','Not used'),"</span>"))+
  guides(shape=guide_legend(title='Ultrafast Bootstrap',order=1,keywidth=0.5,keyheight=0.5),
         fill=guide_legend(title='Site Concordance Factor',order=2,keywidth=0.5,keyheight=0.5),  #,override.aes=list(fill=c('white','green'),colour=c('black','black'))
         colour='none')+ 
  theme(text=element_text(family="work",size=10),plot.margin=margin(0,2,0,0),
        legend.position=c(0.2,0.8),legend.margin=margin(0,0,0,1),
        legend.text=element_text(size=16,margin=margin(0,-5,0,-5)), 
        legend.title=element_text(size=18,margin=margin(b=-6)),
        legend.box.margin=margin(3,3,3,3),legend.spacing=unit(3,'mm'),
        legend.box.background=element_rect(colour='black',linewidth=0.5))
g_tree
