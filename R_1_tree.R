# Script info ----
#Script to plot trees with genomes used for mapping UCEs
#file type: .tree or .treefile, (optional) metadata in .csv
#Written with R version 4.3.3 (2024-02-29)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: April 2024

## Load necessary libraries ----
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
#packageDescription("ggtree", fields="Version")  #Check version
library(ggtext)
library(showtext)
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE) #Change to FALSE to turn it off

# Data prep ----
## Read files ----
tree_default <- read.tree("data/siteconcord.i50.gen.map.215loci.cf.tree")  #IQ-TREE output with site concordance factors
tree_id70_i50 <- read.tree("data/siteconcord.i50.map.id70.1761loci.cf.tree")

#Run below once
#system('mkdir export')

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ----

# Default (id92.5, i50.215loci)----
## Root tree ----
tree_default_rooted <- tree_default
MRCA(tree_default_rooted,c("Actinia_tenebrosa_GCF_009602425_1","Nematostella_vectensis_GCF_932526225_1"))
tree_default_rooted <- ape::root(tree_default_rooted,node=41,resolve.root=TRUE)   #Actiniaria
#Add info to the plot
treedf_default <- as_tibble(tree_default_rooted)
treedf_default = treedf_default %>% #left_join(treedf_default, info, by='label') %>% 
  mutate(SHalrt=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub("\\/.*","",label)),
         UfBS=case_when(grepl("C",label) ~ "NA",
                       TRUE ~ gsub("\\/.*","",str_extract(label,"([^/]+)(?:/[^/]+){1}$"))),
         sCF=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub(".*\\/","",label)),
         SHalrt=as.numeric(SHalrt),
         UfBS=as.numeric(UfBS),
         sCF=as.numeric(sCF),
         SHalrt_cat=factor(case_when(SHalrt>=85 ~ 'SHalrt ≥ 85',  SHalrt<85 ~ 'SHalrt < 85',TRUE ~ NA)),
         UfBS_cat=factor(case_when(UfBS>=95 ~ 'UfBS ≥ 95',  UfBS<95 ~ 'UfBS < 95', TRUE ~ NA)),
         sCF_cat=factor(case_when(sCF>=34 & sCF<50 ~ 'sCF ≥ 34',  sCF>=50 ~ 'sCF ≥ 50',  sCF<34 ~ 'sCF < 34', TRUE ~ NA)),
         inProbe=factor(case_when(grepl("digitifera|Nematostella",label) ~ 'Both',
                           grepl("Discosoma|Exaiptasia|Montastraea",label) ~ 'Genome',
                           grepl("Orbicella|damicornis",label) ~ 'Transcriptome',
                           TRUE ~ 'Not used'),levels=c('Genome','Transcriptome','Both','Not used'))) #IGNORE WARNING
tree_default_rooted <- treeio::as.treedata(treedf_default) #convert it back from tibble to tree object
#Plot with node suppport values and PSH labels
MRCA(tree_default_rooted,c("Acropora_millepora_GCF_013753865_1","Desmophyllum_pertusum_GCA_029204205_1")) #all Scleractinia
MRCA(tree_default_rooted,c("Pocillopora_meandrina_Cyanophora_v1","Desmophyllum_pertusum_GCA_029204205_1")) #Robust
MRCA(tree_default_rooted,c("Acropora_millepora_GCF_013753865_1","Porites_evermanni_GCA_942486025_1")) #Complex
MRCA(tree_default_rooted,c("Discosoma_sp_ReefGenomics_v1","Amplexidiscus_fenestrafer_ReefGenomics_v1")) #Corallimorpharia
MRCA(tree_default_rooted,c("Actinia_tenebrosa_GCF_009602425_1","Nematostella_vectensis_GCF_932526225_1")) #Actiniaria

###Gene dist ----
gene_dist_default <- ape::cophenetic.phylo(tree_default_rooted@phylo) %>% as.data.frame() %>% 
  mutate(genome1=rownames(.),
         genus1=case_when(is.na(str_extract(as.character(genome1),'.*([A-Z].*)( .*)',group=1)==T) ~ str_extract(as.character(genome1),'.*([A-Z].*)(\".*)',group=1),
                          TRUE ~ str_extract(as.character(genome1),'.*([A-Z].*)( .*)',group=1))) %>% 
  pivot_longer(cols=c(1:30),names_to='genome2',values_to='dist') %>% 
  mutate(genus2=case_when(is.na(str_extract(as.character(genome2),'.*([A-Z].*)( .*)',group=1)==T) ~ str_extract(as.character(genome2),'.*([A-Z].*)(\".*)',group=1),
                          TRUE ~ str_extract(as.character(genome2),'.*([A-Z].*)( .*)',group=1))) %>% 
  group_by(genome1,genus1,genus2) %>% summarise(dist_mean=mean(dist)) %>% ungroup()
#  filter(genus1=='Acropora'|genus1=='Nematostella'|genus1=='Discosoma'|genus1=='Exaiptasia'|genus1=='Montastraea')
write.csv(gene_dist_default,row.names=F,file='export/genetic_distances_mean_215loci.csv')
  #Edit as needed

#gene_dist_default <- readxl::read_xlsx('/Users/jc313394/Library/CloudStorage/OneDrive-JamesCookUniversity/PhD/Manuscripts/Chp1_Map-UCE-to-Genome/4_forthesis/Mera_MappingUCE_v4_tables.xlsx',sheet='S3.genDist')
#gene_dist_default = gene_dist_default %>% mutate(remove=case_when(genus2==nearest_genome ~ 'no',TRUE ~ 'yes')) %>% filter(remove=='no')

#Edit tip label again so that latin names are italicised
tree_default_rooted@phylo$tip.label <- 
  case_when(grepl('Pocillopora_meandrina_GCA',tree_default_rooted@phylo$tip.label) ~ 'italic("Pocillopora meandrina")~-NCBI',
            grepl('sp_',tree_default_rooted@phylo$tip.label) ~ 
              paste0('italic("',str_extract(tree_default_rooted@phylo$tip.label,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
            TRUE ~ paste0('italic("',gsub("_"," ",str_extract(tree_default_rooted@phylo$tip.label,"([A-Z].*)_([A-Z].*)",group=1)),'")'))

## Plot ----
### Standalone Plot ----
g_default <- ggtree(tree_default_rooted)+
  geom_tiplab(aes(colour=inProbe,family='work'),size=8,linetype=NULL,parse=T)+ #parse=T needed for partial italic of tip labels
  xlim(-0.07,2.1)+ geom_treescale(x=0,y=15,fontsize=10,offset=0.5,family='work')+
  geom_cladelab(node=44,label=" Scleractinia - Robust",offset=0.63,angle=90,vjust=1.2,hjust=0.5,barcolour='navy',barsize=2,fontsize=10,fontface=2)+
  geom_cladelab(node=38,label=" Scleractinia - Complex",offset=0.415,angle=90,vjust=1.2,hjust=0.5,barcolour='navy',barsize=2,fontsize=10,fontface=2)+
  geom_cladelab(node=41,label=" Corallimorpharia",offset=0.5,angle=90,vjust=1.2,hjust=0.5,barcolour='darkgreen',barsize=2,fontsize=10,fontface=2)+
  geom_cladelab(node=59,label=" Actiniaria",offset=0.4,angle=90,vjust=1.2,hjust=0.5,barcolour='salmon',barsize=2,fontsize=10,fontface=2)+
  geom_nodepoint(aes(shape=UfBS_cat, subset=UfBS_cat == 'UfBS ≥ 95'),size=2,fill='white')+ 
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 34'),size=1,shape=21)+  
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 50'),size=1,shape=21)+
  scale_shape_manual(values=21,labels=c('UfBS ≥ 95 & SHalrt ≥ 85'))+
  scale_fill_manual(values=c('white','green'))+
  scale_colour_manual(values=c('blue','red','purple','black'))+
  guides(colour='none',  shape=guide_legend(title='Node Support',order=1,keywidth=0.5,keyheight=0.5),
         fill=guide_legend(title='Site Concordance Factor',order=2,keywidth=0.5,keyheight=0.5))+
  theme(text=element_text(family="work",size=10), 
        legend.position='inside',legend.position.inside=c(0.2,0.8),legend.margin=margin(0,0,0,0),
        legend.text=element_text(size=16,margin=margin(l=1)), legend.title=element_text(size=18,margin=margin(b=1)),
        legend.box.margin=margin(3,3,3,3),legend.spacing=unit(3,'mm'),
        legend.box.background=element_rect(colour='black',linewidth=0.5))
#Add truncate line for outgroup rooting
g_default$data[g_default$data$node %in% c(31),"x"] <- -0.05 #outgroup
g_default = g_default+annotate('text',x=-0.03,y=1.9,label='~~',angle=90,size=7)+annotate('text',x=-0.035,y=1.9,label='~~',angle=90,size=7)
g_default
ggsave("export/mapgenome_tree_default_i50_215loci.png",plot=g_default,device="png",width=6,height=5,units="in",dpi=300)
tree_default_tiporder <- ggtree::get_taxa_name(g_default) 
summary(g_default$data$UfBS,na.rm=T)
#g_default$data$UfBS[which(is.na(g_default$data$UfBS)==F)] 
summary(g_default$data$UfBS_cat[which(is.na(g_default$data$UfBS_cat)==F)])[2]/sum(summary(g_default$data$UfBS_cat[which(is.na(g_default$data$UfBS_cat)==F)]))
summary(g_default$data$SHalrt_cat[which(is.na(g_default$data$SHalrt_cat)==F)])
sum(summary(g_default$data$sCF_cat[which(is.na(g_default$data$sCF_cat)==F)])[2:3])/sum(summary(g_default$data$sCF_cat[which(is.na(g_default$data$sCF_cat)==F)]))

### Tree plot for composite ----
g_default_tree <- ggtree(tree_default_rooted,linewidth=0.3)+
  geom_tiplab(aes(colour=inProbe,family='work'),size=5,linetype=NULL,parse=T,offset=0.02)+ #parse=T needed for partial italic of tip labels
  xlim(-0.07,2.1)+ geom_treescale(x=0,y=15,fontsize=10,linesize=0.3,offset=0.5,family='work')+
  geom_cladelab(node=44,label=" Scleractinia - Robust",offset=0.645,fontface=2,angle=90,vjust=1.4,hjust=0.5,barcolour='navy',barsize=2,fontsize=7)+
  geom_cladelab(node=38,label=" Scleractinia - Complex",offset=0.43,fontface=2,angle=90,vjust=1.4,hjust=0.5,barcolour='navy',barsize=2,fontsize=7)+
  geom_cladelab(node=41,label=" Corallimorpharia",offset=0.5,fontface=2,angle=90,vjust=1.4,hjust=0.5,barcolour='darkgreen',barsize=2,fontsize=7)+
  geom_cladelab(node=59,label=" Actiniaria",offset=0.4,fontface=2,angle=90,vjust=1.4,hjust=0.5,barcolour='salmon',barsize=2,fontsize=7)+
  geom_nodepoint(aes(shape=UfBS_cat, subset=UfBS_cat == 'UfBS ≥ 95'),size=2,fill='white')+ 
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 34'),size=1,shape=21)+  
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 50'),size=1,shape=21)+
  scale_shape_manual(values=21)+
  scale_fill_manual(values=c('white','green'))+
  scale_colour_manual(values=c('blue','red','purple','black'))+  #,labels=paste0("<span style='margin-top:0px; color:",c('blue','red','purple','black'),"'>",c('Genome','Transcriptome','Both','Not used'),"</span>"))+
  guides(colour='none',  shape=guide_legend(title='Node Support',order=1,keywidth=0.5,keyheight=0.5),
         fill=guide_legend(title='Site Concordance Factor',order=2,keywidth=0.5,keyheight=0.5))+ #,override.aes=list(fill=c('white','green'),colour=c('black','black'))
  theme(text=element_text(family="work",size=10),plot.margin=margin(0,2,0,0),
        legend.position='inside',legend.position.inside=c(0.2,0.8),legend.margin=margin(0,0,0,1),
        legend.text=element_text(size=16,margin=margin(0,-5,0,1)), legend.title=element_text(size=18,margin=margin(b=1)),
        legend.box.margin=margin(3,3,3,3),legend.spacing=unit(3,'mm'),legend.box.background=element_rect(colour='black',linewidth=0.3))
#Add truncate line for outgroup rooting
g_default_tree$data[g_default_tree$data$node %in% c(31),"x"] <- -0.06 #outgroup
g_default_tree = g_default_tree+annotate('text',x=-0.043,y=1.9,label='~~',angle=90,size=7)+annotate('text',x=-0.047,y=1.9,label='~~',angle=90,size=7)
g_default_tree
tree_default_tiporder <- ggtree::get_taxa_name(g_default_tree) 


# id70 i50 ----
## Root tree ----
tree_id70_i50_rooted <- tree_id70_i50
MRCA(tree_id70_i50_rooted,c("Renilla_muelleri_ReefGenomics_v1_id70uce","Dendronephthya_gigantea_GCF_004324835_1_id70uce"))
tree_id70_i50_rooted <- ape::root(tree_id70_i50_rooted,node=48,resolve.root=TRUE)   #Octocorallia
#Add info to the plot
treedf_id70_i50 <- as_tibble(tree_id70_i50_rooted)
treedf_id70_i50 = treedf_id70_i50 %>% #left_join(treedf_id70_i50, info, by='label') %>% 
  mutate(SHalrt=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub("\\/.*","",label)),
         UfBS=case_when(grepl("C",label) ~ "NA",
                        TRUE ~ gsub("\\/.*","",str_extract(label,"([^/]+)(?:/[^/]+){1}$"))),
         sCF=case_when(grepl("C",label) ~ "NA",  TRUE ~ gsub(".*\\/","",label)),
         SHalrt=as.numeric(SHalrt),
         UfBS=as.numeric(UfBS),
         sCF=as.numeric(sCF),
         SHalrt_cat=factor(case_when(SHalrt>=85 ~ 'SHalrt ≥ 85',  SHalrt<85 ~ 'SHalrt < 85',TRUE ~ NA)),
         UfBS_cat=factor(case_when(UfBS>=95 ~ 'UfBS ≥ 95',  UfBS<95 ~ 'UfBS < 95', TRUE ~ NA)),
         sCF_cat=factor(case_when(sCF>=34 & sCF<50 ~ 'sCF ≥ 34',  sCF>=50 ~ 'sCF ≥ 50',  sCF<34 ~ 'sCF < 34', TRUE ~ NA)),
         inProbe=factor(case_when(grepl("digitifera|Nematostella",label) ~ 'Both',
                                  grepl("Discosoma|Exaiptasia|Montastraea",label) ~ 'Genome',
                                  grepl("Orbicella|damicornis",label) ~ 'Transcriptome',
                                  TRUE ~ 'Not used'),levels=c('Genome','Transcriptome','Both','Not used'))) #IGNORE WARNING
tree_id70_i50_rooted <- treeio::as.treedata(treedf_id70_i50) #convert it back from tibble to tree object
#Plot with node suppport values and PSH labels
MRCA(tree_id70_i50_rooted,c("Acropora_millepora_GCF_013753865_1_id70uce","Desmophyllum_pertusum_GCA_029204205_1_id70uce")) #all Scleractinia
MRCA(tree_id70_i50_rooted,c("Pocillopora_meandrina_Cyanophora_v1_id70uce","Desmophyllum_pertusum_GCA_029204205_1_id70uce")) #Robust
MRCA(tree_id70_i50_rooted,c("Acropora_millepora_GCF_013753865_1_id70uce","Porites_evermanni_GCA_942486025_1_id70uce")) #Complex
MRCA(tree_id70_i50_rooted,c("Discosoma_sp_ReefGenomics_v1_id70uce","Amplexidiscus_fenestrafer_ReefGenomics_v1_id70uce")) #Corallimorpharia
MRCA(tree_id70_i50_rooted,c("Actinia_tenebrosa_GCF_009602425_1_id70uce","Nematostella_vectensis_GCF_932526225_1_id70uce")) #Actiniaria
MRCA(tree_id70_i50_rooted,c("Renilla_muelleri_ReefGenomics_v1_id70uce","Dendronephthya_gigantea_GCF_004324835_1_id70uce")) #Octo

#Edit tip label again so that latin names are italicised
tree_id70_i50_rooted@phylo$tip.label <- 
  case_when(grepl('Pocillopora_meandrina_GCA',tree_id70_i50_rooted@phylo$tip.label) ~ 'italic("Pocillopora meandrina")~-NCBI',
            grepl('sp_',tree_id70_i50_rooted@phylo$tip.label) ~ 
              paste0('italic("',str_extract(tree_id70_i50_rooted@phylo$tip.label,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
            TRUE ~ paste0('italic("',gsub("_"," ",str_extract(tree_id70_i50_rooted@phylo$tip.label,"([A-Z].*)_([A-Z].*)",group=1)),'")'))

### Gene dist export----
gene_dist_id70_i50 <- ape::cophenetic.phylo(tree_id70_i50_rooted@phylo) %>% as.data.frame() %>% 
  mutate(genome1=rownames(.),
         genus1=case_when(is.na(str_extract(as.character(genome1),'.*([A-Z].*)( .*)',group=1)==T) ~ str_extract(as.character(genome1),'.*([A-Z].*)(\".*)',group=1),
                          TRUE ~ str_extract(as.character(genome1),'.*([A-Z].*)( .*)',group=1))) %>% 
  pivot_longer(cols=c(1:34),names_to='genome2',values_to='dist') %>% 
  mutate(genus2=case_when(is.na(str_extract(as.character(genome2),'.*([A-Z].*)( .*)',group=1)==T) ~ str_extract(as.character(genome2),'.*([A-Z].*)(\".*)',group=1),
                          TRUE ~ str_extract(as.character(genome2),'.*([A-Z].*)( .*)',group=1))) %>% 
  group_by(genome1,genus1,genus2) %>% summarise(dist_mean=mean(dist)) %>% ungroup()
#  filter(genus1=='Acropora'|genus1=='Nematostella'|genus1=='Discosoma'|genus1=='Exaiptasia'|genus1=='Montastraea')
write.csv(gene_dist_id70_i50,row.names=F,file='export/genetic_distances_mean_id70_1761loci.csv')
#Edit as needed


## Plot ----
### Standalone Plot ----
g_id70_i50 <- ggtree(tree_id70_i50_rooted)+
  geom_tiplab(size=8,linetype=NULL,parse=T)+ #parse=T needed for partial italic of tip labels
  xlim(-0.05,1.8)+ geom_treescale(x=0,y=12,fontsize=10,offset=0.5,family='work')+
  geom_cladelab(node=52,label="Scleractinia - Robust",offset=0.5,offset.text=0.01,fontface=2,angle=90,vjust=1.2,hjust=0.5,barcolour='navy',barsize=2,fontsize=10,align=T)+
  geom_cladelab(node=42,label="Scleractinia - Complex",offset=0.5,offset.text=0.01,fontface=2,angle=90,vjust=1.2,hjust=0.5,barcolour='navy',barsize=2,fontsize=10,align=T)+
  geom_cladelab(node=51,label="Corallimorpharia",offset=0.43,offset.text=0.02,fontface=2,barcolour='darkgreen',barsize=2,fontsize=10)+
  geom_cladelab(node=46,label="Actiniaria",offset=0.5,offset.text=0.01,fontface=2,angle=90,vjust=1.2,hjust=0.5,barcolour='salmon',barsize=2,fontsize=10,align=T)+
  geom_cladelab(node=67,label="Octocorallia",offset=0.4,offset.text=0.02,fontface=2,barcolour='#FFD400',barsize=2,fontsize=10)+
  geom_nodepoint(aes(fill=UfBS_cat, subset=UfBS_cat == 'UfBS ≥ 95'),size=2,shape=21)+ 
  geom_nodepoint(aes(colour=sCF_cat, subset=sCF_cat == 'sCF ≥ 34'),size=1,shape=21,fill='white')+  
  geom_nodepoint(aes(colour=sCF_cat, subset=sCF_cat == 'sCF ≥ 50'),size=1,shape=21,fill='green')+
  #  geom_nodepoint(aes(subset=SHalrt_cat == 'SHalrt ≥ 90'),size=0.5,shape=21,fill='black')+ 
  scale_fill_manual(values=c('white','white','green'),labels=c('UfBS ≥ 95 & SHalrt ≥ 85','sCF ≥ 34','sCF ≥ 50'))+
  scale_colour_manual(values=c('black','black'))+
  guides(fill=guide_legend(title='Branch support',order=1,keywidth=0.5,keyheight=0.5),
         colour=guide_legend(title='Site Concordance Factor',order=2,keywidth=0.5,keyheight=0.5,
                             override.aes=list(fill=c('white','green'),colour=c('black','black'))))+
  theme(text=element_text(family="work",size=10), 
        legend.position='inside',legend.position.inside=c(0.15,0.8),legend.margin=margin(0,0,0,0),
        legend.text=element_text(size=16,margin=margin(l=1)), legend.title=element_text(size=18,margin=margin(b=1)),
        legend.box.margin=margin(3,3,3,3),legend.spacing=unit(3,'mm'),
        legend.box.background=element_rect(colour='black',linewidth=0.5))
#Add truncate line for outgroup rooting
g_id70_i50$data[g_id70_i50$data$node %in% c(35),"x"] <- -0.05 #Root
g_id70_i50 = g_id70_i50+annotate('text',x=-0.03,y=1.9,label='~~',angle=90,size=7)+annotate('text',x=-0.035,y=1.9,label='~~',angle=90,size=7)
g_id70_i50
ggsave("export/mapgenome_tree_id70.i50.1761loci.png",plot=g_id70_i50,device="png",width=6,height=5,units="in",dpi=300)
tree_id70_i50_tiporder <- ggtree::get_taxa_name(g_id70_i50) 
summary(g_id70_i50$data$UfBS_cat[which(is.na(g_id70_i50$data$UfBS_cat)==F)])
summary(g_id70_i50$data$SHalrt_cat[which(is.na(g_id70_i50$data$SHalrt_cat)==F)])
sum(summary(g_id70_i50$data$sCF_cat[which(is.na(g_id70_i50$data$sCF_cat)==F)])[2:3])/sum(summary(g_id70_i50$data$sCF_cat[which(is.na(g_id70_i50$data$sCF_cat)==F)]))


### Tree plot for composite ----
g_id70_i50_tree <- ggtree(tree_id70_i50_rooted,linewidth=0.3)+
  geom_tiplab(aes(colour=inProbe,family='work'),size=5,linetype=NULL,parse=T,offset=0.02)+ #parse=T needed for partial italic of tip labels
  xlim(-0.05,1.85)+ geom_treescale(x=0,y=12,fontsize=10,linesize=0.3,offset=0.5,family='work')+
  geom_cladelab(node=52,label="Scleractinia - Robust",offset=0.45,offset.text=0.01,fontface=2,angle=90,vjust=1.4,hjust=0.5,barcolour='navy',barsize=2,fontsize=7,align=T)+
  geom_cladelab(node=42,label="Scleractinia - Complex",offset=0.45,offset.text=0.01,fontface=2,angle=90,vjust=1.2,hjust=0.5,barcolour='navy',barsize=2,fontsize=7,align=T)+
  geom_cladelab(node=51,label="Corallimorpharia",offset=0.45,offset.text=0.02,fontface=2,barcolour='darkgreen',barsize=2,fontsize=7)+
  geom_cladelab(node=46,label="Actiniaria",offset=0.45,offset.text=0.01,fontface=2,angle=90,vjust=1.2,hjust=0.5,barcolour='salmon',barsize=2,fontsize=7,align=T)+
  geom_cladelab(node=67,label="Octocorallia",offset=0.4,offset.text=0.02,fontface=2,barcolour='#FFD400',barsize=2,fontsize=7)+
  geom_nodepoint(aes(shape=UfBS_cat, subset=UfBS_cat == 'UfBS ≥ 95'),size=2,fill='white')+ 
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 34'),size=1,shape=21)+  
  geom_nodepoint(aes(fill=sCF_cat, subset=sCF_cat == 'sCF ≥ 50'),size=1,shape=21)+
  scale_shape_manual(values=21,labels=c('UfBS ≥ 95 & SHalrt ≥ 85'))+
  scale_fill_manual(values=c('white','green'))+
  scale_colour_manual(values=c('blue','red','purple','black'))+  #,labels=paste0("<span style='margin-top:0px; color:",c('blue','red','purple','black'),"'>",c('Genome','Transcriptome','Both','Not used'),"</span>"))+
  guides(colour='none',  shape=guide_legend(title='Branch support',order=1,keywidth=0.5,keyheight=0.5),
         fill=guide_legend(title='Site Concordance Factor',order=2,keywidth=0.5,keyheight=0.5))+ #,override.aes=list(fill=c('white','green'),colour=c('black','black'))
  theme(text=element_text(family="work",size=10),plot.margin=margin(-1,0,0,0),
        legend.position='inside',legend.position.inside=c(0.15,0.9),legend.margin=margin(0,0,0,1),
        legend.text=element_text(size=16,margin=margin(0,-5,0,1)), legend.title=element_text(size=18,margin=margin(b=1)),
        legend.box.margin=margin(3,3,3,3),legend.spacing=unit(3,'mm'),legend.box.background=element_rect(colour='black',linewidth=0.3))
#Add truncate line for outgroup rooting
g_id70_i50_tree$data[g_id70_i50_tree$data$node %in% c(35),"x"] <- -0.05 #Root
g_id70_i50_tree = g_id70_i50_tree+annotate('text',x=-0.035,y=1.85,label='~~',angle=90,size=7)+annotate('text',x=-0.04,y=1.85,label='~~',angle=90,size=7)
g_id70_i50_tree
tree_id70_i50_tiporder <- ggtree::get_taxa_name(g_id70_i50_tree) 



