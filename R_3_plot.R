# Script info ----
#Script to plot data wrangled outputs from mapping UCEs to genome
#Requirements: datawrangle.R and outputs
#Written with R version 4.3.3 (2024-02-29)
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: April 2024

## Load prerequisite ----
source(file="R_2_datawrangle.R") #Ignore warnings
library(grid)
library(gridExtra)
library(patchwork)
library(scales)
library(showtext)
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE) #Change to FALSE to turn it off
hex_values <- c("#525252","#238b45","#4292c6","#6a51a3","#ff7f00")

N50 <- readxl::read_xlsx(path='/path-to-MappingUCE_tables.xlsx',sheet='1.genomes',skip=1,n_max=35) %>% 
  select(Taxon,`NCBI RefSeq / other identifiers`,N50) %>% 
  mutate(Taxon=paste0(sub(' ','_',sub('\\.','',Taxon)),'_',str_extract(`NCBI RefSeq / other identifiers`,'([A-Z]).*',group=1))) %>% 
  select(-`NCBI RefSeq / other identifiers`) %>% 
  mutate(genome=case_when(grepl('Pocillopora_meandrina_G',Taxon) ~ 'italic("Pocillopora meandrina")~-NCBI',
                          grepl('sp_',Taxon) ~ paste0('italic("',str_extract(Taxon,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
                          TRUE ~ paste0('italic("',gsub("_"," ",str_extract(Taxon,"([A-Z].*)_([A-Z].*)",group=1)),'")'))) %>% 
  mutate(genome=factor(genome,levels=tree_id70_i50_tiporder))

# Relaxed search only hereinafter----

## loci prop ----
all_relaxed_type_flipped = all_relaxed_type %>% 
  mutate(genome=as.factor(genome),
         genome=factor(genome,levels=rev(levels(genome))),
         Type=factor(Type,levels=rev(levels(Type))),
         percent=percent*100)
## loci prop: transcriptome separated ----
all_relaxed_type_trans_flipped = all_relaxed_type_trans %>% 
  mutate(genome=as.factor(genome),
         genome=factor(genome,levels=rev(levels(genome))),
         Type=factor(Type,levels=rev(levels(Type))),
         percent=percent*100,
         origin=factor(origin,levels=c('transcriptome','genome')))

## Distances ----
summary(all_relaxed_scaf$corrected_Distance,na.rm=T) #mean:255556 median:70614 range:90-26468394
sd(all_relaxed_scaf$corrected_Distance,na.rm=T) 
all_relaxed_scaf %>% summarise(n=n()) 

g_relaxed_hist1 <- 
  ggplot(subset(all_relaxed_scaf,corrected_Distance!='NA' & corrected_Distance<1000000),aes(x=corrected_Distance))+
  geom_histogram(colour='black',fill='lightgrey',binwidth=10000,linewidth=0.1)+
  annotate(geom='text',x=300000,y=1700,size=4,hjust=0,family='work',
           label=paste0('Mean = ',comma_format()(round(summary(all_relaxed_scaf$corrected_Distance,na.rm=T)[[4]],0)),' ± ',comma_format()(sd(all_relaxed_scaf$corrected_Distance,na.rm=T)),' SD'))+
  annotate(geom='text',x=300000,y=1600,size=4,hjust=0,family='work',label=paste0('Median = ',comma_format()(summary(all_relaxed_scaf$corrected_Distance,na.rm=T)[[3]])))+
  annotate(geom='text',x=300000,y=1500,size=4,hjust=0,family='work',
           label=paste0('Range = ',summary(all_relaxed_scaf$corrected_Distance,na.rm=T)[[1]],' - ',comma_format()(summary(all_relaxed_scaf$corrected_Distance,na.rm=T)[[6]])))+
  scale_x_continuous(labels=comma)+ labs(y='Frequency',x='',title='All mapped UCE loci')+ theme_bw()+
    theme(text=element_text(family="work",size=10),plot.margin=unit(c(0,0,0,-5),"mm"),plot.title=element_text(size=12),
          axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.5,"mm"),axis.text=element_text(size=8),
          axis.title.y=element_text(size=11,margin=margin(0,5,0,-5)),axis.title.x=element_blank(),
          panel.border=element_rect(linewidth=0.2),panel.grid=element_line(linewidth=0.2))
g_relaxed_hist2 <- 
  ggplot(subset(all_relaxed_scaf,corrected_Distance < 1100 & corrected_Distance > 0),aes(x=corrected_Distance))+
  geom_histogram(colour='black',fill='lightgrey',binwidth=20,linewidth=0.1)+
  annotate(geom='text',x=200,y=30,size=4,hjust=0,family='work',
           label=paste0('n = ',count(subset(all_relaxed_scaf,corrected_Distance < 1100 & corrected_Distance > 0))[[1]],
                        ' (',round(count(subset(all_relaxed_scaf,corrected_Distance < 1100 & corrected_Distance > 0))/count(all_relaxed_scaf),digits=3)*100,' %)'))+
  scale_y_continuous(n.breaks=4)+ theme_bw()+
  labs(title='Distances < 1,100 bp',x='',y=NULL)+
  theme(text=element_text(family="work",size=10),plot.margin=unit(c(0,0,0,1),"mm"),plot.title=element_text(size=12),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.5,"mm"),axis.text=element_text(size=8),
        axis.title=element_blank(),panel.border=element_rect(linewidth=0.2),panel.grid=element_line(linewidth=0.2))
g_relaxed_xlab <- cowplot::get_plot_component(
  ggplot()+labs(x='Distance (bp) from previous UCE locus on same scaffold')+
    theme(text=element_text(family='work',size=10),axis.title=element_text(size=11,margin=margin(0,0,0,0))),
  'xlab-b')
g_relaxed_hist1 + g_relaxed_hist2 + g_relaxed_xlab + plot_layout(widths=c(2,1),heights=c(1000,1),design='AB
              CC') 
ggsave("export/mappedUCE_distances_relaxed.png",device="png",width=4,height=2,units="in",dpi=300)


## Common UCEs: on Ahya ----
common_uce_relaxed_onAhya <- common_uce_chromo_relaxed %>% 
  pivot_longer(names_to="genome",values_to="Chr",cols=c(Acropora_millepora_G:Nematostella_vectensis_G)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Acropora_hyacinthus_G=factor(Acropora_hyacinthus_G,levels=paste0("Chr",1:14))) %>% 
  filter(Acropora_hyacinthus_G!="NA")
g_common_uce_relaxed_onAhya <- ggplot(common_uce_relaxed_onAhya, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+ facet_wrap(~Acropora_hyacinthus_G)
g_common_uce_relaxed_onAhya <- g_common_uce_relaxed_onAhya + geom_point(data=common_uce_relaxed_onAhya,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ scale_y_continuous(breaks=c(1:15,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),axis.text.x=element_blank(),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-5)))
ggsave("export/mappedUCE_common_genomeWchromosome_relaxed_onAhya.png",plot=g_common_uce_relaxed_onAhya,device="png",width=8,height=8,units="in",dpi=300)

## Common UCEs: on Nvec ----
common_uce_relaxed_onNvec <- common_uce_chromo_relaxed %>% #select(-c(n,Montipora_capitata_C)) %>% 
  pivot_longer(names_to="genome",values_to="Chr",cols=c(Acropora_hyacinthus_G:Montipora_capitata_C)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Nematostella_vectensis_G=factor(Nematostella_vectensis_G,levels=paste0("Chr",1:15))) %>% 
  filter(Nematostella_vectensis_G!="NA")
g_common_uce_relaxed_onNvec <- ggplot(common_uce_relaxed_onNvec, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+ facet_wrap(~Nematostella_vectensis_G)
g_common_uce_relaxed_onNvec <- g_common_uce_relaxed_onNvec + geom_point(data=common_uce_relaxed_onNvec,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ scale_y_continuous(breaks=c(1:15,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),axis.text.x=element_blank(),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-5)))
ggsave("export/mappedUCE_common_genomeWchromosome_relaxed_onNvec.png",plot=g_common_uce_relaxed_onNvec,device="png",width=8,height=8,units="in",dpi=300)

## UCEs in phylogeny ----
uce_inphylogeny_relaxed_toplot <- uce_inphylogeny_relaxed %>% select(-c(Montipora_capitata_C)) %>% 
  pivot_longer(names_to="genome",values_to="Chr",cols=c(Acropora_millepora_G:Nematostella_vectensis_G)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Acropora_hyacinthus_G=factor(Acropora_hyacinthus_G,levels=paste0("Chr",1:15))) %>% filter(Acropora_hyacinthus_G!="NA")

g_uce_inphylo_relaxed <- ggplot(uce_inphylogeny_relaxed_toplot, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+ facet_wrap(~Acropora_hyacinthus_G)
g_uce_inphylo_relaxed <- g_uce_inphylo_relaxed + 
  geom_point(data=uce_inphylogeny_relaxed_toplot,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ scale_y_continuous(breaks=c(1:14,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),axis.text.x=element_blank(),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-8)))
ggsave("export/mappedUCE_inphylogeny_genomeWchromosome_relaxed.png",plot=g_uce_inphylo_relaxed,device="png",width=8,height=8,units="in",dpi=300)

## Align with tree ----
all_relaxed_type_flipped_Wtree <- all_relaxed_type_flipped %>% 
  mutate(genome=case_when(grepl('Pocillopora_meandrina_G',genome) ~ 'italic("Pocillopora meandrina")~-NCBI',
                          grepl('sp_',genome) ~ paste0('italic("',str_extract(genome,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
                          TRUE ~ paste0('italic("',gsub("_"," ",str_extract(genome,"([A-Z].*)_([A-Z].*)",group=1)),'")'))) %>% 
  mutate(genome=factor(genome,levels=tree_id70_i50_tiporder)) 
all_relaxed_type_trans_flipped_Wtree <- all_relaxed_type_trans_flipped %>% 
  mutate(genome=case_when(grepl('Pocillopora_meandrina_G',genome) ~ 'italic("Pocillopora meandrina")~-NCBI',
                          grepl('sp_',genome) ~ paste0('italic("',str_extract(genome,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
                          TRUE ~ paste0('italic("',gsub("_"," ",str_extract(genome,"([A-Z].*)_([A-Z].*)",group=1)),'")'))) %>% 
  mutate(genome=factor(genome,levels=tree_id70_i50_tiporder)) 

all_relaxed_trans_Wtree_loci <- all_relaxed_type_trans_flipped_Wtree %>% 
  group_by(genome) %>% summarise(Mapped=sum(n))
phyluce_harvest_loci <- readxl::read_xlsx(path='/path-to-MappingUCE_tables.xlsx',sheet='S2.phyluce',skip=4,col_names=F) %>% 
  select(`...3`,`...9`,`...15`) %>%  #label and uce written default, id70
  mutate(genome=case_when(grepl('Pocillopora_meandrina_G',`...3`) ~ 'italic("Pocillopora meandrina")~-NCBI',
                          grepl('sp_',`...3`) ~ paste0('italic("',str_extract(`...3`,"([A-Z].*)_(sp_.*)",group=1),'")~sp.'),
                          TRUE ~ paste0('italic("',gsub("_"," ",str_extract(`...3`,"([A-Z].*)_([A-Z].*)",group=1)),'")'))) %>% 
  mutate(genome=factor(genome,levels=tree_id70_i50_tiporder)) %>% select(-`...3`) %>% rename(Harvested_default=`...9`) %>% rename(Harvested_id70=`...15`) %>% 
  filter(!is.na(genome)) %>% relocate(genome,.before='Harvested_default')
all_relaxed_trans_Wtree_sum <- left_join(phyluce_harvest_loci,all_relaxed_trans_Wtree_loci) %>% 
  pivot_longer(cols=Harvested_default:Mapped,names_to='UCE',values_to='loci')
all_relaxed_type_trans_flipped_Wtree = rbind(all_relaxed_type_trans_flipped_Wtree,
                                             data.frame(Type=rep('Unclassified',4),genome=c("italic(\"Dendronephthya gigantea\")","italic(\"Paramuricea clavata\")","italic(\"Renilla muelleri\")","italic(\"Heliopora coerulea\")"),origin=rep('genome',4),Order=rep('Octocorallia',4),n=rep(NA,4),percent=rep(NA,4)))

### FigB ----
g_horiz_type_relaxed_sum_Wtree <- 
  ggplot(subset(all_relaxed_trans_Wtree_sum,UCE!='Harvested_default'),aes(y=genome,x=loci))+
  geom_col(aes(fill=UCE),colour='black',linewidth=0.3,position=position_dodge(width=-1),width=0.8)+
  #geom_vline(xintercept=100,linewidth=0.3,colour='red')+
  geom_text(data=subset(all_relaxed_trans_Wtree_sum,loci<100 & UCE=='Mapped'),aes(label=loci,family='work'),x=270,size=5,nudge_y=-0.25)+
  geom_text(data=subset(all_relaxed_trans_Wtree_sum,loci<100 & UCE=='Harvested'),aes(label=loci,family='work'),x=270,size=5,nudge_y=0.25)+
  scale_fill_manual(values=c('#4292c6','lightgrey'),labels=c('Harvested','Mapped'))+
  scale_x_continuous(position='top',name=NULL,breaks=c(seq(0,1500,500)))+
  scale_y_discrete(limits=rev,name=NULL)+ theme_bw()+
  guides(fill=guide_legend(title=NULL,keywidth=0.3,keyheight=0.3,nrow=2))+
  theme(text=element_text(family="work",size=10),plot.margin=margin(-5,2,-5,-20),legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5),legend.box.spacing=unit(-5,'mm'),legend.position='top',#legend.position='inside',legend.position.inside=c(0.7,0.95),
        legend.text=element_text(size=13,margin=margin(-5,0,-5,1)),legend.key.spacing=unit(0.1,'mm'),
        axis.text.x=element_text(size=12),axis.text.y=element_blank(),axis.title=element_blank(),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.4,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.grid.major=element_line(linewidth=0.2))
### FigC ----
g_horiz_type_trans_relaxed_Wtree <- 
  ggplot(all_relaxed_type_trans_flipped_Wtree,aes(y=genome,x=percent))+
  geom_col(aes(fill=Type,alpha=origin),position='fill',linewidth=0.2,colour='white')+
  geom_text(aes(label=n,family="work",group=Type:origin),colour='black',size=5,position=position_fill(vjust=0.5),check_overlap=TRUE)+
  scale_alpha_manual(values=c(0.5,1),labels=c('Transcriptome','Genome'))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0),position='top',name='Proportion of mapped UCE locus type')+
  scale_y_discrete(limits=rev)+ theme_bw()+
  guides(fill=guide_legend(reverse=TRUE,keywidth=0.22,keyheight=0.4,order=1),
         alpha=guide_legend(title='Probe origin',reverse=TRUE,keywidth=0.4,keyheight=0.4,order=2,
                            override.aes=list(fill=c("#ff7f00","#ff7f00"),alpha=c(1,0.5))))+ 
  theme(text=element_text(family="work",size=12),plot.margin=margin(-1,2,0,0),legend.margin=margin(0,-5,0,-8),
        legend.title=element_text(size=16,margin=margin(0,0,1,0)),legend.text=element_text(size=14,margin=margin(0,0,0,1)),
        axis.title.x=element_text(size=15),axis.text.x=element_text(size=12),
        axis.title.y=element_blank(),axis.text.y=element_blank(),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.4,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.grid.major=element_line(linewidth=0.2))

g_id70_i50_tree + g_horiz_type_relaxed_sum_Wtree + g_horiz_type_trans_relaxed_Wtree  + 
  plot_layout(widths=c(6,1,4)) +
  plot_annotation(tag_levels='A') & theme(plot.tag=element_text(family="work",size=28))
ggsave("export/mapgenome_composite_id70.png",device="png",width=8,height=5,units="in",dpi=300)

## Regression ----
loci_recovered <- left_join(all_relaxed_trans_Wtree_loci,N50) %>% select(-Taxon) %>% left_join(.,phyluce_harvest_loci)

g_loci_recovered1 <- ggplot(loci_recovered)+
  geom_point(aes(y=N50,x=Harvested_default),size=1,shape=21,colour='black',fill='#4292c6')+
  geom_smooth(aes(y=N50,x=Harvested_default),method='lm',linewidth=0.5,colour='#4292c6',fill='#4292c6',alpha=0.3)+
  annotate(geom="text",x=400,y=13000000,size=4,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=loci_recovered,N50~Harvested_default))$coef[1],digits=3),' + ',
                        round(summary(lm(data=loci_recovered,N50~Harvested_default))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=400,y=12000000,size=4,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=loci_recovered,N50~Harvested_default))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=400,y=11000000,size=4,hjust=0,family='work',parse=T,
           label=paste0("italic('p') ==",round(summary(lm(data=loci_recovered,N50~Harvested_default))$coef[8],digits=3)))+
  geom_point(aes(y=N50,x=Harvested_id70),size=1,shape=21,colour='black',fill='#0BFCFF')+
  geom_smooth(aes(y=N50,x=Harvested_id70),method='lm',linewidth=0.5,colour='#0BFCFF',fill='#0BFCFF',alpha=0.3)+
  annotate(geom="text",x=700,y=4000000,size=4,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=loci_recovered,N50~Harvested_id70))$coef[1],digits=3),' + ',
                        round(summary(lm(data=loci_recovered,N50~Harvested_id70))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=700,y=3000000,size=4,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=loci_recovered,N50~Harvested_id70))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=700,y=2000000,size=4,hjust=0,family='work',parse=T,
           label=paste0("italic('p') ==",round(summary(lm(data=loci_recovered,N50~Harvested_id70))$coef[8],digits=3)))+
  labs(y='N50 of genome assembly (Mbp)',x='Number of UCE loci harvested')+  theme_bw()+ 
  scale_y_continuous(limits=c(0,3e7),labels=(seq(0,30,10)))+
  theme(text=element_text(family="work",size=16),plot.margin=margin(0,0,0,0),axis.ticks=element_line(linewidth=0.2))
g_loci_recovered2 <- ggplot(loci_recovered)+
  geom_point(aes(y=N50,x=Mapped),size=1,shape=21,fill='darkgrey')+
  geom_smooth(aes(y=N50,x=Mapped),method='lm',linewidth=0.5,colour='darkgrey',fill='darkgrey',alpha=0.3)+
  annotate(geom="text",x=500,y=18000000,size=4,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=loci_recovered,N50~Mapped))$coef[1],digits=3),' + ',
                        round(summary(lm(data=loci_recovered,N50~Mapped))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=500,y=16000000,size=4,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=loci_recovered,N50~Mapped))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=500,y=14000000,size=4,hjust=0,family='work',parse=T,
           label=paste0("italic('p') ==",round(summary(lm(data=loci_recovered,N50~Mapped))$coef[8],digits=3)))+
  xlab('Number of UCE loci mapped')+ scale_y_continuous(limits=c(0,5e7),labels=(seq(0,50,10)))+ theme_bw()+ 
  theme(text=element_text(family="work",size=16),plot.margin=margin(0,0,0,2),
        axis.title.y=element_blank(),axis.ticks=element_line(linewidth=0.2))
g_loci_recovered3 <- ggplot(loci_recovered)+
  geom_point(aes(y=Harvested_default,x=Mapped),size=1,shape=21,fill='#4292c6')+
  geom_smooth(aes(y=Harvested_default,x=Mapped),method='lm',linewidth=0.5,colour='#4292c6',fill='#4292c6',alpha=0.3)+
  annotate(geom="text",x=1000,y=900,size=4,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=loci_recovered,Harvested_default~Mapped))$coef[1],digits=3),' + ',
                        round(summary(lm(data=loci_recovered,Harvested_default~Mapped))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=1000,y=830,size=4,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=loci_recovered,Harvested_default~Mapped))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=1000,y=760,size=4,hjust=0,family='work',parse=T,
           label=paste0("italic('p') < 0.001"))+ #paste0("italic('p') == ",round(summary(lm(data=loci_recovered,Harvested_default~Mapped))$coef[8],digits=3)
  geom_point(aes(y=Harvested_id70,x=Mapped),size=1,shape=21,fill='#0BFCFF')+
  geom_smooth(aes(y=Harvested_id70,x=Mapped),method='lm',linewidth=0.5,colour='#0BFCFF',fill='#0BFCFF',alpha=0.3)+
  annotate(geom="text",x=500,y=1950,size=4,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=loci_recovered,Harvested_id70~Mapped))$coef[1],digits=3),' + ',
                        round(summary(lm(data=loci_recovered,Harvested_id70~Mapped))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=500,y=1880,size=4,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=loci_recovered,Harvested_id70~Mapped))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=500,y=1810,size=4,hjust=0,family='work',parse=T,
           label=paste0("italic('p') ==",round(summary(lm(data=loci_recovered,Harvested_id70~Mapped))$coef[8],digits=3)))+
  labs(y='Number of UCE loci harvested',x='Number of UCE loci mapped')+ theme_bw()+ 
  theme(text=element_text(family="work",size=16),plot.margin=margin(0,0,0,3),axis.ticks=element_line(linewidth=0.2))

g_loci_recovered1 + g_loci_recovered2 + g_loci_recovered3 + 
  plot_annotation(tag_levels='A',theme=theme(text=element_text(family='work',size=20)))
ggsave("export/mapgenome_supplementary_regression.png",device="png",width=6,height=3,units="in",dpi=300)
  
## Genetic distance ----
gene_dist <- readxl::read_xlsx('/path-to-MappingUCE_tables.xlsx', sheet='S3.genDist',skip=1) %>% select(c(1:6)) %>% 
  mutate(genome=case_when(grepl('Pocillopora meandrina -NCBI',genome) ~ 'italic("Pocillopora meandrina")~-NCBI',
                          grepl('sp.',genome) ~ paste0('italic("',str_extract(genome,"([A-Z].*) (sp.*)",group=1),'")~sp.'),
                          TRUE ~ paste0('italic("',str_extract(genome,"([A-Z].* [a-z].*)",group=1),'")')))
gene_dist_toplot <- left_join(gene_dist,loci_recovered, by="genome") %>% 
  filter(!is.na(inProbe_sister_genome)) %>% filter(genus_tocompare==inProbe_sister_genome) %>% 
  mutate(dist_mean_default=as.numeric(dist_mean_default))

g_gene_dist <- ggplot(gene_dist_toplot)+
  geom_point(aes(y=Harvested_default,x=dist_mean_default),size=1,shape=21,fill='#4292c6')+
  geom_smooth(aes(y=Harvested_default,x=dist_mean_default),method='lm',linewidth=0.5,colour='#4292c6',fill='#4292c6',alpha=0.3)+ 
  annotate(geom='text',x=0.5,y=950,size=6,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=gene_dist_toplot,Harvested_default~dist_mean_default))$coef[1],digits=3),
                        ' ',round(summary(lm(data=gene_dist_toplot,Harvested_default~dist_mean_default))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=0.5,y=880,size=6,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=gene_dist_toplot,Harvested_default~dist_mean_default))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=0.5,y=810,size=6,hjust=0,family='work',parse=T,
           label="italic('p') < 0.001")+ #paste0("italic('p') == ",round(summary(lm(data=loci_recovered,Harvested_default~Mapped))$coef[8],digits=3)))+
  geom_point(aes(y=Harvested_id70,x=dist_mean_id70),size=1,shape=21,fill='#0BFCFF')+
  geom_smooth(aes(y=Harvested_id70,x=dist_mean_id70),method='lm',linewidth=0.5,colour='#0BFCFF',fill='#0BFCFF',alpha=0.3)+
  annotate(geom="text",x=0.3,y=1650,size=6,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=gene_dist_toplot,Harvested_id70~dist_mean_id70))$coef[1],digits=3),
                        ' ',round(summary(lm(data=gene_dist_toplot,Harvested_id70~dist_mean_id70))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=0.3,y=1580,size=6,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=gene_dist_toplot,Harvested_id70~dist_mean_id70))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=0.3,y=1510,size=6,hjust=0,family='work',parse=T,
           label=paste0("italic('p') ==",round(summary(lm(data=gene_dist_toplot,Harvested_id70~dist_mean_id70))$coef[8],digits=3)))+
  geom_point(aes(y=Mapped,x=dist_mean_default),size=1,shape=21,fill='darkgrey')+
  geom_smooth(aes(y=Mapped,x=dist_mean_default),method='lm',linewidth=0.5,colour='darkgrey',fill='darkgrey',alpha=0.3)+ 
  annotate(geom='text',x=0.1,y=750,size=6,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$coef[1],digits=3),
                        ' ',round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=0.1,y=680,size=6,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=0.1,y=610,size=6,hjust=0,family='work',parse=T,
           label="italic('p') < 0.001")+ #paste0("italic('p') == ",round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$coef[8],digits=3)))+ 
  theme_bw()+ labs(y='Number of UCE loci',x='Mean genetic distance to the nearest genome used in probe set')+
  theme(text=element_text(family="work",size=20),axis.ticks=element_line(linewidth=0.2),axis.title=element_text(size=18))
g_gene_dist
ggsave("export/mapgenome_geneticdistances.png",plot=g_gene_dist,device="png",width=4,height=4,units="in",dpi=300)

#if keeping mapped separate...
g_gene_dist_2 <- ggplot(gene_dist_toplot)+
  geom_point(aes(y=Mapped,x=dist_mean_default),size=1,shape=21,fill='darkgrey')+
  geom_smooth(aes(y=Mapped,x=dist_mean_default),method='lm',linewidth=0.5,colour='darkgrey',fill='darkgrey',alpha=0.3)+ 
  annotate(geom='text',x=0.5,y=950,size=5,hjust=0,family='work',
           label=paste0("y = ",round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$coef[1],digits=3),
                        ' ',round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$coef[2],digits=3),'x'))+
  annotate(geom='text',x=0.5,y=880,size=5,hjust=0,family='work',parse=T,
           label=paste0("Adj~R^2==",round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$adj.r.squared,digits=3)))+
  annotate(geom='text',x=0.5,y=810,size=5,hjust=0,family='work',parse=T,
           label=paste0("italic('p') == ",round(summary(lm(data=gene_dist_toplot,Mapped~dist_mean_default))$coef[8],digits=3)))+ #"italic('p') < 0.001")+ 
  theme_bw()+ labs(y='Number of UCE loci mapped',x=NULL)+
  theme(text=element_text(family="work",size=16),axis.ticks=element_line(linewidth=0.2),plot.margin=margin(0,0,0,5))
g_gene_xlab <- cowplot::get_plot_component(ggplot()+labs(x='Mean genetic distance to the nearest genome used in probe set')+
                                             theme(axis.title=element_text(size=14,family='work',margin=margin(0,0,0,0))), "xlab-b")
g_gene_dist_1 + g_gene_dist_2 + g_gene_xlab + plot_layout(design='AB
                                                          CC',heights=c(40,1))
#wrap_plots(g_gene_ylab) + (g_gene_dist_1 & labs(title='Harvested') & theme(axis.title=element_blank())) + (g_gene_dist_2 & theme(axis.title=element_blank())) + g_gene_xlab + 
ggsave("export/mapgenome_geneticdistances.png",device="png",width=6,height=3,units="in",dpi=300)

