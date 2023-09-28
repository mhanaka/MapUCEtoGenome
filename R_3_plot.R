# Script info ----
#Script to plot data wrangled outputs from mapping UCEs to genome
#Requirements: datawrangle.R and outputs
#Written with R-4.2.2
#If on Rstudio, I recommend  toggling 'document outline' so you can see the headings
#Author: Hanaka Mera
#Email: hanaka.mera[@]my.jcu.edu.au
#Last modified: May 2023

## Load prerequisite ----
source(file="datawrangle.R")
library(grid)
library(gridExtra)
library(showtext)
font_add_google(name="Raleway", family="raleway") 
font_add_google(name="Roboto", family="roboto") 
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE) #Change to FALSE to turn it off
hex_values <- c("#525252","#238b45","#4292c6","#6a51a3","#ff7f00")

# default ----
## loci prop: bar plot ----
ggplot(all_type,aes(y=percent,x=Type,fill=Type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=percent*100),vjust=-0.2,size=2)+
  labs(x=NULL,y=NULL)+
  scale_fill_brewer(palette="Blues",direction=-1)+
  #scale_y_continuous(breaks=cumsum(dt_type$percent)-dt_type$percent/2,labels=dt_type$percent*100)+
  facet_wrap(~genome)+theme_bw()+
  theme(legend.position="none",strip.text=element_text(size=7,face="italic"),
        axis.text=element_text(size=7,angle=90,hjust=1))

## loci prop: Stacked ----
all_type_flipped = all_type %>% 
  mutate(genome=as.factor(genome),
         genome=factor(genome,levels=rev(levels(genome))),
         Type=factor(Type,levels=rev(levels(Type))),
         percent=percent*100)
g_horiz_type <- ggplot(all_type_flipped,aes(y=genome,x=percent,fill=Type))+
  geom_col(position="fill",alpha=0.7)+
  geom_text(aes(y=genome,x=percent,label=n,family="work"),size=3,position=position_fill(vjust=0.5))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0))+
  guides(fill=guide_legend(reverse=TRUE,nrow=2))+
  ggforce::facet_col(vars(Order),scales="free_y",space="free")+
  theme_bw()+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position="bottom",legend.text=element_text(size=9,margin=margin(l=-4)),
        legend.key.size=unit(1,"mm"),legend.margin=margin(-10,0,-5,-10),
        strip.text=element_text(size=11,angle=0,face='bold',hjust=0,margin=margin(0,0,0.1,0)),
        strip.background=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(face='bold.italic',size=9),axis.text.x=element_text(size=9),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.3,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.spacing=unit(0.5,"mm"),panel.grid.major=element_line(linewidth=0.2))
ggsave("export/mappedUCE_Types_allgenome.png",plot=g_horiz_type,device="png",width=6,height=5,units="cm",dpi=500)

## loci prop: transcriptome separated ----
all_type_trans_flipped = all_type_trans %>% 
  mutate(genome=as.factor(genome),
         genome=factor(genome,levels=rev(levels(genome))),
         Type=factor(Type,levels=rev(levels(Type))),
         percent=percent*100)
g_horiz_type_trans1 <- ggplot(subset(all_type_trans_flipped,origin=="genome"),
                                    aes(y=genome,x=percent,fill=Type))+
  geom_col(position="fill",alpha=0.7)+
  geom_text(aes(y=genome,x=percent,label=n,family="work"),size=2,position=position_fill(vjust=0.5))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0))+
  guides(fill=guide_legend(reverse=TRUE,nrow=2))+
  ggforce::facet_col(vars(Order),scales="free_y",space="free")+ theme_bw()+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position="bottom",legend.text=element_text(size=9,margin=margin(l=-4)),
        legend.key.size=unit(1,"mm"),legend.margin=margin(-10,0,-5,-10),
        strip.text=element_text(size=11,angle=0,face='bold',hjust=0,margin=margin(0,0,0.1,0)),
        strip.background=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(face='bold.italic',size=8),axis.text.x=element_text(size=9),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.3,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.spacing=unit(0.5,"mm"),panel.grid.major=element_line(linewidth=0.2))
g_horiz_type_trans2 <- ggplot(subset(all_type_trans_flipped,origin=="transcriptome"),
                                    aes(y=genome,x=percent,fill=Type))+
  geom_col(position="fill",alpha=0.7)+
  geom_text(aes(y=genome,x=percent,label=n,family="work"),size=2,position=position_fill(vjust=0.5))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0))+
  guides(fill=guide_legend(reverse=TRUE,nrow=2))+
  ggforce::facet_col(vars(Order),scales="free_y",space="free")+ theme_bw()+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position="bottom",legend.text=element_text(size=9,margin=margin(l=-4)),
        legend.key.size=unit(1,"mm"),legend.margin=margin(-10,0,-5,-10),
        strip.text=element_text(size=11,angle=0,face='bold',hjust=0,margin=margin(0,0,0.1,0)),
        strip.background=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(face='bold.italic',size=8),axis.text.x=element_text(size=9),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.3,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.spacing=unit(0.5,"mm"),panel.grid.major=element_line(linewidth=0.2))
dev.off()
png("export/mappedUCE_Types_allgenome_transcriptomeAware.png",width=14,height=5,units="in",res=300)
g_horiz_type_trans <- grid.arrange(textGrob("(A) Genome only",hjust=2.4,gp=gpar(fontsize=10,fontfamily="work",fontface="bold")),
                                   textGrob("(B) Transcriptome only",hjust=1.8,gp=gpar(fontsize=10,fontfamily="work",fontface="bold")),
                                   g_horiz_type_trans1,g_horiz_type_trans2,
                                   nrow=2,heights=c(1,13),widths=c(1,1),padding=unit(-1,"mm"))
dev.off()

## Distances from previous UCE ----
g_hist1 <- ggplot(subset(all_scaf,corrected_Distance <= 1000 & corrected_Distance > 0),aes(x=corrected_Distance))+
  geom_histogram(colour='black',fill='lightgrey',binwidth=20,linewidth=0.1)+
  scale_y_continuous(n.breaks=4)+
  labs(y="Frequency",x="",fill="")+theme_bw()+ #Distances between UCE loci
  theme(text=element_text(family="work",size=10),plot.margin=unit(c(2,0,-3,2),"mm"),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.5,"mm"),
        panel.border=element_rect(linewidth=0.2),panel.grid=element_line(linewidth=0.2))
g_hist2 <- ggplot(subset(all_scaf,corrected_Distance > 1000 & corrected_Distance < 1000000),aes(x=corrected_Distance))+ 
  geom_histogram(colour='black',fill='lightgrey',binwidth=20000,linewidth=0.1)+
  scale_x_continuous(labels=scales::comma)+
  labs(y="",x="",fill="")+theme_bw()+  #Distances between UCE loci
  theme(text=element_text(family="work",size=10),plot.margin=unit(c(2,2,-3,0),"mm"),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.5,"mm"),
        panel.border=element_rect(linewidth=0.2),panel.grid=element_line(linewidth=0.2))
g_hist <- grid.arrange(g_hist1,g_hist2,nrow=1,widths=c(1,2))
dev.off()
png("export/mappedUCE_distances.png",width=8,height=5,units="in",res=300)
g_hist <- grid.arrange(g_hist,textGrob("Distance (bp) from previous UCE locus on same scaffold",
                                       gp=gpar(fontsize=10,fontfamily="work")),
                       nrow=2,heights=c(12,1),padding=unit(-2,"mm"))
dev.off()

## Common UCEs ----
common_uce <- common_uce_chromo %>% select(-c(n,Montipora_capitata_C)) %>% 
  pivot_longer(names_to="genome",values_to="Chr",
               cols=c(Acropora_millepora_G:Nematostella_vectensis_G)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Acropora_hyacinthus_G=factor(Acropora_hyacinthus_G,levels=paste0("Chr",1:15))) %>% 
  filter(Acropora_hyacinthus_G!="NA")
g_common_uce <- ggplot(common_uce, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+
  facet_wrap(~Acropora_hyacinthus_G)
g_common_uce <- g_common_uce + geom_point(data=common_uce,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ #c("#4292c6","#6a51a3","#ff7f00","#41ab5d")
  scale_y_continuous(breaks=c(1:14,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-8)),
        axis.text.x=element_blank())
ggsave("export/mappedUCE_common_genomeWchromosome.png",plot=g_common_uce,
       device="png",width=8,height=8,units="in",dpi=300)

## UCEs in phylogeny ----
uce_inphylogeny_toplot <- uce_inphylogeny %>% select(-c(Montipora_capitata_C)) %>% 
  pivot_longer(names_to="genome",values_to="Chr",
               cols=c(Acropora_millepora_G:Nematostella_vectensis_G)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Acropora_hyacinthus_G=factor(Acropora_hyacinthus_G,levels=paste0("Chr",1:15))) %>% 
  filter(Acropora_hyacinthus_G!="NA")

g_uce_inphylo <- ggplot(uce_inphylogeny_toplot, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+
  facet_wrap(~Acropora_hyacinthus_G)
g_uce_inphylo <- g_uce_inphylo + 
  geom_point(data=uce_inphylogeny_toplot,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ #c("#4292c6","#6a51a3","#ff7f00","#41ab5d")
  scale_y_continuous(breaks=c(1:14,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-8)),
        axis.text.x=element_blank())
ggsave("export/mappedUCE_inphylogeny_genomeWchromosome.png",plot=g_uce_inphylo,
       device="png",width=8,height=8,units="in",dpi=300)

# loose ----
## loci prop: Stacked ----
all_loose_type_flipped = all_loose_type %>% 
  mutate(genome=as.factor(genome),
         genome=factor(genome,levels=rev(levels(genome))),
         Type=factor(Type,levels=rev(levels(Type))),
         percent=percent*100)
g_horiz_type_loose <- ggplot(all_loose_type_flipped,aes(y=genome,x=percent,fill=Type))+
  geom_col(position="fill",alpha=0.7)+
  geom_text(aes(y=genome,x=percent,label=n,family="work"),size=3,position=position_fill(vjust=0.5))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0))+
  guides(fill=guide_legend(reverse=TRUE,nrow=2))+
  ggforce::facet_col(vars(Order),scales="free_y",space="free")+
  theme_bw()+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position="bottom",legend.text=element_text(size=9,margin=margin(l=-4)),
        legend.key.size=unit(1,"mm"),legend.margin=margin(-10,0,-5,-10),
        strip.text=element_text(size=11,angle=0,face='bold',hjust=0,margin=margin(0,0,0.1,0)),
        strip.background=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(face='bold.italic',size=9),axis.text.x=element_text(size=9),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.3,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.spacing=unit(0.5,"mm"),panel.grid.major=element_line(linewidth=0.2))
ggsave("export/mappedUCE_Types_allgenome_loose.png",plot=g_horiz_type_loose,device="png",width=6,height=5,units="cm",dpi=500)

## loci prop: transcriptome separated ----
all_loose_type_trans_flipped = all_loose_type_trans %>% 
  mutate(genome=as.factor(genome),
         genome=factor(genome,levels=rev(levels(genome))),
         Type=factor(Type,levels=rev(levels(Type))),
         percent=percent*100)
g_horiz_type_trans_loose1 <- ggplot(subset(all_loose_type_trans_flipped,origin=="genome"),
                                    aes(y=genome,x=percent,fill=Type))+
  geom_col(position="fill",alpha=0.7)+
  geom_text(aes(y=genome,x=percent,label=n,family="work"),size=2,position=position_fill(vjust=0.5))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0))+
  guides(fill=guide_legend(reverse=TRUE,nrow=2))+
  ggforce::facet_col(vars(Order),scales="free_y",space="free")+ theme_bw()+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position="bottom",legend.text=element_text(size=9,margin=margin(l=-4)),
        legend.key.size=unit(1,"mm"),legend.margin=margin(-10,0,-5,-10),
        strip.text=element_text(size=11,angle=0,face='bold',hjust=0,margin=margin(0,0,0.1,0)),
        strip.background=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(face='bold.italic',size=8),axis.text.x=element_text(size=9),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.3,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.spacing=unit(0.5,"mm"),panel.grid.major=element_line(linewidth=0.2))
g_horiz_type_trans_loose2 <- ggplot(subset(all_loose_type_trans_flipped,origin=="transcriptome"),
                                    aes(y=genome,x=percent,fill=Type))+
  geom_col(position="fill",alpha=0.7)+
  geom_text(aes(y=genome,x=percent,label=n,family="work"),size=2,position=position_fill(vjust=0.5))+
  scale_fill_manual(values=hex_values,na.translate = FALSE)+
  scale_x_continuous(expand=c(0.01,0))+
  guides(fill=guide_legend(reverse=TRUE,nrow=2))+
  ggforce::facet_col(vars(Order),scales="free_y",space="free")+ theme_bw()+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position="bottom",legend.text=element_text(size=9,margin=margin(l=-4)),
        legend.key.size=unit(1,"mm"),legend.margin=margin(-10,0,-5,-10),
        strip.text=element_text(size=11,angle=0,face='bold',hjust=0,margin=margin(0,0,0.1,0)),
        strip.background=element_blank(), axis.title=element_blank(),
        axis.text.y=element_text(face='bold.italic',size=8),axis.text.x=element_text(size=9),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.3,"mm"),
        panel.border=element_rect(linewidth=0.1),panel.spacing=unit(0.5,"mm"),panel.grid.major=element_line(linewidth=0.2))
dev.off()
png("export/mappedUCE_Types_allgenome_loose_transcriptomeAware.png",width=14,height=5,units="in",res=300)
g_horiz_type_trans_loose <- 
  grid.arrange(textGrob("(A) Genome only",hjust=2.4,gp=gpar(fontsize=10,fontfamily="work",fontface="bold")),
               textGrob("(B) Transcriptome only",hjust=1.8,gp=gpar(fontsize=10,fontfamily="work",fontface="bold")),
               g_horiz_type_trans_loose1,g_horiz_type_trans_loose2,
               nrow=2,heights=c(1,13),widths=c(1,1),padding=unit(-1,"mm"))

dev.off()

## Distances ----
g_loose_hist1 <- ggplot(subset(all_loose_scaf,corrected_Distance <= 1000 & corrected_Distance > 0),aes(x=corrected_Distance))+
  geom_histogram(colour='black',fill='lightgrey',binwidth=20,linewidth=0.1)+
  scale_y_continuous(n.breaks=4)+
  labs(y="Frequency",x="",fill="")+theme_bw()+ #Distances between UCE loci
  theme(text=element_text(family="work",size=10),plot.margin=unit(c(2,0,-3,2),"mm"),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.5,"mm"),
        panel.border=element_rect(linewidth=0.2),panel.grid=element_line(linewidth=0.2))
g_loose_hist2 <- ggplot(subset(all_loose_scaf,corrected_Distance > 1000 & corrected_Distance < 1000000),aes(x=corrected_Distance))+ 
  geom_histogram(colour='black',fill='lightgrey',binwidth=20000,linewidth=0.1)+
  scale_x_continuous(labels=scales::comma)+
  labs(y="",x="",fill="")+theme_bw()+  #Distances between UCE loci
  theme(text=element_text(family="work",size=10),plot.margin=unit(c(2,2,-3,0),"mm"),
        axis.ticks=element_line(linewidth=0.1),axis.ticks.length=unit(0.5,"mm"),
        panel.border=element_rect(linewidth=0.2),panel.grid=element_line(linewidth=0.2))
g_loose_hist <- grid.arrange(g_loose_hist1,g_loose_hist2,nrow=1,widths=c(1,2))
dev.off()
png("export/mappedUCE_distances_loose.png",width=8,height=5,units="in",res=300)
g_loose_hist <- grid.arrange(g_loose_hist,textGrob("Distance (bp) from previous UCE locus on same scaffold",
                                                   gp=gpar(fontsize=10,fontfamily="work")),
                       nrow=2,heights=c(12,1),padding=unit(-2,"mm"))
dev.off()

## Common UCEs ----
common_uce_loose <- common_uce_chromo_loose %>% select(-c(n,Montipora_capitata_C)) %>% 
  pivot_longer(names_to="genome",values_to="Chr",
               cols=c(Acropora_millepora_G:Nematostella_vectensis_G)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Acropora_hyacinthus_G=factor(Acropora_hyacinthus_G,levels=paste0("Chr",1:14))) %>% 
  filter(Acropora_hyacinthus_G!="NA")
g_common_uce_loose <- ggplot(common_uce_loose, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+
  facet_wrap(~Acropora_hyacinthus_G)
g_common_uce_loose <- g_common_uce_loose + geom_point(data=common_uce_loose,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ #c("#4292c6","#6a51a3","#ff7f00","#41ab5d")
  scale_y_continuous(breaks=c(1:15,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-8)),
        axis.text.x=element_blank())
ggsave("export/mappedUCE_common_genomeWchromosome_loose.png",plot=g_common_uce_loose,
       device="png",width=8,height=8,units="in",dpi=300)

## UCEs in phylogeny ----
uce_inphylogeny_loose_toplot <- uce_inphylogeny_loose %>% select(-c(Montipora_capitata_C)) %>% 
  pivot_longer(names_to="genome",values_to="Chr",
               cols=c(Acropora_millepora_G:Nematostella_vectensis_G)) %>% 
  mutate(UCE=as.numeric(UCE),
         Chr=as.numeric(sub("Chr","",Chr)),
         genome=as.factor(genome),
         Acropora_hyacinthus_G=factor(Acropora_hyacinthus_G,levels=paste0("Chr",1:15))) %>% 
  filter(Acropora_hyacinthus_G!="NA")

g_uce_inphylo_loose <- ggplot(uce_inphylogeny_loose_toplot, aes(y=0,x=log(UCE)))+
  geom_point(size=1)+
  facet_wrap(~Acropora_hyacinthus_G)
g_uce_inphylo_loose <- g_uce_inphylo_loose + 
  geom_point(data=uce_inphylogeny_loose_toplot,aes(y=Chr,x=log(UCE),colour=genome))+
  scale_colour_manual(values=rev(hex_values))+ #c("#4292c6","#6a51a3","#ff7f00","#41ab5d")
  scale_y_continuous(breaks=c(1:14,1))+
  labs(x="UCE numbers (log) *Not actual position on the chromosome",y="",colour="")+ theme_bw()+
  theme(text=element_text(family="work",size=20),
        legend.position="bottom",legend.text=element_text(margin=margin(l=-8)),
        axis.text.x=element_blank())
ggsave("export/mappedUCE_inphylogeny_genomeWchromosome_loose.png",plot=g_uce_inphylo_loose,
       device="png",width=8,height=8,units="in",dpi=300)
