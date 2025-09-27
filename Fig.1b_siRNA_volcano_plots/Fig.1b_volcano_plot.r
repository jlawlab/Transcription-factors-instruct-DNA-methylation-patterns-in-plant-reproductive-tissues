#!/usr/bin/env Rscript
# Author: Guanghui Xu, Law Lab, Salk Institute
# Read me: This code generates volcano plot showing the differentially expressed siRNA clusters in rim22. Clusters that are downregulated, unaffected, and upregulated are shown as blue, black, and red circles, respectively. The number of clusters in each category is indicated in the correspondingly colored boxes. Clusters overlap with methyl-cutting assay targets are highlighted.

################################# Load required packages
#################################
library(ggplot2)
library(reshape2)

################################# set the parameters
#################################
DEseq_table = "DESeq_24ntsiRNA_cluster_rim22.txt"
highlight_clsy12 = "methyl-cutting_clsy12-dep.bed"
highlight_clsy3 = "methyl-cutting_clsy3-dep.bed"
genotype_id = "rim22"

################################# generate a list that contains the DE cluster ID
#################################
diff.id=list()
diff.id.up=list()
diff.id.none=list()
data.rbind=NULL
for (i in 1){
 print(i)
 genotype.one=genotype_id[i]
 # input data
 data.one=read.csv(DEseq_table, sep="\t", header=T)
 data.one$geno=genotype_id[i]
 data.rbind=rbind(data.rbind,data.one)
 
 # diff id
 down.index = (data.one$log2FoldChange < -1) & (data.one$padj<=0.01)
 up.index = (data.one$log2FoldChange > 1) & (data.one$padj<=0.01)
 down.index[is.na(down.index)]="FALSE"
 up.index[is.na(up.index)]="FALSE"
 down.index = as.logical(down.index)
 up.index = as.logical(up.index)
 none.index=(!down.index)&(!up.index) 

 diff.id[[genotype.one]]=as.character(subset(data.one,down.index)$Cluster)
 diff.id.up[[genotype.one]]=as.character(subset(data.one,up.index)$Cluster)
 diff.id.none[[genotype.one]]=as.character(subset(data.one,none.index)$Cluster)
}

################################# classification: ClusterUp ClusterDown None
#################################
up.index=(data.rbind$log2FoldChange > 1) & (data.rbind$padj <= 0.01)
down.index=(data.rbind$log2FoldChange < (-1)) & (data.rbind$padj <= 0.01)
down.index[is.na(down.index)]="FALSE"
up.index[is.na(up.index)]="FALSE"
down.index = as.logical(down.index)
up.index = as.logical(up.index)
none.index=(!up.index) & (!down.index)
class.DE=NULL
class.DE[up.index]="ClusterUp"
class.DE[down.index]="ClusterDown"
class.DE[none.index]="None"
data.rbind$Class=class.DE
data.rbind=data.rbind[!is.na(data.rbind$padj),]
# specify the genotype order that will present in the figure grid
data.rbind$geno=factor(data.rbind$geno,levels=as.character(genotype_id))

################################# stat
################################# calculate number of Down/Up/None smRNA cluster, fold change range and p value range
res.all=NULL
for(k in 1:length(genotype_id)){
  sample.one=subset(data.rbind,geno==genotype_id[k])
  
  # down & up index
  down.index = (sample.one$log2FoldChange < -1) & (sample.one$padj<=0.01)
  up.index = (sample.one$log2FoldChange > 1) & (sample.one$padj<=0.01)
  down.index[is.na(down.index)]="FALSE"
  up.index[is.na(up.index)]="FALSE"
  down.index = as.logical(down.index)
  up.index = as.logical(up.index)
  none.index=(!down.index)&(!up.index) 
  
  # down & up number
  down.num=sum(down.index,na.rm=T)
  up.num=sum(up.index,na.rm=T)
  none.num=nrow(data.one)-down.num-up.num
   
  # FC range
  maxFC=max(sample.one$log2FoldChange)
  minFC=min(sample.one$log2FoldChange)
  maxP=max(-log10(sample.one[sample.one$padj>0,]$padj))
  
  # 
  res.one=c(genotype_id[k],down.num, up.num, none.num, maxFC, minFC, maxP)
  res.all=rbind(res.all,res.one)
  }

  res.all=data.frame(res.all)
  names(res.all)=c("geno", "ClusterDown", "ClusterUp", "None", "maxFC", "minFC", "maxP")

stat.melt=melt(res.all,id=c("geno", "maxFC", "minFC", "maxP"))
names(stat.melt)=c("geno","maxFC","minFC","maxP","Class","nDE")

stat.melt$maxFC = as.numeric(as.character(stat.melt$maxFC))
stat.melt$minFC = as.numeric(as.character(stat.melt$minFC))
stat.melt$maxP = as.numeric(as.character(stat.melt$maxP))
stat.melt$nDE = as.numeric(as.character(stat.melt$nDE))

################################# x,y coordinate of geom_text
#################################
# ClusterDown
  stat.melt.ClusterDown=subset(stat.melt,Class=="ClusterDown")
  stat.melt.ClusterDown$x=min(stat.melt.ClusterDown$minFC)*0.8
  stat.melt.ClusterDown$y=round(stat.melt.ClusterDown$maxP)*0.75

# ClusterUp
  stat.melt.ClusterUp=subset(stat.melt,Class=="ClusterUp")
  stat.melt.ClusterUp$x=max(stat.melt.ClusterUp$maxFC)*0.8
  stat.melt.ClusterUp$y=round(stat.melt.ClusterUp$maxP)*0.75

# ClusterNone
  stat.melt.None=subset(stat.melt,Class=="None")
  stat.melt.None$x=round(stat.melt.None$minFC)*0.1
  stat.melt.None$y=round(stat.melt.None$maxP)*0.75
  stat.melt.None$x=0

stat.melt.res=rbind(stat.melt.ClusterDown,stat.melt.ClusterUp,stat.melt.None)

################################# plot met-cut assay associated Clusters
#################################
highlight1.dat = read.table(highlight_clsy12, sep="\t", header=F)
highlight1_ID=as.character(highlight1.dat$V4)
data.rbind.hlt1=subset(data.rbind, Cluster %in% highlight1_ID)
colnames(data.rbind.hlt1)=c("chr", "start", "end", "Cluster", "baseMean", "FC", "lfcSE", "stat", "pvalue", "adj", "DE", "geno", "Class")

highlight2.dat = read.table(highlight_clsy3, sep="\t", header=F)
highlight2_ID=as.character(highlight2.dat$V4)
data.rbind.hlt2=subset(data.rbind, Cluster %in% highlight2_ID)
colnames(data.rbind.hlt2)=c("chr", "start", "end", "Cluster", "baseMean", "FC", "lfcSE", "stat", "pvalue", "adj", "DE", "geno", "Class")

################################# volcano plot
#################################
p1=ggplot(data.rbind,aes(log2FoldChange,-log10(padj),col=Class)) +
    geom_point(alpha=0.2) +
    #geom_point(aes(x=2, y=30), color="red") +
    geom_point(data=data.rbind.hlt1, aes(FC,-log10(adj)), color="orange", shape=17, size=1) +
    geom_point(data=data.rbind.hlt2, aes(FC,-log10(adj)), color="yellow", shape=17, size=1) +
    labs(title=paste("siRNA clusters (n=", nrow(data.one), ")", sep="")) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=c(-1,1),linetype="dotted",size=1) +
    geom_hline(yintercept=-log10(0.01),linetype="dotted",size=1) +
    geom_label(data=stat.melt.res,fill="white",show.legend=FALSE,aes(x,y,label=nDE),size=3)+
    #facet_wrap(~geno, scales="free_y",ncol = as.numeric(1)) +
    scale_color_manual(values=c("blue","red","black")) +
    xlab("log2(FC)")+
    ylab("-log10(p)")

pdf("Fig.1b_volcano_plot.pdf",width=3.4, height=2.2)
  plot(p1)
dev.off()


