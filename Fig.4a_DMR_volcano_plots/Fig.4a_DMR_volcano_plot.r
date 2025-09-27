#!/usr/bin/env Rscript
# Author: Guanghui Xu, Law Lab, Salk Institute
# Read me: This code generates volcano plot of DMRs. DMRs that overlap HyperTE loci, siren loci, both, or neither are colored in blue, pink, purple, or black, respectively. The total number of hypo and hyper DMRs is indicated in black boxes in the center. The number of DMRs overlap with the above categories are shown in the correspondingly colored boxes. 

################################# Load required packages
################################# 
library(ggplot2)
library(reshape2)
source("Fig.4a_DMR_volcano_plot_color_settings.r")

################################# set the parameters
#################################
fileName = data.frame(csv.file.name=c("Ovules_rim-cr_CHH_allDMR.bed", "Ovules_rim16_CHH_allDMR.bed", "Ovules_rim22_CHH_allDMR.bed", "Anthers_rim-cr_CHH_allDMR.bed", "Anthers_rim16_CHH_allDMR.bed", "Anthers_rim22_CHH_allDMR.bed"), genotype_id=c("Ov rim-cr DMRs", "Ov rim16 DMRs", "Ov rim22 DMRs", "An rim-cr DMRs", "An rim16 DMRs", "An rim22 DMRs"))
genotype_id=unique(as.character(fileName$genotype_id))

################################# 
#################################
data.rbind = NULL
stat.all = NULL
for (i in 1:length(genotype_id)){
#for (i in 1:2){
print(i)
 genotype.one=genotype_id[i]
 # input data
 data.one=read.csv(as.character(fileName$csv.file.name)[i],sep="\t", header=F)
 names(data.one) = c("chr", "start", "end", "ID", "padj", "WT", "MUT", "DIFF", "context", "DMR", "locus_type","location")
 data.one$geno=genotype_id[i]
 
 # UpDownNum, FC, P, RPKM range
  UpDownNum = table(data.one$DMR)
  UpDownNum2=c(UpDownNum["hyperDMR"], UpDownNum["hypoDMR"])
  names(UpDownNum2) = c("hyperDMR", "hypoDMR")
  UpDownNum2[is.na(UpDownNum2)]=0

  maxC=max(data.one[,c("WT","MUT")], na.rm=T)
  minC=min(data.one[,c("WT","MUT")], na.rm=T)
  maxDiff=max(data.one[,c("DIFF")], na.rm=T)
  minDiff=min(data.one[,c("DIFF")], na.rm=T)
  maxP=max(-log10(data.one[data.one$padj>0,]$padj), na.rm=T)
 
 stat.one = c(UpDownNum2, maxC=maxC, minC=minC, maxDiff=maxDiff, minDiff=minDiff, maxP=maxP, geno=genotype.one)
 stat.all = rbind(stat.all, stat.one)
 data.rbind=rbind(data.rbind,data.one)
 }

data.rbind$ID = as.character(data.rbind$ID)
data.rbind$context = as.character(data.rbind$context)
data.rbind$DMR = as.character(data.rbind$DMR)
data.rbind$locus_type = as.character(data.rbind$locus_type)
data.rbind$location = as.character(data.rbind$location)
data.rbind$locus_type[data.rbind$locus_type=="."]="other"

stat.all = data.frame(stat.all)
stat.all$hyperDMR = as.numeric(as.character(stat.all$hyperDMR))
stat.all$hypoDMR = as.numeric(as.character(stat.all$hypoDMR))
stat.all$maxC = as.numeric(as.character(stat.all$maxC))
stat.all$minC = as.numeric(as.character(stat.all$minC))
stat.all$maxDiff = as.numeric(as.character(stat.all$maxDiff))
stat.all$minDiff = as.numeric(as.character(stat.all$minDiff))
stat.all$maxP = as.numeric(as.character(stat.all$maxP))
stat.all$geno = as.character(stat.all$geno)

stat.melt=melt(stat.all,id=c("geno", "maxC", "minC", "maxDiff", "minDiff", "maxP"))
names(stat.melt)=c("geno","maxC","minC","maxDiff", "minDiff","maxP", "DMR","nDMR")

stat.melt$maxC = as.numeric(as.character(stat.melt$maxC))
stat.melt$minC = as.numeric(as.character(stat.melt$minC))
stat.melt$maxDiff = as.numeric(as.character(stat.melt$maxDiff))
stat.melt$minDiff = as.numeric(as.character(stat.melt$minDiff))
stat.melt$maxP = as.numeric(as.character(stat.melt$maxP))
stat.melt$geno = as.character(stat.melt$geno)
stat.melt$DMR = as.character(stat.melt$DMR)
stat.melt$nDMR = as.numeric(as.character(stat.melt$nDMR))

################################# x,y coordinate of geom_text
#################################
xrange = max(stat.melt$maxDiff, abs(stat.melt$minDiff))
yrange = max(stat.melt$maxP)

# hypoDMR
  stat.melt.hypoDMR=subset(stat.melt,DMR=="hypoDMR")
  stat.melt.hypoDMR$x = stat.melt.hypoDMR$minDiff*0.75
  stat.melt.hypoDMR$y = stat.melt.hypoDMR$maxP*0.75
  stat.melt.hypoDMR$x2 = -xrange*0.75
  stat.melt.hypoDMR$y2 = yrange*0.75

# hyperDMR
  stat.melt.hyperDMR=subset(stat.melt,DMR=="hyperDMR")
  stat.melt.hyperDMR$x = stat.melt.hyperDMR$maxDiff*0.75
  stat.melt.hyperDMR$y = stat.melt.hyperDMR$maxP*0.75
  stat.melt.hyperDMR$x2 = xrange*0.75
  stat.melt.hyperDMR$y2 = yrange*0.75
  #stat.melt.ClusterDown$x=5

stat.melt.res=rbind(stat.melt.hypoDMR,stat.melt.hyperDMR)
stat.melt.res$ID = paste(stat.melt.res$geno, stat.melt.res$DMR, sep="_")

################################# add hyperTE and siren info
#################################
DE.num = data.frame(table(data.rbind[, names(data.rbind) %in% c("geno", "DMR","locus_type")]))
DE.num = subset(DE.num, DMR!="negDMR")
DE.num$ID = paste(DE.num$geno, DE.num$DMR, sep="_")
DE.num2 = cbind(DE.num, stat.melt.res[match(DE.num$ID, stat.melt.res$ID),c("maxC", "minC", "maxDiff", "minDiff", "maxP")])

DE.num.hypo = subset(DE.num2, DMR=="hypoDMR")
DE.num.hyper = subset(DE.num2, DMR=="hyperDMR")

DE.num.hypo$x = DE.num.hypo$minDiff*0.95
DE.num.hyper$x = DE.num.hyper$maxDiff*0.95
DE.num.hypo$x2 = -xrange*0.95
DE.num.hyper$x2 = xrange*0.95
DE.num3 = rbind(DE.num.hypo, DE.num.hyper)

DE.num.y = NULL
DE.num.y[DE.num3$locus_type=="hyperTE"] =subset(DE.num3, locus_type=="siren")$maxP*0.9
DE.num.y[DE.num3$locus_type=="siren"] = subset(DE.num3, locus_type=="hyperTE")$maxP*0.7
DE.num.y[DE.num3$locus_type=="hyperTE_w_siren"] = subset(DE.num3, locus_type=="hyperTE_w_siren")$maxP*0.5
DE.num.y[DE.num3$locus_type=="otherclsy34"] = subset(DE.num3, locus_type=="otherclsy34")$maxP*0.3
DE.num.y[DE.num3$locus_type=="other"] = subset(DE.num3, locus_type=="other")$maxP*0.1
DE.num3$y = DE.num.y

DE.num.y2 = NULL
DE.num.y2[DE.num3$locus_type=="hyperTE"] = yrange*0.9
DE.num.y2[DE.num3$locus_type=="siren"] = yrange*0.7
DE.num.y2[DE.num3$locus_type=="hyperTE_w_siren"] = yrange*0.5
DE.num.y2[DE.num3$locus_type=="otherclsy34"] = yrange*0.3
DE.num.y2[DE.num3$locus_type=="other"] = yrange*0.1
DE.num3$y2 = DE.num.y2

data.rbind.negDMR = subset(data.rbind, DMR=="negDMR")
data.rbind.posDMR = subset(data.rbind, DMR!="negDMR")
data.rbind.posDMR = data.rbind.posDMR[sample(1:nrow(data.rbind.posDMR), nrow(data.rbind.posDMR), replace=F),]

data.rbind.negDMR$geno=factor(data.rbind.negDMR$geno,levels=as.character(fileName$genotype_id))
data.rbind.posDMR$geno=factor(data.rbind.posDMR$geno,levels=as.character(fileName$genotype_id))
data.rbind.posDMR$locus_type=factor(data.rbind.posDMR$locus_type,levels=c("hyperTE", "siren", "hyperTE_w_siren", "otherclsy34", "other"))

stat.melt.res$geno = factor(stat.melt.res$geno,levels=as.character(fileName$genotype_id))
DE.num3$geno = factor(DE.num3$geno,levels=as.character(fileName$genotype_id))
DE.num3$locus_type = factor(DE.num3$locus_type,levels=c("hyperTE", "siren", "hyperTE_w_siren", "otherclsy34", "other"))
DE.num3$DMR = as.character(DE.num3$DMR)
DE.num3$DMR = factor(DE.num3$DMR, levels=c("hypoDMR", "hyperDMR"))

################################# volcano plot
#################################
# volcano plot with fixed y axis range
p1=ggplot(data.rbind.negDMR,aes(DIFF, -log(padj,10))) +
    geom_point(alpha=0.2,size=0.4, color="darkgrey") +
    geom_hline(yintercept= -log(0.01, 10),linetype="dotted",size=0.25,col="blue") +
    geom_vline(xintercept= c(-0.1, 0.1),linetype="dotted",size=0.25,col="blue") +
    geom_point(data=data.rbind.posDMR, aes(DIFF,-log(padj,10), color=locus_type), shape=16, size=0.4, alpha=0.75) + 
    labs(title="DMR volcano plot") + 
    theme_bw(base_size = 7) +
    coord_cartesian(xlim=c(-xrange, xrange)) +
    geom_label(data=stat.melt.res,fill="white",show.legend=FALSE,aes(x2,y2,label=nDMR),size=1.5)+
    geom_label(data=DE.num3, fill="white",show.legend=FALSE,aes(x2,y2,label=Freq, col=locus_type),size=1.5,label.padding = unit(0.1, "lines"))+
    facet_wrap(~geno, scales="fixed", dir="v", ncol = 2) +
    scale_color_manual(values=c("#5CC9E2","#EA0A8C","purple", "gray20")) +
    xlab("mCHH difference")+
    ylab("-log10(p)") +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        #axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(size = 0.25),
        panel.grid = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25),
        panel.grid.minor = element_line(size = 0.25),
        strip.background = element_rect(size = 0.25))

pdf("Fig.4a_DMR_volcano_plot.pdf", width=3, height=4)
 plot(p1)
dev.off()
            
