#!/usr/bin/env Rscript
# Author: Guanghui Xu, Law Lab, Salk Institute

################################# Load required packages
#################################
library(gplots)
library(RColorBrewer)

################################# set parameters
################################# 
fileName = data.frame(csv.file.name=c("DESeq_24ntsiRNA_cluster_Anthers_pol-iv.txt",
"DESeq_24ntsiRNA_cluster_Anthers_clsy3.txt",
"DESeq_24ntsiRNA_cluster_Anthers_clsy34.txt",
"DESeq_24ntsiRNA_cluster_Anthers_clsy12.txt",
"DESeq_24ntsiRNA_cluster_Anthers_rim16.txt",
"DESeq_24ntsiRNA_cluster_Anthers_rim22.txt",
"DESeq_24ntsiRNA_cluster_Anthers_rim-cr.txt",
"DESeq_24ntsiRNA_cluster_Anthers_rim-q.txt"), genotype_id=c("pol-iv", "clsy3", "clsy3,4", "clsy1,2", "rim16", "rim22", "rim-cr", "rim-q"))

subgroups = c(
"siRNA_clusters_rim16_wrim22_n241.bed",
"siRNA_clusters_rim22_notrim16_n233.bed",
"siRNA_cluster_rim16_notrim22_n28.bed",
"siRNA_clusters_rim-q_only_n103.bed",
"siRNA_clusters_HyperTE_only_n110.bed")

# set the column, based on which the ranking algrithom will be conducted
within_group_ranking_method = "strength" # "strength" or "clustering"
if (within_group_ranking_method=="strength"){
    clustering_samples = list()
    clustering_samples[[1]] = "rim-q"
    clustering_samples[[2]] = "rim22"
    clustering_samples[[3]] = "rim16"
    clustering_samples[[4]] = "rim-q"
    clustering_samples[[5]] = "clsy3,4"
} else if (within_group_ranking_method=="clustering") {
    clustering_samples = list()
    clustering_samples[[1]] = c("rim16","rim22")
    clustering_samples[[2]] = c("rim22", "rim-q")
    clustering_samples[[3]] = c("rim16","rim22")
    clustering_samples[[4]] = c("rim-q")
    clustering_samples[[5]] = c("clsy3", "clsy3,4")
}

################################# rank the Clusters within the subgroups; then combine results from different subgroups
#################################
ClusterID.list = list()
res.all = NULL

for (j in 1:length(subgroups)){
    print(j)
    subset = read.table(subgroups[j], sep="\t", header=F)
    ClusterID = as.character(subset$V4)
    ClusterID.list[[j]] = ClusterID

    # combine logFC info
    log2FC.info = NULL
    for (i in 1:nrow(fileName)){
        print(i)
        genotype.one = as.character(fileName$genotype_id)[i]
        DE.one = read.table(as.character(fileName$csv.file.name)[i], sep="\t", header=T)
        DE.one = DE.one[match(ClusterID, as.character(DE.one$Cluster)),]
        log2FC.info = cbind(log2FC.info, DE.one$log2FoldChange)
    }
    log2FC.info = data.frame(log2FC.info)
    names(log2FC.info) = fileName$genotype_id
    rownames(log2FC.info)=ClusterID
    log2FC.info=data.matrix(log2FC.info)

    # fake heatmap to get the ranking info
    if (within_group_ranking_method=="clustering"){
        p1=heatmap.2(log2FC.info[,colnames(log2FC.info) %in% clustering_samples[[j]]],col=bluered, 
            Rowv=T, 
            Colv=F, 
            scale="none", 
            margins = c(12, 2)) # margin of the x axis (bottom); margin of the y axis (right)
    }

    # rank log2FC.info
    if (within_group_ranking_method=="strength"){
        log2FC.rank = log2FC.info[order(log2FC.info[,clustering_samples[[j]]]), ]
    } else if (within_group_ranking_method=="clustering") {log2FC.rank = log2FC.info[p1$rowInd, ]}

    res.all =rbind(res.all, log2FC.rank)

}

################################# RowSideColors for subgroups
################################# 
col =  rep(c("#B2DF8A","#FB9A99","#1F78B4","#33A02C","#A6CEE3"),10)
ClusterID.all = rownames(res.all)
RowSideColors = NULL

subgroups.size = NULL
for (k in 1:length(subgroups)){
    RowSideColors[ClusterID.all %in% ClusterID.list[[k]]] = col[k]
    subgroups.size = c(subgroups.size, length(ClusterID.list[[k]]))
}

################################# RowSideColors.all
################################# 
RowSideColors.all = RowSideColors

################################# rowsep
################################# draw horizontal lines to distinguish different subgroups
subgroups.size2 = c(0, subgroups.size)

rowsep = NULL
rowsep[1]=0
for (n in 2:length(subgroups.size2)){
    dat.one=subgroups.size2[n] + rowsep[n-1]
    rowsep = c(rowsep, dat.one)
}

################################# heatmap
#################################
pdf ("Fig.5b_heatmap.pdf" ,width=5+nrow(fileName)/4,height=10)
    par(cex.main=1, font.main=1)
    p2=heatmap.2(res.all, 
        #col=colorpanel(20, "blue", "white", "red"),
        col=bluered,
        Rowv=NA, 
        Colv=NA, 
        scale="none", 
        trace="none",
        srtCol=60,
        cexRow=0.01,
        cexCol=1.5,
        RowSideColors = as.matrix(RowSideColors.all),
        key=TRUE,
        key.xlab="log2FC",
        colsep=1:nrow(fileName),
        #rowsep=rowsep,
        sepcolor=rep("white"),
        sepwidth=c(0.01,0.01),
        breaks = seq(-8,8,1),
        lhei=c(2,10), # layout of the column dengrogram and the heatmap
        lwid=c(4,10), # layout of the row dengrogram and the heatmap
        margins = c(12, 2) # margin of the x axis (bottom); margin of the y axis (right)
        )
    legend("left", col=col, pch=15, inset=c(-0, 0), c("rim16 & rim22 (241)", "rim22 not rim16 (233)", "rim16 not rim22 (28)", "rim-q only (103)", "HyperTE only (110)"), cex=0.6)
dev.off()


