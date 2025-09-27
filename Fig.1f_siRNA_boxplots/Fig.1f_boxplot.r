#!/usr/bin/env Rscript
# Author: Guanghui Xu, Law Lab, Salk Institute

# Usage: Rscript Fig.1F_boxplot.r <rpkm_file> <sample_sheet> <bed_file> <output_prefix> <zoom_fraction>
# Arguments:
#  	<rpkm_file>      : Tab-delimited file containing RPKM values
#  	<sample_sheet>   : Metadata file specifying sample layout and grouping
#  	<bed_file>     	 : bed files for each subset
#  	<output_prefix>  : Prefix for output PDF and other results
#  	<zoom_factor>  	 : Numeric value (0 < x ≤ 1), e.g., 0.25 to zoom into 25% of max y-axis
#  	<withe_factor>   : adjust the output pdf width

################################# Load required packages
#################################
library(ggplot2)
library(reshape2)
library(ggh4x)
library(patchwork)
library(grid)
source("Fig.1f_boxplot_color_settings.r")

################################# set the parameters
#################################
#args = commandArgs(T)
args=c("Fig.1f_siRNA_rpkm_matrix.txt", "Fig.1f_sample_condition.txt", "Fig.1f_boxplot_panel_layout.txt", "Fig.1f", "0.1", "0.9")
rpkm_file = args[1]
sample_sheet_file = args[2]
bedfile_name = args[3]
prefix = args[4] # prefix of output file name
zoomY.index = as.numeric(as.character(args[5])) # 
widthFactor = as.numeric(as.character(args[6])) 

################################# input data
#################################
dat=read.table(rpkm_file, header=T)
sample_sheet = read.table(sample_sheet_file, sep="\t", header=T)
fileName = read.table(bedfile_name, sep="\t", header=T)
dat = dat[, names(dat) %in% c("Cluster", as.character(sample_sheet$Sample))]
genotype_ID = unique(as.character(sample_sheet$Condition))
col=as.character(my.color)[match(genotype_ID, names(my.color))]

################################# average
#################################
toBeAvg_index = sample_sheet$avg == "yes"
toBeAvg_sample = unique(sample_sheet$Sample[toBeAvg_index])

dat.mean = dat[, names(dat) %in% toBeAvg_sample]
dat.other = dat[, !(names(dat) %in% toBeAvg_sample)]
sample_sheet.mean = subset(sample_sheet, avg=="yes")
sample_sheet.other = subset(sample_sheet, avg=="no")

# mean
conditions = unique(sample_sheet.mean$Condition)
exp_mean = data.frame(Cluster = dat$Cluster)
for (cond in conditions) {
  samples = sample_sheet.mean$Sample[sample_sheet.mean$Condition == cond]
  exp_mean[[paste(cond,"mean",sep="_")]] = rowMeans(dat.mean[, samples, drop = FALSE])
}

# new sample_sheet
exp_mean.dataonly = exp_mean[, -1, drop=FALSE]
dat.new = cbind(dat.other, exp_mean.dataonly)
dat = dat.new

# new sample_sheet
sample_sheet.mean.new = data.frame(Sample=paste(conditions, "mean", sep="_"), Condition=conditions, batch="av", avg="yes")
sample_sheet.new = rbind(sample_sheet.other, sample_sheet.mean.new)[,1:3]
conditions.all = unique(sample_sheet$Condition)
sample_sheet.new$Condition = factor(sample_sheet.new$Condition, levels=conditions.all)
sample_sheet.new = sample_sheet.new[order(sample_sheet.new$Condition),]
sample_sheet = sample_sheet.new

################################# calculate whisker
#################################
whisker.all = NULL
for (j in 1:nrow(fileName)){
	print(j)
	print(as.character(fileName$sample)[j])
	sub = try(read.table(as.character(fileName$fileName)[j], sep="\t", header=F))
	if(class(sub) != "try-error"){
	print(paste("n=", nrow(sub)))
	if(nrow(sub)>0){
		dat.sub = subset(dat, Cluster %in% sub$V4)
		if(nrow(dat.sub)>0){		
		# calculate whisker
		whisker=dat.sub[,-which(names(dat.sub)=="Cluster")]
		whisker.max=max(boxplot(whisker)$stat)
		}
	}
	} else {whisker.max=NA}
	whisker.all= c(whisker.all, whisker.max)
}

################################# plot objects
#################################
plot.list.auto = list()
plot.list.manual = list()
plot.list.manual.reduced = list()
for (j in 1:nrow(fileName)){
	print(j)
	print(as.character(fileName$sample)[j])
	rowID = as.character(fileName$rowID)[j]
	columnID = as.character(fileName$columnID)[j]
	sub = try(read.table(as.character(fileName$fileName)[j], sep="\t", header=F))
	if(class(sub) != "try-error"){
		print(paste("n=", nrow(sub)))
		if(nrow(sub)>0){
			dat.sub = subset(dat, Cluster %in% sub$V4)
			if(nrow(dat.sub)>0){
				# calculate whisker
				whisker=dat.sub[,-which(names(dat.sub)=="Cluster")]
				whisker.max=max(boxplot(whisker)$stat)

				dat.melt = melt(dat.sub, id=c("Cluster"))
				names(dat.melt) = c("Cluster", "Sample", "fpkm")
				sample.pool = unique(as.character(dat.melt$Sample))
				sample_sheet2 = subset(sample_sheet, Sample %in% sample.pool)
				sample_sheet3 = sample_sheet2[match(as.character(dat.melt$Sample), as.character(sample_sheet2$Sample)),]
				dat.melt$genotype = sample_sheet3$Condition
				dat.melt$genotype = factor(dat.melt$genotype, levels=unique(as.character(sample_sheet2$Condition)))
				dat.melt$batch= sample_sheet3$batch
				dat.melt$batch = factor(dat.melt$batch, levels=unique(as.character(sample_sheet2$batch)))
				dat.melt$log2fpkm = log(dat.melt$fpkm +1, 2)
				dat.melt$Sample = factor(dat.melt$Sample, levels=unique(as.character(sample_sheet2$Sample)))
				dat.melt$rowID = rowID
				dat.melt$columnID = columnID

				p1 = ggplot(dat.melt) +
					geom_boxplot(notch = TRUE, notchwidth = 0.5, aes(interaction(batch, genotype), y=fpkm, fill=genotype), width = 0.7, outlier.shape=NA, linewidth=0.25) +
					#geom_label(data=fileName,fill="white",show.legend=FALSE,aes(label=sample),size=3, x=7, y=whisker.max*0.9)+
					#facet_grid(rowID~columnID) +
					coord_cartesian(ylim=c(0, whisker.max)) +
					scale_x_discrete(guide = guide_axis_nested()) +			
					scale_fill_manual(values=col)+
					labs(y = 'siRNAs (rpkm)') +
					labs(x ="") +
					labs(title=paste(as.character(fileName$sample)[j], " (n=", nrow(sub), ")",sep="")) + 
					theme_bw(base_size = 7) +
					theme(
					    legend.position = "none",
					    plot.title = element_text(hjust = 0.5),
					    plot.margin = margin(1, 1, 1, 1, "cm"),
					   # axis.line = element_line(linewidth = 0.5),
					    axis.ticks = element_line(linewidth = 0.25),
					    panel.border = element_rect(linewidth = 0.5),
					    panel.grid = element_line(linewidth = 0.25),
					    panel.grid.major = element_line(linewidth = 0.25),
					    panel.grid.minor = element_line(linewidth = 0.25),
					    strip.background = element_rect(linewidth = 0.25)
					  )

				plot.list.manual[[j]] = p1
			}
		}
	}
}

################################# specify layout
#################################
num.row = length(unique(fileName$rowID))
num.col = length(unique(fileName$columnID))

combined_plot.manual <- Reduce(`+`, plot.list.manual) + 
	  patchwork::plot_layout(ncol = num.col) &
	  theme(plot.margin = margin(15, 15, 15, 15), panel.spacing = unit(5, "lines"))

################################# output
#################################
pdf(paste(prefix, "_boxplot.pdf", sep=""), width=(3+(nrow(sample_sheet)-1)/8)*(num.col-1)*widthFactor, height=2.75)
	print(combined_plot.manual, newpage = FALSE)
dev.off()


