### README
Code and example datasets used in the publication **Xu et al. Transcription factors instruct DNA methylation patterns in plant reproductive tissues (2025) Nature Cell Biology.** Full datasets are available as Supplementary Data.

### Transcription factors instruct DNA methylation patterns in plant reproductive tissues

##### Guanghui Xu<sup>1</sup>, Yuhan Chen<sup>1</sup>, Laura M. Martins<sup>1</sup>, En Li<sup>1</sup>, Fuxi Wang<sup>1</sup>, Tulio Magana<sup>1</sup>, Junlin Ruan<sup>1</sup>, and Julie A. Law<sup>1,2*</sup>

<sup>1</sup> Plant Molecular and Cellular Biology Laboratory, Salk Institute for Biological Studies, La Jolla, CA, 92037, USA

<sup>2</sup> Division of Biological Sciences, University of California, San Diego, La Jolla, CA 92093

*Corresponding Author. Email: jlaw@salk.edu

**Abstract**

DNA methylation is maintained by forming self-reinforcing connections with other repressive chromatin modifications, resulting in stably silenced genes and transposons. However, these mechanisms fail to explain how new methylation patterns are generated. In Arabidopsis, CLASSY3 (CLSY3) targets the RNA-directed DNA methylation (RdDM) machinery to different loci in reproductive tissues, generating distinct epigenomes via unknown mechanism(s). Here, we discovered that several different REPRODUCTIVE MERISTEM (REM) transcription factors are required for methylation at CLSY3 targets specific to anther or ovule tissues. We designate these factors as REM INSTRUCT METHYLATION (RIMs) and demonstrate that disruption of their DNA binding domains, or the motifs they recognize, blocks RdDM. Furthermore, we demonstrate that mis-expression of RIM12 is sufficient to initiate siRNA production at ovule targets in anthers. These findings reveal a critical role for genetic information in targeting DNA methylation in reproductive tissues, expanding our understanding of how methylation is regulated to include inputs from both genetic and epigenetic information.

**Keywords:** DNA methylation, siRNAs, REM transcription factors, RNA-directed DNA methylation, CLSY3 

### Available scripts
#### Fig. 1b siRNA volcano plot
- Fig.1b_volcano_plot.r
- Description: This code generates volcano plot showing the differentially expressed siRNA clusters in rim22. Clusters that are downregulated, unaffected, and upregulated are shown as blue, black, and red circles, respectively. The number of clusters in each category is indicated in the correspondingly colored boxes. Clusters overlap with methyl-cutting assay targets are highlighted. This is the base script used to generate volcano plots for Figures 1b, 1c, 2c, 3b, 3f, 6c, 7g, and Extended Data Figures 1g, 2b, 2c, 2d, 2e, 4a, 5b, 5d, 6f, 6g, 6j.
- Example dataset: Data from Figure 1b.
- Author: Guanghui Xu (Law lab, Salk Institute)

#### Fig. 1f siRNA boxplot
- Fig.1f_boxplot.r
- Fig.1f_boxplot_color_settings.r
- Description: This code generates boxplots comparing siRNA levels across multiple genotypes. This is the base script used to generate boxplots for Figures 1f, 2d, 2e, 3d, 3h, 7e, and Extended Data Figures 3a, 3d, 4c, 4d, 5c.
- Example dataset: Data from Figure 1f.
- Example dataset: Data from Figure 1f.
- Author: Guanghui Xu (Law lab, Salk Institute)

#### Fig. 4a DMR volcano plot
- Fig.4a_DMR_volcano_plot.r
- Fig.4a_DMR_volcano_plot_color_settings.r
- Description: This code generates volcano plot of DMRs. DMRs that overlap HyperTE loci, siren loci, both, or neither are colored in blue, pink, purple, or black, respectively. The total number of hypo and hyper DMRs is indicated in black boxes in the center. The number of DMRs overlap with the above categories are shown in the correspondingly colored boxes. This is the base script used to generate DMR volcano plots for Figures 4a, 4b, and Extended Data Figures 4a.
- Example dataset: Data from Figure 4a.
- Author: Guanghui Xu (Law lab, Salk Institute)

#### Fig. 5b siRNA heatmap
- Fig.5b_heatmap.r
- Description: This script generates heatmaps based on Fold Change valuse got from siRNA DEseq analysis. This is the base script used to generate siRNA heatmaps for Figures 5b and 5d
- Example dataset: Data from Figure 5b.
- Author: Guanghui Xu (Law lab, Salk Institute)

#### Extended Data Fig. 2a REM phylogenetic tree
- Extended_Data_Fig.2a_Phylo_Tree_Visualization_by_L.M..r
- Extended_Data_Fig.2a_REM_phylogenetic_tree_by_L.M..bash
- Description: This analysis reconstructs a phylogeny for 46 REM genes using two custom scripts.
- Author: Laura M. Martins (Law lab, Salk Institute)
