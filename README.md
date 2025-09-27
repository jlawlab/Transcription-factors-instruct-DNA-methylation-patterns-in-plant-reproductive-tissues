### README
Code and example datasets used in the publication Xu et al. Transcription factors instruct DNA methylation patterns in plant reproductive tissues (2025) Nature Cell Biology. Full datasets are available as Supplementary Data.

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
- Description: This code generates volcano plot showing the differentially expressed siRNA clusters in rim22. Clusters that are downregulated, unaffected, and upregulated are shown as blue, black, and red circles, respectively. The number of clusters in each category is indicated in the correspondingly colored boxes. Clusters overlap with methyl-cutting assay targets are highlighted. This is the base script to generate volcano plots in 
- author: Guanghui Xu
