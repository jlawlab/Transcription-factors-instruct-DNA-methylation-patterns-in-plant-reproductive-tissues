#Written by Laura M. Martins — last update: 2025-03-18.

#This analysis reconstructs a phylogeny for 46 REM genes using two custom scripts:
#	LM_REM_phylogenetic_tree.bash — alignment and tree inference
#	LM_Phylo_Tree_Visualization.r — visualization


#	1. REM sequence alignment using MAFFT:

#		a. Downloading of fasta sequences:
#I downloaded a FASTA for the entire B3 TF family and then extracted only the REM TFs.

#conda install -c bioconda seqkit

cut -f1 REM.txt | tail -n +2 > rem_gene_ids.txt #seqkit only expect one column, so I decided to remove the second one.
sed -i 's/ .*//' rem_gene_ids.txt #I ran the two next lines to make sure I remove hidden characters.
tr -d '\r' < rem_gene_ids.txt > temp.txt && mv temp.txt rem_gene_ids.txt
sed 's/$/.*/' rem_gene_ids.txt > rem_gene_ids_regex.txt #I did not have the "." in my ID file so I had to add it with a * to match with my fasta file.
seqkit grep -r -f rem_gene_ids_regex.txt B3_family_seq.fas -o REM_seq.fas


#		b. Alignment using MAFFT:

#conda install -c bioconda mafft

# --auto: auto strategy selection
mafft --auto REM_seq.fas > REM_aligned.fas



#	2. Phylogenetic tree construction using IQ-TREE2 with ultrafast bootstrap:

#conda install -c bioconda iqtree

# -m MFP for ModelFinder Plus. This feature automatically tests a range of substitution models and selects the best-fit model for my data based on statistical criteria (such as BIC or AIC).
# -B 1000 runs ultrafast bootstrapping with 1,000 replicates.
# -nt AUTO for number of threads to use to be the most efficient. 

mkdir -p IQ_Tree2

iqtree2 -s REM_aligned.fas -m MFP -B 1000 -nt AUTO --prefix IQ_Tree2/REM_aligned.fas



#	3. Tree visualization using the following R packages: ape, ggtree, and ggplot2:

#conda-forge::r-ape
#conda install bioconda::bioconductor-ggtree
#conda install conda-forge::r-dplyr
#conda install conda-forge::r-ggplot2
#conda install conda-forge::r-patchwork

Rscript LM_Phylo_Tree_Visualization.r
