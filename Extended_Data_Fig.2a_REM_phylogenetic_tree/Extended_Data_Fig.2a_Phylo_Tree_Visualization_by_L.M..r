#Written by Laura M. Martins — last update: 2025-03-21.


library(ape)
library(ggtree)
library(ggplot2)
library(svglite)


#Here it assumes that my file is a tab-delimited with columns "geneID" and "genesymbol".
mapping <- read.table("REM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


process_tree <- function(tree, mapping) {
  processed_labels <- sapply(tree$tip.label, function(label) {

    gene_part <- strsplit(label, "\\|")[[1]][1]
    gene_id <- sub("\\.[0-9]+$", "", gene_part)
    return(gene_id)
  })
  

  new_labels <- sapply(processed_labels, function(gene) {
    symbol <- mapping$genesymbol[mapping$geneID == gene]
    if (length(symbol) > 0) {
      paste0(gene, " (", symbol, ")")
    } else {
      gene
    }
  })
  tree$tip.label <- new_labels
  return(tree)
}


tree_iq <- read.tree("IQ_Tree2/REM_aligned.fas.treefile")
tree_iq <- process_tree(tree_iq, mapping)


#Below are two ways of representing the tree.
#First way:
plot_iq <- ggtree(tree_iq) +
  geom_tiplab(offset = 0.05) +
  geom_text2(aes(subset = !isTip, label = paste0(label)),
             hjust = -0.3, size = 2.5, color = "black") +
  ggtitle("IQ-TREE2 Phylogeny of REM Genes") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = unit(c(1, 3, 1, 1), "cm"))

ggsave("IQ_Tree2.png", plot = plot_iq, width = 16, height = 10, dpi = 300)
ggsave("IQ_Tree2.pdf", plot = plot_iq, width = 16, height = 10)
ggsave("IQ_Tree2.svg", plot = plot_iq, width = 16, height = 10)


#Second way:
plot_iq <- ggtree(tree_iq) +
  geom_tiplab(align = TRUE, offset = 0.05) +
  geom_text2(aes(subset = !isTip, label = paste0(label)),
             hjust = -0.3, size = 2.5, color = "black") +
  ggtitle("IQ-TREE2 Phylogeny of REM Genes") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  coord_cartesian(clip = "off")

ggsave("IQ_Tree2_2.png", plot = plot_iq, width = 12, height = 8, dpi = 300)
ggsave("IQ_Tree2_2.pdf", plot = plot_iq, width = 12, height = 8)
ggsave("IQ_Tree2_2.svg", plot = plot_iq, width = 12, height = 8)
