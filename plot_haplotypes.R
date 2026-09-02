#!/usr/bin/env Rscript
# Draw the gene structure of each haplotype next to the phylogenetic tree of
# the tagged sequences, as in Figure 4b.
#
# Usage:
#   Rscript plot_haplotypes.R <gene_tags.tsv> <tree.treefile> <output.pdf>
#
#   gene_tags.tsv   long table written by amy_haplotype_tags.py, with the
#                   columns haplotype, order, node, orientation, gene, length
#   tree.treefile   Newick tree written by IQ-TREE 2

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(gggenes)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_haplotypes.R <gene_tags.tsv> <tree.treefile> <output.pdf>")
}

tag_file <- args[1]
tree_file <- args[2]
out_file <- args[3]

tags <- read.delim(tag_file, stringsAsFactors = FALSE)
tree <- read.tree(tree_file)

# Order the haplotypes along the tree so that the two panels line up.
tip_order <- tree$tip.label[tree$edge[tree$edge[, 2] <= Ntip(tree), 2]]
tags <- tags[tags$haplotype %in% tip_order, ]
tags$haplotype <- factor(tags$haplotype, levels = rev(tip_order))

# Lay the gene copies end to end along each haplotype.
tags <- tags %>%
  arrange(haplotype, order) %>%
  group_by(haplotype) %>%
  mutate(
    end   = cumsum(length),
    start = end - length,
    forward = orientation == "+"
  ) %>%
  ungroup()

gene_colours <- c(
  AMY1  = "#E8413C",
  AMY2A = "#C4CC33",
  AMY2B = "#F5A623",
  AMYP1 = "#4CAF50"
)

structure_plot <- ggplot(
  tags,
  aes(xmin = start, xmax = end, y = haplotype, fill = gene, forward = forward)
) +
  geom_gene_arrow() +
  scale_fill_manual(values = gene_colours, na.value = "grey70") +
  labs(x = "Length (bp)", y = NULL, fill = NULL) +
  theme_genes() +
  theme(legend.position = "top")

pdf(out_file, width = 11, height = max(4, 0.18 * length(tip_order)))

layout(matrix(c(1, 2), nrow = 1), widths = c(1, 2.4))
par(mar = c(4, 1, 2, 0))
plot(tree, show.tip.label = FALSE, direction = "rightwards")

print(structure_plot, newpage = FALSE)

invisible(dev.off())
cat("wrote", out_file, "\n")
