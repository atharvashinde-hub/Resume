# --------------------------- Load Libraries ---------------------------
library(DESeq2)
library(fgsea)
library(tidyverse)
library(msigdbr)
library(ggplot2)
library(gridExtra)
library(HGNChelper)  # if needed for gene symbol corrections

# --------------------------- Step 1: Load Data ---------------------------
df <- read.csv("/Users/atharva/Downloads/rna_seq_with_gene_names.txt", sep="\t", header=TRUE)
length(df)
rownames(df) <- make.unique(df$gene_name)  # ensure unique gene names
length(df)
# Extract counts matrix (columns 3-7 correspond to samples; adjust if needed)
counts <- as.matrix(df[, 3:7])
length(counts)
# --------------------------- Step 2: DESeq2 Analysis ---------------------------
sample_info <- data.frame(
  condition = factor(c(rep("KO", 3), rep("WT", 2))),
  row.names = colnames(counts)
)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
length(dds)
res <- results(dds, contrast = c("condition", "KO", "WT"))
length(res)
res$gene <- rownames(res)
res <- res[!is.na(res$padj), ]
length(res)
# --------------------------- Step 3: Prepare Ranked Genes ---------------------------
# Create a ranking metric: -log10(pvalue) * sign(log2FoldChange)
ranked_genes <- setNames(-log10(res$pvalue) * sign(res$log2FoldChange), res$gene)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Add a small jitter to break ties
set.seed(123)
ranked_genes <- ranked_genes + runif(length(ranked_genes), min = -1e-6, max = 1e-6)
length(ranked_genes)
# --------------------------- Step 4: Load Mouse Gene Sets ---------------------------

library(fgsea)

# Replace with the correct path to your downloaded file
gmt_file <- "/Users/atharva/Downloads/m3.all.v2024.1.Mm.symbols.gmt"
gene_sets <- gmtPathways(gmt_file)
length(gene_sets)
# Filter gene sets: keep only those with 10-5000 genes for focused analysis
gene_sets <- gene_sets[sapply(gene_sets, function(x) length(x) >= 10 & length(x) <= 500)]
length(gene_sets)
# --------------------------- Step 5: Run FGSEA ---------------------------
set.seed(123)
fgsea_res <- fgsea(
  pathways = gene_sets,
  stats = ranked_genes,
  minSize = 10,
  maxSize = 500,
  nperm = 10000
)
print(head(fgsea_res[order(fgsea_res$pval), ], 3))
length(fgsea_res)
# Report total number of pathways in FGSEA results
cat("Total number of pathways in FGSEA results:", nrow(fgsea_res), "\n")

# --------------------------- Step 6: Explore NES Distribution ---------------------------
cat("NES summary:\n")
print(summary(fgsea_res$NES))
ggplot(fgsea_res, aes(x = NES)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Normalized Enrichment Scores (NES)",
       x = "NES", y = "Frequency") +
  theme_minimal()

# --------------------------- Step 7: Select Top and Bottom Gene Sets ---------------------------
if(nrow(fgsea_res) >= 10){
  # Sort FGSEA results by NES in descending order
  fgsea_sorted <- fgsea_res[order(fgsea_res$NES, decreasing = TRUE), ]
  top20 <- head(fgsea_sorted, 15)
  bottom20 <- head(fgsea_res[order(fgsea_res$NES, decreasing = FALSE), ], 15)
  top10 <- head(fgsea_sorted, 10)
  bottom10 <- head(fgsea_res[order(fgsea_res$NES, decreasing = FALSE), ], 10)
  selected_gene_sets <- rbind(top20, bottom20)
} else {
  cat("FGSEA results contain fewer than 10 pathways; using all available pathways.\n")
  selected_gene_sets <- fgsea_res
}
print(selected_gene_sets[, c("pathway", "NES", "pval", "padj")])

# --------------------------- Step 8: Create a Summary (Lollipop) Plot ---------------------------
# Filter significant gene sets for top and bottom data (padj < 0.25)
top20_signif <- subset(top20, padj < 0.25)
bottom20_signif <- subset(bottom20, padj < 0.25)

# Plot for Top Enriched Gene Sets
nes_plot_top <- ggplot(top20_signif, aes(x = NES, y = reorder(pathway, NES), 
                                         color = NES > 0, size = -log10(padj))) +
  geom_segment(aes(xend = 0, yend = pathway), linewidth = 0.8) +
  geom_point() +
  scale_color_manual(values = c( "blue"), labels = c("Up")) +
  labs(
    title = "Top 20 Enriched Gene Sets (padj < 0.25)",
    x = "Normalized Enrichment Score (NES)",
    y = "Gene Set Pathway",
    color = "Direction",
    size = "-log10(padj)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
print(nes_plot_top)

# Plot for Bottom Enriched Gene Sets
nes_plot_bottom <- ggplot(bottom20_signif, aes(x = NES, y = reorder(pathway, NES), 
                                               color = NES > 0, size = -log10(padj))) +
  geom_segment(aes(xend = 0, yend = pathway), linewidth = 0.8) +
  geom_point() +
  scale_color_manual(values = c("red", "blue"), labels = c("Down", "Up")) +
  labs(
    title = "Bottom 20 Enriched Gene Sets (padj < 0.25)",
    x = "Normalized Enrichment Score (NES)",
    y = "Gene Set Pathway",
    color = "Direction",
    size = "-log10(padj)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))
print(nes_plot_bottom)

# --------------------------- Step 9: Define Function for Individual Enrichment Plots with Leading Edge Labels ---------------------------
plot_with_leading_edge <- function(pathway_name) {
  gene_set <- gene_sets[[pathway_name]]
  overlap_genes <- intersect(gene_set, names(ranked_genes))
  if (length(overlap_genes) == 0) {
    message("Pathway ", pathway_name, " has no overlap with the ranked gene list. Skipping annotation.")
    return(plotEnrichment(gene_set, ranked_genes) +
             labs(title = pathway_name, x = "Rank in Ordered Dataset", y = "Enrichment Score") +
             theme_minimal() +
             theme(plot.title = element_text(size = 8)))
  }
  p <- plotEnrichment(gene_set, ranked_genes) +
    labs(title = pathway_name, x = "Rank in Ordered Dataset", y = "Enrichment Score") +
    theme_minimal() +
    theme(plot.title = element_text(size = 8))
  running_sum <- tryCatch({
    fgsea:::calcGseaStat(ranked_genes, gene_set)
  }, error = function(e) {
    message("Error calculating running sum for ", pathway_name, ": ", e$message)
    return(NULL)
  })
  if (is.null(running_sum)) {
    return(p)
  }
  pathway_row <- fgsea_res[fgsea_res$pathway == pathway_name, ]
  if (nrow(pathway_row) > 0) {
    le_genes <- pathway_row$leadingEdge[[1]]
    le_positions <- which(names(ranked_genes) %in% le_genes)
    if (length(le_positions) > 0) {
      le_df <- data.frame(
        rank = le_positions,
        gene = names(ranked_genes)[le_positions],
        ES = running_sum[le_positions]
      )
      p <- p + geom_text(data = le_df, 
                         aes(x = rank, y = ES, label = gene),
                         size = 2, angle = 45, hjust = 0, vjust = 0,
                         check_overlap = TRUE)
    }
  }
  return(p)
}

# --------------------------- Step 10: Generate and Arrange Individual Enrichment Plots ---------------------------
selected_pathways <- top20$pathway
plot_list <- lapply(selected_pathways, plot_with_leading_edge)
plot_list <- plot_list[!sapply(plot_list, is.null)]
grid.arrange(grobs = plot_list, ncol = 3)

selected_pathways <- bottom20$pathway
plot_list <- lapply(selected_pathways, plot_with_leading_edge)
plot_list <- plot_list[!sapply(plot_list, is.null)]
grid.arrange(grobs = plot_list, ncol = 3)

# --------------------------- Step 11: Print Leading Edge Genes for Each Selected Gene Set ---------------------------
for (pathway_name in top20$pathway) {
  pathway_row <- fgsea_res[fgsea_res$pathway == pathway_name, ]
  if (nrow(pathway_row) > 0) {
    le_genes <- pathway_row$leadingEdge[[1]]
    if (!is.null(le_genes) && length(le_genes) > 0) {
      cat("\nPathway:", pathway_name, "\nLeading Edge Genes:", paste(le_genes, collapse = ", "), "\n")
    } else {
      cat("\nPathway:", pathway_name, "\nNo leading edge genes found.\n")
    }
  } else {
    cat("\nPathway:", pathway_name, "\nNot found in FGSEA results.\n")
  }
}
