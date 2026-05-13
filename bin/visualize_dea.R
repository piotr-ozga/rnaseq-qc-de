#!/usr/bin/env Rscript
# Visualize DEA results: Volcano Plot, PCA, and Heatmap

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(tibble)
    library(DESeq2)
})

# --- Reproducibility ---
set.seed(42)

# --- Argument Handling ---
args            <- commandArgs(trailingOnly = TRUE)
annotated_tsv   <- args[1]
vst_rds         <- args[2]
samplesheet     <- args[3]
ref_level       <- args[4]
outdir          <- args[5]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Data Loading & Cleaning ---
message("Loading data...")
res      <- read_tsv(annotated_tsv, show_col_types = FALSE)
vst_obj  <- read_rds(vst_rds)
meta     <- read_csv(samplesheet, show_col_types = FALSE)

# Extract matrix from DESeqTransform object
vst_matrix <- assay(vst_obj)

# Set colors for control and treated sample
colors_tc <- setNames(
    c("#4DAF4A", "#E41A1C"), 
    c(ref_level, unique(meta$condition)[unique(meta$condition) != ref_level])
)


# Remove version suffixes from Gene IDs for consistency
res <- res |> mutate(gene_id_clean = gsub("\\..*$", "", gene_id))
rownames(vst_matrix) <- gsub("\\..*$", "", rownames(vst_matrix))

# --- Volcano Plot ---
message("Generating Volcano Plot...")

# Handle padj = 0 by replacing with very small value for -log10 calculation
res_volc <- res |>
    mutate(log_p = -log10(ifelse(padj == 0 | is.na(padj), 1e-300, padj)))

# Determine Y-axis ceiling based on data density (99.5th percentile)
upper_limit <- quantile(res_volc$log_p[res_volc$log_p < 300], 0.995, na.rm = TRUE)
y_limit <- upper_limit * 1.1

res_volc <- res_volc |>
    mutate(
        log_p_capped = ifelse(log_p > y_limit, y_limit, log_p),
        # Jitter points on the ceiling to prevent overlap
        log_p_plot = ifelse(log_p_capped >= y_limit, y_limit + runif(n(), -0.5, 0.5), log_p_capped),
        significance = case_when(
            padj < 0.05 & log2FoldChange > 1  ~ "Up",
            padj < 0.05 & log2FoldChange < -1 ~ "Down",
            TRUE                              ~ "Not Significant"
        ),
        label = coalesce(gene_name, gene_id_clean)
    )

# Select top 15 most significant genes for annotation
top_genes <- res_volc |>
    filter(significance != "Not Significant") |>
    arrange(padj) |>
    slice_head(n = 15)

# Create dynamic title based on experimental groups
conditions <- unique(meta$condition)
comp_title <- paste(conditions[conditions != ref_level], "vs", ref_level)

p_volcano <- ggplot(res_volc, aes(x = log2FoldChange, y = log_p_plot, color = significance)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "Not Significant" = "#D1D1D1")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30") +
    geom_label_repel(
        data = top_genes, aes(label = label),
        size = 3, fontface = "bold", fill = alpha("white", 0.9),
        box.padding = 0.2, point.padding = 0.3, force = 10,
        segment.size = 0.3, min.segment.length = 0, max.overlaps = Inf,
        show.legend = FALSE
    ) + theme_classic() +
    labs(
        title = comp_title,
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"Adjusted P-value")
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, face="bold")
    )

ggsave(file.path(outdir, "volcano_plot.pdf"), plot = p_volcano, width = 10, height = 8)

# --- PCA Plot ---
message("Generating PCA plot...")

pca <- plotPCA(vst_obj, intgroup = "condition") +
    geom_point(size = 5) +
    scale_color_manual(values = colors_tc) +
    geom_text_repel(aes(label = name), size = 3, fontface = "bold", show.legend = FALSE) +
    theme_bw() +
    labs(title = "Sample Clustering (PCA)", color = "condition") +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
    )

ggsave(file.path(outdir, "pca_plot.pdf"), plot = pca, width = 8, height = 7)

# --- Annotated Heatmap ---
message("Generating Heatmap...")

# Extract IDs for top 40 genes
top_heatmap_ids <- res |>
    filter(!is.na(padj)) |>
    arrange(padj) |>
    slice_head(n = 40) |>
    pull(gene_id_clean)

heatmap_data <- vst_matrix[rownames(vst_matrix) %in% top_heatmap_ids, , drop = FALSE]

if (nrow(heatmap_data) > 5) {
    id_to_name <- res |>
    filter(gene_id_clean %in% rownames(heatmap_data)) |>
    distinct(gene_id_clean, .keep_all = TRUE)

    # Sort and rename rows
    heatmap_data <- heatmap_data[match(id_to_name$gene_id_clean, rownames(heatmap_data)), ]
    rownames(heatmap_data) <- coalesce(id_to_name$gene_name, id_to_name$gene_id_clean)

    # Metadata for legends
    ann_color <- data.frame(Group = meta$condition)
    rownames(ann_color) <- meta$sample

    heatmap <- pheatmap(
        heatmap_data,
        annotation_col            = ann_color,
        # annotation_colors must be a list
        annotation_colors         = list(Group = colors_tc),
        scale                     = "row", # Row-wise Z-score scaling
        clustering_distance_rows  = "correlation",
        main                      = "Top 40 Differentially Expressed Genes",
        color                     = colorRampPalette(c("#313695", "white", "#A50026"))(100),
        fontsize_row              = 8,
        fontsize_col              = 10,
        cellheight                = 15,
        annotation_legend         = TRUE,
        legend                    = TRUE,
        border_color              = NA,
        filename                  = file.path(outdir, "heatmap.pdf")
    )
}

message("Visualization finished.")