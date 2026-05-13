#!/usr/bin/env Rscript
# Differential Expression Analysis using DESeq2 + tximport

suppressPackageStartupMessages({
    library(tximport)
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(tibble)
})

# --- Argument Handling ---
args            <- commandArgs(trailingOnly = TRUE)
salmon_dir      <- args[1]
samplesheet     <- args[2]
outdir          <- args[3]
reference_level <- args[4]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load metadata ---
message("Loading metadata and setting reference level to: ", reference_level)
col_data <- read_csv(samplesheet, show_col_types = FALSE) |>
    select(sample, condition) |>
    mutate(condition = as.factor(condition)) |>
    mutate(condition = relevel(condition, ref = reference_level))

# --- File Validation ---
# Check Salmon output to make sure specific 'quant.genes.sf' files exists
# before feeding tximport, since Nextflow passes only directories
message("Verifying Salmon input files...")
quant_files <- file.path(salmon_dir, col_data$sample, "quant.genes.sf")
names(quant_files) <- col_data$sample

missing_files <- quant_files[!file.exists(quant_files)]
if (length(missing_files) > 0) {
    stop(sprintf("CRITICAL ERROR: Missing Salmon files:\n%s",
                paste(missing_files, collapse = "\n")))
}

# --- DESeq2 Processing ---
message("Import data with tximport...")
txi <- tximport(quant_files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(
    txi,
    colData = col_data |> column_to_rownames("sample"),
    design = ~ condition
)

# Filtering: keep genes with total count >= 10 in at least 50% of samples
# Genes expressed only in one group (zero in others) are retained,
# while globally lowly expressed genes are removed to improve dispersion
keep <- rowSums(counts(dds) >= 10) >= (ncol(dds) / 2)
dds <- dds[keep, ]
message(sprintf("Retained %d genes after filtering.", sum(keep)))

dds <- DESeq(dds)

# --- Exporting Results ---
res <- results(dds)

results_df <- as.data.frame(res) |>
    rownames_to_column("gene_id") |>
    arrange(padj)

write_tsv(results_df, file.path(outdir, "results.tsv"))

# Use min of 1000 or total genes for VST nsub to handle small test datasets
vst_nsub <- min(1000, nrow(dds))
vst_vals <- vst(dds, blind = FALSE, nsub = vst_nsub)

# Export VST normalized counts for PCA and Heatmap
vst_df <- as.data.frame(assay(vst_vals)) |> tibble::rownames_to_column("gene_id")
write_tsv(vst_df, file.path(outdir, "vst_counts.tsv"))

# Save RDS for reporting
saveRDS(dds, file.path(outdir, "dds.rds"))
saveRDS(vst_vals, file.path(outdir, "vst.rds"))

message(sprintf("Found %d significant DEGs (padj < 0.05).",
                sum(results_df$padj < 0.05, na.rm = TRUE)
))

message("Differential expression analysis finished.")