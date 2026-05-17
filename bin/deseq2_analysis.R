#!/usr/bin/env Rscript
# Differential Expression Analysis using DESeq2 + tximport

suppressPackageStartupMessages({
    library(tximport)
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(tibble)
    library(apeglm)
})

# --- Argument Handling ---
args            <- commandArgs(trailingOnly = TRUE)
salmon_dir      <- args[1]
samplesheet     <- args[2]
reference_level <- args[3]
lfc_threshold   <- as.numeric(args[4])
padj_threshold  <- as.numeric(args[5])
outdir          <- args[6]

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

# Filtering: keep genes with at least 10 counts in at least 50% of samples
keep <- rowSums(counts(dds) >= 10) >= (ncol(dds) / 2)
dds <- dds[keep, ]
message(sprintf("Retained %d genes after filtering.", sum(keep)))

dds <- DESeq(dds)

# --- Exporting Results ---
message("Shrinking LFC estimates...")
message("Using coefficient: ", resultsNames(dds)[2])
res <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")

results_df <- as.data.frame(res) |>
    rownames_to_column("gene_id") |>
    arrange(padj)
write_tsv(results_df, file.path(outdir, "results.tsv"))

vst_nsub <- min(1000, nrow(dds))
vst_vals <- vst(dds, blind = FALSE, nsub = vst_nsub)

saveRDS(dds, file.path(outdir, "dds.rds"))
saveRDS(vst_vals, file.path(outdir, "vst.rds"))

message(sprintf("Found %d significant DEGs (padj < %.2f, |LFC| > %.1f).",
                sum(results_df$padj < padj_threshold & abs(results_df$log2FoldChange) > lfc_threshold, na.rm = TRUE),
                padj_threshold,
                lfc_threshold
))

message("Differential expression analysis finished.")
