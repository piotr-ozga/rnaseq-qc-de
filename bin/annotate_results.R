#!/usr/bin/env Rscript
# Annotate DESeq2 results using a GTF file
# Extracts: gene_name, gene_biotype, chromosome, start, end, strand

suppressPackageStartupMessages({
    library(rtracklayer)
    library(readr)
    library(dplyr)
})

# --- Argument Handling ---
args            <- commandArgs(trailingOnly = TRUE)
results_tsv     <- args[1]
gtf_file        <- args[2]
outdir          <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load DESeq2 results ---
message("Loading DESeq2 results: ", results_tsv)
results_df <- read_tsv(results_tsv, show_col_types = FALSE)

# --- Parse GTF ---
message("Parsing GTF: ", gtf_file)

# Skip dozens of unused columns to significantly reduce RAM usage
gtf <- rtracklayer::import(
    gtf_file,
    feature.type = "gene",
    colnames = c("gene_id", "gene_name", "gene_biotype", "gene_type")
)

# Convert GRanges to data frame, keep only columns that are reliably present
gtf_df <- as.data.frame(gtf) |>
    # 'gene_id' is mandatory in GTF spec, others may be absent
    select(
        gene_id,
        any_of(c("gene_name", "gene_biotype", "gene_type")),
        seqnames, start, end, strand
    ) |>
    # Normalise the biotype column name across GTF sources:
    # - Ensembl GTFs use 'gene_biotype'
    # - GENCODE GTFs use 'gene_type'
    # - Some builds export both as synonyms
    # - gene_type only -> rename to gene_biotype
    # - both present -> drop gene_type, keep gene_biotype
    # - gene_biotype only -> no changes
    (\(df) if ("gene_type" %in% names(df) && !"gene_biotype" %in% names(df))
                rename(df, gene_biotype = gene_type)
            else if ("gene_type" %in% names(df))
                select(df, -gene_type)
            else df
    )()|>
   rename(chromosome = seqnames)

if (!"gene_name" %in% colnames(gtf_df)) {
    gtf_df <- gtf_df |> mutate(gene_name = gene_id)
}

gtf_df <- gtf_df |>
    mutate(
        # Strip Ensembl version suffix if present (e.g. ENSG00000001.5 -> ENSG00000001)
        gene_id = sub("\\.\\d+$", "", gene_id),
        # If gene_name is NA, use gene_id so the output column is never empty
        gene_name = coalesce(gene_name, gene_id),
        # Flag primary chromosomes to exclude scaffolds and PAR-region duplicates in next step
        # that usually contains underscore or dot in the name
        is_primary_chrom = !grepl("[_\\.]", chromosome)
    ) |>
    arrange(gene_id, desc(is_primary_chrom)) |>
    distinct(gene_id, .keep_all = TRUE) |>
    select(-is_primary_chrom)

message(sprintf("GTF parsed: %d gene records loaded.", nrow(gtf_df)))

# --- Strip version suffix from results gene_id as well ---
results_df <- results_df |>
    mutate(gene_id = sub("\\.\\d+$", "", gene_id))

# --- Merge ---
n_before <- nrow(results_df)

annotated_df <- results_df |>
    left_join(gtf_df, by = "gene_id") |>
    relocate(
        gene_id,
        any_of(c("gene_name", "gene_biotype", "chromosome", "start", "end", "strand")),
        baseMean, log2FoldChange, lfcSE,
        any_of(c("stat")),
        pvalue, padj
    )

# Use chromosome presence as the true "was found in GTF" signal since annotation rate 
# computed on gene_name is never NA because of coalesce
n_annotated <- sum(!is.na(annotated_df$chromosome))
message(sprintf(
    "Annotation complete: %d / %d genes matched in GTF (%.1f%%).",
    n_annotated, n_before, 100 * n_annotated / n_before
))

# Warn if annotation rate is poor
if (n_annotated / n_before < 0.5) {
    warning(
        "Less than 50% of genes were annotated.",
        "Check that the GTF matches the genome build used for Salmon indexing."
    )
}

# --- Write outputs ---
out_tsv <- file.path(outdir, "results_annotated.tsv")
write_tsv(annotated_df, out_tsv)
message("Annotated results written to: ", out_tsv)

# --- Summary report ---
biotype_summary <- annotated_df |>
    filter(!is.na(gene_biotype)) |>
    count(gene_biotype, sort = TRUE)

summary_lines <- c(
    sprintf("total_genes\t%d",          n_before),
    sprintf("genes_annotated\t%d",      n_annotated),
    sprintf("annotation_rate\t%.4f",    n_annotated / n_before),
    "",
    "# gene_biotype counts",
    sprintf("%s\t%d", biotype_summary$gene_biotype, biotype_summary$n)
)

writeLines(summary_lines, file.path(outdir, "annotation_summary.txt"))
message("Annotation summary written.")

message("Annotation finished.")