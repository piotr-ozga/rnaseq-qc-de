# rnaseq-qc-de

A Nextflow pipeline for processing RNA-seq data - from raw reads or SRA accessions to differential expression results, annotation, and plots.

**Disclaimer:** Please note that this pipeline is a work in progress. Some functionalities may not work correctly or are still under active development.

---

## Motivation

This project was born out of a desire to move beyond simple tutorials and build a modular Nextflow pipeline from scratch. My primary goals were to solidify DSL2
skills and experience the full development lifecycle from raw data acquisition to biological insights. By architecting this tool, I gained hands-on experience in managing complex data flows and addressing the practical challenges of building bioinformatics pipelines.

---

## Pipeline Steps

The pipeline is organized into five subworkflows that run sequentially:

### `INPUT_HANDLER` - Input Handling & Validation

Reads the samplesheet CSV and routes each sample to the appropriate source. Samples with local FASTQ paths are loaded directly; samples with an SRR accession are downloaded from NCBI SRA using `sra-tools` (`prefetch` + `fasterq-dump` + `pigz`). The subworkflow validates every row for required fields (`sample`, `condition` and either `srr` or `read1`) and emits a unified reads channel regardless of origin. An optional `NCBI_API_KEY` Nextflow secret can be set to avoid rate-limiting.

### `QC_AND_TRIMMING` - Quality Control & Adapter Trimming

Runs `FastQC` on raw reads, then trims adapters and low-quality bases with `Trim Galore`, and runs `FastQC` again on the trimmed output. Both raw and post-trim QC results are forwarded to MultiQC via the `multiqc_files` channel topic. The module handles both paired-end and single-end data automatically. The `--fastqc` flag is explicitly filtered from `ext.args` to prevent duplicate QC runs.

### `SALMON_ALIGN` - Indexing & Quantification

Builds a decoy-aware `Salmon` index from the reference genome and transcriptome FASTAs (decoys are extracted automatically from chromosome headers). Each sample is then quantified with `salmon quant` using `--gcBias` and `--seqBias` correction flags and the `--geneMap` option to produce both transcript-level (`quant.sf`) and gene-level (`quant.genes.sf`) output. The index is built once and broadcast to all samples via `.collect()`.

### `RNASEQ_DEA` - Differential Expression Analysis, Annotation & Visualization

A dedicated subworkflow that chains three processes:

- **`DESEQ2_ANALYSIS`**: Imports gene-level counts from all Salmon output directories using `tximport`, constructs a `DESeqDataSet`, applies a count filter (≥ 10 counts in ≥50% of samples), and runs `DESeq2`. LFC estimates are shrunk using `lfcShrink()` with `type = "normal"` to reduce noise from low-count genes. Exports raw results (`results.tsv`), and serialized R objects (`dds.rds`, `vst.rds`) for downstream use.

- **`ANNOTATE_RESULTS`**: Parses the provided GTF file with `rtracklayer` (loading only `gene`-type features and minimal set of columns to reduce memory usage). Strips Ensembl version suffixes from gene IDs, normalises the biotype column across GTF sources (Ensembl `gene_biotype` vs. GENCODE `gene_type`), and performs a left join onto the DESeq2 results table. Outputs an annotated TSV and a plain-text summary with annotation rate and biotype counts. Warns if fewer than 50% of genes are matched, which typically indicates a GTF/genome build mismatch.

- **`VISUALIZE_DEA`**: Produces three PDF figures using `ggplot2`, `ggrepel`, `DESeq2` and `pheatmap`:
  - **Volcano plot**: colour-coded by significance (padj < `--padj_threshold`, |log₂FC| > `--lfc_threshold`), with the top `--volcano_labels` significant genes labelled; Y-axis is capped at the 99.5th percentile to prevent extreme p-values from compressing the plot.
  - **PCA plot**: sample-level clustering using VST counts from `DESeq2::plotPCA`.
  - **Heatmap**: row-wise Z-score scaled expression of the top `--heatmap_genes` most significant genes, with hierarchical clustering by correlation distance and sample group annotations.

### `REPORTING` - MultiQC Aggregation

Collects all files pushed to the `multiqc_files` channel topic (FastQC zips, Trim Galore logs, Salmon quant directories) and assembles a software versions table from the `versions` channel topic. Both are passed to `MultiQC` to produce a single HTML report summarising QC metrics across all samples.

## Prerequisites

To run this pipeline, you need:

* **Nextflow** (version >= 25.04.0)
* A container engine (**Docker** or **Singularity**) or **Conda/Mamba**.
* NCBI API Key (Optional): To avoid rate-limiting during SRA downloads, it is recommended to set your API key as a Nextflow secret:
`nextflow secrets set NCBI_API_KEY 'your_api_key_here'`

## Usage

### 1. Prepare the Samplesheet

Create a `.csv` file containing your sample metadata. The pipeline accepts local FASTQ paths or SRA accession numbers.

**Example `samplesheet.csv`:**

```
sample,condition,srr,read1,read2
WT_REP1,control,,data/WT_1.fastq.gz,data/WT_2.fastq.gz
KO_REP1,treatment,SRR1234567,,
```

*Note: For single-end reads, leave the `read2` column blank. For SRA accessions, leave `read1` and `read2` blank.*

### 2. Run the Pipeline

```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --genome_fasta reference/genome.fa \
    --transcriptome_fasta reference/transcriptome.fa \
    --gtf reference/genes.gtf \
    --ref_level control \
    -profile docker
```

**Required Parameters:**
* `--samplesheet`: Path to the CSV samplesheet.
* `--genome_fasta`: Path to the reference genome FASTA.
* `--transcriptome_fasta`: Path to the reference transcriptome FASTA.
* `--gtf`: Path to the annotation GTF.
* `--ref_level`: The baseline condition for DESeq2 (e.g. `control`).

**Optional Parameters**
* `--lfc_threshold`: Log2 fold change threshold (default: 1.0) .
* `--padj_threshold`: Adjusted p-value threshold (default: 0.05).
* `--heatmap_genes`: Number of top genes shown on heatmap (default: 40).
* `--volcano_labels`: Number of gene labels on volcano plot (default: 15)

**Profiles:**
Use `-profile docker`, `-profile singularity`, or `-profile conda`.

## Output Structure

Results are written to `--outdir` (default: `./results`):

```
results/
├── raw_fastq/              # FASTQ files downloaded from SRA
├── fastqc/
│   ├── raw/                # FastQC reports for raw reads
│   └── trimmed/            # FastQC reports for trimmed reads
├── trimgalore/
│   ├── logs/               # Trimming reports
│   └── *.fq.gz             # Trimmed reads (if save_trimmed = true)
├── salmon/
│   ├── index/              # Salmon index
│   └── <sample>/           # Per-sample quant.sf and quant.genes.sf
├── dea/
│   ├── deseq2/             # results.tsv, dds.rds, vst.rds
│   ├── annotation/         # results_annotated.tsv, annotation_summary.txt
│   └── visualization/      # volcano_plot.pdf, pca_plot.pdf, heatmap.pdf
├── reports/
│   └── multiqc/            # multiqc_report.html
└── pipeline_info/          # Execution timeline, report, DAG, trace
```

## CI/CD

Tests are run automatically on every push and pull request to `main` and `dev` using **GitHub Actions**. The workflow detects which modules were changed and runs only the relevant `nf-test` unit tests, plus a full integration test on every run.

## To Do

* **Downstream Analysis Tests**: Add `nf-test` coverage for R-based modules.
* **Functional Analysis**: Implement additional processes such as Gene Ontology (GO) or Pathway Enrichment (GSEA).

---

## AI Usage Statement

The core data flow, architecture, and logic of this pipeline were designed and implemented manually by the author. AI tools (LLMs) were used as assistants to improve development efficiency, especially for:

* Debugging errors and resolving Nextflow/Groovy syntax issues.
* Organizing and formatting Git commit messages (using AI to structure history based on provided skeletons).
* Suggesting code improvements, R plotting adjustments, and general syntax extensions.
* Assistance in creating README.md.