# RNAmeta

```text
 ____  _   _    _    __  __      _
|  _ \| \ | |  / \  |  \/  | ___| |_ __ _
| |_) |  \| | / _ \ | |\/| |/ _ \ __/ _` |
|  _ <| |\  |/ ___ \| |  | |  __/ || (_| |
|_| \_\_| \_/_/   \_\_|  |_|\___|\__\__,_|
```

RNAmeta is an R package for genomic interval annotation, transcript-region
profiling, statistical summaries, and interactive visualization. Figures can be
viewed interactively or exported as PNG and PDF files.

- [RNAmeta-online](https://rnainformatics.cn/RNAmeta-online/)
- [GitHub repository](https://github.com/Rwtkc/RNAmeta_R)
- [rwtkc@foxmail.com](mailto:rwtkc@foxmail.com)

## Main functions

| Module | Purpose |
| --- | --- |
| Preliminary Analysis | Prepare shared peak annotations for downstream analyses |
| Gene Type | Summarize peak-associated gene types |
| Gene Matrix | Compare shared and sample-specific genes across 2–5 samples |
| Peak Summary | Plot peak distribution, exon type, exon number, exon size, and gene size |
| Transcription | Profile signals around transcription start and end sites |
| Translation | Profile signals around translation start and stop sites |
| Splice Site | Profile signals around donor and acceptor sites |
| Meta Plot | Display peak distributions across transcript regions |

## Installation

Install the current version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("Rwtkc/RNAmeta_R")
```

Install a downloaded source package:

```r
install.packages("RNAmeta_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Input files

RNAmeta reads BED files and species annotation files supplied by the user:

- `bed_file`: one or more BED files;
- `txdb_file`: a TxDb object or TxDb SQLite file;
- `tx_length_file`: a transcript-length `.rda` file;
- `gff_file`: a processed GFF `GRanges` `.rda` file.

## Basic workflow

```r
library(RNAmeta)

beds <- c("sample_a.bed", "sample_b.bed")
labels <- c("Control", "Treatment")

preliminary <- RNAmeta_get_preliminary_analysis_data(
  bed_file = beds,
  txdb_file = "species.txdb.sqlite",
  tx_length_file = "species.txlens.rda",
  gff_file = "species.gff.rda"
)
```

The `preliminary` object can be reused by Gene Matrix, Peak Summary,
Transcription, Translation, and Splice Site modules.

## Gene Type

```r
gene_type <- RNAmeta_gene_type_data(
  bed_file = beds,
  txdb_file = "species.txdb.sqlite",
  gff_file = "species.gff.rda",
  set_group_name = labels
)

RNAmeta_gene_type_plot(gene_type)
RNAmeta_gene_type_export(gene_type, "gene_type.png")
```

## Gene Matrix

```r
gene_matrix <- RNAmeta_gene_matrix_data(preliminary, labels)
RNAmeta_get_gene_matrix_plot(gene_matrix)
RNAmeta_gene_matrix_export(gene_matrix, "gene_matrix.pdf")
```

Gene Matrix accepts two to five samples.

## Peak Summary

```r
RNAmeta_peak_distribution_plot(preliminary, labels)
RNAmeta_peak_exon_type_plot(preliminary, labels)
RNAmeta_peak_exon_number_plot(preliminary, labels)
RNAmeta_peak_exon_size_plot(preliminary, labels)
RNAmeta_peak_gene_size_plot(preliminary, labels)
```

Export the corresponding plots with:

```r
RNAmeta_peak_distribution_export(preliminary, "peak_distribution.png", labels)
RNAmeta_peak_exon_type_export(preliminary, "peak_exon_type.pdf", labels)
RNAmeta_peak_exon_number_export(preliminary, "peak_exon_number.png", labels)
RNAmeta_peak_exon_size_export(preliminary, "peak_exon_size.pdf", labels)
RNAmeta_peak_gene_size_export(preliminary, "peak_gene_size.png", labels)
```

## Transcription, Translation, and Splice Site

Each function returns a named list of individual plots. Use `panel` to display or
export one plot, or keep `panel = "all"` to export all plots.

```r
transcription <- get_transcription_data(preliminary, seed = 1)
transcription_plots <- RNAmeta_get_transcription_plot(transcription, labels)
transcription_plots[["Transcription_start_site"]]

RNAmeta_transcription_export(
  transcription,
  "transcription.png",
  labels,
  panel = "Transcription_start_site"
)

translation <- get_translation_data(preliminary)
translation_plots <- RNAmeta_get_translation_plot(translation, labels)
RNAmeta_translation_export(translation, "translation.pdf", labels)

splice_plots <- RNAmeta_get_splicesite_plot(preliminary, labels)
RNAmeta_splicesite_export(preliminary, "splice_site.png", labels)
```

## Meta Plot

```r
meta <- RNAmeta_meta_plot(
  bed_file = beds,
  txdb_file = "species.txdb.sqlite",
  plot_tx_type = c("mRNA", "ncRNA"),
  set_group_name = labels,
  sample_model = "Equidistance",
  enable_confidence_interval = TRUE,
  seed = 1
)

meta$mRNA
meta$ncRNA

RNAmeta_meta_plot_export(meta, "meta_mrna.png", tx_type = "mRNA")
RNAmeta_meta_plot_export(meta, "meta_ncrna.pdf", tx_type = "ncRNA")
```

## Image export

PNG and PDF output dimensions can be controlled with `width`, `height`, and
`aspect_ratio`. PNG output also supports `scale`.

```r
RNAmeta_peak_distribution_export(
  preliminary,
  "peak_16x9.png",
  set_group_name = labels,
  width = 1600,
  aspect_ratio = c(16, 9),
  scale = 2
)

RNAmeta_peak_distribution_export(
  preliminary,
  "peak_exact.pdf",
  set_group_name = labels,
  width = 1400,
  height = 700
)
```

## Cache

Repeated analyses reuse cached annotation and profile results. Clear a selected
cache or all cached results when input files change:

```r
RNAmeta_clear_cache("annotation")
RNAmeta_clear_cache("site_profile")
RNAmeta_clear_cache("meta")
RNAmeta_clear_cache("all")
```

## License

RNAmeta is distributed under the GNU General Public License version 3 or later.
