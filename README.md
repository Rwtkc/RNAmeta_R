# RNAmeta

RNAmeta is an R package for annotating genomic intervals, summarizing peak
features, profiling transcript landmarks, and generating interactive figures
with headless PNG/PDF export.

- **Online server:** [RNAmeta-online](https://rnainformatics.cn/RNAmeta-online/)
- **Source code:** [GitHub](https://github.com/Rwtkc/RNAmeta_R)
- **Contact:** [rwtkc@foxmail.com](mailto:rwtkc@foxmail.com)

## Features

- Shared peak annotation and in-memory caching across analysis modules.
- Preliminary Analysis, Gene Type, and Gene Matrix summaries.
- Peak distribution, exon type, exon number, exon size, and gene size charts.
- Transcription, Translation, and Splice Site profiles with independent panels.
- Meta Plot analysis for transcript, mRNA, and ncRNA regions.
- Reproducible sampling through explicit random seeds.
- PNG and PDF export on Windows, macOS, and headless Linux systems.

Motif discovery, Genome Browser/IGV, and RNA Structure are intentionally outside
the scope of this R package. These features remain available only where the
required external services and runtime environments are provided separately.

## Installation

### From GitHub

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

options(repos = BiocManager::repositories())
remotes::install_github("Rwtkc/RNAmeta_R")
```

### From a source package

```r
install.packages("RNAmeta_0.2.0.tar.gz", repos = NULL, type = "source")
```

RNAmeta uses Bioconductor annotation infrastructure. On Linux, installation of
the `rsvg` dependency may additionally require the system package
`librsvg2-dev` (Debian/Ubuntu) or `librsvg2-devel` (Fedora/RHEL).

## Annotation inputs

Species annotation files are supplied by path and loaded only when needed. The
main workflow accepts:

- `txdb_file`: a TxDb object or TxDb SQLite file;
- `tx_length_file`: a transcript-length `.rda` file;
- `gff_file`: a processed GFF `GRanges` `.rda` file;
- `bed_file`: one or more BED paths, or supported `GRanges` objects where
  documented.

Species-specific annotation packages do not need to be installed.

## Quick start

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

The result is a named list shared by downstream modules. Named access is
recommended; legacy positional access remains supported.

```r
names(preliminary)
# cct_list, peak_gr_list, gff, peak_gr_for_motifscan, motif_parameter,
# flanking_size, txdb_features, gene_max_tx_lengths, short_flanking_size,
# tx_lengths
```

## Gene Type and Gene Matrix

```r
gene_type <- RNAmeta_gene_type_data(
  bed_file = beds,
  txdb_file = "species.txdb.sqlite",
  gff_file = "species.gff.rda",
  set_group_name = labels
)

RNAmeta_gene_type_plot(gene_type)
RNAmeta_gene_type_export(gene_type, "gene_type.png", width = 1600,
                         aspect_ratio = c(16, 9))

gene_matrix <- RNAmeta_gene_matrix_data(preliminary, labels)
RNAmeta_get_gene_matrix_plot(gene_matrix)
RNAmeta_gene_matrix_export(gene_matrix, "gene_matrix.pdf")
```

Gene Matrix supports two to five samples and uses the same center-site and
transcript-to-gene mapping workflow as RNAmeta-online.

## Peak summary charts

```r
RNAmeta_peak_distribution_plot(preliminary, labels)
RNAmeta_peak_exon_type_plot(preliminary, labels)
RNAmeta_peak_exon_number_plot(preliminary, labels)
RNAmeta_peak_exon_size_plot(preliminary, labels)
RNAmeta_peak_gene_size_plot(preliminary, labels)

RNAmeta_peak_distribution_export(preliminary, "peak_distribution.png", labels)
RNAmeta_peak_exon_type_export(preliminary, "peak_exon_type.pdf", labels)
RNAmeta_peak_exon_number_export(preliminary, "peak_exon_number.png", labels)
RNAmeta_peak_exon_size_export(preliminary, "peak_exon_size.pdf", labels)
RNAmeta_peak_gene_size_export(preliminary, "peak_gene_size.png", labels)
```

## Transcription, Translation, and Splice Site

Site-profile functions return a named list of independent figures by default.
Use `panel` to select one figure. When `panel = "all"`, export functions append
the panel name to each output filename automatically.

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

# Multiple transcript types return a named RNAmeta_meta_plots list.
meta$mRNA
meta$ncRNA

RNAmeta_meta_plot_export(meta, "meta_mrna.png", tx_type = "mRNA")
RNAmeta_meta_plot_export(meta, "meta_ncrna.pdf", tx_type = "ncRNA")
```

The `seed` controls random site sampling and confidence-interval resampling
without changing the caller's global random state.

## Export dimensions

All PNG/PDF exporters support explicit dimensions or an aspect ratio. PNG
export also supports `scale`. Exported figures use a white background.

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

Static export is generated through R and `rsvg`; it does not require RStudio
Viewer, X11, or an active display.

## Cache management

RNAmeta caches annotation bundles, per-sample annotation results, site profiles,
and Meta Plot calculations in memory.

```r
RNAmeta_clear_cache("annotation")
RNAmeta_clear_cache("site_profile")
RNAmeta_clear_cache("meta")
RNAmeta_clear_cache("all")

# Backward-compatible annotation-only function
RNAmeta_clear_annotation_cache()
```

Each cache keeps at most 128 entries by default and uses least-recently-used
eviction. The limit can be changed for the current R session:

```r
options(RNAmeta.cache.max_entries = 64L)
```

## License

RNAmeta is distributed under the GNU General Public License version 3 or later.
