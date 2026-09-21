# RNAmeta 0.2.0

## Analysis

- Added a shared annotation bundle and bounded in-memory cache layer.
- Reworked Preliminary Analysis to reuse one consistent peak annotation result.
- Aligned Gene Matrix center-site annotation and transcript-to-gene mapping with
  RNAmeta-online.
- Added per-sample caching for Transcription, Translation, and Splice Site
  profiles.
- Improved Meta Plot mapping, sampling, reproducibility, and density caching.

## Visualization and export

- Rebuilt all maintained modules with a consistent interactive chart system.
- Added headless PNG and PDF export with configurable dimensions, aspect ratio,
  and PNG scale.
- Added independent Site Profile panels with single-panel or batch export.
- Restored chart titles, axes, labels, region guides, and white export
  backgrounds.
- Updated Gene Matrix and transcript-region layouts to match RNAmeta-online.

## Package scope

- The maintained package includes Preliminary Analysis, Gene Type, Gene Matrix,
  Peak summaries, Transcription, Translation, Splice Site, and Meta Plot.
- Motif, Genome Browser/IGV, and RNA Structure are not included because they
  require external engines or browser runtimes outside the package scope.
