# Reproducibility of the glioblastoma super-enhancer landscape

A fully computational re-analysis of four published H3K27ac ChIP-seq series,
asking two questions that single-cohort enhancer catalogues cannot answer:

1. **How much of the glioblastoma super-enhancer landscape reproduces between
   independent cohorts**, and how does cross-cohort agreement compare with the
   ceiling set by sampling noise within a single cohort?
2. **How often is the nearest gene the gene an enhancer actually contacts?**
   H3K27ac HiChIP loops make this measurable rather than assumed.

No new data are generated. Every input is a published, processed file from GEO
or the UCSC Xena hub.

## Data

| Series | Content used here | Build | Role |
|---|---|---|---|
| GSE119776 | H3K27ac signal (MACS2 pileup bedGraph), 89 glioma stem-cell lines + 9 NSC | hg19 | adult model-system cohort |
| GSE145646 | H3K27ac signal (bigWig), 41 patient tumours + 3 NSC | hg38 | adult patient cohort |
| GSE162976 | H3K27ac, H3K27M, H3K27me3 signal **and** published MACS2 peaks, 21 H3K27ac samples; RNA-seq TPM | hg38 | peak-caller calibration + pediatric external test |
| GSE262089 | H3K27ac HiChIP loops (FitHiChIP, 10 kb anchors), H3K27ac and ATAC peaks, RNA-seq, 4 GSC lines + NSC | hg38 | 3D enhancer-target assignment |
| TCGA GBM / LGG (Xena) | bulk RNA-seq, n=172 / n=530 | — | independent expression validation |

Total: 657 samples catalogued, 142 adult and 21 pediatric H3K27ac tracks used.

DepMap CRISPR dependency data is **not** included: depmap.org serves a
bot-verification page to non-browser clients, so that file cannot be fetched
programmatically.

## Design decisions worth knowing

- **Two cohorts published signal but no peaks.** Rather than pick a threshold by
  eye, the signal-to-peak caller is *calibrated* against GSE162976, the one
  series that published both a signal track and MACS2 peaks for the same GSMs
  (`11_calibrate_caller.py`). The selected parameters are then applied unchanged
  to every cohort, so cross-cohort differences cannot be an artefact of the
  caller.
- **Published peaks are not used for the comparison itself**, even where they
  exist, so that all cohorts are processed identically.
- **Native-build processing.** Peak calling and stitching happen in each
  cohort's own genome build; only the final enhancer coordinates are lifted to
  hg38. Whole signal tracks are never lifted.
- **Fixed 5 kb bins, not a union merge.** A union of stitched enhancer domains
  across 160+ samples collapses into megabase blocks, so recurrence is computed
  per bin.
- **Every agreement statistic is reported against a split-half noise ceiling**
  computed within each cohort, because between-cohort agreement is
  uninterpretable without knowing what two halves of one cohort achieve.
- Super-enhancers follow ROSE: promoter-proximal peaks dropped, 12.5 kb
  stitching, ranking by total H3K27ac signal, cutoff at the inflection point of
  the ranked curve (found numerically, not at a fixed rank).

## Pipeline

```
00_install_env.sh          conda environment (samtools, bedtools, MACS3,
                           deepTools, bedops, UCSC tools, MEME, pyBigWig)
01_geo_metadata.py         per-sample metadata from GEO SOFT files
02_filelist_sizes.py       download volume from GEO's own manifest
03_download.py             resumable download, calibration cohort first
03b_verify_downloads.py    byte-size verification and repair
04_annotation.sh           GENCODE v44 genes/TSS/promoters, blacklists, chains
04c_hg19_promoters.sh      hg19 promoter windows
05_build_manifest.py       one row per usable assay file
06_bedgraph_to_bigwig.sh   hg19 bedGraph cohort -> bigWig
07_download_motif_ref.sh   hg38 sequence + JASPAR motifs
08_download_validation.sh  TCGA GBM/LGG expression (Xena)
10_signal_peaks.py         signal -> enriched regions (calibratable)
11_calibrate_caller.py     tune the caller against published MACS2 peaks
12_stitch_superenhancers.py  ROSE-style stitching and SE ranking
13_call_all.py             call peaks + SEs for every sample, lift to hg38
20_loop_targets.py         3D contact vs nearest-TSS assignment
21_loop_controls.py        single-TSS subset, distance-matched null,
                           cross-line reproducibility
30_cross_cohort.py         recurrence, split-half ceilings, core domains
40_motifs.py               motif enrichment, reproducible vs sporadic
50_expression.py           recurrence vs independent expression
60_figures.py              600 dpi PNG + LZW TIFF figures
run_pipeline.sh            unattended orchestration of the stages above
```

Run `bash code/run_pipeline.sh` after `00_install_env.sh`; it waits for
downloads, calibrates, converts, and calls every sample. Stages are resumable.

## Results so far

Enhancer-to-gene assignment (GSE262089, 4 GSC lines + NSC, 33k-96k loops each):

- The nearest gene is the contacted gene for only **15-19%** of
  enhancer-promoter contacts; **22-31%** when restricted to promoter anchors
  containing a single TSS; a distance-matched decoy gene achieves **8-10%**.
- For only **40-53%** of linked enhancers is the nearest gene even among the
  genes that enhancer contacts.
- Median contact distance 77-104 kb against a median nearest-TSS distance of
  14-20 kb; the median link skips 2-3 intervening genes.
- Genes contacted by the most distinct enhancers in all four GSC lines: MYC
  (with CASC11), CDK6, JUN, JAG1, TRIB2, ZNF608, SPRY1, ZFP36L2, DUSP4, ATXN1.
- Contact-target agreement between GSC lines (mean Jaccard 0.52-0.59) is
  **not higher than GSC-versus-NSC agreement** (0.48-0.55), so at 10 kb anchor
  resolution the enhancer-target wiring looks largely cell-type invariant; what
  differs between tumour and normal is which enhancers are active, not which
  genes they reach.

Cross-cohort super-enhancer reproducibility, motif enrichment and expression
validation are computed by stages that run after the full download completes.

## Caveats

- HiChIP anchors are 10 kb bins, so enhancer and promoter positions are coarse
  and fine-scale rewiring cannot be detected.
- Loop depth differs several-fold between samples (33k-96k loops), so any
  "present in GSC, absent in NSC" statement needs depth matching before it can
  be read as biology.
- The two adult cohorts published no processed RNA-seq, so cohort-internal
  expression validation is possible only in the pediatric series and the four
  HiChIP GSC lines; TCGA provides the adult expression comparison instead.
- GSE119776 is hg19 and is lifted after calling; 48 of 38,899 TSSs do not map
  between builds.

## Layout

```
code/      analysis scripts (this repository)
data/      manifests and reference annotation; raw signal is not committed
results/   per-sample calls, agreement tables, gene assignments
figures/   600 dpi PNG and LZW TIFF
logs/      stage logs
```
