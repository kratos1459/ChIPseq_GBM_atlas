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
  stitching, ranking by total H3K27ac signal. The ROSE inflection cutoff turned
  out to track data provenance (median 280 SEs per GSC line vs 753 per
  pediatric sample), so the **primary definition is the top 500 domains per
  sample**; inflection results are kept alongside, and top-250 / top-1000 are
  run as robustness checks.
- **Split-half ceilings use half-size groups**, so comparing them with a
  full-group cross-cohort statistic flatters the cross-cohort side. The headline
  comparison is therefore **size-matched** (k samples vs k samples against a
  shared reference set; `SIZEMATCH=1`).
- Average precision with its chance level replaces ROC AUC, which exceeded 0.9
  for every comparison at ~1 positive bin in 40.

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
22_loop_depth_matched.py   subsample loops to equal depth, 20 repeats
30_cross_cohort.py         recurrence, split-half ceilings, core domains
                           env: TOP_N (default 500), BOOT (bootstrap reps),
                           SIZEMATCH=1 (size-matched agreement), CC_OUT
40_motifs.py               SEA on constituent peaks, vs sporadic and shuffled
50_expression.py           recurrence vs independent expression
60_figures.py              600 dpi PNG + LZW TIFF figures
70_write_report.py         short Methods & Results report (docx)
80_manuscript.py           full manuscript + supplementary (docx), every
                           number read from results/ at build time
90_sync_github.sh          mirror code, tables, figures, manuscript for git
run_pipeline.sh            unattended orchestration of the stages above
```

Rebuild the manuscript after any change upstream:

```
TOP_N=500 BOOT=200 SIZEMATCH=1 python code/30_cross_cohort.py
python code/60_figures.py
python3 code/80_manuscript.py      # needs python-docx, pandas
```

Run `bash code/run_pipeline.sh` after `00_install_env.sh`; it waits for
downloads, calibrates, converts, and calls every sample. Stages are resumable.

## Results (fixed top-500 definition; manuscript in `manuscript/`)

**Cross-cohort reproducibility** (89 GSC lines, 41 tumours, 12 NSC, 21 pediatric)

- GSC-line recurrence predicts recurrent tumour bins with average precision
  **0.83** (chance 0.03). Against the half-size split-half ceiling (0.87) that
  looks like 96%, stable at top-250/500/1000.
- **Size-matched**, GSC lines recover **87% (95% range 80-94%)** of the
  agreement between two independent sets of 20 tumours, and fall short in
  50/50 draws; ~85% at k = 6 and 10. A real but modest cell-line deficit.
- Adult cohorts recover only **~60%** of within-pediatric agreement (k = 10).
  Agreement with NSC cannot be separated from NSC's own noise (k = 6).
- 985 bins are recurrent in both adult cohorts; 87 are rare in NSC, giving
  **30 core domains (0.44 Mb)** incl. SOX9, LFNG, ID2-AS1, CCND2-AS1, WBP2,
  C9orf50, CDC25B. Core size: 18 / 30 / 76 domains at top-250 / 500 / 1000.

**Enhancer-to-gene assignment** (GSE262089, 4 GSC lines + NSC)

- Nearest gene = contacted gene for **15-19%** of links; **22-31%** for
  single-TSS anchors; distance-matched decoy **8-10%**.
- Median contact distance 77-104 kb vs nearest TSS 14-20 kb.
- "GSC-specific" contacted genes: **1,397 as called -> 436 at matched depth
  -> 107 robust** (>=80% of 20 subsamples), incl. BRD4, JUN, FOSL2, CEBPB,
  FOXG1, GOLIM4, ADAM17, LAMC1, KIF5B.

**Motifs and expression**

- 52 GC-box motifs enriched vs sporadic peaks; **none** vs a
  dinucleotide-shuffled control (core 56.9% GC vs 48.9%): composition, not TF
  signal; 79 peaks, underpowered.
- Genes near core domains are more highly expressed in TCGA GBM (Cliff's
  delta 0.25, p = 0.008), pediatric HGG and LGG; genome-wide, NSC recurrence
  correlates with expression nearly as well as GSC recurrence (0.20 vs 0.25).

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
manuscript/ manuscript and supplementary information (docx)
results/   per-sample calls, agreement tables, gene assignments
figures/   600 dpi PNG and LZW TIFF
logs/      stage logs
```
