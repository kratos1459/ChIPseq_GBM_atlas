# Glioblastoma enhancer atlas: public ChIP-seq inventory (started 2026-11-12)
Independent of GBM Visium, KXA and Moringa projects. Fully computational.

| GEO | Content (from GEO/papers) | Use |
|---|---|---|
| GSE119776 | H3K27ac ChIP-seq + RNA-seq, 43 glioma stem cell (GSC) lines and 9 NSCs (Mack et al.) | core GSC super-enhancer set |
| GSE145646 | H3K27ac ChIP-seq of GBM patient tumours (Xu et al.; 34 patients aggregated in a 2025 Cell Rep paper) | patient-tumour enhancer landscape |
| GSE262089 | H3K27ac HiChIP, H3K27ac ChIP-seq, ATAC-seq, RNA-seq, scRNA-seq, 4 patients + NSC | 3D enhancer-promoter links |
| GSE162976 | pHGG ChIP-seq (H3K27ac, H3K27me3, H3K27M), Hi-C, ATAC | H3-mutant vs wild-type contrast |

Status: GEO is reachable from this machine. No ChIP tools installed yet (macs3, bowtie2, samtools, bedtools missing); to be installed through conda before any processing.
Next: per-sample metadata table (GSM, mark, tissue, subtype, read length) from GEO SOFT files; decide raw-FASTQ vs processed peaks per series.
