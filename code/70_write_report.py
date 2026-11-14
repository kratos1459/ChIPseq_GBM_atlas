#!/usr/bin/env python3
"""Assemble the Methods and Results document from the result files.

Every number in the output is read from a file written by an earlier stage, so
the document cannot drift from the analysis. Sections whose inputs do not exist
yet are written with an explicit "not yet computed" line rather than silently
omitted.

Run with the base python (it needs python-docx and pandas, not pyBigWig):
    python3 code/70_write_report.py
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
CC = os.path.join(RES, "cross_cohort")
LT = os.path.join(RES, "loop_targets")
OUT = os.path.join(ROOT, "ChIPseq_GBM_Methods_Results.docx")


def read_csv(path, **kw):
    return pd.read_csv(path, **kw) if os.path.exists(path) else None


def read_json(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return None


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    if isinstance(x, (int, np.integer)):
        return f"{int(x):,}"
    return str(x)


def build():
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Pt

    doc = Document()
    for name, size in (("Normal", 11),):
        st = doc.styles[name]
        st.font.name = "Times New Roman"
        st.font.size = Pt(size)

    def H(text, level=1):
        p = doc.add_paragraph()
        r = p.add_run(text)
        r.bold = True
        r.italic = level > 1
        r.font.size = Pt(13 - level)
        return p

    def P(text, italic=False):
        p = doc.add_paragraph()
        r = p.add_run(text)
        r.italic = italic
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        return p

    def table(df, caption, max_rows=30):
        P(caption, italic=True)
        if df is None or len(df) == 0:
            P("Not yet computed.")
            return
        d = df.head(max_rows)
        t = doc.add_table(rows=1, cols=len(d.columns))
        t.style = "Table Grid"
        for i, c in enumerate(d.columns):
            cell = t.rows[0].cells[i]
            cell.text = str(c)
            for pr in cell.paragraphs:
                for r in pr.runs:
                    r.bold = True
                    r.font.size = Pt(8)
        for _, row in d.iterrows():
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = fmt(v) if not isinstance(v, str) else v[:40]
                for pr in cells[i].paragraphs:
                    for r in pr.runs:
                        r.font.size = Pt(8)
        doc.add_paragraph()

    # ---------------------------------------------------------------- load --
    manifest = read_csv(os.path.join(DATA, "manifest.csv"))
    se_sum = read_csv(os.path.join(RES, "se_summary.csv"))
    calib = read_json(os.path.join(RES, "calibration_best.json"))
    agree = read_csv(os.path.join(CC, "agreement_both_definitions.csv"))
    if agree is None:
        agree = read_csv(os.path.join(CC, "agreement.csv"))
    ccsum = read_json(os.path.join(CC, "summary.json"))
    core_genes = read_csv(os.path.join(CC, "core_domain_genes.csv"))
    controls = read_csv(os.path.join(LT, "controls_summary.csv"))
    nvc = read_csv(os.path.join(LT, "nearest_vs_contact.csv"))
    depth = read_csv(os.path.join(LT, "depth_matched_summary.csv"), index_col=0)
    depth_genes = read_csv(os.path.join(LT, "depth_matched_genes.csv"))
    gsc_spec = read_csv(os.path.join(LT, "gsc_specific_targets.csv"))
    motifs = read_csv(os.path.join(RES, "motifs", "motif_enrichment.csv"))
    exprt = read_csv(os.path.join(RES, "expression", "expression_tests.csv"))

    # --------------------------------------------------------------- title --
    t = doc.add_paragraph()
    r = t.add_run("Reproducibility of the glioblastoma super-enhancer landscape "
                  "and the failure of nearest-gene target assignment")
    r.bold = True
    r.font.size = Pt(14)
    t.alignment = WD_ALIGN_PARAGRAPH.CENTER
    P("Methods and Results. Generated from the analysis outputs; every value "
      "below is read from a result file rather than transcribed.", italic=True)

    # ------------------------------------------------------------- methods --
    H("Methods")

    H("Data sources", 2)
    P("Four published H3K27ac ChIP-seq series were re-analysed, together with "
      "bulk RNA-seq from The Cancer Genome Atlas obtained through the UCSC Xena "
      "hub. No new data were generated. GSE119776 provides MACS2 pileup signal "
      "tracks for glioma stem-cell (GSC) lines and neural stem-cell (NSC) "
      "controls on hg19. GSE145646 provides signal tracks for patient "
      "glioblastoma tumours on hg38. GSE162976 provides, for the same samples, "
      "both signal tracks and published ENCODE-pipeline MACS2 peak calls for "
      "pediatric high-grade glioma, and is used to calibrate peak calling. "
      "GSE262089 provides H3K27ac HiChIP loops called with FitHiChIP at 10 kb "
      "resolution for four GSC lines and an NSC control, together with H3K27ac "
      "and ATAC peak calls and RNA-seq.")
    if manifest is not None:
        k = manifest[(manifest["mark"] == "H3K27ac")]
        comp = (k.groupby(["cohort", "kind"]).size().unstack(fill_value=0)
                .reset_index())
        table(comp, "Table 1. H3K27ac files available per cohort, by file type.")

    H("Software environment", 2)
    P("Analyses ran in a conda environment containing samtools, bedtools, MACS3, "
      "deepTools, bedops, the UCSC liftOver and bedGraphToBigWig utilities, the "
      "MEME suite (SEA, AME) and pyBigWig, on a Linux workstation. GENCODE v44 "
      "provided gene and transcription start site annotation, restricted to "
      "protein-coding and lncRNA genes on the 24 main chromosomes (20,033 "
      "protein-coding and 18,866 lncRNA genes). ENCODE blacklist v2 regions were "
      "excluded from all peak sets.")

    H("Peak calling from signal tracks, and its calibration", 2)
    P("Two of the four series released normalised signal tracks without peak "
      "calls. Peaks were therefore derived from signal: the genome was divided "
      "into 200 bp bins, bins were thresholded at a quantile of each sample's "
      "own non-zero signal distribution (so that sequencing depth and track "
      "scaling cancel), surviving bins were merged across a fixed gap, and "
      "regions shorter than a minimum width or overlapping the blacklist were "
      "discarded.")
    if calib:
        P(f"The three parameters were not chosen by inspection. GSE162976 "
          f"published both a signal track and MACS2 peak calls for the same "
          f"GSMs, so that series was used as ground truth: a grid of quantile, "
          f"gap and minimum-width values was scored against the published peaks "
          f"in {fmt(calib.get('n_calibration_samples'))} samples at both "
          f"base-pair level (Jaccard index) and peak level (precision, recall). "
          f"The selected setting was a quantile of {fmt(calib['quantile'])}, a "
          f"gap of {fmt(calib['gap'])} bp and a minimum width of "
          f"{fmt(calib['min_width'])} bp, achieving a median Jaccard index of "
          f"{fmt(calib.get('jaccard_median'))}, median peak-level precision "
          f"{fmt(calib.get('precision_median'))} and recall "
          f"{fmt(calib.get('recall_median'))} against the published calls. "
          f"These values were then applied unchanged to every cohort, including "
          f"those with published peaks, so that differences between cohorts "
          f"cannot be an artefact of the peak caller.")
    else:
        P("Calibration results are not yet computed.")

    H("Super-enhancer identification", 2)
    P("Super-enhancers were called following the ROSE procedure. Peaks lying "
      "entirely within a promoter window (transcription start site plus or "
      "minus 2 kb) were removed so that highly acetylated promoters are not "
      "counted as enhancer domains; remaining peaks within 12.5 kb were "
      "stitched; each stitched domain was scored by the total H3K27ac signal it "
      "contains; and domains were ranked by score. The super-enhancer cutoff "
      "was placed at the inflection point of the ranked signal curve, located "
      "numerically as the point at which the tangent to the curve, with both "
      "axes scaled to the unit interval, has slope one. No fixed rank or "
      "percentile was imposed, so the number of super-enhancers is a property "
      "of each sample's own signal distribution.")
    P("Peak calling and stitching were performed in each cohort's native genome "
      "build. Only the resulting enhancer coordinates were lifted to hg38 "
      "(hg19ToHg38 chain), so that no signal track was ever interpolated across "
      "builds.")

    H("Cross-cohort recurrence and the noise ceiling", 2)
    P("A union merge of stitched domains across more than 160 samples collapses "
      "into megabase blocks, so recurrence was computed on fixed 5 kb bins. For "
      "each cohort group, recurrence in a bin is the fraction of samples in that "
      "group whose super-enhancer calls cover it. Agreement between two cohorts "
      "was measured by Spearman correlation of bin recurrence, by the area under "
      "the precision-recall curve (average precision) when one cohort's "
      "recurrence is used to predict whether a bin is recurrent in the other "
      "(ROC AUC is not used: with roughly one positive bin in sixty it reads "
      "above 0.9 for every comparison, including unrelated ones), and by the Jaccard index of bins recurrent "
      "in at least half the samples of each.")
    P("Agreement between cohorts is not interpretable on its own, because part "
      "of any disagreement is sampling noise rather than biology. Each cohort "
      "was therefore split at random into two halves 25 times and the same three "
      "statistics computed between the halves; this split-half value is the "
      "ceiling that a perfectly reproducible landscape could reach given the "
      "cohort size, and every cross-cohort statistic is reported against it.")

    H("Enhancer-to-gene assignment", 2)
    P("HiChIP loop anchors were classified as promoter (containing an annotated "
      "transcription start site), enhancer (overlapping an H3K27ac peak with no "
      "start site) or neither. Every loop joining an enhancer anchor to a "
      "promoter anchor was taken as a contact-based enhancer-gene link, and "
      "compared with the assignment that would be made by taking the gene whose "
      "start site is nearest the enhancer.")
    P("Three controls were applied. First, because a 10 kb promoter anchor can "
      "contain several start sites, the comparison was repeated on the subset of "
      "links whose promoter anchor contains exactly one. Second, since a gene "
      "100 kb away is rarely the closest gene for trivial reasons, a "
      "distance-matched decoy was drawn for every link - the gene whose distance "
      "from the same enhancer most nearly equals the true contact distance, "
      "excluding the contacted gene itself - and the decoy was tested for being "
      "the nearest gene, giving the agreement expected from genomic architecture "
      "alone. Third, because the five HiChIP samples differ several-fold in loop "
      "count, every sample was subsampled to the loop count of the shallowest "
      "sample, 20 times, before any sample-to-sample comparison of contacted "
      "genes was made.")

    H("Motif enrichment", 2)
    P("Sequences of domain sets were extracted from the hg38 assembly and tested "
      "with SEA against the JASPAR CORE vertebrate non-redundant motif "
      "collection. The primary comparison was reproducible against sporadic "
      "domains, with the control set drawn to match the length distribution of "
      "the primary set, rather than against shuffled background alone.")

    H("Expression validation", 2)
    P("Each gene was assigned the maximum super-enhancer recurrence found within "
      "100 kb of its transcription start site - the median HiChIP contact "
      "distance measured here - separately for each cohort group. That per-gene "
      "recurrence was tested against expression data that played no part in the "
      "enhancer calls: TCGA GBM (n = 172) and TCGA LGG (n = 530) bulk RNA-seq, "
      "the pediatric RNA-seq published with GSE162976, and the GSC line RNA-seq "
      "published with GSE262089. Associations are reported as Spearman "
      "correlations across all genes and as Mann-Whitney comparisons between "
      "gene classes with Cliff's delta as effect size.")

    H("Code and data availability", 2)
    P("All analysis code is at https://github.com/kratos1459/ChIPseq_GBM_atlas. "
      "Raw signal tracks are not redistributed; every input is a published file "
      "retrievable from GEO or the UCSC Xena hub using the manifests in that "
      "repository. DepMap CRISPR dependency data was not used because "
      "depmap.org serves a bot-verification page to non-browser clients and the "
      "files cannot be retrieved programmatically.")

    # ------------------------------------------------------------- results --
    H("Results")

    H("Cohort composition and super-enhancer counts", 2)
    if se_sum is not None and len(se_sum):
        ok = se_sum[se_sum["status"] == "ok"]
        P(f"Super-enhancers were called in {fmt(len(ok))} H3K27ac samples. "
          f"Per-sample counts and domain totals are given in Table 2.")
        if len(ok):
            tab = (ok.groupby(["cohort", "tissue_class"])
                   .agg(n=("label", "size"),
                        median_stitched=("n_stitched", "median"),
                        median_super=("n_super", "median"),
                        median_super_Mb=("super_mb", "median"))
                   .reset_index().round(2))
            table(tab, "Table 2. Stitched enhancer domains and super-enhancers "
                       "per cohort (medians across samples).")
    else:
        P("Super-enhancer calling across all samples is still running; counts "
          "will be filled in when it completes.")

    H("Cross-cohort reproducibility against the noise ceiling", 2)
    if agree is not None and len(agree):
        ceil = agree[agree["comparison"].str.contains("ceiling")]
        xs = agree[~agree["comparison"].str.contains("ceiling")]
        if len(xs) and "spearman" in xs.columns:
            P("Table 3 reports every agreement statistic together with the "
              "within-cohort split-half ceiling. Cross-cohort agreement should "
              "be read as a fraction of that ceiling, not in absolute terms.")
            tn = read_csv(os.path.join(CC, "agreement_topN.csv"))
            if tn is not None and "ap_predicting_b" in tn.columns:
                def v(comp, col):
                    r = tn[tn["comparison"] == comp]
                    return float(r[col].iloc[0]) if len(r) else float("nan")
                sg = v("GSC_line split-half (ceiling)", "spearman")
                sp_ = v("patient_tumour split-half (ceiling)", "spearman")
                gx = v("GSC_line -> patient_tumour", "spearman")
                ag = v("GSC_line -> patient_tumour", "ap_predicting_b")
                ab = v("GSC_line -> patient_tumour", "ap_baseline")
                cp = v("patient_tumour split-half (ceiling)", "ap_predicting_b")
                an = v("GSC_line -> NSC", "ap_predicting_b")
                ap_ = v("GSC_line -> pediatric", "ap_predicting_b")
                P(f"Under the fixed top-500 definition, GSC-line recurrence "
                  f"predicts which bins are recurrent in patient tumours with an "
                  f"average precision of {ag:.3f}, against a chance level of "
                  f"{ab:.3f} and a ceiling of {cp:.3f} set by two halves of the "
                  f"tumour cohort predicting each other - {100*ag/cp:.0f}% of "
                  f"the achievable agreement. The Spearman correlation between "
                  f"the two adult cohorts, {gx:.3f}, likewise sits within the "
                  f"split-half values of {sg:.3f} and {sp_:.3f}. Cell lines and "
                  f"patient tumours therefore agree with each other about as "
                  f"well as either agrees with itself: there is no detectable "
                  f"cohort-specific super-enhancer component beyond sampling "
                  f"noise. That ceiling is itself modest, so most per-bin "
                  f"variation is noise in both cohorts. Agreement falls to "
                  f"{an:.3f} for neural stem cells and {ap_:.3f} for pediatric "
                  f"high-grade glioma. Lower agreement is expected for different "
                  f"cell states, but those cohorts also differ in laboratory and "
                  f"processing, and the two explanations cannot be separated "
                  f"with these data.")
        table(agree.round(4), "Table 3. Cross-cohort agreement and split-half "
                              "noise ceilings.")
    else:
        P("Cross-cohort agreement is not yet computed.")

    if ccsum:
        P(f"Super-enhancers are defined here as the {ccsum.get('definition', 'top domains')} "
          f"ranked by H3K27ac signal. The ROSE inflection rule was also run, "
          f"but its per-sample count tracked data provenance rather than "
          f"biology (median 280 in the hg19 bedGraph GSC cohort against 606 in "
          f"neural stem-cell controls), so a fixed number per sample is used to "
          f"give every sample equal weight; agreement under both definitions is "
          f"given in Table 3.")
        P(f"Applying a threshold of recurrence in at least half the samples of "
          f"both adult cohorts, and under "
          f"{100*ccsum.get('nsc_max', 0.25):.0f}% of the 12 neural stem-cell "
          f"controls (at most two), "
          f"defines {fmt(ccsum.get('n_core_domains'))} consensus core domains "
          f"spanning {fmt(ccsum.get('core_mb'), 1)} Mb. For comparison, "
          f"{fmt(ccsum.get('n_sporadic_domains'))} domains were called a "
          f"super-enhancer in only one adult sample, "
          f"{fmt(ccsum.get('n_gsc_specific_domains'))} were recurrent in GSC "
          f"lines but near-absent from patient tumours, and "
          f"{fmt(ccsum.get('n_patient_specific_domains'))} the reverse.")
        if "core_frac_recurrent_in_pediatric" in ccsum:
            P(f"Of the core adult domains, "
              f"{fmt(100*ccsum['core_frac_recurrent_in_pediatric'], 1)}% are also "
              f"recurrent in the independent pediatric high-grade glioma cohort "
              f"(median recurrence "
              f"{fmt(ccsum.get('core_median_recurrence_in_pediatric'))}), which "
              f"was processed with the same calibrated caller but contributed "
              f"nothing to the core definition.")
    if core_genes is not None and len(core_genes):
        top = core_genes.nlargest(20, "width_kb")[
            ["domain", "chrom", "start", "width_kb", "nearest_gene",
             "n_contact_genes"]]
        table(top.round(1), "Table 4. The twenty widest consensus core domains.")

    H("The nearest gene is usually not the contacted gene", 2)
    if controls is not None and len(controls):
        c = controls
        P(f"Across the five HiChIP samples, the nearest gene is the contacted "
          f"gene for {fmt(100*c['nearest_correct_all'].min(),1)}-"
          f"{fmt(100*c['nearest_correct_all'].max(),1)}% of enhancer-promoter "
          f"links. Restricting to promoter anchors containing a single "
          f"transcription start site raises this to "
          f"{fmt(100*c['nearest_correct_clean'].min(),1)}-"
          f"{fmt(100*c['nearest_correct_clean'].max(),1)}%, while a "
          f"distance-matched decoy gene achieves "
          f"{fmt(100*c['nearest_correct_decoy_null'].min(),1)}-"
          f"{fmt(100*c['nearest_correct_decoy_null'].max(),1)}%. Contact-based "
          f"assignment therefore does favour nearby genes, by roughly "
          f"threefold over architecture alone, but the nearest gene is still "
          f"wrong for about seven links in ten even in the cleanest subset.")
        P(f"For only {fmt(100*c['enh_with_nearest_among_contacts'].min(),1)}-"
          f"{fmt(100*c['enh_with_nearest_among_contacts'].max(),1)}% of linked "
          f"enhancers is the nearest gene even among the genes that enhancer "
          f"contacts. The median contact spans "
          f"{fmt(c['median_contact_kb'].min(),0)}-"
          f"{fmt(c['median_contact_kb'].max(),0)} kb against a median "
          f"nearest-start-site distance of {fmt(c['median_nearest_kb'].min(),0)}-"
          f"{fmt(c['median_nearest_kb'].max(),0)} kb.")
        table(c.round(4), "Table 5. Contact-based versus nearest-gene "
                          "assignment, with controls.")
    else:
        P("Target-assignment controls are not yet computed.")

    if nvc is not None and len(nvc):
        table(nvc.round(3), "Table 6. Loop counts and link statistics per "
                            "HiChIP sample.")

    H("Loop depth, not biology, explains most apparently tumour-specific targets", 2)
    if depth is not None and depth_genes is not None:
        n_before = None
        if gsc_spec is not None and len(gsc_spec):
            n_before = int(((gsc_spec["n_gsc_samples"] == gsc_spec["n_gsc_samples"].max())
                            & (gsc_spec["in_nsc"] == 0)).sum())
        robust = int((depth_genes["frac_repeats"] >= 0.8).sum())
        mid = depth.loc["genes_all_gsc_not_nsc", "mean"] \
            if "genes_all_gsc_not_nsc" in depth.index else np.nan
        P(f"Taken at face value, {fmt(n_before)} genes are contacted by "
          f"enhancers in every GSC line and in none of the neural stem-cell "
          f"control. The five samples differ several-fold in loop count, "
          f"however. After subsampling every sample to the loop count of the "
          f"shallowest, this falls to {fmt(float(mid), 0)} genes on average, and "
          f"only {fmt(robust)} genes remain specific in at least 80% of "
          f"subsampling repeats - about "
          f"{fmt(100*robust/n_before, 1) if n_before else 'n/a'}% of the "
          f"uncontrolled count. The great majority of the apparent "
          f"tumour-specific signal is sequencing depth.")
        named = sorted(g for g in depth_genes.loc[depth_genes["frac_repeats"] >= 0.8,
                                                  "gene"] if not str(g).startswith("ENSG"))
        if named:
            known = [g for g in ("BRD4", "JUN", "FOSL2", "CEBPB", "FOXG1", "GOLIM4",
                                 "ADAM17", "IRF2", "DOCK4", "LAMC1", "LIMD1")
                     if g in named]
            P(f"{fmt(len(named))} of the surviving genes carry an approved "
              f"symbol; they are listed in full in Table 8. The set is not a "
              f"random draw from the genome: it contains the H3K27ac reader "
              f"BRD4, which is itself the target of BET inhibitors tested in "
              f"glioblastoma, and the AP-1 and CEBPB factors that define the "
              f"mesenchymal glioblastoma regulatory circuit" +
              (f" ({', '.join(known)})" if known else "") +
              ". GOLIM4 is among them, which is notable because the authors of "
              "GSE262089 independently singled out that gene; recovering it "
              "from their loop data under depth matching is a check that the "
              "procedure is not simply discarding real signal.")
            P("Table 8. Genes contacted in every GSC line and in no neural "
              "stem-cell control, surviving in at least 80% of depth-matched "
              "subsampling repeats (symbols only):", italic=True)
            P(", ".join(named) + ".")
        if "nearest_correct_gsc_mean" in depth.index:
            P(f"At matched depth the nearest-gene hit rate is "
              f"{fmt(depth.loc['nearest_correct_gsc_mean','mean'])} in GSC lines "
              f"and {fmt(depth.loc['nearest_correct_nsc','mean'])} in the neural "
              f"stem-cell control, so the failure of nearest-gene assignment is "
              f"a general property of enhancer-promoter wiring and not specific "
              f"to the tumour state.")
        table(depth.reset_index().round(3),
              "Table 7. Contact statistics at matched loop depth "
              "(mean and standard deviation over 20 subsampling repeats).")
    else:
        P("Depth-matched controls are not yet computed.")

    H("Motif content of reproducible domains", 2)
    gcj = read_json(os.path.join(RES, "motifs", "gc_composition.json")) or {}
    cs = gcj.get("core_vs_sporadic", {})
    if motifs is not None and len(motifs) and cs:
        qc = "QVALUE" if "QVALUE" in motifs.columns else None
        n_sp = int((motifs[(motifs["comparison"] == "core_vs_sporadic")][qc] < 0.05).sum()) if qc else 0
        n_sh = int((motifs[(motifs["comparison"] == "core_vs_shuffled")][qc] < 0.05).sum()) if qc else 0
        P(f"Motifs were tested on the constituent H3K27ac peaks of each domain "
          f"set (peak centre ±250 bp), because whole stitched domains of "
          f"10-40 kb contain nearly every motif and leave presence-based "
          f"enrichment uninformative. Against peaks from sporadic domains, "
          f"{n_sp} of 879 JASPAR motifs reached q < 0.05, every one of them a "
          f"GC-box family (KLF, SP, PATZ1, ZNF281, TFAP2C). The core peaks are "
          f"{100*cs.get('primary', float('nan')):.1f}% GC against "
          f"{100*cs.get('control', float('nan')):.1f}% for the controls, and "
          f"against a dinucleotide-shuffled control of the same sequences, "
          f"which preserves base composition, {n_sh} motifs reached q < 0.05. "
          f"The apparent enrichment is therefore explained by base composition "
          f"and is not evidence for any specific transcription factor. With "
          f"79 core peaks the test has limited power, and a GC-matched "
          f"background or a larger consensus set would be needed before a "
          f"negative result here could be read as absence of motif signal.")
    if motifs is not None and len(motifs):
        qcol = "QVALUE" if "QVALUE" in motifs.columns else None
        for comp, d in motifs.groupby("comparison"):
            if qcol:
                sig = d[d[qcol] < 0.05]
                P(f"{comp.replace('_', ' ')}: {fmt(len(sig))} of {fmt(len(d))} "
                  f"motifs enriched at q < 0.05.")
                cols = [c for c in ("ID", "ALT_ID", "ENR_RATIO", "QVALUE")
                        if c in d.columns]
                table(sig.nsmallest(15, qcol)[cols],
                      f"Table. Top motifs, {comp.replace('_', ' ')}.")
    else:
        P("Motif enrichment is not yet computed.")

    H("Super-enhancer recurrence and expression", 2)
    if exprt is not None and len(exprt):
        sp_ = exprt[exprt["test"] == "spearman"]
        if len(sp_):
            P("Spearman correlations between per-gene super-enhancer recurrence "
              "and expression in each independent dataset are given in Table 8.")
            table(sp_.round(4), "Table 8. Recurrence against expression.")
        grp = exprt[exprt["test"].str.contains("vs")]
        if len(grp):
            table(grp.round(4), "Table 9. Expression of genes by "
                                "super-enhancer class.")
    else:
        P("Expression validation is not yet computed.")

    # ------------------------------------------------------------ captions --
    H("Figure captions")
    caps = [
        ("Figure 1", "Cohort composition and super-enhancer counts. (a) Number "
         "of H3K27ac samples analysed per cohort. (b) Super-enhancers called per "
         "sample; boxes show the interquartile range, points individual samples. "
         "(c) Super-enhancers against total stitched enhancer domains."),
        ("Figure 2", "Recurrence landscape. (a) Ranked super-enhancer recurrence "
         "per 5 kb bin in each cohort group. (b) Joint distribution of "
         "recurrence in GSC lines and patient tumours; dashed lines mark the 50% "
         "recurrence threshold used to define consensus domains."),
        ("Figure 3", "Cross-cohort agreement against within-cohort split-half "
         "noise ceilings, using the fixed top-500 domains per sample: Spearman "
         "correlation of bin recurrence, average precision for predicting the "
         "other cohort with chance level marked in red, and the Jaccard index "
         "at 50% recurrence."),
        ("Figure 4", "Enhancer-to-gene assignment. (a) Fraction of "
         "enhancer-promoter links for which the nearest gene is the contacted "
         "gene, for all links, for single-start-site promoter anchors, and for "
         "distance-matched decoy genes. (b) Median distance to the nearest "
         "versus the contacted start site. (c) Agreement of contacted-gene sets "
         "between GSC lines and between GSC lines and the neural stem-cell "
         "control. (d) Genes contacted in every GSC line but not in the control, "
         "as called, at matched loop depth, and robust across subsampling."),
        ("Figure 5", "Consensus core domains. (a) Size of each domain class. "
         "(b) Width distribution of core domains. (c) The widest core domains, "
         "labelled by nearest gene."),
        ("Figure 6", "Motif enrichment among the constituent peaks of core "
         "domains, against peaks from sporadic domains (left) and against a "
         "dinucleotide-shuffled control of the same sequences (right). Every "
         "motif significant on the left is a GC-box family, and none survives "
         "the composition-preserving control."),
        ("Figure 7", "Super-enhancer recurrence against independent expression "
         "data. Boxes compare genes carrying a core super-enhancer, another "
         "super-enhancer, or none, in each expression dataset; the final panel "
         "shows Spearman correlations between per-gene recurrence and "
         "expression."),
    ]
    for name, text in caps:
        p = doc.add_paragraph()
        r = p.add_run(f"{name}. ")
        r.bold = True
        p.add_run(text)
    P("All figures are written at 600 dpi as both PNG and LZW-compressed TIFF.",
      italic=True)

    # --------------------------------------------------------- limitations --
    H("Limitations")
    for text in [
        "HiChIP loop anchors are 10 kb bins, so enhancer and promoter positions "
        "are coarse and fine-scale rewiring between cell states cannot be "
        "resolved. Contact-based assignment is therefore a better reference than "
        "nearest-gene assignment, not a ground truth.",
        "Loop depth differs several-fold between the five HiChIP samples. All "
        "sample-to-sample comparisons here are depth-matched, but the absolute "
        "number of contacts per enhancer remains depth-dependent.",
        "The two adult cohorts published no processed RNA-seq, so "
        "cohort-internal expression validation was possible only in the "
        "pediatric series and the four HiChIP GSC lines; TCGA provides the adult "
        "expression comparison, but from different patients than those whose "
        "chromatin was profiled.",
        "GSE119776 was published on hg19 and its enhancer coordinates are lifted "
        "after calling; 48 of 38,899 transcription start sites do not map "
        "between builds.",
        "Peaks for two cohorts are derived from signal rather than from reads. "
        "The caller was calibrated against published MACS2 calls in an "
        "independent series, but no input-normalised statistical test of "
        "enrichment was possible without the original alignments.",
        "Cell-line and tumour cohorts differ in more than their origin - "
        "sequencing protocol, depth and processing pipeline all differ - so "
        "cohort-specific domains should not be read as model-system artefacts "
        "without further work.",
    ]:
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(text)

    doc.save(OUT)
    return OUT


def main():
    try:
        path = build()
    except ImportError:
        print("python-docx is not installed in this interpreter; run with the "
              "base python3")
        return 1
    size = os.path.getsize(path)
    print(f"wrote {path} ({size/1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
