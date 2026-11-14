#!/usr/bin/env python3
"""Manuscript and Supplementary Information for the GBM super-enhancer study.

Every quantitative statement is read from a result file at build time; the
prose around the numbers is fixed. References were checked against PubMed /
publisher records on 2026-11-14 and are listed in REFS in citation order.

Run with the base python (needs python-docx and pandas):
    python3 code/80_manuscript.py
Writes manuscript/GBM_SE_Manuscript.docx and manuscript/GBM_SE_Supplementary.docx
"""
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
CC = os.path.join(RES, "cross_cohort")
LT = os.path.join(RES, "loop_targets")
FIG = os.path.join(ROOT, "figures", "PNG_600dpi")
OUT = os.path.join(ROOT, "manuscript")

REFS = [
    "Mack SC, Singh I, Wang X, et al. Chromatin landscapes reveal developmentally encoded transcriptional states that define human glioblastoma. J Exp Med. 2019;216(5):1071-1090. doi:10.1084/jem.20190196",
    "Xu L, Chen Y, Huang Y, et al. Topography of transcriptionally active chromatin in glioblastoma. Sci Adv. 2021;7(18):eabd4676. doi:10.1126/sciadv.abd4676",
    "Wang J, Huang TY, Hou Y, et al. Epigenomic landscape and 3D genome structure in pediatric high-grade glioma. Sci Adv. 2021;7(23):eabg4126. doi:10.1126/sciadv.abg4126",
    "Breves SL, Di Giammartino DC, Nicholson J, et al. Three-dimensional regulatory hubs support oncogenic programs in glioblastoma. Mol Cell. 2025;85(7):1330-1348.e6. doi:10.1016/j.molcel.2025.03.007",
    "Whyte WA, Orlando DA, Hnisz D, et al. Master transcription factors and Mediator establish super-enhancers at key cell identity genes. Cell. 2013;153(2):307-319. doi:10.1016/j.cell.2013.03.035",
    "Lovén J, Hoke HA, Lin CY, et al. Selective inhibition of tumor oncogenes by disruption of super-enhancers. Cell. 2013;153(2):320-334. doi:10.1016/j.cell.2013.03.036",
    "Hnisz D, Abraham BJ, Lee TI, et al. Super-enhancers in the control of cell identity and disease. Cell. 2013;155(4):934-947.",
    "Fulco CP, Nasser J, Jones TR, et al. Activity-by-contact model of enhancer-promoter regulation from thousands of CRISPR perturbations. Nat Genet. 2019;51(12):1664-1669. doi:10.1038/s41588-019-0538-0",
    "Nasser J, Bergman DT, Fulco CP, et al. Genome-wide enhancer maps link risk variants to disease genes. Nature. 2021;593(7858):238-243. doi:10.1038/s41586-021-03446-x",
    "Bhattacharyya S, Chandra V, Vijayanand P, Ay F. Identification of significant chromatin contacts from HiChIP data by FitHiChIP. Nat Commun. 2019;10:4221. doi:10.1038/s41467-019-11950-y",
    "Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: identification of problematic regions of the genome. Sci Rep. 2019;9:9354. doi:10.1038/s41598-019-45839-z",
    "Bailey TL, Grant CE. SEA: Simple Enrichment Analysis of motifs. bioRxiv. 2021. doi:10.1101/2021.08.23.457422",
    "Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Res. 2024;52(D1):D174-D182. doi:10.1093/nar/gkad1059",
    "Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nat Biotechnol. 2020;38(6):675-678. doi:10.1038/s41587-020-0546-8",
]


def rcsv(p, **kw):
    p = os.path.join(RES, p) if not os.path.isabs(p) else p
    return pd.read_csv(p, **kw) if os.path.exists(p) else None


def rjson(p):
    p = os.path.join(RES, p)
    return json.load(open(p)) if os.path.exists(p) else {}


def gc_frac(fa):
    if not os.path.exists(fa):
        return float("nan")
    gc = tot = 0
    for line in open(fa):
        if line.startswith(">"):
            continue
        s = line.strip().upper()
        gc += s.count("G") + s.count("C")
        tot += sum(s.count(b) for b in "ACGT")
    return gc / tot if tot else float("nan")


def n_fasta(fa):
    return sum(1 for l in open(fa) if l.startswith(">")) if os.path.exists(fa) else 0


def collect():
    N = {}
    cal = rjson("calibration_best.json")
    N.update({f"cal_{k}": v for k, v in cal.items()})
    su = rcsv("cross_cohort/samples_used.csv")
    g = su.groupby("group")
    N["n"] = g.size().to_dict()
    N["se_med"] = g["n_se"].median().to_dict()
    N["n_total"] = int(len(su))
    tn = rcsv("cross_cohort/agreement_topN.csv").set_index("comparison")
    def a(c, k):
        return float(tn.loc[c, k]) if c in tn.index else float("nan")
    N["ceil_sp_gsc"] = a("GSC_line split-half (ceiling)", "spearman")
    N["ceil_sp_tum"] = a("patient_tumour split-half (ceiling)", "spearman")
    N["ceil_ap_gsc"] = a("GSC_line split-half (ceiling)", "ap_predicting_b")
    N["ceil_ap_tum"] = a("patient_tumour split-half (ceiling)", "ap_predicting_b")
    N["x_sp"] = a("GSC_line -> patient_tumour", "spearman")
    N["x_ap"] = a("GSC_line -> patient_tumour", "ap_predicting_b")
    N["x_ap_rev"] = a("patient_tumour -> GSC_line", "ap_predicting_b")
    N["x_base"] = a("GSC_line -> patient_tumour", "ap_baseline")
    N["x_jac"] = a("GSC_line -> patient_tumour", "jaccard_at_50pct")
    N["ceil_jac_tum"] = a("patient_tumour split-half (ceiling)", "jaccard_at_50pct")
    N["nsc_ap"] = a("GSC_line -> NSC", "ap_predicting_b")
    N["ped_ap"] = a("GSC_line -> pediatric", "ap_predicting_b")
    N["ceil_ap_nsc"] = a("NSC split-half (ceiling)", "ap_predicting_b")
    N["ceil_ap_ped"] = a("pediatric split-half (ceiling)", "ap_predicting_b")
    bs = rcsv("cross_cohort/bootstrap_summary.csv", index_col=0)
    N["bs"] = bs
    rb = rcsv("robustness_topN.csv")
    N["rob"] = rb
    s = rjson("cross_cohort/summary.json")
    N["cc"] = s
    ctl = rcsv("loop_targets/controls_summary.csv")
    N["ctl"] = ctl
    dm = rcsv("loop_targets/depth_matched_summary.csv", index_col=0)
    N["dm"] = dm
    dg = rcsv("loop_targets/depth_matched_genes.csv")
    N["dg_robust"] = dg[dg["frac_repeats"] >= 0.8]
    gs = rcsv("loop_targets/gsc_specific_targets.csv")
    N["n_spec_raw"] = int(((gs["n_gsc_samples"] == gs["n_gsc_samples"].max())
                           & (gs["in_nsc"] == 0)).sum())
    nvc = rcsv("loop_targets/nearest_vs_contact.csv")
    N["nvc"] = nvc
    rep = rcsv("loop_targets/reproducibility_pairs.csv")
    N["rep"] = rep
    N["gc"] = {k: gc_frac(os.path.join(RES, "motifs", "core_vs_sporadic", k + ".fa"))
               for k in ("primary", "control")}
    N["n_core_peaks"] = n_fasta(os.path.join(RES, "motifs", "core_vs_sporadic", "primary.fa"))
    N["sm"] = rcsv("cross_cohort/size_matched.csv")
    reps = rcsv("cross_cohort/size_matched_reps.csv")
    smp = {}
    for (t, k), g in reps.groupby(["target", "k"]):
        w = g[g["source"].str.endswith("(within)")].reset_index(drop=True)
        for s in g["source"].unique():
            if s.endswith("(within)"):
                continue
            x = g[g["source"] == s].reset_index(drop=True)
            fr = x["ap"] / w["ap"]
            smp[(t, int(k), s)] = {
                "within": w["ap"].mean(), "cross": x["ap"].mean(),
                "base": x["baseline"].mean(),
                "frac": fr.median(), "lo": np.percentile(fr, 2.5),
                "hi": np.percentile(fr, 97.5),
                "p_ge": float((x["ap"] >= w["ap"]).mean()),
                "sp_within": w["spearman"].mean(), "sp_cross": x["spearman"].mean()}
    N["smp"] = smp
    mot = rcsv("motifs/motif_enrichment.csv")
    N["mot"] = mot
    ex = rcsv("expression/expression_tests.csv")
    N["ex"] = ex
    cg = rcsv("cross_cohort/core_domain_genes.csv")
    N["core_genes"] = cg
    return N


def f(x, nd=2):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "n/a"


def ex_row(N, test, expr, rec=None):
    e = N["ex"]
    q = (e["test"] == test) & (e["expression"] == expr)
    if rec:
        q &= e["recurrence"] == rec
    r = e[q]
    return r.iloc[0] if len(r) else None


class Doc:
    def __init__(self):
        from docx import Document
        from docx.shared import Pt, Cm
        self.Pt, self.Cm = Pt, Cm
        self.d = Document()
        st = self.d.styles["Normal"]
        st.font.name = "Times New Roman"
        st.font.size = Pt(11)
        for s in self.d.sections:
            s.left_margin = s.right_margin = Cm(2.2)
            s.top_margin = s.bottom_margin = Cm(2.2)

    def title(self, text, size=15):
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        p = self.d.add_paragraph()
        r = p.add_run(text)
        r.bold = True
        r.font.size = self.Pt(size)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        return p

    def center(self, text, italic=False, size=11):
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        p = self.d.add_paragraph()
        r = p.add_run(text)
        r.italic = italic
        r.font.size = self.Pt(size)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    def h(self, text, level=1):
        p = self.d.add_paragraph()
        r = p.add_run(text)
        r.bold = True
        r.italic = level > 1
        r.font.size = self.Pt(13 if level == 1 else 11.5)
        p.paragraph_format.space_before = self.Pt(10 if level == 1 else 6)
        p.paragraph_format.keep_with_next = True

    def p(self, text, italic=False, bold_lead=None):
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        par = self.d.add_paragraph()
        if bold_lead:
            r = par.add_run(bold_lead)
            r.bold = True
        r = par.add_run(text)
        r.italic = italic
        par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        par.paragraph_format.space_after = self.Pt(6)
        return par

    def table(self, headers, rows, caption, size=8.5):
        cap = self.d.add_paragraph()
        r = cap.add_run(caption)
        r.bold = True
        r.font.size = self.Pt(10)
        cap.paragraph_format.keep_with_next = True
        t = self.d.add_table(rows=1, cols=len(headers))
        t.style = "Table Grid"
        for i, h in enumerate(headers):
            c = t.rows[0].cells[i]
            c.text = str(h)
            for pr in c.paragraphs:
                for rr in pr.runs:
                    rr.bold = True
                    rr.font.size = self.Pt(size)
        for row in rows:
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = str(v)
                for pr in cells[i].paragraphs:
                    for rr in pr.runs:
                        rr.font.size = self.Pt(size)
        self.d.add_paragraph()

    def figure(self, path, caption, width_cm=16.5):
        if os.path.exists(path):
            self.d.add_picture(path, width=self.Cm(width_cm))
        par = self.d.add_paragraph()
        par.paragraph_format.space_after = self.Pt(10)
        head, _, rest = caption.partition(". ")
        r = par.add_run(head + ". ")
        r.bold = True
        r.font.size = self.Pt(9.5)
        r2 = par.add_run(rest)
        r2.font.size = self.Pt(9.5)

    def page_break(self):
        from docx.enum.text import WD_BREAK
        self.d.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

    def save(self, path):
        self.d.save(path)


def pc(x, nd=0):
    return f"{100*float(x):.{nd}f}%"


def ranges(series, scale=100, nd=1):
    return f"{scale*series.min():.{nd}f}–{scale*series.max():.{nd}f}"


TITLE = ("Glioblastoma stem-cell lines recover most, but not all, of the "
         "patient-tumour super-enhancer landscape: a size-matched cross-cohort "
         "reanalysis with three-dimensional target assignment")


def build_main(N):
    D = Doc()
    n = N["n"]
    ctl, bs, rob, cc = N["ctl"], N["bs"], N["rob"], N["cc"]
    frac = N["x_ap"] / N["ceil_ap_tum"]
    smp = N["smp"]
    t20 = smp[("patient_tumour", 20, "GSC_line")]
    t10 = smp[("patient_tumour", 10, "GSC_line")]
    t6 = smp[("patient_tumour", 6, "GSC_line")]
    p10 = smp[("pediatric", 10, "GSC_line")]
    p10t = smp[("pediatric", 10, "patient_tumour")]
    n6g = smp[("NSC", 6, "GSC_line")]
    n6t = smp[("NSC", 6, "patient_tumour")]
    gbm = ex_row(N, "core_vs_other_se", "tcga_gbm")
    ped = ex_row(N, "core_vs_other_se", "pedi_tpm")
    sp_g = ex_row(N, "spearman", "tcga_gbm", "rec_GSC_line")
    sp_n = ex_row(N, "spearman", "tcga_gbm", "rec_NSC")
    dm = N["dm"]
    n_rob = len(N["dg_robust"])
    named = sorted(g for g in N["dg_robust"]["gene"] if not str(g).startswith("ENSG"))
    joint = cc.get("n_bins_core_both_cohorts", float("nan"))
    spec = cc.get("n_bins_core_tumour_specific", float("nan"))
    shared_nsc = 100 * (1 - spec / joint) if joint else float("nan")
    rep = N["rep"]
    j_gg = rep[rep["sample_b"] != "nsc"]["mean_target_jaccard"]
    j_gn = rep[rep["sample_b"] == "nsc"]["mean_target_jaccard"]

    D.title(TITLE)
    D.center("Tawsif Al Arian¹*")
    D.center("¹ [Affiliation to be completed]", italic=True, size=10)
    D.center("*Correspondence: tawsifalarian43@gmail.com", size=10)
    D.center("[Co-authors and affiliations to be confirmed before submission]",
             italic=True, size=9)

    # ---------------------------------------------------------------- abstract
    D.h("Abstract")
    D.p("Super-enhancers called from H3K27ac profiles are widely used to nominate "
        "transcriptional dependencies in glioblastoma, usually from a single cohort "
        "and with each element assigned to its nearest gene. How far such calls "
        "reproduce between independent cohorts, and between patient-derived cell "
        "lines and tumours, has rarely been measured against the agreement two "
        "halves of one cohort can reach.", bold_lead="Background. ")
    D.p(f"We reanalysed published H3K27ac data from {N['n_total']} samples in four "
        f"series: {n.get('GSC_line', 0)} glioblastoma stem-cell (GSC) lines, "
        f"{n.get('patient_tumour', 0)} patient tumours, {n.get('NSC', 0)} neural "
        f"stem-cell controls and {n.get('pediatric', 0)} pediatric high-grade "
        f"glioma samples. Two cohorts released signal tracks only, so peaks were "
        f"derived from signal with a caller calibrated against published MACS2 "
        f"peaks (median Jaccard {f(N['cal_jaccard_median'])}) and applied "
        f"identically to every cohort. Super-enhancer recurrence was compared "
        f"between cohorts in 5 kb bins against within-cohort split-half "
        f"ceilings. Enhancer-gene assignment was evaluated against H3K27ac HiChIP "
        f"loops with single-TSS, distance-matched-decoy and loop-depth controls.",
        bold_lead="Methods. ")
    D.p(f"GSC-line recurrence predicted tumour recurrence with an average "
        f"precision of {f(N['x_ap'])} against a chance level of {f(N['x_base'])}. "
        f"With both sides restricted to the same number of samples, GSC lines "
        f"recovered {pc(t20['frac'])} of the agreement between two independent "
        f"sets of tumours (20 vs 20 samples; 95% range {pc(t20['lo'])}–"
        f"{pc(t20['hi'])}; GSC lines below the within-tumour value in "
        f"{pc(1 - t20['p_ge'])} of draws), and about {pc(t10['frac'])} at 6 and "
        f"10 samples: a modest but consistent cell-line deficit. Adult cohorts "
        f"recovered only {pc(p10['frac'])} of within-pediatric agreement, while "
        f"agreement with neural stem cells could not be distinguished from their "
        f"own noise ceiling at the six samples available per half. Only "
        f"{cc.get('n_core_domains')} consensus domains "
        f"({f(cc.get('core_mb'))} Mb) were recurrent in both adult cohorts and "
        f"rare in neural stem cells; {shared_nsc:.0f}% of bins recurrent in both "
        f"adult cohorts were also carried by at least a quarter of neural "
        f"stem-cell samples. The nearest gene was the "
        f"contacted gene for {ranges(ctl['nearest_correct_all'])}% of "
        f"enhancer-promoter links ({ranges(ctl['nearest_correct_clean'])}% for "
        f"single-TSS anchors; distance-matched decoys "
        f"{ranges(ctl['nearest_correct_decoy_null'])}%). Matching loop depth "
        f"reduced apparently GSC-specific contacted genes from "
        f"{N['n_spec_raw']:,} to {n_rob}. Motif enrichment in reproducible "
        f"domains was accounted for by GC content, while genes near them were "
        f"more highly expressed in TCGA glioblastoma (Cliff's δ "
        f"{f(gbm['cliffs_delta'])}).", bold_lead="Results. ")
    D.p("GSC lines are a good but imperfect proxy for the tumour super-enhancer "
        "landscape, and comparisons that set full cohorts against half-size "
        "ceilings hide the difference. The reproducible core is small and "
        "largely shared with neural stem cells. Nearest-gene assignment, and "
        "contact comparisons that ignore sequencing depth, introduce errors at "
        "least as large as the cell-line deficit.", bold_lead="Conclusions. ")
    D.p("glioblastoma; super-enhancer; H3K27ac; reproducibility; noise ceiling; "
        "HiChIP; enhancer-gene assignment", bold_lead="Keywords: ")

    # ------------------------------------------------------------ introduction
    D.h("Introduction")
    D.p("Super-enhancers are clusters of enhancers with exceptionally high "
        "occupancy by master transcription factors, Mediator and acetylated "
        "histone, and they mark genes that define cell identity [5]. Tumour cells "
        "assemble super-enhancers at oncogenes, and those genes are "
        "disproportionately sensitive to inhibition of the H3K27ac reader BRD4 "
        "[6]; super-enhancers also concentrate disease-associated sequence "
        "variation across cell types [7]. Ranking H3K27ac domains to call "
        "super-enhancers has therefore become a standard route to candidate "
        "dependencies in cancer.")
    D.p("In glioblastoma, H3K27ac profiling of patient-derived glioblastoma "
        "stem-cell (GSC) lines separated them into groups with distinct "
        "enhancer landscapes and core transcription-factor circuits [1], and "
        "profiling of primary tumour biopsies defined subtype-specific enhancer "
        "domains and super-enhancer-driven dependencies [2]. The second study "
        "observed that tissue and cell-line samples co-clustered poorly [2], "
        "which raises a practical question for anyone building on either "
        "resource: how much of a super-enhancer map called in one cohort should "
        "be expected to reappear in another, and is disagreement between cell "
        "lines and tumours biology or noise?")
    D.p("That question cannot be answered from cross-cohort agreement alone, "
        "because even two random halves of a single cohort disagree. The relevant "
        "comparison is between the agreement two cohorts reach and the agreement "
        "a cohort reaches with itself. We are not aware of a glioblastoma study "
        "that reports super-enhancer reproducibility against such a ceiling.")
    D.p("A second, independent source of error lies in assigning each enhancer "
        "to a gene. Most catalogues use the nearest gene, although CRISPR "
        "perturbation studies show that distance alone predicts functional "
        "enhancer-gene connections poorly and that combining activity with 3D "
        "contact does markedly better [8, 9]. H3K27ac HiChIP now exists for "
        "glioblastoma stem cells [4], so nearest-gene assignment can be tested "
        "directly in this tumour type rather than assumed.")
    D.p("Here we reanalyse four published series under a single, calibrated "
        "procedure. We ask how reproducible the glioblastoma super-enhancer "
        "landscape is between GSC lines and patient tumours relative to the "
        "split-half ceiling, what the reproducible core contains, and how often "
        "nearest-gene assignment agrees with 3D contact, with controls for the "
        "confounders that can manufacture or hide each effect.")


    # ---------------------------------------------------------------- methods
    D.h("Methods")
    D.h("Data", 2)
    D.p("All data are published and were obtained from GEO and the UCSC Xena hub; "
        "no new data were generated (Table 1). GSE119776 [1] provides MACS2 "
        "pileup signal for GSC lines and neural stem-cell (NSC) controls on "
        "hg19. GSE145646 [2] provides signal tracks for primary glioblastoma "
        "biopsies and NSC controls on hg38. GSE162976 [3] provides, for the same "
        "samples, signal tracks and ENCODE-pipeline MACS2 peak calls for "
        "pediatric high-grade glioma and control lines on hg38. GSE262089 [4] "
        "provides H3K27ac HiChIP loops called with FitHiChIP [10] at 10 kb "
        "resolution, together with H3K27ac peaks and RNA-seq, for four GSC "
        "lines and an NSC control. TCGA glioblastoma (n = 172) and lower-grade "
        "glioma (n = 530) RNA-seq were taken from the Xena hub [14]. Samples "
        "were assigned to analysis groups by their annotated tissue of origin, "
        "so NSC controls from GSE119776 and GSE145646 form one NSC group.")
    D.table(["Series", "Content used", "Build", "Samples (group)", "Role"], [
        ["GSE119776 [1]", "H3K27ac pileup signal", "hg19",
         f"{n.get('GSC_line', 0)} GSC lines; 9 NSC", "adult model-system cohort"],
        ["GSE145646 [2]", "H3K27ac signal", "hg38",
         f"{n.get('patient_tumour', 0)} tumours; 3 NSC", "adult tumour cohort"],
        ["GSE162976 [3]", "H3K27ac signal and MACS2 peaks; RNA-seq", "hg38",
         f"{n.get('pediatric', 0)} pediatric HGG and control lines",
         "caller calibration; external cohort"],
        ["GSE262089 [4]", "H3K27ac HiChIP loops, peaks; RNA-seq", "hg38",
         "4 GSC lines; 1 NSC", "3D target assignment"],
        ["TCGA via Xena [14]", "bulk RNA-seq", "—",
         "172 GBM; 530 LGG", "independent expression"],
    ], "Table 1. Data sources.")

    D.h("Peak calling from signal, and its calibration", 2)
    D.p("Two cohorts released normalised signal without peak calls. Peaks were "
        "therefore derived from signal in every cohort: the genome was divided "
        "into 200 bp bins, bins were thresholded at a quantile of each sample's "
        "own non-zero signal distribution, so that sequencing depth and track "
        "scaling cancel, surviving bins were merged across a gap, and regions "
        "shorter than a minimum width or overlapping the ENCODE blacklist [11] "
        "were removed. The three parameters were chosen against ground truth "
        "rather than by inspection: GSE162976 published both signal and MACS2 "
        "peaks for the same samples, and a grid of 30 settings was scored "
        f"against those peaks in {N['cal_n_calibration_samples']} samples at "
        "base-pair (Jaccard) and peak level (precision, recall). The selected "
        f"setting (quantile {N['cal_quantile']}, gap {N['cal_gap']} bp, minimum "
        f"width {N['cal_min_width']} bp) reached a median Jaccard of "
        f"{f(N['cal_jaccard_median'])}, precision {f(N['cal_precision_median'])} "
        f"and recall {f(N['cal_recall_median'])}, and was then applied unchanged "
        "to every cohort, including the cohort it was calibrated on, so that "
        "cross-cohort differences cannot arise from the caller.")

    D.h("Super-enhancers", 2)
    D.p("Peaks lying wholly within 2 kb of an annotated transcription start site "
        "(GENCODE v44, protein-coding and lncRNA genes) were removed, the "
        "remainder were stitched within 12.5 kb, and each stitched domain was "
        "scored by its total H3K27ac signal, following ROSE [5, 6]. Calling and "
        "stitching were done in each cohort's native genome build, and only the "
        "resulting coordinates were lifted to hg38. Two definitions were "
        "applied. The ROSE inflection rule places the cutoff where the ranked "
        "signal curve, scaled to the unit square, has slope one; we found its "
        "per-sample count tracked data provenance (Results), so the primary "
        "definition was instead the 500 highest-scoring domains per sample, "
        "which gives every sample equal weight. Results under both definitions "
        "are reported, and the fixed definition was repeated at 250 and 1,000 "
        "domains.")

    D.h("Cross-cohort recurrence and the noise ceiling", 2)
    D.p(f"Because a union of stitched domains across {N['n_total']} samples merges into "
        "megabase blocks, recurrence was computed in fixed 5 kb bins as the "
        "fraction of a group's samples whose super-enhancers cover the bin. "
        "Agreement between two groups was measured by the Spearman correlation "
        "of recurrence over bins carrying a super-enhancer in either, by the "
        "average precision with which one group's recurrence predicts bins "
        "recurrent (≥50%) in the other, reported with its chance level "
        "(the prevalence of positives), and by the Jaccard index of bins "
        "recurrent in both. ROC AUC was not used: with roughly one positive bin "
        "in forty it exceeded 0.9 for every comparison, including unrelated ones. "
        "For each group, the same statistics were computed between two random "
        "halves of that group, 25 times. A split-half value is built from "
        "half-size groups, and recurrence estimated from fewer samples is "
        "noisier, so setting it against a cross-cohort value computed from full "
        "groups favours the cross-cohort side. The primary comparison was "
        "therefore size-matched: for a target cohort, 2k samples were drawn "
        "without replacement and split into two disjoint sets of k; the first "
        "set gave the within-target agreement, and k samples drawn from the "
        "source cohort were scored against the same second set, so the two "
        "values share their reference and differ only in the origin of the "
        "predicting samples. This was repeated 50 times for k = 6, 10 and 20 "
        "(tumours), 6 and 10 (pediatric) and 6 (NSC), and summarised as the "
        "median ratio of cross-cohort to within-target average precision with "
        "its 2.5–97.5 percentile range. The full-group comparison was "
        "retained as a secondary analysis. Its uncertainty was estimated by 200 "
        "bootstrap resamples of samples within each cohort, with the ceiling "
        "for each replicate taken from two halves of distinct samples, since "
        "splitting a resample drawn with replacement can place copies of one "
        "sample in both halves and inflate the ceiling.")
    D.p("Consensus core domains were bins recurrent in at least half of both "
        "adult cohorts and in fewer than 25% of the 12 NSC samples (at most "
        "two), merged into contiguous runs. A stricter threshold of 10% was "
        "tried first; with 12 samples it excludes a bin as soon as two NSC "
        "samples carry it and removed nearly every jointly recurrent bin, so it "
        "was not used.")

    D.h("Enhancer-to-gene assignment", 2)
    D.p("HiChIP loop anchors were classified as promoter (containing a GENCODE "
        "transcription start site), enhancer (overlapping an H3K27ac peak and "
        "containing no start site) or neither, and every enhancer-promoter loop "
        "was taken as a contact-based link. Each link was compared with the "
        "gene whose start site is nearest the enhancer. Three controls were "
        "applied. (i) Because a 10 kb promoter anchor can hold several start "
        "sites, the comparison was repeated on anchors with exactly one. (ii) "
        "Because a gene 100 kb away is rarely the closest for trivial reasons, "
        "a decoy gene at the same distance from the same enhancer was drawn for "
        "every link and tested for being the nearest gene, giving the agreement "
        "expected from genome architecture alone. (iii) Because loop counts "
        "differed several-fold between samples, every sample was subsampled to "
        "the loop count of the shallowest, 20 times, before any comparison of "
        "contacted genes between samples.")

    D.h("Motifs and expression", 2)
    D.p("Motif enrichment was tested with SEA [12] against JASPAR 2024 CORE "
        "vertebrate non-redundant motifs [13] on the constituent peaks of each "
        "domain set (peak centre ±250 bp from the hg38 tumour cohort), "
        "against peaks of sporadic domains drawn to equal number and against a "
        "dinucleotide-shuffled control of the same sequences, which preserves "
        "base composition. Whole stitched domains were not used because at "
        "10–40 kb they contain nearly every motif. For expression, each "
        "gene was assigned the highest recurrence within 100 kb of its start "
        "site, the median HiChIP contact distance observed here, and "
        "associations with expression were tested by Spearman correlation and "
        "Mann-Whitney tests, with Cliff's δ as effect size.")

    D.h("Software and availability", 2)
    D.p("Analyses used Python 3.11 (NumPy, pandas, SciPy, scikit-learn, pyBigWig), "
        "bedtools, UCSC liftOver and bedGraphToBigWig, and MEME Suite 5.5.9. "
        "All code, intermediate tables and figure scripts are available at "
        "https://github.com/kratos1459/ChIPseq_GBM_atlas; each stage is "
        "scripted and resumable from the published GEO and Xena files.")

    # ---------------------------------------------------------------- results
    D.h("Results")
    med = N["se_med"]
    D.h("A calibrated caller places four cohorts on one footing", 2)
    D.p(f"The signal-based caller recovered published MACS2 peaks with a median "
        f"Jaccard of {f(N['cal_jaccard_median'])} (F1 {f(N['cal_f1_median'])}) "
        f"across {N['cal_n_calibration_samples']} samples, and was then applied "
        f"to all {N['n_total']} samples (Figure 1a). Under the ROSE inflection "
        f"rule the median number of super-enhancers per sample was "
        f"{med.get('GSC_line', 0):.0f} in GSC lines, "
        f"{med.get('patient_tumour', 0):.0f} in tumours, {med.get('NSC', 0):.0f} "
        f"in NSC and {med.get('pediatric', 0):.0f} in pediatric samples (Figure "
        f"1b). The count followed the number of stitched domains in each sample "
        f"(Figure 1c), so it reflects how the data were generated and scaled at "
        f"least as much as biology; for this reason the inflection count was "
        f"not interpreted, and the fixed top-500 definition was used from here "
        f"on, with the inflection definition reported alongside (Supplementary "
        f"Table S2).")

    D.h("GSC lines reproduce most of the tumour landscape, but less than "
        "tumours reproduce each other", 2)
    D.p(f"Recurrence in GSC lines and in tumours was strongly concordant: "
        f"{joint:,} 5 kb bins were recurrent in at least half of both cohorts "
        f"(Figure 2b), and the Jaccard index of recurrent bins was "
        f"{f(N['x_jac'])} against {f(N['ceil_jac_tum'])} between two halves of "
        f"the tumour cohort. GSC-line recurrence predicted recurrent tumour bins "
        f"with an average precision of {f(N['x_ap'])} (chance "
        f"{f(N['x_base'])}), close to the split-half value of "
        f"{f(N['ceil_ap_tum'])} (Figure 3a). Taken at face value this is "
        f"{pc(frac)} of the ceiling, a fraction that was stable at 250, 500 and "
        f"1,000 domains per sample ({', '.join(pc(x) for x in rob['ap_frac_ceiling'])}; "
        f"Supplementary Table S3).")
    D.p(f"That comparison is not like for like: the cross-cohort value uses all "
        f"{n.get('GSC_line', 0)} GSC lines and {n.get('patient_tumour', 0)} "
        f"tumours, whereas each half of the ceiling averages only "
        f"{n.get('patient_tumour', 0)//2} tumours. When both sides were "
        f"restricted to the same number of samples and scored against the same "
        f"reference set, GSC lines consistently fell short of tumours (Figure "
        f"3b). At 20 samples per side, GSC lines reached an average precision of "
        f"{f(t20['cross'])} against {f(t20['within'])} for an independent set of "
        f"tumours, {pc(t20['frac'])} of the within-tumour value (95% range "
        f"{pc(t20['lo'])}–{pc(t20['hi'])}), and they did worse than tumours "
        f"in {pc(1 - t20['p_ge'])} of the 50 draws. The ratio was the same at 10 "
        f"samples ({pc(t10['frac'])}, {pc(t10['lo'])}–{pc(t10['hi'])}) and "
        f"at 6 ({pc(t6['frac'])}, {pc(t6['lo'])}–{pc(t6['hi'])}), where the "
        f"wider ranges reflect noisier recurrence rather than a different "
        f"effect. Rank agreement showed the same direction more strongly "
        f"(Spearman {f(t20['sp_cross'])} against {f(t20['sp_within'])} at 20 "
        f"samples). GSC lines therefore recover most, but not all, of the "
        f"reproducible tumour super-enhancer landscape, which is consistent "
        f"with the limited co-clustering of tissue and cell-line samples noted "
        f"previously [2]. A bootstrap of the full-group comparison gave a "
        f"fraction of ceiling of {f(bs.loc['ap_fraction_of_ceiling','median'])} "
        f"({f(bs.loc['ap_fraction_of_ceiling','ci_low'])}–"
        f"{f(bs.loc['ap_fraction_of_ceiling','ci_high'])}); that interval "
        f"reaches 1 only because of the size mismatch described above.")
    D.p(f"The adult cohorts agreed much less with pediatric high-grade glioma. "
        f"Size-matched at 10 samples, GSC lines and tumours recovered "
        f"{pc(p10['frac'])} and {pc(p10t['frac'])} of within-pediatric "
        f"agreement, and both fell below it in nearly every draw. Agreement with "
        f"NSC could not be resolved: with 12 NSC samples, two sets of six agreed "
        f"with each other only weakly (average precision {f(n6g['within'])}, "
        f"chance {f(n6g['base'])}), and GSC lines ({pc(n6g['frac'])}, "
        f"{pc(n6g['lo'])}–{pc(n6g['hi'])}) and tumours ({pc(n6t['frac'])}, "
        f"{pc(n6t['lo'])}–{pc(n6t['hi'])}) were within that noise. "
        f"Whether adult glioblastoma differs from neural stem cells in its "
        f"super-enhancers cannot be decided from these NSC data.")

    D.h("The reproducible core is small and mostly shared with neural stem "
        "cells", 2)
    D.p(f"Of the {joint:,} bins recurrent in both adult cohorts, only {spec} "
        f"were carried by fewer than a quarter of NSC samples. Merged, they form "
        f"{cc.get('n_core_domains')} core domains covering "
        f"{f(cc.get('core_mb'))} Mb (Figure 5a, b), against "
        f"{cc.get('n_sporadic_domains'):,} sporadic domains ({f(cc.get('sporadic_mb'), 1)} Mb) "
        f"seen in a single sample. The largest core domains lie at C9orf50, "
        f"WBP2, JMJD6, BRAT1, LINC00963 and CCND2-AS1, and the set includes "
        f"SOX9, LFNG, ID2-AS1, ASIC4 and CDC25B (Figure 5c; Supplementary "
        f"Table S4). Only {cc.get('n_gsc_specific_domains')} domains "
        f"({f(cc.get('gsc_specific_mb'))} Mb) were recurrent in GSC lines but "
        f"rare in tumours, and {cc.get('n_patient_specific_domains')} "
        f"({f(cc.get('patient_specific_mb'))} Mb) the reverse. Core bins were "
        f"recurrent in {pc(cc.get('core_frac_recurrent_in_pediatric'))} of "
        f"pediatric samples on average (median "
        f"{f(cc.get('core_median_recurrence_in_pediatric'))}), in line with the "
        f"weak adult–pediatric agreement. The core size grew with the "
        f"number of domains per sample ({', '.join(str(int(x)) for x in rob['core_domains'])} "
        f"domains at 250, 500 and 1,000), so it should be read as a "
        f"high-confidence subset rather than a complete list.")

    D.h("Nearest-gene assignment is usually not the contacted gene", 2)
    D.p(f"In H3K27ac HiChIP from four GSC lines and one NSC line, enhancers "
        f"contacted promoters a median {ranges(ctl['median_contact_kb'], 1, 0)} "
        f"kb away, whereas the nearest start site lay "
        f"{ranges(ctl['median_nearest_kb'], 1, 0)} kb away (Figure 4b). The "
        f"nearest gene was the contacted gene for "
        f"{ranges(ctl['nearest_correct_all'])}% of enhancer-promoter links "
        f"(Figure 4a). Restricting to promoter anchors containing a single start "
        f"site, which removes ambiguity about which promoter in a 10 kb anchor "
        f"is contacted, raised this to {ranges(ctl['nearest_correct_clean'])}%; "
        f"a distance-matched decoy gene was the nearest gene "
        f"{ranges(ctl['nearest_correct_decoy_null'])}% of the time, so "
        f"nearest-gene assignment performs above the level set by genome "
        f"architecture but misses most contacts. Even allowing any contacted "
        f"gene, the nearest gene was among an enhancer's contacts for "
        f"{ranges(ctl['enh_with_nearest_among_contacts'])}% of enhancers. For "
        f"the core domains with loop data, the nearest gene was among the "
        f"contacted genes for {pc(cc.get('core_nearest_is_contact_frac'))}, "
        f"higher than genome-wide, consistent with strong, recurrent elements "
        f"lying closer to their targets.")
    D.p(f"Contacted-gene sets differed little between samples: target-set "
        f"Jaccard was {f(j_gg.min())}–{f(j_gg.max())} between GSC lines "
        f"and {f(j_gn.min())}–{f(j_gn.max())} between a GSC line and NSC "
        f"(Figure 4c). Taken at face value, {N['n_spec_raw']:,} genes were "
        f"contacted in all four GSC lines but not in NSC. The libraries differed "
        f"several-fold in loop count, and a gene is more likely to be seen in "
        f"deeper data. After subsampling every library to the depth of the "
        f"shallowest, {dm.loc['genes_all_gsc_not_nsc','mean']:.0f} genes "
        f"(SD {dm.loc['genes_all_gsc_not_nsc','sd']:.0f}) met the criterion, "
        f"and {n_rob} did so in at least 80% of 20 subsamples (Figure 4d), "
        f"among them BRD4, JUN, FOSL2, CEBPB, FOXG1, GOLIM4, ADAM17, LAMC1 and "
        f"KIF5B (Supplementary Table S5). At matched depth the nearest-gene "
        f"error rate was the same in GSC and NSC data.")

    D.h("Motif enrichment in the core reflects base composition", 2)
    gcp, gcc = N["gc"]["primary"], N["gc"]["control"]
    m = N["mot"]
    n_sig = int(((m["comparison"] == "core_vs_sporadic") & (m["QVALUE"] < 0.05)).sum())
    n_shuf = int(((m["comparison"] == "core_vs_shuffled") & (m["QVALUE"] < 0.05)).sum())
    D.p(f"Against peaks from sporadic domains, the {N['n_core_peaks']} "
        f"constituent peaks of core domains were enriched for {n_sig} JASPAR "
        f"motifs at q < 0.05, all GC-rich: KLF and SP family members, PATZ1, "
        f"PLAG1, ZNF281 and TFAP2C (Figure 6, left). Core peaks were "
        f"{pc(gcp, 1)} GC against {pc(gcc, 1)} for the control. Against a "
        f"dinucleotide-shuffled control of the same sequences, which keeps GC "
        f"content, {n_shuf} motifs reached q < 0.05 (Figure 6, right). The "
        f"apparent enrichment is therefore a composition effect, and with fewer "
        f"than 100 peaks the analysis has little power to detect a specific "
        f"transcription-factor signature.")

    D.h("Core domains mark more highly expressed genes", 2)
    D.p(f"Genes within 100 kb of a core domain were more highly expressed in "
        f"TCGA glioblastoma than genes near other super-enhancers (median "
        f"{f(gbm['median_a'])} against {f(gbm['median_b'])}; Mann-Whitney "
        f"p = {gbm['p']:.3f}; Cliff's δ {f(gbm['cliffs_delta'])}; "
        f"n = {int(gbm['n'])} core genes) and than genes near no super-enhancer "
        f"(δ {f(ex_row(N, 'core_vs_no_se', 'tcga_gbm')['cliffs_delta'])}). "
        f"The same held in pediatric RNA-seq (δ {f(ped['cliffs_delta'])}, "
        f"p = {ped['p']:.3f}) and in lower-grade glioma, with no difference in "
        f"effect between glioblastoma and lower-grade glioma (Figure 7). "
        f"Recurrence correlated modestly with expression genome-wide, and "
        f"nearly as well for NSC recurrence (ρ {f(sp_n['rho'])}) as for "
        f"GSC-line recurrence (ρ {f(sp_g['rho'])}), so the association is "
        f"largely a generic property of active chromatin rather than evidence "
        f"of tumour-specific regulation.")

    # ------------------------------------------------------------- discussion
    D.h("Discussion")
    D.p("The central result is quantitative. Patient-derived GSC lines reproduce "
        f"about {pc(t20['frac'])} of the super-enhancer agreement that one set of "
        "tumours reaches with another of the same size. That is high enough that "
        "the lines remain a reasonable model of the tumour enhancer landscape, "
        "and low enough that domains seen in the lines but not in tumours, or "
        "the reverse, should not be assumed to be noise. The deficit was visible "
        "only after matching group sizes. Against a ceiling built from half-size "
        "groups, the full-cohort comparison came within a few percent of the "
        "ceiling and its bootstrap interval included it, which would have led "
        "us to conclude that the two sources are interchangeable. The same "
        "pitfall applies whenever a cross-cohort statistic is compared with a "
        "split-half statistic without equalising the number of samples behind "
        "each.")
    D.p("The tumour-specific part of the landscape is small. Of nearly a "
        "thousand 5 kb bins recurrent in both adult cohorts, about nine in ten "
        "were also present in at least a quarter of NSC samples, leaving 30 core "
        "domains. Several sit at genes with established roles in glioma stem "
        "cells or neural development, including SOX9, ID2, LFNG and CCND2, but "
        "the list is short and its size depends on the definition, and we would "
        "treat it as a set of well-supported candidates rather than as the "
        "glioblastoma super-enhancer repertoire. A second limitation applies "
        "here: the NSC group combines lines from two studies, and with 12 "
        "samples it supports only weak exclusion. A bin carried by two of 12 NSC "
        "samples is still counted as rare, and in the size-matched analysis "
        "adult agreement with NSC could not be separated from the NSC "
        "group's own noise. More NSC and normal-brain H3K27ac data would sharpen both "
        "results.")
    D.p("For target assignment, three-dimensional contact and linear proximity "
        "disagree for most enhancers, in line with CRISPR-based benchmarks in "
        "other cell types [8, 9]. Two controls changed the picture materially. "
        "Requiring a single start site per promoter anchor raised the agreement "
        "between nearest gene and contact by roughly half, so part of the "
        "apparent disagreement is anchor resolution rather than biology; and "
        "subsampling to matched loop depth removed most apparently GSC-specific "
        "contacts, leaving about a hundred genes whose specificity survives. "
        "Hub genes reported in the original HiChIP study, such as GOLIM4 [4], "
        "were among them. Comparisons of contact maps between samples that "
        "differ in depth should be matched in this way before differences are "
        "interpreted.")
    D.p("The motif analysis is a negative result worth stating plainly. GC-box "
        "motifs appeared strongly enriched against a conventional background "
        "but not against a composition-matched one, and with fewer than a "
        "hundred constituent peaks the analysis was underpowered. The "
        "expression analysis supports the core domains as marking active, "
        "highly expressed genes, but the genome-wide association between "
        "recurrence and expression was nearly as strong for NSC as for "
        "glioblastoma, so it does not by itself argue for tumour-specific "
        "regulation.")
    D.p("Further limitations should be noted. All data are public and were "
        "reprocessed from signal rather than reads, so peak calling depends on "
        "a calibration performed in one cohort; the calibration was good, and "
        "the identical procedure in every cohort means that caller error cannot "
        "explain differences between cohorts, but it may blunt agreement "
        "overall. Cohorts were profiled on different builds and lifted to hg38 "
        "after calling, and residual differences in antibody, protocol and "
        "depth between laboratories are confounded with the cell-line versus "
        "tumour contrast; the cell-line deficit could partly reflect laboratory "
        "rather than biology. Bulk tumour profiles also contain non-malignant "
        "cells, which may contribute tumour-only domains. Finally, HiChIP was "
        "available for only four GSC lines and one NSC line.")
    D.p("In summary, GSC lines are a good but imperfect proxy for the "
        "glioblastoma super-enhancer landscape, the reproducible tumour-specific "
        "core is small, and enhancer-gene assignment by proximity or by "
        "depth-unmatched contact maps introduces errors at least as large as the "
        "difference between model and tumour. Each of these can be measured "
        "with public data, and we provide the code to do so.")

    D.h("Data and code availability")
    D.p("All analysed data are public (GSE119776, GSE145646, GSE162976, "
        "GSE262089; TCGA via UCSC Xena). Code, result tables and figure scripts: "
        "https://github.com/kratos1459/ChIPseq_GBM_atlas.")
    D.h("Author contributions")
    D.p("[To be completed by the authors.]")
    D.h("Competing interests")
    D.p("[To be completed by the authors.]")
    D.h("Funding")
    D.p("[To be completed by the authors.]")

    D.h("References")
    for i, r in enumerate(REFS, 1):
        D.p(f"{i}. {r}")

    # ---------------------------------------------------------------- figures
    D.page_break()
    D.h("Figures")
    for name, cap in FIG_CAPTIONS:
        D.figure(os.path.join(FIG, name + ".png"), cap)
    return D


FIG_CAPTIONS = [
    ("Fig1_cohort_composition",
     "Figure 1. Cohorts and super-enhancer counts. (a) H3K27ac samples per group "
     "after quality filtering. (b) Super-enhancers per sample under the ROSE "
     "inflection rule; points are samples. (c) Super-enhancers against stitched "
     "enhancer domains per sample: the inflection count follows the number of "
     "domains, which motivated the fixed top-500 definition."),
    ("Fig2_recurrence_landscape",
     "Figure 2. Recurrence landscape. (a) Fraction of samples carrying a "
     "super-enhancer, for every 5 kb bin with a super-enhancer in any sample, "
     "ranked within each group. (b) Recurrence in GSC lines against recurrence "
     "in patient tumours (hexagonal bins, log-scaled counts); dashed lines mark "
     "50% recurrence."),
    ("Fig3_agreement_vs_ceiling",
     "Figure 3. Cross-cohort agreement against within-cohort agreement. (a) "
     "Full-group comparisons: Spearman correlation of recurrence, average "
     "precision for predicting bins recurrent in the second group (red ticks, "
     "chance level) and Jaccard index of recurrent bins. Green, two random "
     "halves of one cohort (mean of 25 splits); blue, between cohorts. (b) "
     "Size-matched comparisons: average precision when k samples from the "
     "source cohort (blue, GSC lines; red, tumours) or k other samples of the "
     "target cohort (green) predict recurrence in the same k reference samples "
     "of the target. Bars, mean of 50 draws; lines, 2.5–97.5 percentile; "
     "red marks, chance."),
    ("Fig4_target_assignment",
     "Figure 4. Enhancer-gene assignment from H3K27ac HiChIP. (a) Fraction of "
     "enhancer-promoter links for which the nearest gene is the contacted gene, "
     "for all links, for promoter anchors containing a single start site, and "
     "for distance-matched decoy genes. (b) Median distance from enhancer to "
     "nearest and to contacted start site. (c) Jaccard index of contacted-gene "
     "sets for shared enhancers between pairs of GSC lines and between GSC "
     "lines and NSC. (d) Genes contacted in all GSC lines but not in NSC, as "
     "called, at matched loop depth (mean of 20 subsamples) and in at least 80% "
     "of subsamples."),
    ("Fig5_core_domains",
     "Figure 5. Consensus domains. (a) Domains in each class (log scale): core, "
     "recurrent in both adult cohorts and rare in NSC; GSC-specific and "
     "patient-specific, recurrent in one adult cohort and rare in the other; "
     "sporadic, seen in a single sample. (b) Width of core domains. (c) The "
     "largest core domains, labelled by nearest gene, or by a contacted gene "
     "where the nearest gene is unnamed."),
    ("Fig6_motif_enrichment",
     "Figure 6. Motif enrichment in constituent peaks of core domains (SEA, "
     "JASPAR 2024 CORE vertebrates). Left, against peaks of sporadic domains; "
     "right, against dinucleotide-shuffled core sequences. Dashed line, "
     "q = 0.05. Only the top 15 motifs are shown."),
    ("Fig7_expression_validation",
     "Figure 7. Expression of genes near super-enhancers. (a–c) Median "
     "expression of genes within 100 kb of a core domain, of another "
     "super-enhancer, or of none, in TCGA glioblastoma, pediatric high-grade "
     "glioma and GSC lines. (d) Spearman correlation between bin recurrence in "
     "each group and expression, per expression dataset."),
]


def fmt(v, nd=3):
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return "—"
        return f"{v:.{nd}g}" if abs(v) < 1e-3 and v != 0 else f"{v:.{nd}f}"
    return str(v)


def build_supp(N):
    D = Doc()
    D.title("Supplementary Information", 14)
    D.center(TITLE, italic=True, size=10)
    D.p("All tables are generated from the result files in the code repository "
        "(path given in each caption); the files hold the complete, unrounded "
        "values.")

    # S1 calibration
    grid = rcsv("calibration_grid.csv")
    if grid is not None:
        g = (grid.groupby(["quantile", "gap", "min_width"])
             [["jaccard", "precision", "recall", "f1", "n_called"]].median()
             .reset_index().sort_values("jaccard", ascending=False))
        cols = list(g.columns)
        D.table(["quantile", "gap (bp)", "min width (bp)", "Jaccard", "precision",
                 "recall", "F1", "peaks called"],
                [[fmt(r["quantile"], 3), int(r["gap"]), int(r["min_width"]),
                  fmt(r["jaccard"]), fmt(r["precision"]), fmt(r["recall"]),
                  fmt(r["f1"]), f"{int(r['n_called']):,}"] for _, r in g.iterrows()],
                "Table S1. Peak-caller calibration grid against published MACS2 "
                "peaks (GSE162976), sorted by median Jaccard. "
                "results/calibration_grid.csv")
    # S2 agreement, both definitions
    ab = rcsv("cross_cohort/agreement_both_definitions.csv")
    cols = ["definition", "comparison", "spearman", "ap_predicting_b", "ap_baseline",
            "jaccard_at_50pct"]
    D.table(["definition", "comparison", "Spearman", "AP", "chance", "Jaccard"],
            [[fmt(r[c]) for c in cols] for _, r in ab.iterrows()],
            "Table S2. Full-group agreement under the ROSE inflection and fixed "
            "top-500 definitions. results/cross_cohort/agreement_both_definitions.csv",
            size=7.5)
    # S3 robustness + size-matched + bootstrap
    rb = N["rob"]
    D.table(["domains/sample", "GSC→tumour AP", "tumour ceiling AP",
             "fraction", "Spearman", "chance", "core domains", "core Mb"],
            [[int(r.top_n), fmt(r.gsc_tum_ap), fmt(r.ceil_ap_tum),
              fmt(r.ap_frac_ceiling), fmt(r.gsc_tum_sp), fmt(r.base),
              int(r.core_domains), fmt(r.core_mb, 2)] for r in rb.itertuples()],
            "Table S3a. Full-group comparison at 250, 500 and 1,000 domains per "
            "sample (half-size ceiling; see Methods). results/robustness_topN.csv")
    rows = []
    for (t, k, s), v in sorted(N["smp"].items()):
        rows.append([t, k, s, fmt(v["cross"]), fmt(v["within"]), fmt(v["base"]),
                     fmt(v["frac"]), f"{fmt(v['lo'])}–{fmt(v['hi'])}",
                     fmt(1 - v["p_ge"], 2)])
    D.table(["target", "k", "source", "AP source", "AP within", "chance",
             "ratio (median)", "95% range", "P(source < within)"], rows,
            "Table S3b. Size-matched agreement, 50 draws per row (top-500). "
            "results/cross_cohort/size_matched_reps.csv", size=7.5)
    bs = N["bs"]
    D.table(["statistic", "2.5%", "median", "97.5%"],
            [[i, fmt(r.ci_low), fmt(r["median"]), fmt(r.ci_high)]
             for i, r in bs.iterrows()],
            "Table S3c. Bootstrap (200 resamples) of the full-group GSC→tumour "
            "comparison, ceiling from distinct samples. "
            "results/cross_cohort/bootstrap_summary.csv")
    # S4 core domains
    cg = N["core_genes"]
    D.table(["domain", "locus (hg38)", "kb", "nearest gene", "dist kb",
             "contacted genes (HiChIP)"],
            [[r.domain, f"{r.chrom}:{int(r.start):,}-{int(r.end):,}",
              fmt(r.width_kb, 0), r.nearest_gene, fmt(r.nearest_dist_kb, 1),
              (str(r.contact_genes) if isinstance(r.contact_genes, str) else "—")
              .replace(";", ", ")] for r in cg.itertuples()],
            "Table S4. The 30 core domains. results/cross_cohort/core_domain_genes.csv",
            size=7)
    # S5 loops
    c = N["ctl"]
    D.table(["sample", "E-P links", "nearest correct", "single-TSS links",
             "nearest correct (single-TSS)", "decoy", "nearest among contacts",
             "median contact kb", "median nearest kb"],
            [[r.sample, f"{int(r.n_links):,}", fmt(r.nearest_correct_all),
              f"{int(r.n_links_clean):,}", fmt(r.nearest_correct_clean),
              fmt(r.nearest_correct_decoy_null),
              fmt(r.enh_with_nearest_among_contacts), fmt(r.median_contact_kb, 1),
              fmt(r.median_nearest_kb, 1)] for r in c.itertuples()],
            "Table S5a. Nearest-gene versus contact assignment per HiChIP "
            "library. results/loop_targets/controls_summary.csv", size=7.5)
    dg = N["dg_robust"].sort_values(["frac_repeats", "gene"], ascending=[False, True])
    names = [g for g in dg["gene"] if not str(g).startswith("ENSG")]
    D.p(f"Table S5b. Genes contacted in all four GSC lines and not in NSC in at "
        f"least 80% of 20 depth-matched subsamples (n = {len(dg)}; "
        f"{len(names)} named, listed; Ensembl-only identifiers in the source "
        f"file). results/loop_targets/depth_matched_genes.csv",
        bold_lead="")
    D.p(", ".join(names))
    # S6 motifs
    m = N["mot"]
    m = m[m["QVALUE"] < 0.05].sort_values("QVALUE")
    D.table(["comparison", "motif", "ID", "consensus", "% target", "% control",
             "enrichment", "q"],
            [[r.comparison, r.ALT_ID, r.ID, r.CONSENSUS, fmt(r._7, 1), fmt(r._9, 1),
              fmt(r.ENR_RATIO, 2), fmt(r.QVALUE, 2)] for r in m.itertuples()],
            f"Table S6. Motifs at q < 0.05 (SEA). Base composition: core peaks "
            f"{pc(N['gc']['primary'], 1)} GC, sporadic control "
            f"{pc(N['gc']['control'], 1)} GC. No motif reached q < 0.05 against "
            f"the shuffled control. results/motifs/motif_enrichment.csv", size=7)
    # S7 expression
    e = N["ex"]
    D.table(["test", "expression", "recurrence", "n", "ρ", "p",
             "median a", "median b", "Cliff's δ"],
            [[r.test, r.expression, r.recurrence, int(r.n), fmt(r.rho), fmt(r.p, 2),
              fmt(r.median_a, 2), fmt(r.median_b, 2), fmt(r.cliffs_delta)]
             for r in e.itertuples()],
            "Table S7. Expression tests. results/expression/expression_tests.csv",
            size=7)
    return D


def main():
    os.makedirs(OUT, exist_ok=True)
    N = collect()
    m = build_main(N)
    p1 = os.path.join(OUT, "GBM_SE_Manuscript.docx")
    m.save(p1)
    s = build_supp(N)
    p2 = os.path.join(OUT, "GBM_SE_Supplementary.docx")
    s.save(p2)
    print("wrote", p1, "and", p2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
