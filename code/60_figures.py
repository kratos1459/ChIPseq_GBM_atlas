#!/usr/bin/env python3
"""Publication figures for the GBM enhancer-reproducibility atlas.

Every panel is written at 600 dpi as both PNG and LZW-compressed TIFF. Each
figure is independent and wrapped so that a missing input skips that figure
rather than killing the run.

  Fig 1  cohort composition and per-sample super-enhancer counts
  Fig 2  recurrence distributions and the GSC-versus-tumour recurrence plane
  Fig 3  cross-cohort agreement against split-half noise ceilings
  Fig 4  nearest-gene versus 3D-contact assignment, with the decoy null
  Fig 5  core domain catalogue and its pediatric generalisation
  Fig 6  motif enrichment, reproducible versus sporadic domains
  Fig 7  super-enhancer recurrence against independent expression data
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
CC = os.path.join(RES, "cross_cohort")
LT = os.path.join(RES, "loop_targets")
FIG = os.path.join(ROOT, "figures")
PNG = os.path.join(FIG, "PNG_600dpi")
TIF = os.path.join(FIG, "TIFF_600dpi")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelsize": 8, "axes.titlesize": 9, "axes.linewidth": 0.6,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
GREY, BLUE, RED, GREEN, ORANGE = "#555555", "#1f4fd8", "#d62728", "#2e8b2e", "#e6a800"
PALETTE = {"GSC_line": BLUE, "patient_tumour": RED, "NSC": GREY,
           "pediatric": GREEN, "other": ORANGE}


def save(fig, name):
    os.makedirs(PNG, exist_ok=True)
    os.makedirs(TIF, exist_ok=True)
    fig.savefig(os.path.join(PNG, f"{name}.png"), dpi=600)
    fig.savefig(os.path.join(TIF, f"{name}.tif"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw", "dpi": (600, 600)})
    plt.close(fig)
    print(f"  wrote {name}.png / .tif")


def letter(ax, s):
    ax.text(-0.14, 1.06, s, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left")


def fig1():
    p = os.path.join(CC, "samples_used.csv")
    if not os.path.exists(p):
        return print("  fig1: no samples_used.csv")
    d = pd.read_csv(p)
    order = [g for g in ("GSC_line", "patient_tumour", "NSC", "pediatric")
             if g in set(d["group"])]
    names = {"GSC_line": "GSC lines", "patient_tumour": "patient tumours",
             "NSC": "neural stem cells", "pediatric": "pediatric HGG"}
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    n = d.groupby("group").size().reindex(order)
    axes[0].barh(range(len(n)), n.values, color=[PALETTE.get(g, ORANGE) for g in order],
                 height=0.6)
    for i, v in enumerate(n.values):
        axes[0].text(v + 1, i, str(int(v)), va="center", fontsize=6)
    axes[0].set_yticks(range(len(n)))
    axes[0].set_yticklabels([names[g] for g in order])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("H3K27ac samples")
    letter(axes[0], "a")

    data = [d.loc[d["group"] == g, "n_se"].values for g in order]
    bp = axes[1].boxplot(data, widths=0.55, patch_artist=True, showfliers=False)
    for patch, g in zip(bp["boxes"], order):
        patch.set_facecolor(PALETTE.get(g, ORANGE))
        patch.set_alpha(0.45)
        patch.set_linewidth(0.6)
    rng = np.random.default_rng(1)
    for i, vals in enumerate(data, 1):
        axes[1].plot(rng.normal(i, 0.055, len(vals)), vals, ".", ms=2, color="k", alpha=0.45)
    axes[1].set_xticks(range(1, len(order) + 1))
    axes[1].set_xticklabels([names[g].split()[0] for g in order], rotation=20)
    axes[1].set_ylabel("super-enhancers (ROSE inflection)")
    letter(axes[1], "b")

    for g in order:
        s = d[d["group"] == g]
        axes[2].scatter(s["n_stitched"], s["n_se"], s=6, alpha=0.7, linewidths=0,
                        color=PALETTE.get(g, ORANGE), label=names[g])
    axes[2].set_xlabel("stitched enhancer domains")
    axes[2].set_ylabel("super-enhancers")
    axes[2].legend(frameon=False, fontsize=5.5, markerscale=2)
    letter(axes[2], "c")
    fig.tight_layout()
    save(fig, "Fig1_cohort_composition")


def fig2():
    p = os.path.join(CC, "recurrence.csv.gz")
    if not os.path.exists(p):
        return print("  fig2: no recurrence table")
    r = pd.read_csv(p)
    cols = ([c for c in r.columns if c.startswith("top_rec_")]
            or [c for c in r.columns if c.startswith("se_rec_")])
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for c in cols:
        g = c.split("rec_", 1)[1]
        v = np.sort(r[c].to_numpy())[::-1]
        axes[0].plot(np.arange(1, len(v) + 1) / 1e3, v, lw=1.1,
                     color=PALETTE.get(g, ORANGE), label=g.replace("_", " "))
    axes[0].set_xlabel("super-enhancer bins (thousands, ranked)")
    axes[0].set_ylabel("fraction of samples")
    axes[0].legend(frameon=False)
    letter(axes[0], "a")

    if ("top_rec_GSC_line" in r.columns or "se_rec_GSC_line" in r.columns):
        p = "top_rec_" if "top_rec_GSC_line" in r.columns else "se_rec_"
        x, y = r[p + "GSC_line"], r[p + "patient_tumour"]
        axes[1].hexbin(x, y, gridsize=38, bins="log", cmap="Greys", linewidths=0)
        axes[1].axhline(0.5, color=RED, lw=0.6, ls="--")
        axes[1].axvline(0.5, color=RED, lw=0.6, ls="--")
        axes[1].set_xlabel("recurrence in GSC lines")
        axes[1].set_ylabel("recurrence in patient tumours")
        both = ((x >= .5) & (y >= .5)).sum()
        axes[1].set_title(f"{both:,} bins recurrent in both", fontsize=8)
    letter(axes[1], "b")
    fig.tight_layout()
    save(fig, "Fig2_recurrence_landscape")


def fig3():
    p = os.path.join(CC, "agreement_topN.csv")
    if not os.path.exists(p):
        p = os.path.join(CC, "agreement.csv")
    if not os.path.exists(p):
        return print("  fig3: no agreement table")
    a = pd.read_csv(p)
    order = ["GSC_line split-half (ceiling)", "patient_tumour split-half (ceiling)",
             "GSC_line -> patient_tumour", "patient_tumour -> GSC_line",
             "GSC_line -> NSC", "patient_tumour -> NSC",
             "GSC_line -> pediatric", "patient_tumour -> pediatric",
             "NSC split-half (ceiling)", "pediatric split-half (ceiling)"]
    a = a.set_index("comparison").reindex([o for o in order if o in set(a["comparison"])])
    labels = [s.replace("_", " ").replace(" split-half (ceiling)", ", split-half")
              .replace("patient tumour", "tumour").replace("GSC line", "GSC")
              .replace(" -> ", " → ") for s in a.index]
    metrics = [("spearman", "Spearman ρ", (-0.1, 0.5)),
               ("ap_predicting_b", "average precision", (0, 1)),
               ("jaccard_at_50pct", "Jaccard, ≥50% recurrence", (0, 1))]
    sm_path = os.path.join(CC, "size_matched_reps.csv")
    has_sm = os.path.exists(sm_path)
    fig = plt.figure(figsize=(7.2, 6.0 if has_sm else 3.3), layout="constrained")
    top, bot = (fig.subfigures(2, 1, height_ratios=[1.1, 1.0], hspace=0.04)
                if has_sm else (fig, None))
    axes = top.subplots(1, 3, sharey=True)
    y = np.arange(len(a))
    cols = [GREEN if "split-half" in s else BLUE for s in a.index]
    for ax, (col, lab, lim) in zip(axes, metrics):
        if col not in a.columns:
            ax.set_visible(False)
            continue
        v = a[col].to_numpy(dtype=float)
        ax.barh(y, np.nan_to_num(v), color=cols, height=0.62)
        for yi, vi in zip(y, v):
            if np.isfinite(vi):
                off = 0.01 * (lim[1] - lim[0])
                ax.text(max(vi, 0) + off, yi, f"{vi:.2f}",
                        va="center", ha="left", fontsize=5.5)
        if col == "ap_predicting_b" and "ap_baseline" in a.columns:
            base = a["ap_baseline"].to_numpy(dtype=float)
            ok = np.isfinite(base)
            ax.plot(base[ok], y[ok], "|", color=RED, ms=7, mew=1.2,
                    label="chance (prevalence)")
            ax.legend(frameon=False, fontsize=5.5, loc="lower right")
        ax.set_xlim(*lim)
        ax.set_xlabel(lab)
        ax.axvline(0, color="k", lw=0.4)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels, fontsize=6)
    axes[0].invert_yaxis()
    letter(axes[0], "a")
    top.text(0.01, 0.0, "green: two random halves of one cohort; blue: between "
             "cohorts, full groups. Fixed top-500 domains per sample.", fontsize=5.5)
    if has_sm:
        sm = pd.read_csv(sm_path)
        bx = bot.subplots(1, 3, gridspec_kw={"width_ratios": [3, 2, 1]})
        tnames = {"patient_tumour": "target: patient tumours",
                  "pediatric": "target: pediatric HGG", "NSC": "target: NSC"}
        snames = {"GSC_line": "GSC lines", "patient_tumour": "tumours"}
        for ax, tgt in zip(bx, ["patient_tumour", "pediatric", "NSC"]):
            s = sm[sm["target"] == tgt]
            ks = sorted(s["k"].unique())
            srcs = [x for x in s["source"].unique() if not x.endswith("(within)")]
            series = [f"{tgt} (within)"] + srcs
            w = 0.8 / len(series)
            for j, src in enumerate(series):
                col = GREEN if src.endswith("(within)") else PALETTE.get(src, ORANGE)
                lab = "within target" if src.endswith("(within)") else snames.get(src, src)
                for i, k in enumerate(ks):
                    v = s[(s["source"] == src) & (s["k"] == k)]["ap"].dropna()
                    xpos = i + (j - (len(series) - 1) / 2) * w
                    ax.bar(xpos, v.mean(), w * 0.9, color=col, alpha=0.8,
                           label=lab if i == 0 else None)
                    lo, hi = np.percentile(v, [2.5, 97.5])
                    ax.plot([xpos, xpos], [lo, hi], color="k", lw=0.6)
            base = s.groupby("k")["baseline"].mean().reindex(ks).to_numpy()
            ax.plot(np.arange(len(ks)), base, "_", color=RED, ms=14, mew=1.3,
                    label="chance")
            ax.set_xticks(range(len(ks)))
            ax.set_xticklabels([f"k = {k}" for k in ks])
            ax.set_ylim(0, 1)
            ax.set_title(tnames[tgt], fontsize=7.5)
            ax.legend(frameon=False, fontsize=5.5, loc="upper left")
        bx[0].set_ylabel("average precision\n(k samples vs k samples)")
        letter(bx[0], "b")
        bot.text(0.01, 0.0, "Both sides average exactly k distinct samples; bars, "
                 "mean of 50 draws; lines, 2.5–97.5 percentile.", fontsize=5.5)
    save(fig, "Fig3_agreement_vs_ceiling")


def fig4():
    p = os.path.join(LT, "controls_summary.csv")
    if not os.path.exists(p):
        return print("  fig4: no loop controls")
    c = pd.read_csv(p)
    fig, axes = plt.subplots(1, 4, figsize=(9.2, 2.5))
    x = np.arange(len(c))
    w = 0.26
    axes[0].bar(x - w, c["nearest_correct_all"], w, color=GREY, label="all links")
    axes[0].bar(x, c["nearest_correct_clean"], w, color=BLUE,
                label="single-TSS anchors")
    axes[0].bar(x + w, c["nearest_correct_decoy_null"], w, color=RED,
                label="distance-matched null")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(c["sample"], rotation=30)
    axes[0].set_ylabel("nearest gene is the contacted gene")
    axes[0].legend(frameon=False, fontsize=6)
    letter(axes[0], "a")

    axes[1].bar(x - w / 2, c["median_nearest_kb"], w, color=GREY, label="nearest TSS")
    axes[1].bar(x + w / 2, c["median_contact_kb"], w, color=BLUE, label="contacted TSS")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(c["sample"], rotation=30)
    axes[1].set_ylabel("median distance (kb)")
    axes[1].legend(frameon=False, fontsize=6)
    letter(axes[1], "b")

    rp = os.path.join(LT, "reproducibility_pairs.csv")
    if os.path.exists(rp):
        r = pd.read_csv(rp)
        gsc = r[r["sample_b"] != "nsc"]["mean_target_jaccard"]
        nsc = r[r["sample_b"] == "nsc"]["mean_target_jaccard"]
        axes[2].boxplot([gsc.dropna(), nsc.dropna()], widths=0.5,
                        patch_artist=True, showfliers=False)
        for patch, col in zip(axes[2].patches if hasattr(axes[2], "patches") else [],
                              [BLUE, GREY]):
            patch.set_facecolor(col)
        axes[2].set_xticklabels(["GSC vs GSC", "GSC vs NSC"], rotation=15)
        axes[2].set_ylabel("target-set Jaccard")
    letter(axes[2], "c")

    # depth matching: how much of the "GSC-specific" claim survives
    dm = os.path.join(LT, "depth_matched_summary.csv")
    gs = os.path.join(LT, "gsc_specific_targets.csv")
    dg = os.path.join(LT, "depth_matched_genes.csv")
    if os.path.exists(dm) and os.path.exists(gs) and os.path.exists(dg):
        b = pd.read_csv(gs)
        n_before = int(((b["n_gsc_samples"] == b["n_gsc_samples"].max())
                        & (b["in_nsc"] == 0)).sum())
        g = pd.read_csv(dg)
        n_after = int((g["frac_repeats"] >= 0.8).sum())
        d = pd.read_csv(dm, index_col=0)
        n_mid = float(d.loc["genes_all_gsc_not_nsc", "mean"]) \
            if "genes_all_gsc_not_nsc" in d.index else np.nan
        vals = [n_before, n_mid, n_after]
        axes[3].bar(range(3), vals, color=[GREY, ORANGE, GREEN], width=0.6)
        for i, v in enumerate(vals):
            if np.isfinite(v):
                axes[3].text(i, v, f"{v:,.0f}", ha="center", va="bottom", fontsize=6)
        axes[3].set_xticks(range(3))
        axes[3].set_xticklabels(["as called", "depth-\nmatched", "robust\n(>=80%)"],
                                fontsize=6)
        axes[3].set_ylabel("genes contacted in all GSC, not NSC")
    letter(axes[3], "d")
    fig.tight_layout()
    save(fig, "Fig4_target_assignment")


def fig5():
    sp = os.path.join(CC, "summary.json")
    gp = os.path.join(CC, "core_domain_genes.csv")
    if not (os.path.exists(sp) and os.path.exists(gp)):
        return print("  fig5: no core domain outputs")
    s = json.load(open(sp))
    g = pd.read_csv(gp)
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    keys = [k for k in ("n_core_domains", "n_gsc_specific_domains",
                        "n_patient_specific_domains", "n_sporadic_domains")
            if k in s]
    vals5 = [s[k] for k in keys]
    axes[0].bar(range(len(keys)), vals5, color=[GREEN, BLUE, RED, GREY][: len(keys)])
    axes[0].set_yscale("log")
    for i, v in enumerate(vals5):
        axes[0].text(i, v * 1.15, f"{v:,}", ha="center", fontsize=6)
    axes[0].set_xticks(range(len(keys)))
    axes[0].set_xticklabels([k.replace("n_", "").replace("_domains", "").replace("_", " ")
                             for k in keys], rotation=25, fontsize=6)
    axes[0].set_ylabel("domains")
    letter(axes[0], "a")

    axes[1].hist(g["width_kb"], bins=40, color=GREEN, alpha=0.85)
    axes[1].set_xlabel("core domain width (kb)")
    axes[1].set_ylabel("domains")
    axes[1].set_yscale("log")
    letter(axes[1], "b")

    top = g.nlargest(14, "width_kb").iloc[::-1]
    def _lab(r):
        n = r["nearest_gene"] if isinstance(r["nearest_gene"], str) else ""
        if n and not n.startswith("ENSG"):
            return n[:16]
        cg = r.get("contact_genes", "")
        if isinstance(cg, str) and cg:
            return (cg.split(";")[0] + " (contact)")[:20]
        return (n or "-")[:16]
    lab = [_lab(r) for _, r in top.iterrows()]
    axes[2].barh(range(len(top)), top["width_kb"], color=GREEN, height=0.65)
    axes[2].set_yticks(range(len(top)))
    axes[2].set_yticklabels(lab, fontsize=6)
    axes[2].set_xlabel("width (kb)")
    letter(axes[2], "c")
    fig.tight_layout()
    save(fig, "Fig5_core_domains")


def fig6():
    p = os.path.join(RES, "motifs", "motif_enrichment.csv")
    if not os.path.exists(p):
        return print("  fig6: no motif table")
    m = pd.read_csv(p)
    qcol = "QVALUE" if "QVALUE" in m.columns else None
    ecol = "ENR_RATIO" if "ENR_RATIO" in m.columns else None
    if not (qcol and ecol):
        return print("  fig6: unexpected SEA columns")
    comps = [c for c in m["comparison"].unique()][:2]
    fig, axes = plt.subplots(1, len(comps), figsize=(3.6 * len(comps), 3.0),
                             squeeze=False)
    name_col = "ALT_ID" if "ALT_ID" in m.columns else "ID"
    for ax, comp in zip(axes[0], comps):
        d = m[(m["comparison"] == comp) & np.isfinite(m[qcol])]
        d = d.nsmallest(15, qcol).iloc[::-1]
        ax.barh(range(len(d)), -np.log10(d[qcol].clip(lower=1e-300)),
                color=BLUE, height=0.65)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels([str(x)[:14] for x in d[name_col]], fontsize=6)
        ax.set_xlabel("-log10 q")
        ax.axvline(-np.log10(0.05), color=RED, lw=0.8, ls="--")
        ax.set_title(comp.replace("_", " "), fontsize=7)
    xmax = max(a.get_xlim()[1] for a in axes[0])
    for a in axes[0]:
        a.set_xlim(0, xmax)
    axes[0][-1].text(-np.log10(0.05), axes[0][-1].get_ylim()[0], " q = 0.05",
                     color=RED, fontsize=6, va="bottom")
    letter(axes[0][0], "a")
    fig.tight_layout()
    save(fig, "Fig6_motif_enrichment")


def fig7():
    p = os.path.join(RES, "expression", "gene_recurrence_expression.csv.gz")
    t = os.path.join(RES, "expression", "expression_tests.csv")
    if not os.path.exists(p):
        return print("  fig7: no expression table")
    g = pd.read_csv(p)
    ecols = [c for c in ("tcga_gbm", "pedi_tpm", "gsc_cpm") if c in g.columns]
    if not ecols or "class" not in g.columns:
        return print("  fig7: expression or class column missing")
    fig, axes = plt.subplots(1, len(ecols) + 1, figsize=(2.3 * (len(ecols) + 1), 2.6),
                             squeeze=False)
    order = ["core", "other_se", "no_se"]
    cols = [GREEN, BLUE, GREY]
    for ax, ec in zip(axes[0], ecols):
        data = [g.loc[g["class"] == o, ec].dropna().values for o in order]
        bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False)
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
            patch.set_linewidth(0.6)
        ax.set_xticks(range(1, len(order) + 1))
        ax.set_xticklabels(["core SE", "other SE", "no SE"], rotation=20, fontsize=6)
        ylab = {"tcga_gbm": "TCGA GBM, median log2(norm. count + 1)",
                "pedi_tpm": "pediatric HGG, median log2(TPM + 1)",
                "gsc_cpm": "GSC lines, median log2(CPM + 1)"}
        ax.set_ylabel(ylab.get(ec, ec), fontsize=6.5)
    for ax, s in zip(axes[0], "abcd"):
        letter(ax, s)
    ax = axes[0][-1]
    if os.path.exists(t):
        tt = pd.read_csv(t)
        d = tt[(tt["test"] == "spearman")].dropna(subset=["rho"])
        d = d[d["recurrence"].isin(["rec_GSC_line", "rec_patient_tumour",
                                    "rec_NSC", "rec_pediatric"])]
        if len(d):
            d = d.sort_values("rho")
            ax.barh(range(len(d)), d["rho"], color=BLUE, height=0.6)
            ax.set_yticks(range(len(d)))
            ax.set_yticklabels([f"{a}/{b.replace('rec_','')}"
                                for a, b in zip(d["expression"], d["recurrence"])],
                               fontsize=5.5)
            ax.set_xlabel("Spearman rho, recurrence vs expression")
    letter(axes[0][0], "a")
    fig.tight_layout()
    save(fig, "Fig7_expression_validation")


def main():
    print("building figures")
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7):
        try:
            fn()
        except Exception as exc:
            print(f"  {fn.__name__} failed: {type(exc).__name__}: {exc}")
    n = len(glob.glob(os.path.join(PNG, "*.png")))
    print(f"{n} figures in {PNG} and {TIF}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
