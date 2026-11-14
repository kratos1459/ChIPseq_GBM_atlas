#!/usr/bin/env python3
"""v2 analysis definitions (2026-11-14), applied as a reviewed patch.

Why, from the v1 run:
  * Super-enhancer counts from the ROSE inflection rule tracked data provenance
    (median 280 in the hg19 bedGraph GSC cohort, 606 in NSC controls), so the
    fixed top-500 domains per sample becomes the primary definition.
  * "<10% recurrence in NSC" with only 12 NSC samples excludes a bin as soon as
    two NSC samples carry it; it removed 528 of 539 jointly recurrent bins.
    Relaxed to <25%, i.e. at most 2 of 12.
  * ROC AUC read 0.93-0.999 for every comparison because ~1 bin in 60 is
    positive; average precision is reported instead.
"""
import re
import sys

C30 = "code/30_cross_cohort.py"
C50 = "code/50_expression.py"
C60 = "code/60_figures.py"
C70 = "code/70_write_report.py"


def sub(text, old, new, count=1, path=""):
    if old not in text:
        sys.exit(f"patch failed in {path}: pattern not found:\n{old[:120]}")
    return text.replace(old, new, count)


# ------------------------------------------------------------ 30_cross ----
t = open(C30).read()
t = sub(t, "from sklearn.metrics import roc_auc_score",
        "from sklearn.metrics import average_precision_score, roc_auc_score", path=C30)
t = sub(t, "TOP_N = 500 ",
        "NSC_MAX = 0.25       # at most 2 of 12 NSC samples (v2; was 0.10)\n"
        "OTHER_MAX = 0.25     # 'absent' from the other adult cohort (v2; was 0.10)\n"
        "TOP_N = 500 ", path=C30)
# average precision alongside AUC, in the ceiling and in cross-cohort stats
t = sub(t, '    out = {"spearman": [], "auc": [], "jaccard": []}',
        '    out = {"spearman": [], "auc": [], "ap": [], "jaccard": []}', path=C30)
t = sub(t, '            out["auc"].append(roc_auc_score(yb, ra[keep]))',
        '            out["auc"].append(roc_auc_score(yb, ra[keep]))\n'
        '            out["ap"].append(average_precision_score(yb, ra[keep]))', path=C30)
t = sub(t, '    ya_t, yb_t = ra >= thresh, rb >= thresh',
        '    res["ap_predicting_b"] = (\n'
        '        float(average_precision_score(yb, ra[keep])) if 0 < yb.sum() < yb.size else np.nan\n'
        '    )\n'
        '    res["ap_baseline"] = float(yb.mean())\n'
        '    ya_t, yb_t = ra >= thresh, rb >= thresh', path=C30)
# ceiling rows: carry AP (both the inflection and the top-N blocks)
t = t.replace('"auc_predicting_b": c.get("auc", (np.nan,))[0],',
              '"auc_predicting_b": c.get("auc", (np.nan,))[0],\n'
              '                "ap_predicting_b": c.get("ap", (np.nan,))[0],')
# recurrence table also carries the top-N recurrence
t = sub(t, "    any_se = np.zeros(bi.nbins, dtype=bool)\n    for r in rec.values():\n        any_se |= r > 0\n",
        "    rec_top_early = {g: np.vstack(v).mean(0) for g, v in top_by_group.items() if len(v) >= 3}\n"
        "    any_se = np.zeros(bi.nbins, dtype=bool)\n    for r in rec.values():\n        any_se |= r > 0\n"
        "    for r in rec_top_early.values():\n        any_se |= r > 0\n", path=C30)
t = sub(t, '    for g, r in rec_enh.items():\n        loci[f"enh_rec_{g}"] = r[idx]\n',
        '    for g, r in rec_enh.items():\n        loci[f"enh_rec_{g}"] = r[idx]\n'
        '    for g, r in rec_top_early.items():\n        loci[f"top_rec_{g}"] = r[idx]\n', path=C30)
# core / specific / pediatric blocks use the top-N recurrence
a = t.index("    # ---- consensus core domains")
b = t.index("        # ---- genes at core domains")
block = t[a:b]
block = block.replace('summary = {"binsize": BINSIZE, "n_bins_any_se": int(len(idx))}',
                      'summary = {"binsize": BINSIZE, "n_bins_any_se": int(len(idx)),\n'
                      '               "definition": f"fixed top-{TOP_N} domains per sample",\n'
                      '               "nsc_max": NSC_MAX, "other_max": OTHER_MAX}\n'
                      '    R = rec_top')
block = re.sub(r'\brec\[', 'R[', block)
block = block.replace('rec.get("NSC"', 'R.get("NSC"')
block = block.replace('if "GSC_line" in rec and "patient_tumour" in rec:',
                      'if "GSC_line" in R and "patient_tumour" in R:')
block = block.replace('if "pediatric" in rec:', 'if "pediatric" in R:')
block = block.replace('core & (nsc < 0.10)', 'core & (nsc < NSC_MAX)')
block = block.replace('< 0.1)', '< OTHER_MAX)')
block = block.replace('<10% of NSC', f'<25% of NSC')
t = t[:a] + block + t[b:]
open(C30, "w").write(t)
print("patched", C30)

# ---------------------------------------------------------- 50_expression --
t = open(C50).read()
t = sub(t, '    cols = [c for c in rec.columns if c.startswith("se_rec_")]',
        '    cols = ([c for c in rec.columns if c.startswith("top_rec_")]\n'
        '            or [c for c in rec.columns if c.startswith("se_rec_")])', path=C50)
t = sub(t, '        genes[c.replace("se_rec_", "rec_")] = out[c]',
        '        genes[c.replace("top_rec_", "rec_").replace("se_rec_", "rec_")] = out[c]', path=C50)
t = sub(t, '(g["rec_adult_min"] >= 0.5) & (nsc < 0.10)',
        '(g["rec_adult_min"] >= 0.5) & (nsc < 0.25)', path=C50)
open(C50, "w").write(t)
print("patched", C50)

# ------------------------------------------------------------ 60_figures --
t = open(C60).read()
t = sub(t, '    cols = [c for c in r.columns if c.startswith("se_rec_")]',
        '    cols = ([c for c in r.columns if c.startswith("top_rec_")]\n'
        '            or [c for c in r.columns if c.startswith("se_rec_")])', path=C60)
t = sub(t, '        g = c[len("se_rec_"):]',
        '        g = c.split("rec_", 1)[1]', path=C60)
t = sub(t, '        x, y = r["se_rec_GSC_line"], r["se_rec_patient_tumour"]',
        '        p = "top_rec_" if "top_rec_GSC_line" in r.columns else "se_rec_"\n'
        '        x, y = r[p + "GSC_line"], r[p + "patient_tumour"]', path=C60)
t = sub(t, '    if "se_rec_GSC_line" in r.columns and "se_rec_patient_tumour" in r.columns:',
        '    if ("top_rec_GSC_line" in r.columns or "se_rec_GSC_line" in r.columns):', path=C60)
t = sub(t, '               ("auc_predicting_b", "AUC predicting target cohort"),',
        '               ("ap_predicting_b", "average precision, target cohort"),', path=C60)
t = sub(t, '    p = os.path.join(CC, "agreement.csv")\n    if not os.path.exists(p):\n        return print("  fig3: no agreement table")\n    a = pd.read_csv(p)',
        '    p = os.path.join(CC, "agreement_topN.csv")\n'
        '    if not os.path.exists(p):\n        p = os.path.join(CC, "agreement.csv")\n'
        '    if not os.path.exists(p):\n        return print("  fig3: no agreement table")\n    a = pd.read_csv(p)', path=C60)
open(C60, "w").write(t)
print("patched", C60)

# ---------------------------------------------------------- 70_report ------
t = open(C70).read()
t = sub(t, '"the ROC curve when one cohort\'s recurrence is used to predict whether a "\n      "bin is recurrent in the other, ',
        '"the precision-recall curve (average precision) when one cohort\'s "\n      "recurrence is used to predict whether a bin is recurrent in the other "\n      "(ROC AUC is not used: with roughly one positive bin in sixty it reads "\n      "above 0.9 for every comparison, including unrelated ones), ', path=C70)
t = sub(t, '    agree = read_csv(os.path.join(CC, "agreement.csv"))',
        '    agree = read_csv(os.path.join(CC, "agreement_both_definitions.csv"))\n'
        '    if agree is None:\n        agree = read_csv(os.path.join(CC, "agreement.csv"))', path=C70)
open(C70, "w").write(t)
print("patched", C70)
