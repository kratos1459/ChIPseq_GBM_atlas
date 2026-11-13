#!/usr/bin/env python3
"""Build a per-sample metadata table and a processed-file inventory for the
four GEO series in the GBM enhancer atlas.

Reads nothing but GEO FTP. Writes:
  data/geo_soft/<GSE>_family.soft.gz   raw SOFT records (kept for provenance)
  data/sample_metadata.csv             one row per GSM
  data/suppl_inventory.csv             one row per series-level processed file
"""
import gzip
import io
import os
import re
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
SOFT = os.path.join(DATA, "geo_soft")
os.makedirs(SOFT, exist_ok=True)

SERIES = ["GSE119776", "GSE145646", "GSE262089", "GSE162976"]
UA = {"User-Agent": "Mozilla/5.0 (chipseq-atlas metadata fetch)"}


def ftp_prefix(gse: str) -> str:
    stub = gse[:-3] + "nnn" if len(gse) > 6 else gse
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}"


def get(url: str, tries: int = 5, timeout: int = 180) -> bytes:
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return fh.read()
        except Exception as exc:          # transient FTP/HTTP hiccups are common
            last = exc
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def fetch_soft(gse: str) -> str:
    dest = os.path.join(SOFT, f"{gse}_family.soft.gz")
    if not (os.path.exists(dest) and os.path.getsize(dest) > 1000):
        raw = get(f"{ftp_prefix(gse)}/soft/{gse}_family.soft.gz")
        with open(dest, "wb") as fh:
            fh.write(raw)
    with gzip.open(dest, "rt", errors="replace") as fh:
        return fh.read()


def parse_soft(gse: str, text: str):
    """Return (series_fields, [sample_dicts])."""
    series = {"gse": gse}
    samples, cur = [], None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gse": gse, "gsm": line.split("=", 1)[1].strip(), "chars": []}
            continue
        if line.startswith("^SERIES") or line.startswith("^PLATFORM"):
            if cur:
                samples.append(cur)
                cur = None
            continue
        if not line.startswith("!"):
            continue
        try:
            key, val = line[1:].split("=", 1)
        except ValueError:
            continue
        key, val = key.strip(), val.strip()
        if key.startswith("Series_") and cur is None:
            k = key[len("Series_"):]
            if k in ("title", "summary", "overall_design", "pubmed_id"):
                series[k] = (series.get(k, "") + " " + val).strip()
        elif cur is not None and key.startswith("Sample_"):
            k = key[len("Sample_"):]
            if k.startswith("characteristics"):
                cur["chars"].append(val)
            elif k == "relation":
                cur.setdefault("relation", []).append(val)
            elif k == "supplementary_file" or k.startswith("supplementary_file"):
                cur.setdefault("suppl", []).append(val)
            else:
                cur[k] = (cur.get(k, "") + " " + val).strip() if k in cur else val
    if cur:
        samples.append(cur)
    return series, samples


MARKS = [
    ("H3K27ac", r"h3k27ac|k27ac"),
    ("H3K27me3", r"h3k27me3|k27me3"),
    ("H3K4me1", r"h3k4me1|k4me1"),
    ("H3K4me3", r"h3k4me3|k4me3"),
    ("H3K9me3", r"h3k9me3"),
    ("H3K36me3", r"h3k36me3"),
    ("H3K27M", r"h3k27m\b|k27m\b"),
    ("BRD4", r"brd4"),
    ("CTCF", r"ctcf"),
    ("OLIG2", r"olig2"),
    ("SOX2", r"sox2"),
    ("POLR2A", r"pol\s?2|polr2|rna\s?pol"),
    ("EP300", r"ep300|p300"),
    ("input", r"\binput\b|wce|igg|control"),
]


def classify(sample: dict) -> dict:
    blob = " ".join(
        [sample.get("title", ""), sample.get("source_name_ch1", "")]
        + sample.get("chars", [])
    ).lower()
    mark = ""
    for name, pat in MARKS:
        if re.search(pat, blob):
            mark = name
            break
    strat = sample.get("library_strategy", "")
    if not mark and strat.lower() in ("rna-seq", "atac-seq"):
        mark = strat
    # tissue / model class
    if re.search(r"\bnsc\b|neural stem|normal brain|non-?tumou?r|astrocyte", blob):
        tissue = "normal/NSC"
    elif re.search(r"patient|tumou?r|biopsy|surgical|primary gbm|resect", blob):
        tissue = "patient tumour"
    elif re.search(r"\bgsc\b|stem cell line|cell line|xenograft|pdx", blob):
        tissue = "GSC/line"
    else:
        tissue = ""
    srx = ""
    for rel in sample.get("relation", []):
        m = re.search(r"(SRX\d+|SRP\d+)", rel)
        if m:
            srx = m.group(1)
            break
    return {
        "gse": sample["gse"],
        "gsm": sample["gsm"],
        "title": sample.get("title", ""),
        "mark": mark,
        "tissue_class": tissue,
        "library_strategy": strat,
        "instrument": sample.get("instrument_model", ""),
        "organism": sample.get("organism_ch1", ""),
        "source_name": sample.get("source_name_ch1", ""),
        "characteristics": " | ".join(sample.get("chars", [])),
        "sra": srx,
        "n_suppl": len(sample.get("suppl", [])),
        "suppl_files": " | ".join(
            os.path.basename(u) for u in sample.get("suppl", []) if "ftp" in u or "/" in u
        ),
    }


def suppl_listing(gse: str):
    """Series-level processed files, from the FTP suppl/ directory index."""
    try:
        html = get(f"{ftp_prefix(gse)}/suppl/", tries=3, timeout=120).decode(
            "utf-8", "replace"
        )
    except Exception as exc:
        print(f"  ! suppl listing failed for {gse}: {exc}", file=sys.stderr)
        return []
    rows = []
    # NCBI's FTP-over-HTTP index: <a href="file">file</a> ... size
    for m in re.finditer(r'href="([^"?][^"]*)"', html):
        name = m.group(1)
        if name in ("../",) or name.endswith("/"):
            continue
        rows.append({"gse": gse, "file": name})
    # sizes, when the index supplies them
    sizes = dict(
        (a, b)
        for a, b in re.findall(
            r'href="([^"?][^"]*)"[^<]*</a>\s*</td>\s*<td[^>]*>[^<]*</td>\s*<td[^>]*>\s*([0-9.]+\s*[KMGT]?)',
            html,
        )
    )
    for r in rows:
        r["size"] = sizes.get(r["file"], "")
    return rows


def main():
    import pandas as pd

    all_samples, all_suppl, series_rows = [], [], []
    for gse in SERIES:
        print(f"== {gse}")
        text = fetch_soft(gse)
        series, samples = parse_soft(gse, text)
        series_rows.append(series)
        rows = [classify(s) for s in samples]
        all_samples += rows
        print(f"   {len(rows)} samples; title: {series.get('title','')[:90]}")
        sup = suppl_listing(gse)
        all_suppl += sup
        print(f"   {len(sup)} series-level processed files")

    sm = pd.DataFrame(all_samples)
    sm.to_csv(os.path.join(DATA, "sample_metadata.csv"), index=False)
    pd.DataFrame(all_suppl).to_csv(os.path.join(DATA, "suppl_inventory.csv"), index=False)
    pd.DataFrame(series_rows).to_csv(os.path.join(DATA, "series_metadata.csv"), index=False)

    print("\n== samples per series x mark ==")
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(pd.crosstab(sm["gse"], sm["mark"].replace("", "(unclassified)")))
    print("\n== tissue class ==")
    print(pd.crosstab(sm["gse"], sm["tissue_class"].replace("", "(unclassified)")))
    print(f"\nwrote {len(sm)} samples -> data/sample_metadata.csv")


if __name__ == "__main__":
    main()
