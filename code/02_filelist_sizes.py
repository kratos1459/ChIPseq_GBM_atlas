#!/usr/bin/env python3
"""Pull GEO's per-file manifest (suppl/filelist.txt) for each series so the
download volume is known before anything is fetched. Writes data/filelist.csv.

Also re-classifies marks with 'input' tested FIRST, so that samples titled
e.g. GBM1_H3K27AC_INPUT are not counted as H3K27ac ChIP.
"""
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module

m = import_module("01_geo_metadata") if False else None  # avoid numeric-module import

import time
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
SERIES = ["GSE119776", "GSE145646", "GSE262089", "GSE162976"]
UA = {"User-Agent": "Mozilla/5.0 (chipseq-atlas manifest fetch)"}


def ftp_prefix(gse):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}"


def get(url, tries=5, timeout=180):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=timeout
            ) as fh:
                return fh.read()
        except Exception as exc:
            last = exc
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def main():
    import pandas as pd

    rows = []
    for gse in SERIES:
        txt = get(f"{ftp_prefix(gse)}/suppl/filelist.txt").decode("utf-8", "replace")
        lines = [l for l in txt.splitlines() if l.strip()]
        hdr = lines[0].split("\t")
        for line in lines[1:]:
            parts = line.split("\t")
            d = dict(zip(hdr, parts))
            name = d.get("Name") or parts[1] if len(parts) > 1 else ""
            size = d.get("Size") or ""
            try:
                size = int(str(size).replace(",", ""))
            except ValueError:
                size = 0
            rows.append(
                {
                    "gse": gse,
                    "file": name,
                    "bytes": size,
                    "type": d.get("Type", ""),
                }
            )
        print(f"{gse}: {len(lines)-1} files listed")

    fl = pd.DataFrame(rows)
    fl["gsm"] = fl["file"].str.extract(r"^(GSM\d+)")
    fl["ext"] = fl["file"].str.replace(r"^GSM\d+_", "", regex=True).str.extract(
        r"(\.[A-Za-z0-9.]+)$"
    )
    fl.to_csv(os.path.join(DATA, "filelist.csv"), index=False)

    sm = pd.read_csv(os.path.join(DATA, "sample_metadata.csv")).fillna("")
    # input-first re-classification
    def reclass(r):
        blob = f"{r['title']} {r['source_name']} {r['characteristics']}".lower()
        if re.search(r"input|\bwce\b|\bigg\b", blob):
            return "input"
        return r["mark"]

    sm["mark2"] = sm.apply(reclass, axis=1)
    sm.to_csv(os.path.join(DATA, "sample_metadata.csv"), index=False)

    j = fl.merge(sm[["gsm", "gse", "mark2", "tissue_class"]], on=["gsm", "gse"], how="left")
    gb = (
        j.groupby(["gse", "mark2", "ext"])
        .agg(n=("file", "size"), gb=("bytes", lambda s: round(s.sum() / 1e9, 2)))
        .reset_index()
        .sort_values(["gse", "gb"], ascending=[True, False])
    )
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print("\n== files by series / mark / extension ==")
        print(gb[gb["n"] > 0].to_string(index=False))
    print("\n== H3K27ac ChIP only (no input) ==")
    k = j[(j["mark2"] == "H3K27ac")]
    print(
        k.groupby(["gse", "ext"])
        .agg(n=("file", "size"), gb=("bytes", lambda s: round(s.sum() / 1e9, 2)))
        .to_string()
    )
    print(f"\ntotal across all listed files: {j['bytes'].sum()/1e9:.1f} GB")
    j.to_csv(os.path.join(DATA, "filelist_annotated.csv"), index=False)


if __name__ == "__main__":
    main()
