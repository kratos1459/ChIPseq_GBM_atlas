#!/usr/bin/env python3
"""Verify every downloaded GEO file against the byte size in GEO's manifest and
re-fetch any that disagree (a resumed download that overlapped a second writer
produces a file larger than the original, which is silently corrupt).

Run after 03_download.py. Exits non-zero if anything is still wrong.
"""
import concurrent.futures as cf
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RAW = os.path.join(DATA, "raw")


def url_for(gse, gsm, fname):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{fname}"


def refetch(url, dest, expect):
    tmp = dest + ".refetch"
    for attempt in range(3):
        if os.path.exists(tmp):
            os.remove(tmp)
        p = subprocess.run(
            ["wget", "-q", "--tries=6", "--timeout=120", "-O", tmp, url],
            capture_output=True, text=True,
        )
        got = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        if p.returncode == 0 and got == expect:
            os.replace(tmp, dest)
            return dest, "repaired", got
    if os.path.exists(tmp):
        os.remove(tmp)
    return dest, f"STILL-BAD got={got} want={expect}", got


def main():
    fl = pd.read_csv(os.path.join(DATA, "filelist.csv"))
    fl = fl[fl["gsm"].notna()]
    bad, missing, ok = [], 0, 0
    for _, r in fl.iterrows():
        dest = os.path.join(RAW, r["gse"], r["file"])
        if not os.path.exists(dest):
            missing += 1
            continue
        if os.path.getsize(dest) != int(r["bytes"]):
            bad.append((url_for(r["gse"], r["gsm"], r["file"]), dest, int(r["bytes"])))
        else:
            ok += 1
    print(f"present-and-correct={ok}  size-mismatch={len(bad)}  not-downloaded={missing}")
    for _, d, _ in bad:
        print("  bad:", os.path.basename(d))
    if not bad:
        return 0
    fails = 0
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for dest, status, got in ex.map(lambda j: refetch(*j), bad):
            print(f"  {status:<28} {got/1e6:8.1f} MB  {os.path.basename(dest)}")
            if status != "repaired":
                fails += 1
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
