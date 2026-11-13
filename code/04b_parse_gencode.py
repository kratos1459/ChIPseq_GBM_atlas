#!/usr/bin/env python3
"""Read GENCODE gene records on stdin, write gene-body and TSS BEDs.

Keeps protein_coding and lncRNA genes. BED6+2: chrom start end name score
strand gene_type gene_id.
"""
import re
import sys

ATTR = {k: re.compile(k + r' "([^"]+)"') for k in ("gene_type", "gene_name", "gene_id")}


def attr(blob, key):
    m = ATTR[key].search(blob)
    return m.group(1) if m else ""


def main(gene_out, tss_out):
    genes, tss = [], []
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        blob = f[8]
        gtype = attr(blob, "gene_type")
        if gtype not in ("protein_coding", "lncRNA"):
            continue
        name, gid = attr(blob, "gene_name"), attr(blob, "gene_id")
        chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
        t = start if strand == "+" else end - 1
        genes.append(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{gtype}\t{gid}")
        tss.append(f"{chrom}\t{t}\t{t+1}\t{name}\t0\t{strand}\t{gtype}\t{gid}")
    with open(gene_out, "w") as fh:
        fh.write("\n".join(genes) + "\n")
    with open(tss_out, "w") as fh:
        fh.write("\n".join(tss) + "\n")
    print(f"genes={len(genes)} tss={len(tss)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
