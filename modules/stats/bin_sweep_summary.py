#!/usr/bin/env python3
"""Aggregate DIFF_BREAKS summary outputs across bin sizes.

Nextflow will provide a simple 3-column TSV on stdin:
  bin_size\tcontrast_name\tsummary_path

Each summary file is expected to contain one key-value per line, e.g.:
  bins_tested: 4500000
  bins_FDR_le_threshold: 113
  sig_up: 60
  sig_down: 53
  fdr_threshold: 0.05
  regions_bed: /path/to/regions.bed  (optional)

We keep parsing tolerant (':' or '\t' separators) to avoid brittleness.

Outputs a single TSV to --out.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, Tuple


KEY_ALIASES = {
    "bins_tested": ["bins_tested", "tested_bins", "n_bins_tested"],
    "sig_total": [
        "bins_FDR_le_threshold",
        "significant_bins",
        "sig_bins",
        "sig_total",
        "n_sig",
    ],
    "sig_up": ["sig_up", "significant_up", "n_sig_up"],
    "sig_down": ["sig_down", "significant_down", "n_sig_down"],
    "fdr_threshold": ["fdr_threshold", "fdr", "padj_threshold"],
    "regions_bed": ["regions_bed", "bed", "regions"],
}


def parse_summary(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not os.path.exists(path):
        return out

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # key: value or key\tvalue
            if ":" in line:
                key, val = line.split(":", 1)
            elif "\t" in line:
                key, val = line.split("\t", 1)
            else:
                continue
            key = key.strip()
            val = val.strip()
            if key:
                out[key] = val
    return out


def pick(d: Dict[str, str], canonical: str) -> str:
    for k in KEY_ALIASES.get(canonical, [canonical]):
        if k in d:
            return d[k]
    return ""


def to_int(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip()
    if not s:
        return ""
    # remove commas and approximate markers
    s2 = s.replace(",", "")
    s2 = re.sub(r"^~", "", s2)
    # keep leading integer
    m = re.match(r"^-?\d+", s2)
    return m.group(0) if m else s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    rows = []
    for raw in sys.stdin:
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split("\t")
        if len(parts) < 3:
            continue
        bin_size, contrast, summary_path = parts[0], parts[1], "\t".join(parts[2:])
        d = parse_summary(summary_path)

        bins_tested = to_int(pick(d, "bins_tested"))
        sig_total = to_int(pick(d, "sig_total"))
        sig_up = to_int(pick(d, "sig_up"))
        sig_down = to_int(pick(d, "sig_down"))
        fdr_thr = pick(d, "fdr_threshold")
        regions_bed = pick(d, "regions_bed")

        rows.append((
            to_int(bin_size),
            contrast,
            bins_tested,
            sig_total,
            sig_up,
            sig_down,
            fdr_thr,
            regions_bed,
        ))

    # Sort by bin_size then contrast for stable output
    def sort_key(r: Tuple[str, ...]):
        try:
            bs = int(r[0])
        except Exception:
            bs = 10**18
        return (bs, r[1])

    rows.sort(key=sort_key)

    with open(args.out, "w", encoding="utf-8") as out:
        out.write(
            "bin_size\tcontrast\tbins_tested\tsig_total\tsig_up\tsig_down\tfdr_threshold\tregions_bed\n"
        )
        for r in rows:
            out.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
