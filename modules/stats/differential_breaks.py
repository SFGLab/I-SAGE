#!/usr/bin/env python3
"""
Differential break analysis on binned bedGraph tracks.

Method (default): Fisher exact test per bin vs rest-of-genome:
  For each bin:
    [[k_case,  total_case - k_case],
     [k_ctrl,  total_ctrl - k_ctrl]]
  p-value = fisher_exact(two-sided)
  log2FC  = log2( (k_case/total_case + eps) / (k_ctrl/total_ctrl + eps) )
  FDR     = Benjamini-Hochberg

Inputs:
  --samples-tsv: sample_id <tab> plus_bedGraph <tab> minus_bedGraph
  --contrast: name:CASE1,CASE2:CTRL1,CTRL2  (case and control lists by sample_id)
Outputs:
  <out-prefix>.tsv
  <out-prefix>.volcano.png
"""

from __future__ import annotations

import argparse
import math
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # ✅ headless-safe for HPC
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact


def read_bedgraph(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "value"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "value": np.float64},
    )
    return df


def load_sample_total_track(plus_path: str, minus_path: str) -> pd.DataFrame:
    plus = read_bedgraph(plus_path)
    minus = read_bedgraph(minus_path)

    key = ["chrom", "start", "end"]
    merged = plus.merge(minus, on=key, how="outer", suffixes=("_plus", "_minus"))
    merged["value_plus"] = merged["value_plus"].fillna(0.0)
    merged["value_minus"] = merged["value_minus"].fillna(0.0)
    merged["count"] = merged["value_plus"] + merged["value_minus"]
    return merged[key + ["count"]]


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (returns q-values)."""
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty_like(ranked)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev
    out = np.empty_like(q)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def parse_contrast(s: str) -> Tuple[str, List[str], List[str]]:
    parts = s.split(":")
    if len(parts) != 3:
        raise ValueError("Contrast must be 'name:CASE1,CASE2:CTRL1,CTRL2'")
    name = parts[0]
    case = [x for x in parts[1].split(",") if x]
    ctrl = [x for x in parts[2].split(",") if x]
    if not case or not ctrl:
        raise ValueError("Contrast must include at least 1 case and 1 control sample_id")
    return name, case, ctrl


def volcano_plot(df: pd.DataFrame, out_png: str) -> None:
    x = df["log2FC"].to_numpy()
    y = -np.log10(np.clip(df["padj"].to_numpy(), 1e-300, 1.0))
    plt.figure()
    plt.scatter(x, y, s=3)
    plt.xlabel("log2 fold-change (case vs control)")
    plt.ylabel("-log10(FDR)")
    plt.title("Differential breaks volcano")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples-tsv", required=True, help="TSV: sample_id\\tplus_bg\\tminus_bg")
    ap.add_argument("--contrast", required=True, help="name:CASE1,CASE2:CTRL1,CTRL2")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--eps", type=float, default=0.5, help="Pseudocount for log2FC stability")
    args = ap.parse_args()

    cname, case_ids, ctrl_ids = parse_contrast(args.contrast)

    samp = pd.read_csv(
        args.samples_tsv,
        sep="\t",
        header=None,
        names=["sample_id", "plus_bg", "minus_bg"],
    )

    samp_map: Dict[str, Tuple[str, str]] = {
        str(r.sample_id): (str(r.plus_bg), str(r.minus_bg))
        for r in samp.itertuples(index=False)
    }

    missing = [s for s in case_ids + ctrl_ids if s not in samp_map]
    if missing:
        raise SystemExit(f"ERROR: contrast references unknown samples: {missing}")

    def sum_group(sample_ids: List[str]) -> pd.DataFrame:
        acc: pd.DataFrame | None = None
        for sid in sample_ids:
            plus_path, minus_path = samp_map[sid]
            df = load_sample_total_track(str(plus_path), str(minus_path))
            df = df.rename(columns={"count": sid})
            if acc is None:
                acc = df
            else:
                acc = acc.merge(df, on=["chrom", "start", "end"], how="outer")

        assert acc is not None

        for sid in sample_ids:
            if sid in acc.columns:
                acc[sid] = acc[sid].fillna(0.0)

        acc["sum"] = acc[sample_ids].sum(axis=1)
        return acc[["chrom", "start", "end", "sum"]].rename(columns={"sum": "count"})

    case_df = sum_group(case_ids)
    ctrl_df = sum_group(ctrl_ids)

    merged = case_df.merge(
        ctrl_df,
        on=["chrom", "start", "end"],
        how="outer",
        suffixes=("_case", "_ctrl"),
    )
    merged["count_case"] = merged["count_case"].fillna(0.0)
    merged["count_ctrl"] = merged["count_ctrl"].fillna(0.0)

    total_case = float(merged["count_case"].sum())
    total_ctrl = float(merged["count_ctrl"].sum())
    if total_case <= 0 or total_ctrl <= 0:
        raise SystemExit(f"ERROR: total_case or total_ctrl is <= 0 (case={total_case}, ctrl={total_ctrl})")

    pvals = np.empty(len(merged), dtype=float)
    log2fc = np.empty(len(merged), dtype=float)

    eps = float(args.eps)
    case_arr = merged["count_case"].to_numpy()
    ctrl_arr = merged["count_ctrl"].to_numpy()

    for i, (kc, kt) in enumerate(zip(case_arr, ctrl_arr)):
        kc = float(kc)
        kt = float(kt)

        rc = max(total_case - kc, 0.0)
        rt = max(total_ctrl - kt, 0.0)

        table = np.array(
            [[int(round(kc)), int(round(rc))],
             [int(round(kt)), int(round(rt))]],
            dtype=np.int64,
        )

        _, p = fisher_exact(table, alternative="two-sided")
        pvals[i] = float(p)

        pc = (kc + eps) / (total_case + eps)
        pt = (kt + eps) / (total_ctrl + eps)
        log2fc[i] = math.log2(pc / pt)

    merged["pvalue"] = pvals
    merged["padj"] = bh_fdr(pvals)
    merged["log2FC"] = log2fc

    out_tsv = f"{args.out_prefix}.tsv"
    merged_out = merged[["chrom", "start", "end", "count_case", "count_ctrl", "log2FC", "pvalue", "padj"]]
    merged_out.to_csv(out_tsv, sep="\t", index=False)

    out_png = f"{args.out_prefix}.volcano.png"
    volcano_plot(merged_out, out_png)

    sig = int((merged["padj"] <= 0.05).sum())
    with open(f"{args.out_prefix}.summary.txt", "w") as f:
        f.write(f"contrast\t{cname}\n")
        f.write(f"case\t{','.join(case_ids)}\n")
        f.write(f"control\t{','.join(ctrl_ids)}\n")
        f.write(f"total_case\t{total_case}\n")
        f.write(f"total_ctrl\t{total_ctrl}\n")
        f.write(f"bins_tested\t{len(merged)}\n")
        f.write(f"bins_FDR_le_0.05\t{sig}\n")


if __name__ == "__main__":
    main()