#!/usr/bin/env python3
"""
Differential break analysis on binned bedGraph tracks.

Method (default): Fisher exact test per bin vs rest-of-genome (within tested universe):
  For each bin:
    [[k_case,  total_case - k_case],
     [k_ctrl,  total_ctrl - k_ctrl]]
  p-value = fisher_exact(two-sided)
  log2FC  = log2( (k_case/total_case + eps) / (k_ctrl/total_ctrl + eps) )
  FDR     = Benjamini-Hochberg

Replicate handling:
  --replicate-method pooled:
      pool counts across samples then test.
  --replicate-method meta_fisher:
      compute per-case-replicate p-values vs pooled control and combine using Fisher's method;
      report effect as median per-replicate log2FC.

Optional:
  --regions-bed: restrict testing to bins overlapping provided BED regions (totals computed within restricted set)
  --ebv-regex: annotate EBV contigs by chrom name regex and report EBV enrichment among significant bins

Inputs:
  --samples-tsv: sample_id <tab> plus_bedGraph <tab> minus_bedGraph
  --contrast: name:CASE1,CASE2:CTRL1,CTRL2

Outputs (prefix = --out-prefix):
  <prefix>.tsv
  <prefix>.sig.tsv
  <prefix>.sig_up.tsv
  <prefix>.sig_down.tsv
  <prefix>.summary.txt
  <prefix>.volcano.png / .pdf
  <prefix>.ma.png / .pdf
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple, Literal

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # headless-safe (HPC)
import matplotlib.pyplot as plt

from scipy.stats import fisher_exact, chi2


def read_bed_regions(path: str) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Read BED (chrom, start, end) and build per-chrom interval index.

    Returns dict: chrom -> (starts_sorted, ends_sorted, prefix_max_end)
    The prefix max end allows O(log n) overlap queries for many bins.
    """
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64},
    )
    bed = bed.dropna()
    bed = bed[bed["end"] > bed["start"]]

    out: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, sub in bed.groupby("chrom", sort=False):
        sub = sub.sort_values("start", kind="mergesort")
        starts = sub["start"].to_numpy(dtype=np.int64)
        ends = sub["end"].to_numpy(dtype=np.int64)
        prefix_max_end = np.maximum.accumulate(ends)
        out[str(chrom)] = (starts, ends, prefix_max_end)
    return out


def bins_overlap_regions(
    bins: pd.DataFrame,
    regions: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> np.ndarray:
    """Return boolean mask of bins overlapping any region in BED.

    A bin [s,e) overlaps a region if exists interval with start < e and end > s.
    """
    mask = np.zeros(len(bins), dtype=bool)
    for chrom, idx in bins.groupby("chrom").groups.items():
        key = str(chrom)
        if key not in regions:
            continue
        starts, _ends, prefix_max_end = regions[key]
        s = bins.loc[idx, "start"].to_numpy(dtype=np.int64)
        e = bins.loc[idx, "end"].to_numpy(dtype=np.int64)

        j = np.searchsorted(starts, e, side="right") - 1
        ok = j >= 0
        over = np.zeros_like(ok, dtype=bool)
        over[ok] = prefix_max_end[j[ok]] > s[ok]
        mask[idx] = over
    return mask


def read_bedgraph(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "value"],
        dtype={"chrom": str, "start": np.int64, "end": np.int64, "value": np.float64},
    )


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


def volcano_plot(df: pd.DataFrame, out_prefix: str, contrast_name: str, fdr_thr: float, subtitle: str = "") -> None:
    df2 = df.copy()
    df2["neglog10_p"] = -np.log10(df2["pvalue"].clip(lower=1e-300))

    sig = df2["padj"] <= fdr_thr
    sig_up = sig & (df2["log2FC"] > 0)
    sig_down = sig & (df2["log2FC"] < 0)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(df2.loc[~sig, "log2FC"], df2.loc[~sig, "neglog10_p"], s=6, alpha=0.6)
    ax.scatter(df2.loc[sig_up, "log2FC"], df2.loc[sig_up, "neglog10_p"], s=10, marker="^", alpha=0.9)
    ax.scatter(df2.loc[sig_down, "log2FC"], df2.loc[sig_down, "neglog10_p"], s=10, marker="v", alpha=0.9)

    yline = -np.log10(max(fdr_thr, 1e-300))
    ax.axhline(y=yline, linestyle="--", linewidth=1)

    ax.set_xlabel("log2 fold-change (case vs control)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{contrast_name} — Volcano")
    if subtitle:
        ax.text(0.01, 0.99, subtitle, transform=ax.transAxes, ha="left", va="top", fontsize=9)

    fig.tight_layout()
    fig.savefig(f"{out_prefix}.volcano.png", dpi=200)
    fig.savefig(f"{out_prefix}.volcano.pdf")
    plt.close(fig)


def ma_plot(df: pd.DataFrame, out_prefix: str, contrast_name: str, fdr_thr: float, subtitle: str = "") -> None:
    df2 = df.copy()
    df2["A"] = np.log10((df2["count_case"] + df2["count_ctrl"]) / 2.0 + 1.0)

    sig = df2["padj"] <= fdr_thr
    sig_up = sig & (df2["log2FC"] > 0)
    sig_down = sig & (df2["log2FC"] < 0)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(df2.loc[~sig, "A"], df2.loc[~sig, "log2FC"], s=6, alpha=0.6)
    ax.scatter(df2.loc[sig_up, "A"], df2.loc[sig_up, "log2FC"], s=10, marker="^", alpha=0.9)
    ax.scatter(df2.loc[sig_down, "A"], df2.loc[sig_down, "log2FC"], s=10, marker="v", alpha=0.9)

    ax.axhline(y=0.0, linestyle="--", linewidth=1)
    ax.set_xlabel("log10(mean count + 1)")
    ax.set_ylabel("log2 fold-change (case vs control)")
    ax.set_title(f"{contrast_name} — MA")
    if subtitle:
        ax.text(0.01, 0.99, subtitle, transform=ax.transAxes, ha="left", va="top", fontsize=9)

    fig.tight_layout()
    fig.savefig(f"{out_prefix}.ma.png", dpi=200)
    fig.savefig(f"{out_prefix}.ma.pdf")
    plt.close(fig)


def fisher_pvalue_for_bin(k_case: float, k_ctrl: float, total_case: float, total_ctrl: float) -> float:
    rc = max(total_case - k_case, 0.0)
    rt = max(total_ctrl - k_ctrl, 0.0)
    table = np.array(
        [[int(round(k_case)), int(round(rc))],
         [int(round(k_ctrl)), int(round(rt))]],
        dtype=np.int64,
    )
    _, p = fisher_exact(table, alternative="two-sided")
    return float(p)


def log2fc_for_bin(k_case: float, k_ctrl: float, total_case: float, total_ctrl: float, eps: float) -> float:
    pc = (k_case + eps) / (total_case + eps)
    pt = (k_ctrl + eps) / (total_ctrl + eps)
    return float(math.log2(pc / pt))


def combine_pvalues_fisher(pvals: List[float]) -> float:
    clean = [p for p in pvals if 0.0 < p <= 1.0 and not math.isnan(p)]
    if not clean:
        return 1.0
    stat = -2.0 * sum(math.log(max(p, 1e-300)) for p in clean)
    df = 2 * len(clean)
    return float(chi2.sf(stat, df))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples-tsv", required=True, help="TSV: sample_id\tplus_bg\tminus_bg")
    ap.add_argument("--contrast", required=True, help="name:CASE1,CASE2:CTRL1,CTRL2")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument(
        "--regions-bed",
        default=None,
        help="Optional BED file. If provided, only bins overlapping these regions are tested (and totals are computed within these bins).",
    )
    ap.add_argument(
        "--ebv-regex",
        default="",
        help=(
            "Regex to identify EBV contigs by chromosome name for EBV enrichment reporting. "
            "If empty (default), EBV reporting is disabled."
        ),
    )
    ap.add_argument("--eps", type=float, default=0.5, help="Pseudocount for log2FC stability")
    ap.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for significance reporting")
    ap.add_argument(
        "--replicate-method",
        choices=["pooled", "meta_fisher"],
        default="pooled",
        help="pooled sums counts across samples then tests; meta_fisher combines per-replicate p-values (Fisher) and uses median log2FC.",
    )
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

    def load_group_matrix(sample_ids: List[str]) -> pd.DataFrame:
        acc: pd.DataFrame | None = None
        for sid in sample_ids:
            plus_path, minus_path = samp_map[sid]
            df = load_sample_total_track(str(plus_path), str(minus_path)).rename(columns={"count": sid})
            acc = df if acc is None else acc.merge(df, on=["chrom", "start", "end"], how="outer")
        assert acc is not None
        for sid in sample_ids:
            if sid in acc.columns:
                acc[sid] = acc[sid].fillna(0.0)
        return acc

    case_mat = load_group_matrix(case_ids)
    ctrl_mat = load_group_matrix(ctrl_ids)

    case_df = case_mat[["chrom", "start", "end"]].copy()
    ctrl_df = ctrl_mat[["chrom", "start", "end"]].copy()
    case_df["count"] = case_mat[case_ids].sum(axis=1)
    ctrl_df["count"] = ctrl_mat[ctrl_ids].sum(axis=1)

    merged = case_df.merge(
        ctrl_df,
        on=["chrom", "start", "end"],
        how="outer",
        suffixes=("_case", "_ctrl"),
    )
    merged["count_case"] = merged["count_case"].fillna(0.0)
    merged["count_ctrl"] = merged["count_ctrl"].fillna(0.0)

    regions_bed = str(args.regions_bed) if args.regions_bed else None
    if regions_bed:
        regions = read_bed_regions(regions_bed)
        keep = bins_overlap_regions(merged[["chrom", "start", "end"]], regions)
        merged = merged.loc[keep].reset_index(drop=True)
        if len(merged) == 0:
            raise SystemExit(f"ERROR: No bins overlap regions BED: {regions_bed}")

    total_case = float(merged["count_case"].sum())
    total_ctrl = float(merged["count_ctrl"].sum())
    if total_case <= 0 or total_ctrl <= 0:
        raise SystemExit(f"ERROR: total_case or total_ctrl is <= 0 (case={total_case}, ctrl={total_ctrl})")

    pvals = np.empty(len(merged), dtype=float)
    log2fc = np.empty(len(merged), dtype=float)
    eps = float(args.eps)

    rep_method: Literal["pooled", "meta_fisher"] = args.replicate_method

    if rep_method == "pooled" or (len(case_ids) == 1 and rep_method == "meta_fisher"):
        case_arr = merged["count_case"].to_numpy(dtype=float)
        ctrl_arr = merged["count_ctrl"].to_numpy(dtype=float)
        for i, (kc, kt) in enumerate(zip(case_arr, ctrl_arr)):
            pvals[i] = fisher_pvalue_for_bin(float(kc), float(kt), total_case, total_ctrl)
            log2fc[i] = log2fc_for_bin(float(kc), float(kt), total_case, total_ctrl, eps)
        rep_method_used = "pooled"
    else:
        key = ["chrom", "start", "end"]
        case_aligned = merged[key].merge(case_mat[key + case_ids], on=key, how="left").fillna(0.0)
        ctrl_pooled = merged["count_ctrl"].to_numpy(dtype=float)

        case_totals = {sid: float(case_aligned[sid].sum()) for sid in case_ids}

        for i in range(len(merged)):
            kt = float(ctrl_pooled[i])
            per_rep_p: List[float] = []
            per_rep_fc: List[float] = []
            for sid in case_ids:
                tc = case_totals[sid]
                if tc <= 0:
                    continue
                kc = float(case_aligned.at[i, sid])
                per_rep_p.append(fisher_pvalue_for_bin(kc, kt, tc, total_ctrl))
                per_rep_fc.append(log2fc_for_bin(kc, kt, tc, total_ctrl, eps))

            pvals[i] = combine_pvalues_fisher(per_rep_p) if per_rep_p else 1.0
            log2fc[i] = float(np.median(per_rep_fc)) if per_rep_fc else 0.0

        rep_method_used = "meta_fisher"

    merged["pvalue"] = pvals
    merged["padj"] = bh_fdr(pvals)
    merged["log2FC"] = log2fc

    out_prefix = args.out_prefix
    out_tsv = f"{out_prefix}.tsv"

    merged_out = merged[["chrom", "start", "end", "count_case", "count_ctrl", "log2FC", "pvalue", "padj"]].copy()

    ebv_regex = str(args.ebv_regex) if args.ebv_regex is not None else ""
    if ebv_regex.strip():
        try:
            ebv_pat = re.compile(ebv_regex)
            merged_out["is_EBV"] = merged_out["chrom"].astype(str).apply(lambda c: bool(ebv_pat.search(c)))
        except re.error:
            raise SystemExit(f"ERROR: Invalid --ebv-regex pattern: {ebv_regex}")
    else:
        merged_out["is_EBV"] = False

    merged_out.to_csv(out_tsv, sep="\t", index=False)

    fdr_thr = float(args.fdr)
    sig_mask = merged_out["padj"] <= fdr_thr
    sig_total = int(sig_mask.sum())

    sig_out = merged_out.loc[sig_mask, :].copy()
    sig_all_tsv = f"{out_prefix}.sig.tsv"
    sig_up_tsv = f"{out_prefix}.sig_up.tsv"
    sig_down_tsv = f"{out_prefix}.sig_down.tsv"

    sig_out.to_csv(sig_all_tsv, sep="\t", index=False)
    sig_out[sig_out["log2FC"] > 0].to_csv(sig_up_tsv, sep="\t", index=False)
    sig_out[sig_out["log2FC"] < 0].to_csv(sig_down_tsv, sep="\t", index=False)

    sig_up = int((sig_out["log2FC"] > 0).sum())
    sig_down = int((sig_out["log2FC"] < 0).sum())

    # EBV metrics + enrichment
    n_tested = len(merged_out)
    ebv_bins_tested = int(merged_out["is_EBV"].sum())
    ebv_sig_bins = int((sig_mask & merged_out["is_EBV"]).sum())

    pct_ebv_tested = (ebv_bins_tested / n_tested * 100.0) if n_tested else float("nan")
    pct_ebv_sig = (ebv_sig_bins / sig_total * 100.0) if sig_total > 0 else float("nan")

    if ebv_bins_tested > 0 and sig_total > 0:
        # 2x2: EBV vs nonEBV by significant vs not significant
        a = ebv_sig_bins
        b = sig_total - ebv_sig_bins
        c = ebv_bins_tested - ebv_sig_bins
        non_ebv_tested = n_tested - ebv_bins_tested
        d = non_ebv_tested - b
        if d < 0:
            d = 0
        table = np.array([[a, b], [c, d]], dtype=np.int64)
        _, ebv_enrichment_p = fisher_exact(table, alternative="two-sided")
        ebv_enrichment_ratio = (pct_ebv_sig / pct_ebv_tested) if pct_ebv_tested > 0 else float("inf")
    else:
        ebv_enrichment_p = float("nan")
        ebv_enrichment_ratio = float("nan")

    subtitle = f"method={rep_method_used}; n_case={len(case_ids)}; n_ctrl={len(ctrl_ids)}"
    if regions_bed:
        subtitle += f"; regions={Path(regions_bed).name}"

    volcano_plot(merged_out, out_prefix, cname, fdr_thr, subtitle=subtitle)
    ma_plot(merged_out, out_prefix, cname, fdr_thr, subtitle=subtitle)

    with open(f"{out_prefix}.summary.txt", "w") as f:
        f.write(f"contrast\t{cname}\n")
        f.write(f"case\t{','.join(case_ids)}\n")
        f.write(f"control\t{','.join(ctrl_ids)}\n")
        f.write(f"replicate_method\t{rep_method_used}\n")
        f.write(f"n_case\t{len(case_ids)}\n")
        f.write(f"n_control\t{len(ctrl_ids)}\n")
        f.write(f"fdr_threshold\t{fdr_thr}\n")
        f.write(f"regions_bed\t{regions_bed if regions_bed else 'NA'}\n")
        f.write(f"ebv_regex\t{ebv_regex if ebv_regex.strip() else 'NA'}\n")
        f.write(f"ebv_bins_tested\t{ebv_bins_tested}\n")
        f.write(f"ebv_sig_bins\t{ebv_sig_bins}\n")
        f.write(f"ebv_pct_tested\t{pct_ebv_tested}\n")
        f.write(f"ebv_pct_sig\t{pct_ebv_sig}\n")
        f.write(f"ebv_enrichment_ratio\t{ebv_enrichment_ratio}\n")
        f.write(f"ebv_enrichment_pvalue\t{ebv_enrichment_p}\n")
        f.write(f"total_case\t{total_case}\n")
        f.write(f"total_ctrl\t{total_ctrl}\n")
        f.write(f"bins_tested\t{n_tested}\n")
        f.write(f"bins_FDR_le_threshold\t{sig_total}\n")
        f.write(f"sig_up\t{sig_up}\n")
        f.write(f"sig_down\t{sig_down}\n")


if __name__ == "__main__":
    main()
