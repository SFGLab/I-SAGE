#!/usr/bin/env python3
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


def read_bed_regions(path: str) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Read BED (chrom, start, end) and build per-chrom interval index."""
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
    """Return boolean mask of bins overlapping any region in BED."""
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


def load_total_counts(plus_path: str, minus_path: str) -> pd.DataFrame:
    plus = read_bedgraph(plus_path)
    minus = read_bedgraph(minus_path)

    key = ["chrom", "start", "end"]
    merged = plus.merge(minus, on=key, how="outer", suffixes=("_plus", "_minus"))
    merged["value_plus"] = merged["value_plus"].fillna(0.0)
    merged["value_minus"] = merged["value_minus"].fillna(0.0)
    merged["count"] = merged["value_plus"] + merged["value_minus"]
    return merged[key + ["count"]]


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
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


def group_sum(samples: Dict[str, Tuple[str, str]], ids: List[str]) -> pd.DataFrame:
    acc: pd.DataFrame | None = None
    for sid in ids:
        plus, minus = samples[sid]
        df = load_total_counts(str(plus), str(minus)).rename(columns={"count": sid})
        acc = df if acc is None else acc.merge(df, on=["chrom", "start", "end"], how="outer")
    assert acc is not None
    for sid in ids:
        acc[sid] = acc[sid].fillna(0.0)
    acc["count"] = acc[ids].sum(axis=1)
    return acc[["chrom", "start", "end", "count"]]


def fisher_per_bin(case_counts: np.ndarray, ctrl_counts: np.ndarray, eps: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    total_case = float(case_counts.sum())
    total_ctrl = float(ctrl_counts.sum())
    if total_case <= 0 or total_ctrl <= 0:
        raise ValueError("total_case or total_ctrl <= 0")

    pvals = np.empty(case_counts.size, dtype=float)
    log2fc = np.empty(case_counts.size, dtype=float)

    for i, (kc, kt) in enumerate(zip(case_counts, ctrl_counts)):
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

    padj = bh_fdr(pvals)
    return pvals, padj, log2fc


def downsample_counts(counts: np.ndarray, frac: float, rng: np.random.Generator) -> np.ndarray:
    counts_int = np.round(counts).astype(np.int64)
    if frac >= 1.0:
        return counts_int.astype(float)
    return rng.binomial(counts_int, frac).astype(float)


def run_downsample(
    base: pd.DataFrame,
    frac_list: List[float],
    fdr: float,
    reps: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out_rows = []
    case_base = base["count_case"].to_numpy(dtype=float)
    ctrl_base = base["count_ctrl"].to_numpy(dtype=float)

    _, q_full, _ = fisher_per_bin(case_base, ctrl_base)
    sig_full = (q_full <= fdr)

    for frac in frac_list:
        for r in range(reps):
            case_ds = downsample_counts(case_base, frac, rng)
            ctrl_ds = downsample_counts(ctrl_base, frac, rng)
            _, q, _ = fisher_per_bin(case_ds, ctrl_ds)
            sig = (q <= fdr)

            inter = int(np.logical_and(sig, sig_full).sum())
            union = int(np.logical_or(sig, sig_full).sum())
            jacc = (inter / union) if union > 0 else 1.0

            out_rows.append(
                {
                    "fraction": frac,
                    "replicate": r,
                    "sig_bins": int(sig.sum()),
                    "sig_bins_full": int(sig_full.sum()),
                    "jaccard_vs_full": jacc,
                }
            )
    return pd.DataFrame(out_rows)


def run_spikein(
    base: pd.DataFrame,
    n_bins: int,
    effect_mult: float,
    fdr: float,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    ctrl = base["count_ctrl"].to_numpy(dtype=float)

    case = ctrl.copy()

    n_bins = min(n_bins, case.size)
    spike_idx = rng.choice(case.size, size=n_bins, replace=False)

    case[spike_idx] = np.round(case[spike_idx] * effect_mult)

    _, q, _ = fisher_per_bin(case, ctrl)
    called = (q <= fdr)

    tp = int(called[spike_idx].sum())
    fp = int(called.sum() - tp)
    fn = int(n_bins - tp)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 1.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    return pd.DataFrame(
        [
            {
                "spike_bins": n_bins,
                "effect_mult": effect_mult,
                "fdr_thresh": fdr,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": precision,
                "recall": recall,
            }
        ]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples-tsv", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out-prefix", required=True)

    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--downsample-fracs", default="1.0,0.5,0.25,0.1")
    ap.add_argument("--downsample-reps", type=int, default=1)
    ap.add_argument("--spikein-bins", type=int, default=1000)
    ap.add_argument("--spikein-mult", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument(
        "--regions-bed",
        default=None,
        help="Optional BED file. If provided, validation operates only on bins overlapping these regions.",
    )
    args = ap.parse_args()

    cname, case_ids, ctrl_ids = parse_contrast(args.contrast)

    samp = pd.read_csv(args.samples_tsv, sep="\t", header=None, names=["sample_id", "plus_bg", "minus_bg"])
    samp_map: Dict[str, Tuple[str, str]] = {
        str(r.sample_id): (str(r.plus_bg), str(r.minus_bg))
        for r in samp.itertuples(index=False)
    }

    case_df = group_sum(samp_map, case_ids).rename(columns={"count": "count_case"})
    ctrl_df = group_sum(samp_map, ctrl_ids).rename(columns={"count": "count_ctrl"})
    base = case_df.merge(ctrl_df, on=["chrom", "start", "end"], how="outer")
    base["count_case"] = base["count_case"].fillna(0.0)
    base["count_ctrl"] = base["count_ctrl"].fillna(0.0)

    regions_bed = str(args.regions_bed) if args.regions_bed else None
    if regions_bed:
        regions = read_bed_regions(regions_bed)
        keep = bins_overlap_regions(base[["chrom", "start", "end"]], regions)
        base = base.loc[keep].reset_index(drop=True)
        if len(base) == 0:
            raise SystemExit(f"ERROR: No bins overlap regions BED: {regions_bed}")

    fracs = [float(x) for x in args.downsample_fracs.split(",") if x.strip()]

    down_df = run_downsample(base, fracs, args.fdr, args.downsample_reps, args.seed)
    down_df.to_csv(f"{args.out_prefix}.downsample.tsv", sep="\t", index=False)

    plt.figure()
    g = down_df.groupby("fraction")["sig_bins"].mean().reset_index()
    plt.plot(g["fraction"], g["sig_bins"], marker="o")
    plt.xlabel("Downsample fraction")
    plt.ylabel(f"Mean significant bins (FDR ≤ {args.fdr})")
    plt.title(f"Downsample power curve: {cname}")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.downsample.png", dpi=200)
    plt.close()

    spike_df = run_spikein(base, args.spikein_bins, args.spikein_mult, args.fdr, args.seed)
    spike_df.to_csv(f"{args.out_prefix}.spikein.tsv", sep="\t", index=False)

    with open(f"{args.out_prefix}.validation.summary.txt", "w") as f:
        f.write(f"contrast\t{cname}\n")
        f.write(f"regions_bed\t{regions_bed if regions_bed else 'NA'}\n")
        f.write(f"fdr\t{args.fdr}\n")
        f.write(f"downsample_fracs\t{args.downsample_fracs}\n")
        f.write(f"downsample_reps\t{args.downsample_reps}\n")
        f.write(f"spikein_bins\t{args.spikein_bins}\n")
        f.write(f"spikein_mult\t{args.spikein_mult}\n")
        f.write(f"regions_bed\t{regions_bed if regions_bed else 'NA'}\n")


if __name__ == "__main__":
    main()