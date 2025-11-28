#!/usr/bin/env python3

import argparse
from pathlib import Path
from collections import defaultdict

import pysam

# -----------------------------------------------------------
# Known SgrDI recognition sequence (palindromic)
# -----------------------------------------------------------
S_GRDI_SEQ = "CGTCGACG"  # 8 bp

# -----------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="I-BLESS Break Caller: Extract strand-aware break positions from BAM."
    )

    parser.add_argument("--bam", required=True, help="Input BAM file (sorted, indexed).")
    parser.add_argument("--out", required=True, help="Output BED file.")
    parser.add_argument(
        "--mode",
        choices=["per_base", "binned"],
        default="per_base",
        help="Counting mode: per_base or binned.",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=1000,
        help="Bin size when mode == binned.",
    )
    parser.add_argument(
        "--remove-sgrdi",
        action="store_true",
        help="Remove breaks overlapping SgrDI recognition motif.",
    )
    parser.add_argument(
        "--reference-fasta",
        help="Reference genome FASTA (required if --remove-sgrdi).",
    )

    return parser.parse_args()

# -----------------------------------------------------------
# Helpers
# -----------------------------------------------------------

def get_break_position(aln: pysam.AlignedSegment):
    """
    Compute genomic break coordinate and strand for an alignment.

    For I-BLESS:
      - A DSB end is tagged by R1.
      - Forward strand (not reverse): break = reference_start (0-based).
      - Reverse strand: break = reference_end - 1 (0-based, inclusive).
    """
    if aln.is_reverse:
        strand = "-"
        break_pos = aln.reference_end - 1 if aln.reference_end is not None else 0
    else:
        strand = "+"
        break_pos = aln.reference_start if aln.reference_start is not None else 0

    return break_pos, strand


def sgrdi_at_position(ref_fa: pysam.FastaFile, chrom: str, pos: int) -> bool:
    """
    Return True if the SgrDI motif is present at chrom:pos (0-based).
    We assume reference is 0-based coordinate; pysam.fetch uses 0-based, end-exclusive.
    """
    motif_len = len(S_GRDI_SEQ)
    if pos < 0:
        return False

    # Fetch sequence from reference: [pos, pos+motif_len)
    seq = ref_fa.fetch(chrom, pos, pos + motif_len)
    return seq == S_GRDI_SEQ

# -----------------------------------------------------------
# Core logic
# -----------------------------------------------------------

def process_bam(args):
    bam_path = args.bam
    bin_size = args.bin_size

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Reference FASTA (only needed if we remove SgrDI sites)
    ref_fa = None
    if args.remove_sgrdi:
        if not args.reference_fasta:
            raise ValueError(
                "--remove-sgrdi was set but --reference-fasta was not provided."
            )
        ref_fa = pysam.FastaFile(args.reference_fasta)

    # counts[chrom][pos_or_bin][strand] = count
    counts = defaultdict(lambda: defaultdict(lambda: {"+": 0, "-": 0}))

    for aln in bam.fetch(until_eof=True):
        # Skip unwanted reads
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue

        chrom = bam.get_reference_name(aln.reference_id)
        break_pos, strand = get_break_position(aln)

        # SgrDI filtering (if enabled)
        if args.remove_sgrdi and ref_fa is not None:
            if sgrdi_at_position(ref_fa, chrom, break_pos):
                continue

        # Choose coordinate key (per-base vs binned)
        if args.mode == "per_base":
            key_pos = break_pos
        else:
            # Bin start: floor(break_pos / bin_size) * bin_size
            key_pos = (break_pos // bin_size) * bin_size

        counts[chrom][key_pos][strand] += 1

    bam.close()
    if ref_fa is not None:
        ref_fa.close()

    return counts

# -----------------------------------------------------------
# BED writer
# -----------------------------------------------------------

def write_bed(counts, out_path: Path, mode: str, bin_size: int):
    """
    Write results as BED with strand-aware counts.

    Columns:
      chrom, start, end, name, score, strand
    Where:
      - name encodes break_+ or break_-
      - score is the count at that position/bin
    """
    with out_path.open("w") as out:
        for chrom in sorted(counts.keys()):
            for pos in sorted(counts[chrom].keys()):
                strand_counts = counts[chrom][pos]

                for strand in ("+", "-"):
                    count = strand_counts[strand]
                    if count == 0:
                        continue

                    if mode == "per_base":
                        start = pos
                        end = pos + 1
                    else:
                        start = pos
                        end = pos + bin_size

                    name = f"break_{strand}"
                    score = count

                    out.write(
                        f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
                    )

# -----------------------------------------------------------
# Main
# -----------------------------------------------------------

def main():
    args = parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    counts = process_bam(args)
    write_bed(counts, out_path, args.mode, args.bin_size)

    print(f"[I-SAGE] Break calling complete → {out_path}")

if __name__ == "__main__":
    main()
