#!/usr/bin/env python

import argparse
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mode", choices=["per_base", "binned"], default="per_base")
    p.add_argument("--bin-size", type=int, default=1000)
    p.add_argument("--remove-sgrdi", action="store_true")
    return p.parse_args()

def main():
    args = parse_args()
    out = Path(args.out)
    # TODO: implement real logic using pysam
    # For now, just create an empty BED with a header-like comment
    out.write_text("# placeholder BED; break calling not yet implemented\n")

if __name__ == "__main__":
    main()
