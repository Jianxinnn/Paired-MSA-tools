#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from utils.vis import compute_msa_depth, plot_msa_depth


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Visualize per-position MSA abundance (depth). "
            "Supports single-chain A3M and stitched paired A3M (auto-splits chains if query header has _lenN).")
    )
    p.add_argument("a3m", type=str, help="Path to A3M file (single or stitched paired)")
    p.add_argument("--out", type=str, default=None, help="Output image path (.png/.pdf). Default: <a3m>.depth.png")
    p.add_argument("--include_query", action="store_true", help="Include the first query row in depth counts")
    p.add_argument("--normalize", action="store_true", help="Normalize by max depth to [0,1]")
    p.add_argument("--smooth", type=int, default=1, help="Moving-average window for smoothing")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--title", type=str, default=None)
    p.add_argument("--chain_lengths", type=str, default=None, help="Override chain lengths, e.g., '569,471'")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_argparser().parse_args(argv)
    a3m_path = args.a3m
    res = compute_msa_depth(a3m_path, include_query=args.include_query)
    depth = res["depth"]
    chain_lengths = res["chain_lengths"]
    if args.chain_lengths:
        try:
            chain_lengths = [int(x) for x in args.chain_lengths.split(",") if x]
        except Exception as e:
            raise SystemExit(f"Invalid --chain_lengths: {e}")

    out_path = args.out
    if out_path is None:
        out_path = str(Path(a3m_path).with_suffix("") ) + ".depth.png"

    title = args.title
    if not title:
        title = Path(a3m_path).name

    plot_msa_depth(
        depth,
        chain_lengths=chain_lengths,
        title=title,
        out_path=out_path,
        dpi=args.dpi,
        normalize=args.normalize,
        smooth=args.smooth,
        include_query=args.include_query,
    )

    n_total = res["n_rows_total"]
    n_used = res["n_rows_used"]
    n_skip = res["n_rows_skipped"]
    L = res["L"]
    k = 1 if chain_lengths is None else len(chain_lengths)
    print(f"MSA depth plotted -> {out_path}")
    print(f"Aligned length L={L}, chains={k}; rows total={n_total}, used={n_used}, skipped={n_skip}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

