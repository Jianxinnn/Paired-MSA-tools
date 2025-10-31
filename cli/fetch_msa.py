#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional, Sequence

# Make 'fold' importable when running as a plain script
BASE_DIR = Path(__file__).resolve().parent.parent
PARENT_DIR = BASE_DIR.parent
sys.path.insert(0, str(PARENT_DIR))

from fold.utils.seqio import normalize_input_sequences, ensure_dir
from fold.utils.retrieval import run_remote


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Fetch MSA remotely (Protenix or ColabFold) and write per-chain pairing/non_pairing A3M files."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--seq", type=str, help="Sequence string; use ':' to join chains for multimers")
    src.add_argument("--fasta", type=str, help="FASTA file with 1+ entries")

    p.add_argument("--out_dir", type=str, required=True, help="Output root directory")
    p.add_argument("--server", type=str, default="protenix", choices=["protenix", "colabfold"], help="Remote server backend")
    p.add_argument("--host_url", type=str, default=None, help="Override server host URL")
    p.add_argument("--pairing", type=str, default="greedy", choices=["off", "greedy", "complete"], help="Pairing mode for remote backend")
    p.add_argument("--user_agent", type=str, default="paired-msa-tools/0.1", help="HTTP User-Agent header")
    p.add_argument("--email", type=str, default="", help="Contact email for the server")
    p.add_argument("--auth_user", type=str, default=None, help="BasicAuth username (Protenix)")
    p.add_argument("--auth_pass", type=str, default=None, help="BasicAuth password (Protenix)")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    seqs = normalize_input_sequences(args.seq, args.fasta)
    out_dir = ensure_dir(args.out_dir)
    subdirs = run_remote(
        seqs=seqs,
        out_dir=out_dir,
        server=args.server,
        host_url=args.host_url,
        pairing=args.pairing,
        user_agent=args.user_agent,
        email=args.email,
        auth_user=args.auth_user,
        auth_pass=args.auth_pass,
    )
    print("MSA directories (per-chain precomputed_msa_dir):")
    for d in subdirs:
        print(d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
