#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import shutil
import time
from pathlib import Path
from typing import Optional, Sequence

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

from utils.seqio import normalize_input_sequences, ensure_dir
from utils.retrieval import run_remote


def run_single_job(
    seqs: list[str],
    out_dir: str,
    server: str,
    pairing: str,
    host_url: Optional[str],
    auth_user: Optional[str],
    auth_pass: Optional[str],
    email: str,
    user_agent: str,
    sleep_between_steps: float = 3.0,
) -> list[str]:
    server = server.lower()
    out_dir = ensure_dir(out_dir)
    if server == "protenix":
        # One-step
        return run_remote(
            seqs=seqs,
            out_dir=out_dir,
            server=server,
            host_url=host_url,
            pairing=pairing,
            user_agent=user_agent,
            email=email,
            auth_user=auth_user,
            auth_pass=auth_pass,
        )
    else:
        # ColabFold two-step: per-chain non-pairing then complex pairing
        # Step 1: per-chain non-pairing
        tmp_chain_dirs: list[str] = []
        for i, s in enumerate(seqs):
            chain_dir = ensure_dir(os.path.join(out_dir, f"tmp_chain_{i}"))
            run_remote(
                seqs=[s],
                out_dir=chain_dir,
                server="colabfold",
                host_url=host_url,
                pairing="off",
                user_agent=user_agent,
                email=email,
                auth_user=auth_user,
                auth_pass=auth_pass,
            )
            tmp_chain_dirs.append(chain_dir)
            time.sleep(sleep_between_steps)

        # Move non_pairing into final layout
        final_msa_root = ensure_dir(os.path.join(out_dir, "msa"))
        for i, chain_dir in enumerate(tmp_chain_dirs):
            src = os.path.join(chain_dir, "msa", "0")
            dst = os.path.join(final_msa_root, str(i))
            ensure_dir(dst)
            # copy only non_pairing.a3m if exists
            np_src = os.path.join(src, "non_pairing.a3m")
            if os.path.exists(np_src):
                shutil.copy2(np_src, os.path.join(dst, "non_pairing.a3m"))

        # Step 2: complex pairing
        run_remote(
            seqs=seqs,
            out_dir=out_dir,
            server="colabfold",
            host_url=host_url,
            pairing=pairing if pairing != "off" else "greedy",
            user_agent=user_agent,
            email=email,
            auth_user=auth_user,
            auth_pass=auth_pass,
        )
        return [os.path.join(out_dir, "msa", str(i)) for i in range(len(seqs))]


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Run a single MSA job (mono or multimer). For ColabFold, performs the recommended two-step (non-pairing per chain, then pairing for the complex)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--seq", type=str, help="Sequence string; use ':' for multimers")
    src.add_argument("--fasta", type=str, help="FASTA file with 1+ entries")
    p.add_argument("--out_dir", type=str, required=True)
    p.add_argument("--server", choices=["protenix", "colabfold"], default="protenix")
    p.add_argument("--pairing", choices=["off", "greedy", "complete"], default="greedy")
    p.add_argument("--host_url", type=str, default=None)
    p.add_argument("--auth_user", type=str, default=None)
    p.add_argument("--auth_pass", type=str, default=None)
    p.add_argument("--email", type=str, default="")
    p.add_argument("--user_agent", type=str, default="paired-msa-tools/0.1")
    p.add_argument("--sleep_between_steps", type=float, default=3.0)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    seqs = normalize_input_sequences(args.seq, args.fasta)
    out_dirs = run_single_job(
        seqs=seqs,
        out_dir=args.out_dir,
        server=args.server,
        pairing=args.pairing,
        host_url=args.host_url,
        auth_user=args.auth_user,
        auth_pass=args.auth_pass,
        email=args.email,
        user_agent=args.user_agent,
        sleep_between_steps=args.sleep_between_steps,
    )
    print("MSA directories:")
    for d in out_dirs:
        print(d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
