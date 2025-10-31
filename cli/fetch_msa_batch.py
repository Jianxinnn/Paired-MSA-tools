#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

from utils.seqio import ensure_dir, normalize_input_sequences
from utils.retrieval import run_remote


@dataclass
class RateLimiter:
    min_interval_s: float = 5.0  # at most 1 submit per 5 seconds by default
    _last_ts: float = 0.0

    def wait(self):
        now = time.time()
        dt = now - self._last_ts
        if dt < self.min_interval_s:
            time.sleep(self.min_interval_s - dt)
        self._last_ts = time.time()


def process_job(
    name: str,
    seq_or_fasta: str,
    server: str,
    pairing: str,
    out_root: str,
    host_url: Optional[str],
    auth_user: Optional[str],
    auth_pass: Optional[str],
    email: str,
    user_agent: str,
    limiter: RateLimiter,
) -> None:
    out_dir = ensure_dir(os.path.join(out_root, name))
    # Parse sequences (colon-separated string or FASTA path auto-detected)
    is_fasta = os.path.isfile(seq_or_fasta)
    seqs = normalize_input_sequences(None if is_fasta else seq_or_fasta, seq_or_fasta if is_fasta else None)
    limiter.wait()  # rate limit submissions
    run_remote(
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


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Batch MSA retrieval with built-in rate limiting (safe defaults for public servers)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, help="CSV/TSV with columns: name,seq or name,fasta")
    p.add_argument("--delimiter", default=",", help="Field delimiter (',' or '\t')")
    p.add_argument("--out_root", required=True, help="Root output directory for all jobs")
    p.add_argument("--server", choices=["protenix", "colabfold"], default="protenix")
    p.add_argument("--pairing", choices=["off", "greedy", "complete"], default="greedy")
    p.add_argument("--host_url", default=None)
    p.add_argument("--auth_user", default=None)
    p.add_argument("--auth_pass", default=None)
    p.add_argument("--email", default="")
    p.add_argument("--user_agent", default="paired-msa-tools/0.1")
    p.add_argument("--min_interval_s", type=float, default=5.0, help="Min seconds between job submissions")
    p.add_argument("--sleep_between_jobs", type=float, default=2.0, help="Extra sleep between jobs to be gentle")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    limiter = RateLimiter(min_interval_s=args.min_interval_s)
    ensure_dir(args.out_root)
    # Read jobs
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t" if args.delimiter == "\t" else ",")
        jobs: List[dict] = list(reader)
    if not jobs:
        raise RuntimeError("No jobs found in input")
    for job in jobs:
        name = job.get("name")
        seq = job.get("seq")
        fasta = job.get("fasta")
        if not name:
            continue
        if not seq and not fasta:
            continue
        payload = seq if seq else fasta
        print(f"[batch] processing {name}...")
        process_job(
            name=name,
            seq_or_fasta=payload,
            server=args.server,
            pairing=args.pairing,
            out_root=args.out_root,
            host_url=args.host_url,
            auth_user=args.auth_user,
            auth_pass=args.auth_pass,
            email=args.email,
            user_agent=args.user_agent,
            limiter=limiter,
        )
        time.sleep(args.sleep_between_jobs)
    print("[batch] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
