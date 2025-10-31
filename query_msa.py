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

# Make 'fold' importable as a package
BASE_DIR = Path(__file__).resolve().parent
PARENT_DIR = BASE_DIR.parent
sys.path.insert(0, str(PARENT_DIR))

from fold.utils.seqio import normalize_input_sequences, ensure_dir
from fold.utils.retrieval import run_remote


@dataclass
class RateLimiter:
    min_interval_s: float = 5.0
    _last_ts: float = 0.0

    def wait(self):
        now = time.time()
        dt = now - self._last_ts
        if dt < self.min_interval_s:
            time.sleep(self.min_interval_s - dt)
        self._last_ts = time.time()


def run_single_job(
    seqs: List[str],
    out_dir: str,
    server: str,
    pairing: str,
    host_url: Optional[str],
    auth_user: Optional[str],
    auth_pass: Optional[str],
    email: str,
    user_agent: str,
    two_step: bool = False,
    sleep_between_steps: float = 3.0,
) -> List[str]:
    server = server.lower()
    out_dir = ensure_dir(out_dir)
    if server == "colabfold" and two_step and len(seqs) > 1:
        # Step 1: per-chain non-pairing
        tmp_chain_dirs: List[str] = []
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
        for i, chain_dir in enumerate(tmp_chain_dirs):
            src = Path(chain_dir) / "msa" / "0" / "non_pairing.a3m"
            dst = Path(out_dir) / "msa" / str(i)
            ensure_dir(str(dst))
            if src.exists():
                (dst / "non_pairing.a3m").write_text(src.read_text())

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
        return [str(Path(out_dir) / "msa" / str(i)) for i in range(len(seqs))]

    # Default path (Protenix one-step; ColabFold single-step)
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


def run_batch(
    input_path: str,
    out_root: str,
    server: str,
    pairing: str,
    host_url: Optional[str],
    auth_user: Optional[str],
    auth_pass: Optional[str],
    email: str,
    user_agent: str,
    delimiter: str = ",",
    min_interval_s: float = 5.0,
    sleep_between_jobs: float = 2.0,
) -> None:
    limiter = RateLimiter(min_interval_s=min_interval_s)
    ensure_dir(out_root)
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t" if delimiter == "\t" else ",")
        jobs = list(reader)
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
        # Detect sequence list
        is_fasta = os.path.isfile(payload)
        seqs = normalize_input_sequences(None if is_fasta else payload, payload if is_fasta else None)
        out_dir = os.path.join(out_root, name)
        limiter.wait()
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
        time.sleep(sleep_between_jobs)
    print("[batch] done.")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Simple entry point for remote MSA fetching (single or batch).\n"
            "- Single job: provide --seq or --fasta.\n"
            "- Batch job: provide --input CSV/TSV with name,seq or name,fasta.\n"
            "For ColabFold multimers, use --two_step to run per-chain non-pairing then complex pairing."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--seq", type=str, help="Sequence string; use ':' to join chains for multimers")
    mode.add_argument("--fasta", type=str, help="FASTA file with 1+ entries")
    mode.add_argument("--input", type=str, help="CSV/TSV for batch: columns name,seq or name,fasta")

    p.add_argument("--out_dir", type=str, required=True, help="Output directory (single) or root (batch)")
    p.add_argument("--server", choices=["protenix", "colabfold"], default="protenix")
    p.add_argument("--pairing", choices=["off", "greedy", "complete"], default="greedy")
    p.add_argument("--host_url", default=None)
    p.add_argument("--auth_user", default=None)
    p.add_argument("--auth_pass", default=None)
    p.add_argument("--email", default="")
    p.add_argument("--user_agent", default="paired-msa-tools/0.1")
    # ColabFold helper
    p.add_argument("--two_step", action="store_true", help="For ColabFold multimers: run per-chain non-pairing then complex pairing")
    # Batch-only
    p.add_argument("--delimiter", default=",", help="Batch input delimiter (',' or '\t')")
    p.add_argument("--min_interval_s", type=float, default=5.0)
    p.add_argument("--sleep_between_jobs", type=float, default=2.0)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.input:
        run_batch(
            input_path=args.input,
            out_root=args.out_dir,
            server=args.server,
            pairing=args.pairing,
            host_url=args.host_url,
            auth_user=args.auth_user,
            auth_pass=args.auth_pass,
            email=args.email,
            user_agent=args.user_agent,
            delimiter=args.delimiter,
            min_interval_s=args.min_interval_s,
            sleep_between_jobs=args.sleep_between_jobs,
        )
        return 0

    # Single job
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
        two_step=args.two_step,
    )
    print("MSA directories:")
    for d in out_dirs:
        print(d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
