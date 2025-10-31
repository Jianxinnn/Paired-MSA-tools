from __future__ import annotations

import os
import time
from typing import List, Optional, Sequence

from tqdm import tqdm

from .client import RemoteClient
from .seqio import ensure_dir
from .postprocess import (
    read_m8,
    protenix_split_to_pairing_and_nonpairing,
    protenix_nonpairing_only,
    colabfold_post_nonpair,
    colabfold_post_pair,
    extract_tarball_to_dir,
)


def run_remote(
    seqs: Sequence[str],
    out_dir: str,
    server: str,
    host_url: Optional[str],
    pairing: str,
    user_agent: str,
    email: str,
    auth_user: Optional[str],
    auth_pass: Optional[str],
) -> List[str]:
    server = server.lower()
    if host_url is None:
        host_url = (
            "https://protenix-server.com/api/msa" if server == "protenix" else "https://api.colabfold.com"
        )

    client = RemoteClient(
        host_url=host_url, user_agent=user_agent, email=email, auth_user=auth_user, auth_pass=auth_pass
    )

    # Build query
    lines = []
    for i, s in enumerate(seqs):
        lines.append(f">query_{i}\n")
        lines.append(f"{s}\n")
    if lines and lines[-1].endswith("\n"):
        lines[-1] = lines[-1].rstrip("\n")
    query = "".join(lines)

    # Endpoint/mode
    want_pair = pairing in ("greedy", "complete")
    if server == "protenix":
        endpoint = "ticket/msa"
        mode = "env"
        want_pair = False
    else:
        if want_pair and len(seqs) > 1:
            endpoint = "ticket/pair"
            mode = "pairgreedy" if pairing == "greedy" else "paircomplete"
        else:
            endpoint = "ticket/msa"
            mode = "env"

    # Submit/poll
    with tqdm(total=100, bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed}]") as pbar:
        pbar.set_description("SUBMIT")
        out = client.submit(query, endpoint, mode)
        while out.get("status") in {"UNKNOWN", "RATELIMIT"}:
            time.sleep(60)
            out = client.submit(query, endpoint, mode)
        if out.get("status") in {"ERROR", "MAINTENANCE"}:
            raise RuntimeError(f"Server status: {out}")
        job_id = out.get("id")
        pbar.set_description(out.get("status", "PENDING"))
        elapsed = 0
        while out.get("status") in {"UNKNOWN", "RUNNING", "PENDING"}:
            time.sleep(60)
            out = client.status(job_id)
            pbar.set_description(out.get("status", "RUNNING"))
            elapsed += 60
            pbar.n = min(99, int(100 * elapsed / (30 * 60)))
            pbar.refresh()
        pbar.n = 100
        pbar.refresh()

    if out.get("status") != "COMPLETE":
        raise RuntimeError(f"Job did not complete successfully: {out}")

    # Download & extract
    blob = client.download(out["id"])  # .tar.gz
    raw_root = extract_tarball_to_dir(blob, out_dir)

    # Post-process to per-chain files
    final_root = ensure_dir(os.path.join(out_dir, "msa"))
    subdirs: List[str] = []

    if server == "protenix":
        m8_file = os.path.join(raw_root, "uniref_tax.m8")
        has_m8 = os.path.exists(m8_file)
        tax_map = read_m8(m8_file) if has_m8 else {}
        for i, seq in enumerate(seqs):
            seq_dir = ensure_dir(os.path.join(final_root, str(i)))
            a3m_path = os.path.join(raw_root, f"{i}.a3m")
            if os.path.exists(a3m_path):
                if has_m8:
                    protenix_split_to_pairing_and_nonpairing(seq, seq_dir, a3m_path, tax_map)
                else:
                    protenix_nonpairing_only(seq, seq_dir, a3m_path)
                subdirs.append(seq_dir)
            else:
                # write minimal stubs
                open(os.path.join(seq_dir, "non_pairing.a3m"), "w").write(f">query\n{seq}\n")
                open(os.path.join(seq_dir, "pairing.a3m"), "w").write(f">query\n{seq}\n")
                subdirs.append(seq_dir)

    else:
        if endpoint == "ticket/msa":
            for i in range(len(seqs)):
                header = f"query_{i}"
                seq_dir = colabfold_post_nonpair(raw_root, header, final_root)
                subdirs.append(seq_dir)
        else:
            subdirs = colabfold_post_pair(raw_root, final_root)
            for i, seq in enumerate(seqs):
                seq_dir = ensure_dir(os.path.join(final_root, str(i)))
                nonpair = os.path.join(seq_dir, "non_pairing.a3m")
                if not os.path.exists(nonpair):
                    open(nonpair, "w").write(f">query\n{seq}\n")

    return subdirs
