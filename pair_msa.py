#!/usr/bin/env python3
"""
Standalone species-based pairing tool for multimer MSAs.

Given per-chain pairing-candidate A3M files (e.g., query_msa.py produced
`msa/<i>/pairing.a3m` for each chain), this script:
  1) parses headers to extract species IDs (UniRef100_*_<taxid>/ pattern),
  2) computes per-entry similarity to the chain query,
  3) for each species, ranks entries by similarity within each chain, and
  4) pairs them across chains by rank (min count across present chains),
  5) outputs a concatenated, paired A3M (`paired_concat.a3m`).

The concatenated A3M places chain segments back-to-back along columns,
forming one multimer alignment, with one `>query` row equals to the
concatenation of per-chain queries.

No dependency on the Protenix repository.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


HEADER_RE_UNIREF100 = re.compile(r"^UniRef100_[^_]+_(\d+)(?:/|\b)")


@dataclass
class A3MEntry:
    header: str
    raw: str            # with lowercase insertions kept
    aligned: str        # deletion-table applied (lowercase removed)
    similarity: float   # to query
    species: Optional[str]


@dataclass
class ChainA3M:
    query_raw: str
    query_aligned: str
    entries: List[A3MEntry]


def _strip_inserts(seq: str) -> str:
    # Remove lowercase insertions to get the aligned sequence
    return "".join(ch for ch in seq if not ch.islower())


def _species_from_header(header: str) -> Optional[str]:
    # Strip leading '>' and take first token
    name = header[1:].split()[0] if header.startswith(">") else header.split()[0]
    m = HEADER_RE_UNIREF100.match(name)
    if m:
        return m.group(1)  # taxid
    return None


def _parse_a3m(path: str) -> Tuple[str, str, List[Tuple[str, str]]]:
    heads: List[str] = []
    seqs: List[str] = []
    with open(path, "r") as f:
        cur_head = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                heads.append(line)
                cur_head = line
            else:
                if cur_head is None:
                    raise ValueError(f"Sequence line without header in {path}")
                seqs.append(line)
    if not heads:
        raise ValueError(f"Empty A3M: {path}")
    # Query is the first entry whose header starts with '>query'
    q_idx = next((i for i, h in enumerate(heads) if h[1:].lower().startswith("query")), 0)
    query_raw = seqs[q_idx]
    entries = list(zip(heads, seqs))
    return query_raw, _strip_inserts(query_raw), entries


def _build_chain(path: str) -> ChainA3M:
    query_raw, query_aligned, entries = _parse_a3m(path)
    out_entries: List[A3MEntry] = []
    # Precompute aligned for query once
    q = query_aligned
    for h, s in entries:
        aln = _strip_inserts(s)
        # Defensive: lengths should match; if not, skip
        if len(aln) != len(q):
            continue
        # simple similarity: fraction of exact matches over aligned length
        sim = sum(1 for a, b in zip(aln, q) if a == b) / float(len(q))
        sp = _species_from_header(h)
        out_entries.append(A3MEntry(header=h, raw=s, aligned=aln, similarity=sim, species=sp))
    return ChainA3M(query_raw=query_raw, query_aligned=query_aligned, entries=out_entries)


def _concat_rows(rows_per_chain: List[Optional[A3MEntry]], pad_lengths: List[int]) -> Tuple[str, str]:
    # Build combined header and sequence by concatenation in chain order
    headers = []
    seqs = []
    for e, L in zip(rows_per_chain, pad_lengths):
        if e is None:
            headers.append("PAD")
            seqs.append("-" * L)
        else:
            # keep original header token for transparency
            headers.append(e.header[1:].split()[0])
            seqs.append(e.raw)
    combo_header = ">paired|" + "|".join(headers)
    combo_seq = "".join(seqs)
    return combo_header, combo_seq


def pair_by_species(a3m_paths: Sequence[str], out_path: str, min_chains: int = 2, max_rows_per_species: int = 600) -> int:
    # Load all chains
    chains = [ _build_chain(p) for p in a3m_paths ]
    n_chain = len(chains)
    pad_lengths = [len(c.query_aligned) for c in chains]
    # Union of species present across chains
    species_set = set()
    per_chain_species: List[Dict[str, List[A3MEntry]]] = []
    for c in chains:
        d: Dict[str, List[A3MEntry]] = {}
        for e in c.entries:
            if e.species is None:
                continue
            d.setdefault(e.species, []).append(e)
        per_chain_species.append(d)
        species_set.update(d.keys())

    # Prepare output
    lines: List[str] = []
    # First write concatenated query
    lines.append(">query\n")
    lines.append("".join(c.query_raw for c in chains) + "\n")

    total_rows = 0
    for sp in sorted(species_set):
        present = [i for i, d in enumerate(per_chain_species) if sp in d]
        if len(present) < min_chains:
            continue  # require species appears in at least `min_chains` chains

        # Row limit guard: if any chain has > max_rows_per_species for this species, skip
        if any(len(per_chain_species[i][sp]) > max_rows_per_species for i in present):
            continue

        # Build ranked lists by similarity within each chain
        ranked: List[List[A3MEntry]] = []
        for i in range(n_chain):
            lst = per_chain_species[i].get(sp, [])
            lst_sorted = sorted(lst, key=lambda e: e.similarity, reverse=True)
            ranked.append(lst_sorted)

        # Number of pairs = min count among present chains
        k = min(len(ranked[i]) for i in present)
        for r in range(k):
            row_entries: List[Optional[A3MEntry]] = []
            for i in range(n_chain):
                if i in present:
                    row_entries.append(ranked[i][r])
                else:
                    row_entries.append(None)  # pad missing chain with gaps
            h, s = _concat_rows(row_entries, pad_lengths)
            lines.append(h + f"\tsp={sp}\trank={r+1}\n")
            lines.append(s + "\n")
            total_rows += 1

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for ln in lines:
            f.write(ln)
    return total_rows


def _pair_by_index(a3m_paths: Sequence[str], out_path: str) -> int:
    chains = [ _build_chain(p) for p in a3m_paths ]
    pad_lengths = [len(c.query_aligned) for c in chains]
    # write concatenated query
    lines: List[str] = []
    lines.append(">query\n")
    lines.append("".join(c.query_raw for c in chains) + "\n")
    # build non-query lists preserving order of appearance
    non_query_per_chain: List[List[A3MEntry]] = []
    for c in chains:
        lst = [e for e in c.entries if not e.header[1:].lower().startswith("query")]  # keep order
        non_query_per_chain.append(lst)
    k = min(len(lst) for lst in non_query_per_chain)
    total = 0
    for r in range(k):
        row_entries = [ non_query_per_chain[i][r] for i in range(len(chains)) ]
        h, s = _concat_rows(row_entries, pad_lengths)
        lines.append(h + f"\tidx={r+1}\n")
        lines.append(s + "\n")
        total += 1
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for ln in lines:
            f.write(ln)
    return total


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Pair multimer MSAs by species and output a concatenated paired A3M.\n"
            "Input: per-chain pairing A3M files in chain order (e.g., msa/0/pairing.a3m msa/1/pairing.a3m)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("a3m", nargs='+', help="Per-chain pairing A3M files in chain order")
    p.add_argument("--out", type=str, default="paired_concat.a3m", help="Output paired A3M path")
    p.add_argument("--strategy", choices=["species","index"], default="species",
                   help="Pairing strategy: by species (UniRef100 taxid) or by row index (ColabFold pair semantics)")
    p.add_argument("--min_chains", type=int, default=2, help="Minimum chains a species must appear in to be paired")
    p.add_argument("--max_rows_per_species", type=int, default=600, help="Skip species if any chain exceeds this count")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_argparser().parse_args(argv)
    if args.strategy == "species":
        n = pair_by_species(args.a3m, args.out, args.min_chains, args.max_rows_per_species)
    else:
        n = _pair_by_index(args.a3m, args.out)
    print(f"Wrote {n} paired rows to {args.out} using {args.strategy}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
