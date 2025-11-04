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
import logging
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from utils import stitch as stitch_utils


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


def _format_sim_annotation(rows_per_chain: List[Optional[A3MEntry]]) -> str:
    """Format per-chain similarities and their average as a header suffix.

    Missing chains are shown as '-'. Values are rounded to 3 decimals.
    """
    vals: List[str] = []
    sims: List[float] = []
    for e in rows_per_chain:
        if e is None:
            vals.append('-')
        else:
            vals.append(f"{e.similarity:.3f}")
            sims.append(float(e.similarity))
    avg = (sum(sims) / len(sims)) if sims else 0.0
    return f"\tsim={','.join(vals)}\tsim_avg={avg:.3f}"


def pair_by_species(a3m_paths: Sequence[str], out_path: str, min_chains: int = 2, max_rows_per_species: int = 600) -> int:
    # Load all chains
    chains = [ _build_chain(p) for p in a3m_paths ]
    n_chain = len(chains)
    pad_lengths = [len(c.query_aligned) for c in chains]
    logging.getLogger("pair_msa").info(
        "species: loaded %d chains | lengths=%s", n_chain, pad_lengths
    )
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

    logging.getLogger("pair_msa").info(
        "species: candidate species=%d | min_chains=%d | skip if any chain count > %d",
        len(species_set), min_chains, max_rows_per_species,
    )
    # Prepare output
    lines: List[str] = []
    # First write concatenated query with per-chain aligned lengths for downstream tools
    q_header = ">query" + "".join(f"_len{len(c.query_aligned)}" for c in chains)
    lines.append(q_header + "\n")
    lines.append("".join(c.query_raw for c in chains) + "\n")

    total_rows = 0
    processed = 0
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
            sim_note = _format_sim_annotation(row_entries)
            lines.append(h + f"\tsp={sp}\trank={r+1}{sim_note}\n")
            lines.append(s + "\n")
            total_rows += 1
        processed += 1
        if processed % 200 == 0:
            logging.getLogger("pair_msa").info(
                "species: processed %d/%d | rows so far=%d", processed, len(species_set), total_rows
            )

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for ln in lines:
            f.write(ln)
    logging.getLogger("pair_msa").info("species: wrote %d rows -> %s", total_rows, out_path)
    return total_rows


def _pair_by_index(a3m_paths: Sequence[str], out_path: str) -> int:
    chains = [ _build_chain(p) for p in a3m_paths ]
    logging.getLogger("pair_msa").info(
        "index: loaded %d chains", len(chains)
    )
    pad_lengths = [len(c.query_aligned) for c in chains]
    # write concatenated query with per-chain aligned lengths
    lines: List[str] = []
    q_header = ">query" + "".join(f"_len{len(c.query_aligned)}" for c in chains)
    lines.append(q_header + "\n")
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
        sim_note = _format_sim_annotation(row_entries)
        lines.append(h + f"\tidx={r+1}{sim_note}\n")
        lines.append(s + "\n")
        total += 1
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for ln in lines:
            f.write(ln)
    logging.getLogger("pair_msa").info("index: wrote %d rows -> %s", total, out_path)
    return total


def pair_by_rank(
    a3m_paths: Sequence[str],
    out_path: str,
    min_chains: int = 2,
    pad_missing: bool = True,
    max_rows: Optional[int] = None,
) -> int:
    """Non-species pairing: per chain, sort all entries by similarity to query, then zip by rank.

    - If pad_missing is True, rows where some chains lack rank-r entries are padded with gaps,
      and kept only if present chains >= min_chains.
    - If pad_missing is False, only keep ranks r < min(len(chain_i)).
    - max_rows can be used to cap the total number of paired rows.
    """
    chains = [_build_chain(p) for p in a3m_paths]
    n_chain = len(chains)
    pad_lengths = [len(c.query_aligned) for c in chains]
    logging.getLogger("pair_msa").info(
        "rank: loaded %d chains | lengths=%s | pad_missing=%s | min_chains=%d | max_rows=%s",
        n_chain, pad_lengths, str(pad_missing), min_chains, str(max_rows)
    )

    # Prepare output with concatenated query and per-chain aligned lengths
    lines: List[str] = []
    q_header = ">query" + "".join(f"_len{len(c.query_aligned)}" for c in chains)
    lines.append(q_header + "\n")
    lines.append("".join(c.query_raw for c in chains) + "\n")

    # Build ranked non-query lists by similarity per chain
    ranked_per_chain: List[List[A3MEntry]] = []
    per_chain_counts = []
    for c in chains:
        lst = [e for e in c.entries if not e.header[1:].lower().startswith("query")]  # exclude query
        lst_sorted = sorted(lst, key=lambda e: e.similarity, reverse=True)
        ranked_per_chain.append(lst_sorted)
        per_chain_counts.append(len(lst_sorted))
    logging.getLogger("pair_msa").info(
        "rank: per-chain non-query counts=%s", per_chain_counts
    )

    if pad_missing:
        r_max = max((len(lst) for lst in ranked_per_chain), default=0)
    else:
        r_max = min((len(lst) for lst in ranked_per_chain), default=0)

    total_rows = 0
    for r in range(r_max):
        row_entries: List[Optional[A3MEntry]] = []
        present = 0
        for i in range(n_chain):
            lst = ranked_per_chain[i]
            if r < len(lst):
                row_entries.append(lst[r])
                present += 1
            else:
                row_entries.append(None)
        if present >= min_chains:
            h, s = _concat_rows(row_entries, pad_lengths)
            sim_note = _format_sim_annotation(row_entries)
            lines.append(h + f"\trank={r+1}\tpresent={present}/{n_chain}{sim_note}\n")
            lines.append(s + "\n")
            total_rows += 1
            if max_rows is not None and total_rows >= max_rows:
                break

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.writelines(lines)
    logging.getLogger("pair_msa").info("rank: wrote %d rows -> %s", total_rows, out_path)
    return total_rows


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Pair multimer MSAs by species and output a concatenated paired A3M.\n"
            "Input: per-chain pairing A3M files in chain order (e.g., msa/0/pairing.a3m msa/1/pairing.a3m)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("a3m", nargs='*', help="Per-chain pairing A3M files in chain order (ignored when --auto_dir is set)")
    p.add_argument("--out", type=str, default="paired_concat.a3m", help="Output paired A3M path")
    p.add_argument("--strategy", choices=["species","index","rank","stitch"], default="species",
                   help=(
                       "Pairing strategy: \n"
                       "  - species: group by taxid (UniRef100_*_<taxid>/).\n"
                       "  - index: zip rows by position across chains (as-is order).\n"
                       "  - rank: non-species; sort each chain by similarity-to-query and zip by rank, with optional padding.\n"
                       "  - stitch: ColabFold pair.a3m semantics with optional filters.") )
    p.add_argument("--min_chains", type=int, default=2, help="Minimum chains a species must appear in to be paired")
    p.add_argument("--pad_missing", action="store_true", help="In rank mode, allow rows where some chains are missing (pad with gaps)")
    p.add_argument("--max_rows", type=int, default=None, help="Optional cap on number of paired rows (rank mode)")
    p.add_argument("--max_rows_per_species", type=int, default=600, help="Skip species if any chain exceeds this count")
    # Auto-detect inputs
    p.add_argument("--auto_dir", type=str, default=None,
                   help=("Auto-detect inputs from a job directory (output of query_msa/cli fetch).\n"
                         "- For stitch: prefer <auto_dir>/raw/pair.a3m; fallback to per-chain pairing.a3m in <auto_dir>/msa/<i>/.\n"
                         "- For species/index/rank: use per-chain pairing.a3m under <auto_dir>/msa/<i>/ (i=0,1,...)"))
    # ColabFold stitched advanced options
    p.add_argument("--pair_a3m", type=str, default=None, help="Path to ColabFold raw pair.a3m (overrides per-chain inputs in stitch mode)")
    p.add_argument("--min_coverage", type=float, default=0.75)
    p.add_argument("--min_identity", type=float, default=0.15)
    p.add_argument("--max_evalue", type=float, default=None)
    p.add_argument("--min_alnscore", type=float, default=None)
    p.add_argument("--genomic_distance", type=int, default=1, help="Optional Δgene for local filtering (stitch mode)")
    p.add_argument("--use_aligned", action="store_true", help="Strip lowercase insertions and right-pad per-chain blocks to query aligned length when stitching")
    # Logging
    p.add_argument("--verbose", action="store_true", help="Print progress and summary stats during pairing")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_argparser().parse_args(argv)
    # Setup logging
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="[%(levelname)s] %(message)s",
    )
    log = logging.getLogger("pair_msa")
    # Auto-detect inputs from a job directory if provided
    if args.auto_dir is not None:
        root = Path(args.auto_dir)
        raw_dir = root / "raw"
        msa_dir = root / "msa"
        # Collect per-chain pairing files
        chain_files = []
        if msa_dir.exists():
            # find numeric subdirs and pairing.a3m inside
            subdirs = sorted([p for p in msa_dir.iterdir() if p.is_dir() and p.name.isdigit()], key=lambda p: int(p.name))
            for sd in subdirs:
                f = sd / "pairing.a3m"
                if f.exists():
                    chain_files.append(str(f))
        # Stitch preference: raw/pair.a3m
        if args.strategy == "stitch":
            if args.pair_a3m is None and raw_dir.exists():
                p = raw_dir / "pair.a3m"
                if p.exists():
                    args.pair_a3m = str(p)
                    log.info("[auto] using ColabFold container: %s", args.pair_a3m)
            # If still no pair.a3m, fall back to per-chain
            if not chain_files and args.pair_a3m is None:
                raise ValueError(f"--auto_dir given but neither {raw_dir/'pair.a3m'} nor any <msa/i/pairing.a3m> found")
            args.a3m = chain_files
        else:
            if not chain_files:
                raise ValueError(f"--auto_dir given but no per-chain pairing.a3m found under {msa_dir}")
            args.a3m = chain_files
        if chain_files:
            log.info("[auto] detected %d chain pairing files under %s", len(chain_files), msa_dir)

    if args.strategy == "species":
        log.info("strategy=species | inputs=%d chains | min_chains=%d | max_rows_per_species=%s", len(args.a3m), args.min_chains, str(args.max_rows_per_species))
        n = pair_by_species(args.a3m, args.out, args.min_chains, args.max_rows_per_species)
        log.info("done | wrote %d rows -> %s", n, args.out)
    elif args.strategy == "index":
        log.info("strategy=index | inputs=%d chains", len(args.a3m))
        n = _pair_by_index(args.a3m, args.out)
        log.info("done | wrote %d rows -> %s", n, args.out)
    elif args.strategy == "rank":
        if not args.a3m:
            raise ValueError("rank mode requires per-chain pairing A3Ms as positional args (one per chain in order)")
        log.info("strategy=rank | inputs=%d chains | pad_missing=%s | min_chains=%d | max_rows=%s", len(args.a3m), args.pad_missing, args.min_chains, str(args.max_rows))
        n = pair_by_rank(args.a3m, args.out, min_chains=args.min_chains, pad_missing=args.pad_missing, max_rows=args.max_rows)
        log.info("done | wrote %d rows -> %s", n, args.out)
    else:
        # stitch: prefer raw pair.a3m; fallback to per-chain pairing.a3m zipped by index
        entries: List[dict] = []
        # Robustness: if user mistakenly passed a pair.a3m as positional arg, auto-handle it
        if (args.pair_a3m is None and len(args.a3m) == 1 and args.a3m[0].endswith("pair.a3m")
                and Path(args.a3m[0]).exists()):
            args.pair_a3m = args.a3m[0]

        if args.pair_a3m is not None and Path(args.pair_a3m).exists():
            log.info("strategy=stitch | source=pair.a3m | file=%s", args.pair_a3m)
            entries = stitch_utils.parse_paired_a3m(args.pair_a3m)
        elif args.a3m:
            log.info("strategy=stitch | source=per-chain pairing.a3m | chains=%d", len(args.a3m))
            per_chain_lines: List[List[str]] = []
            for pth in args.a3m:
                with open(pth, 'r') as f:
                    per_chain_lines.append(f.read().splitlines())
            parsed = [stitch_utils.parse_msa_lines(lines) for lines in per_chain_lines]
            if not parsed:
                raise ValueError("No inputs provided for stitch mode.")
            k = min(len(p) for p in parsed)
            if len(parsed) == 1:
                raise ValueError(
                    "stitch mode fallback expects per-chain pairing.a3m files; got a single file. "
                    "If you intended to use ColabFold pair container, pass it via --pair_a3m"
                )
            for r in range(k):
                headers: List[str] = []
                sequences: List[str] = []
                covs: List[float] = []
                ids: List[float] = []
                evs: List[float] = []
                als: List[float] = []
                uids: List[str] = []
                upnums: List[int] = []
                has_uniref = True
                for c in range(len(parsed)):
                    e = parsed[c][r]
                    headers.append(e['header'])
                    sequences.append(e['sequence'])
                    covs.append(e['coverage'])
                    ids.append(e['identity'])
                    evs.append(e['evalue'])
                    als.append(e['alnscore'])
                    uids.append(e['uid'])
                    upnums.append(e['uniprot_num'])
                    has_uniref = has_uniref and e['has_uniref']
                entries.append({
                    'headers': headers,
                    'sequences': sequences,
                    'coverages': covs,
                    'identities': ids,
                    'evalues': evs,
                    'alnscores': als,
                    'uids': uids,
                    'uniprot_nums': upnums,
                    'has_uniref': has_uniref,
                    'is_query': (r == 0),
                })
        else:
            raise ValueError("stitch mode requires --pair_a3m or per-chain A3M inputs")

        before = len(entries)
        log.info("stitch: entries(before)=%d | filters: cov=%s id=%s e=%s aln=%s Δgene=%s use_aligned=%s",
                 before,
                 str(args.min_coverage), str(args.min_identity), str(args.max_evalue), str(args.min_alnscore),
                 str(args.genomic_distance), str(args.use_aligned))
        n, filtered = stitch_utils.save_stitched(
            entries,
            args.out,
            min_coverage=args.min_coverage,
            min_identity=args.min_identity,
            max_evalue=args.max_evalue,
            min_alnscore=args.min_alnscore,
            max_genomic_distance=args.genomic_distance,
            use_aligned=args.use_aligned,
        )
        log.info("stitch: entries(after)=%d | wrote -> %s", n, args.out)
        if args.verbose:
            try:
                stats_before = stitch_utils.get_stats(entries)
                stats_after = stitch_utils.get_stats(filtered)
                log.info("stats-before: %s", stats_before)
                log.info("stats-after : %s", stats_after)
            except Exception:
                pass
        print(f"Wrote {n} paired rows to {args.out} using {args.strategy}")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
