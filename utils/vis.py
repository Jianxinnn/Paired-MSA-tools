from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


_LEN_RE = re.compile(r"_len(\d+)")


def _remove_inserts(seq: str) -> str:
    """Remove lowercase insertions from an A3M row and return aligned sequence.

    Keeps only non-lowercase characters; typically 'A-Z' and '-' remain.
    """
    if not seq:
        return seq
    return "".join(ch for ch in seq if not ch.islower())


def _first_entry(lines: Sequence[str]) -> Tuple[str, str]:
    """Return (header, sequence) for the first FASTA entry in lines."""
    header = None
    seq_parts: List[str] = []
    for ln in lines:
        if ln.startswith(">"):
            if header is None:
                header = ln.strip()
            elif header is not None and seq_parts:
                break
        else:
            if header is not None:
                seq_parts.append(ln.strip())
    if header is None:
        raise ValueError("No FASTA header found in A3M file")
    return header, "".join(seq_parts)


def _parse_chain_lengths_from_query_header(header: str) -> Optional[List[int]]:
    """Extract chain lengths from query header like ">query_len569_len471".

    Returns list of lengths or None if not present.
    """
    lens = [int(m) for m in _LEN_RE.findall(header)]
    return lens if lens else None


def compute_msa_depth(
    a3m_path: str,
    include_query: bool = False,
    skip_misaligned: bool = True,
) -> Dict[str, object]:
    """Compute per-position MSA depth (non-gap counts) from an A3M.

    - Removes lowercase insertions in every row before counting.
    - Counts positions with characters other than '-' (and '.') as present.
    - If `include_query` is False, the first entry is excluded from counts.

    Returns a dict with keys:
      depth: np.ndarray shape (L,)
      L: int total aligned length
      n_rows_total: total number of entries in file
      n_rows_used: number of rows included in counts
      n_rows_skipped: rows skipped due to misalignment length (when enabled)
      chain_lengths: Optional[List[int]] from query header
    """
    text = Path(a3m_path).read_text().replace("\x00", "")
    lines = [ln for ln in text.splitlines() if ln and not ln.lstrip().startswith('#')]
    if not lines:
        raise ValueError(f"Empty A3M: {a3m_path}")

    # First pass: query header+sequence and aligned length
    q_header, q_seq_raw = _first_entry(lines)
    q_seq_aln = _remove_inserts(q_seq_raw)
    L = len(q_seq_aln)
    if L == 0:
        raise ValueError("Query aligned length is zero after removing inserts")
    chain_lengths = _parse_chain_lengths_from_query_header(q_header)

    # Second pass: iterate entries and accumulate counts
    depth = np.zeros(L, dtype=np.int32)
    n_total = 0
    n_used = 0
    n_skipped = 0

    cur_header: Optional[str] = None
    cur_seq_parts: List[str] = []

    def _flush_row(header: Optional[str], seq_parts: List[str]):
        nonlocal n_total, n_used, n_skipped
        if header is None:
            return
        n_total += 1
        seq_raw = "".join(seq_parts)
        seq_aln = _remove_inserts(seq_raw)
        if len(seq_aln) != L:
            if skip_misaligned:
                n_skipped += 1
                return
            # pad/truncate to L as a fallback
            if len(seq_aln) < L:
                seq_aln = seq_aln + ("-" * (L - len(seq_aln)))
            else:
                seq_aln = seq_aln[:L]
        # Optionally exclude the first (query) row
        if not include_query and n_total == 1:
            return
        for i, ch in enumerate(seq_aln):
            if ch != '-' and ch != '.':
                depth[i] += 1
        n_used += 1

    for ln in lines:
        if ln.startswith(">"):
            # new entry
            _flush_row(cur_header, cur_seq_parts)
            cur_header = ln.strip()
            cur_seq_parts = []
        else:
            cur_seq_parts.append(ln.strip())
    _flush_row(cur_header, cur_seq_parts)

    return {
        "depth": depth,
        "L": int(L),
        "n_rows_total": int(n_total),
        "n_rows_used": int(n_used),
        "n_rows_skipped": int(n_skipped),
        "chain_lengths": chain_lengths,
    }


def _moving_average(x: np.ndarray, w: int) -> np.ndarray:
    if w is None or w <= 1:
        return x
    w = int(w)
    if w > len(x):
        return x
    c = np.convolve(x, np.ones(w, dtype=float) / float(w), mode="same")
    return c


def plot_msa_depth(
    depth: np.ndarray,
    chain_lengths: Optional[Sequence[int]] = None,
    *,
    title: Optional[str] = None,
    out_path: Optional[str] = None,
    dpi: int = 150,
    normalize: bool = False,
    smooth: int = 1,
    include_query: bool = False,
):
    """Plot per-position MSA depth; split into subplots per chain if lengths provided.

    - depth: integer counts per aligned position
    - chain_lengths: e.g., [569, 471] for two chains; if None, plot as single chain
    - normalize: if True, plot fraction by dividing by max(depth) (or N rows if known)
    - smooth: moving-average window (in positions)
    - include_query: label hint only (legend text)
    """
    import matplotlib.pyplot as plt  # local import to keep base import light

    y = depth.astype(float)
    if normalize:
        denom = float(y.max() if y.size else 1.0)
        if denom > 0:
            y = y / denom
    if smooth and smooth > 1:
        y = _moving_average(y, int(smooth))

    if not chain_lengths:
        chain_lengths = [len(y)]
    if sum(chain_lengths) != len(y):
        chain_lengths = [len(y)]

    # Figure size heuristic
    total_L = len(y)
    n_chain = len(chain_lengths)
    width = min(16.0, 6.0 + total_L / 250.0)
    height = min(2.0 * n_chain, 0.9 + 1.2 * n_chain)
    fig, axes = plt.subplots(n_chain, 1, figsize=(width, height), sharey=True)
    if n_chain == 1:
        axes = [axes]

    start = 0
    y_max = float(np.max(y) if y.size else 1.0)
    for idx, (ax, Lc) in enumerate(zip(axes, chain_lengths), start=1):
        end = start + Lc
        xs = np.arange(1, Lc + 1)
        ax.plot(xs, y[start:end], lw=1.0, color="#1f77b4")
        ax.fill_between(xs, 0, y[start:end], color="#1f77b4", alpha=0.15, linewidth=0)
        ax.set_xlim(1, Lc)
        ax.set_ylabel("depth" + (" (frac)" if normalize else ""))
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
        ax.set_title(f"Chain {idx}/{n_chain} (L={Lc})", fontsize=10)
        start = end
    axes[-1].set_xlabel("position (aligned, inserts removed)")
    if title:
        fig.suptitle(title + ("  [incl query]" if include_query else ""), fontsize=11)
    fig.tight_layout(rect=[0, 0.0, 1, 0.98])
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=dpi)
    return fig, axes
