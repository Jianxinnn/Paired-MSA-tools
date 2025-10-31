from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional


def parse_fasta_string(fasta_string: str) -> Dict[str, str]:
    fasta_dict: Dict[str, str] = {}
    header: Optional[str] = None
    for line in fasta_string.strip().split("\n"):
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            fasta_dict[header] = ""
        else:
            if header is None:
                raise ValueError("FASTA content missing header before sequence lines")
            fasta_dict[header] += line.strip()
    return fasta_dict


def normalize_input_sequences(seq: Optional[str], fasta: Optional[str]) -> List[str]:
    out: List[str] = []
    if seq:
        parts = [p.strip() for p in seq.replace("\n", "").split(":")]
        out = [p for p in parts if p]
    elif fasta:
        d = parse_fasta_string(Path(fasta).read_text())
        out = [v.strip() for v in d.values() if v.strip()]
    else:
        raise ValueError("Either --seq or --fasta must be provided.")
    if not out:
        raise ValueError("No valid sequences found in input.")
    return [s.replace(" ", "").upper() for s in out]


def ensure_dir(path: str) -> str:
    Path(path).mkdir(parents=True, exist_ok=True)
    return os.path.abspath(path)


def write_lines(path: str, lines: List[str]) -> None:
    with open(path, "w") as f:
        for line in lines:
            f.write(line)

