from __future__ import annotations

import io
import os
import tarfile
from pathlib import Path
from typing import Dict, List, Mapping, Tuple

from .seqio import parse_fasta_string, ensure_dir, write_lines


def read_m8(m8_file: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open(m8_file, "r") as infile:
        for line in infile:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            hit_name, ncbi_taxid = cols[1], cols[2]
            mapping[hit_name] = ncbi_taxid
    return mapping


def read_a3m(a3m_file: str) -> Tuple[List[str], List[str], int]:
    heads: List[str] = []
    seqs: List[str] = []
    uniref_index = 0
    with open(a3m_file, "r") as infile:
        query_name = None
        for idx, line in enumerate(infile):
            if line.startswith(">"):
                heads.append(line)
                if idx == 0:
                    query_name = line
                elif idx > 0 and line == query_name:
                    uniref_index = idx
            else:
                seqs.append(line)
    return heads, seqs, uniref_index


def protenix_split_to_pairing_and_nonpairing(
    query_seq: str, seq_dir: str, raw_a3m_path: str, uniref_to_ncbi_taxid: Mapping[str, str]
) -> None:
    heads, msa_seqs, uniref_index = read_a3m(raw_a3m_path)
    uniref100_lines = [">query\n", f"{query_seq}\n"]
    other_lines = [">query\n", f"{query_seq}\n"]

    for idx, (head, msa_seq) in enumerate(zip(heads, msa_seqs)):
        if msa_seq.rstrip("\n") == query_seq:
            continue
        uniref_id = head.split("\t")[0][1:]
        ncbi_taxid = uniref_to_ncbi_taxid.get(uniref_id, None)
        if (ncbi_taxid is not None) and (idx < (uniref_index // 2)):
            if not uniref_id.startswith("UniRef100_"):
                head = head.replace(uniref_id, f"UniRef100_{uniref_id}_{ncbi_taxid}/")
            else:
                head = head.replace(uniref_id, f"{uniref_id}_{ncbi_taxid}/")
            uniref100_lines.extend([head, msa_seq])
        else:
            other_lines.extend([head, msa_seq])

    write_lines(os.path.join(seq_dir, "pairing.a3m"), uniref100_lines)
    write_lines(os.path.join(seq_dir, "non_pairing.a3m"), other_lines)


def protenix_nonpairing_only(query_seq: str, seq_dir: str, raw_a3m_path: str) -> None:
    heads, msa_seqs, _ = read_a3m(raw_a3m_path)
    lines = [">query\n", f"{query_seq}\n"]
    for head, msa_seq in zip(heads, msa_seqs):
        if msa_seq.rstrip("\n") == query_seq:
            continue
        lines.extend([head, msa_seq])
    write_lines(os.path.join(seq_dir, "non_pairing.a3m"), lines)
    write_lines(os.path.join(seq_dir, "pairing.a3m"), [">query\n", f"{query_seq}\n"])


def colabfold_post_nonpair(prefix: str, query_header: str, out_dir: str) -> str:
    env_path = os.path.join(prefix, "bfd.mgnify30.metaeuk30.smag30.a3m")
    uniref_path = os.path.join(prefix, "uniref.a3m")
    if not (os.path.exists(env_path) and os.path.exists(uniref_path)):
        raise FileNotFoundError("ColabFold result missing env/uniref a3m files")
    env = parse_fasta_string(Path(env_path).read_text().replace("\x00", ""))
    uniref = parse_fasta_string(Path(uniref_path).read_text().replace("\x00", ""))
    query_id = str(int(query_header.split("_")[-1])) if "_" in query_header else "0"
    seq_dir = ensure_dir(os.path.join(out_dir, query_id))
    non_pair_f = os.path.join(seq_dir, "non_pairing.a3m")
    query_seq = ""
    for d in (env, uniref):
        for k, v in d.items():
            if k.startswith("query_"):
                query_seq = v
                break
        if query_seq:
            break
    with open(non_pair_f, "w") as f:
        f.write(">query\n" + query_seq + "\n")
        for k, v in env.items():
            if k.startswith("query_"):
                continue
            f.write(f">{k}\n{v}\n")
        for k, v in uniref.items():
            if k.startswith("query_"):
                continue
            f.write(f">{k}\n{v}\n")
    write_lines(os.path.join(seq_dir, "pairing.a3m"), [">query\n", f"{query_seq}\n"])
    return seq_dir


def colabfold_post_pair(prefix: str, out_dir: str) -> List[str]:
    pair_a3m = os.path.join(prefix, "pair.a3m")
    if not os.path.exists(pair_a3m):
        raise FileNotFoundError("ColabFold pairing result missing pair.a3m")
    chunks = Path(pair_a3m).read_text().split("\x00")
    out_subdirs: List[str] = []
    for chunk in chunks[:-1]:
        fasta = parse_fasta_string(chunk)
        query_header = next(iter(fasta.keys()))
        query_id = str(int(query_header.split("_")[-1])) if "_" in query_header else "0"
        seq_dir = ensure_dir(os.path.join(out_dir, query_id))
        with open(os.path.join(seq_dir, "pairing.a3m"), "w") as f:
            for i, (k, v) in enumerate(fasta.items()):
                if k.startswith("query_"):
                    f.write(">query\n" + v + "\n")
                else:
                    f.write(f">{k}\n{v}\n")
        out_subdirs.append(seq_dir)
    return out_subdirs


def extract_tarball_to_dir(blob: bytes, out_root: str) -> str:
    extract_root = ensure_dir(os.path.join(out_root, "raw"))
    with tarfile.open(fileobj=io.BytesIO(blob)) as tar_gz:
        tar_gz.extractall(extract_root)
    return extract_root
