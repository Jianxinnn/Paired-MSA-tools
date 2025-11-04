from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import os
from pathlib import Path


@dataclass
class _UniProtEncoder:
    pa: Dict[str, int]
    ma: List[List[Dict[str, int]]]
    upi_encoding: Dict[str, int]


def _init_uniprot_encoder() -> _UniProtEncoder:
    from string import ascii_uppercase
    pa = {a: 0 for a in ascii_uppercase}
    for a in ["O", "P", "Q"]:
        pa[a] = 1
    ma: List[List[Dict[str, int]]] = [[{} for _ in range(6)], [{} for _ in range(6)]]
    for n, t in enumerate(range(10)):
        for i in [0, 1]:
            for j in [0, 4]:
                ma[i][j][str(t)] = n
    for n, t in enumerate(list(ascii_uppercase) + list(range(10))):
        for i in [0, 1]:
            for j in [1, 2]:
                ma[i][j][str(t)] = n
        ma[1][3][str(t)] = n
    for n, t in enumerate(ascii_uppercase):
        ma[0][3][str(t)] = n
        for i in [0, 1]:
            ma[i][5][str(t)] = n
    upi_encoding: Dict[str, int] = {}
    hex_chars = list(range(10)) + ["A", "B", "C", "D", "E", "F"]
    for n, char in enumerate(hex_chars):
        upi_encoding[str(char)] = n
    return _UniProtEncoder(pa=pa, ma=ma, upi_encoding=upi_encoding)


_ENC = _init_uniprot_encoder()


def _extract_uniprot_id(header: str) -> str:
    """
    Extract UniProt ID that follows the "UniRef" token in the header.

    ColabFold / UniRef headers commonly look like:
      ">UniRef100_Q8ZIN0_511145\t..." or ">UniRef100_UPI002236BB6D\t..."

    We should capture only the core UniProt accession (e.g., "Q8ZIN0"
    or the full UPI "UPI002236BB6D"), stopping at '_' (taxid) or other
    separators. The previous implementation did not stop at '_' which
    could yield "Q8ZIN0_511145" and break numeric encoding downstream,
    disabling genomic-distance filtering.
    """
    pos = header.find("UniRef")
    if pos == -1:
        return ""
    start = header.find('_', pos)
    if start == -1:
        return ""
    start += 1
    end = start
    # Stop on whitespace, second underscore (taxid) or slash
    while end < len(header) and header[end] not in ' _\t/':
        end += 1
    uid = header[start:end]
    # Allow UPI IDs
    if len(uid) >= 3 and uid[:3] == "UPI":
        return uid
    # Regular UniProt accessions are 6 or 10 chars, starting with a letter
    if len(uid) not in [6, 10]:
        return ""
    if not uid[0].isalpha():
        return ""
    return uid


def _uniprot_to_number(uniprot_ids: List[str]) -> List[int]:
    numbers: List[int] = []
    for uni in uniprot_ids:
        if not uni or not uni[0].isalpha():
            numbers.append(0)
            continue
        if uni.startswith("UPI") and len(uni) == 13:
            hex_part = uni[3:]
            num = 0
            tot = 1
            for u in reversed(hex_part):
                if str(u) in _ENC.upi_encoding:
                    num += _ENC.upi_encoding[str(u)] * tot
                    tot *= 16
                else:
                    num = 0
                    break
            numbers.append(num + 10**15)
            continue
        p = _ENC.pa.get(uni[0], 0)
        tot, num = 1, 0
        if len(uni) == 10:
            for n, u in enumerate(reversed(uni[-4:])):
                if str(u) in _ENC.ma[p][n]:
                    num += _ENC.ma[p][n][str(u)] * tot
                    tot *= len(_ENC.ma[p][n].keys())
        for n, u in enumerate(reversed(uni[:6])):
            if n < len(_ENC.ma[p]) and str(u) in _ENC.ma[p][n]:
                num += _ENC.ma[p][n][str(u)] * tot
                tot *= len(_ENC.ma[p][n].keys())
        numbers.append(num)
    return numbers


def parse_msa_lines(lines: List[str]) -> List[Dict]:
    entries: List[Dict] = []
    i = 0
    is_first = True
    while i < len(lines):
        line = lines[i].rstrip()
        if line.startswith('>'):
            header = line
            seq_lines: List[str] = []
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                if lines[i].strip():
                    seq_lines.append(lines[i].rstrip())
                i += 1
            sequence = ''.join(seq_lines)
            header_parts = header.split('\t')
            header_clean = header_parts[0].lstrip('>').replace('UniRef100_', '')
            uid = _extract_uniprot_id(header)
            has_uniref = "UniRef" in header
            uniprot_num = 0
            if uid:
                uniprot_nums = _uniprot_to_number([uid])
                uniprot_num = uniprot_nums[0] if uniprot_nums else 0
            if is_first:
                coverage = 1.0
                identity = 1.0
                evalue = 0.0
                alnscore = float('inf')
                is_first = False
            else:
                coverage = None
                identity = None
                evalue = None
                alnscore = None
                if len(header_parts) >= 10:
                    try:
                        alnscore = float(header_parts[1])
                        identity = float(header_parts[2])
                        evalue = float(header_parts[3])
                        q_start = int(header_parts[4])
                        q_end = int(header_parts[5])
                        q_len = int(header_parts[6])
                        coverage = (q_end - q_start + 1) / q_len
                    except Exception:
                        pass
                if coverage is None:
                    coverage = 0.0
                if identity is None:
                    identity = 0.0
                if evalue is None:
                    evalue = float('inf')
                if alnscore is None:
                    alnscore = 0.0
            entries.append({
                'header': header_clean,
                'sequence': sequence,
                'coverage': coverage,
                'identity': identity,
                'evalue': evalue,
                'alnscore': alnscore,
                'uid': uid,
                'uniprot_num': uniprot_num,
                'has_uniref': has_uniref,
            })
        else:
            i += 1
    return entries


def parse_paired_a3m(path: str) -> List[Dict]:
    # Split by chain blocks (separated by NUL in some outputs)
    raw_msas: Dict[int, List[str]] = {}
    update_M = True
    M: Optional[int] = None
    with open(path, 'r') as f:
        for line in f:
            if "\x00" in line:
                line = line.replace("\x00", "")
                update_M = True
            if line.startswith(">") and update_M:
                M = int(line[1:].rstrip().split('_')[-1])
                update_M = False
                if M not in raw_msas:
                    raw_msas[M] = []
            if M is not None:
                raw_msas[M].append(line.rstrip())
    parsed_msas: Dict[int, List[Dict]] = {}
    for seq_id, lines in raw_msas.items():
        parsed_msas[seq_id] = parse_msa_lines(lines)
    seq_ids = sorted(parsed_msas.keys())
    num_entries_per_chain = [len(parsed_msas[sid]) for sid in seq_ids]
    min_entries = min(num_entries_per_chain) if num_entries_per_chain else 0
    stitched: List[Dict] = []
    for i in range(min_entries):
        headers: List[str] = []
        sequences: List[str] = []
        coverages: List[float] = []
        identities: List[float] = []
        evalues: List[float] = []
        alnscores: List[float] = []
        uids: List[str] = []
        uniprot_nums: List[int] = []
        has_uniref = True
        for sid in seq_ids:
            entry = parsed_msas[sid][i]
            headers.append(entry['header'])
            sequences.append(entry['sequence'])
            coverages.append(entry['coverage'])
            identities.append(entry['identity'])
            evalues.append(entry['evalue'])
            alnscores.append(entry['alnscore'])
            uids.append(entry['uid'])
            uniprot_nums.append(entry['uniprot_num'])
            has_uniref = has_uniref and entry['has_uniref']
        stitched.append({
            'headers': headers,
            'sequences': sequences,
            'coverages': coverages,
            'identities': identities,
            'evalues': evalues,
            'alnscores': alnscores,
            'uids': uids,
            'uniprot_nums': uniprot_nums,
            'has_uniref': has_uniref,
            'is_query': (i == 0),
        })
    return stitched


def parse_single_msas(prefix: str) -> List[Dict]:
    msa_files = ['uniref.a3m', 'bfd.mgnify30.metaeuk30.smag30.a3m']
    all_entries: List[Dict] = []
    seen_sequences = set()
    for msa_file in msa_files:
        msa_path = os.path.join(prefix, msa_file)
        if os.path.exists(msa_path):
            with open(msa_path, 'r') as f:
                lines = f.readlines()
            entries = parse_msa_lines(lines)
            for entry in entries:
                if entry['sequence'] not in seen_sequences:
                    seen_sequences.add(entry['sequence'])
                    all_entries.append({
                        'headers': [entry['header']],
                        'sequences': [entry['sequence']],
                        'coverages': [entry['coverage']],
                        'identities': [entry['identity']],
                        'evalues': [entry['evalue']],
                        'alnscores': [entry['alnscore']],
                        'uids': [entry['uid']],
                        'uniprot_nums': [entry['uniprot_num']],
                        'has_uniref': entry['has_uniref'],
                        'is_query': len(all_entries) == 0,
                    })
    return all_entries


def _calc_genomic_distances(entry: Dict) -> List[int]:
    distances: List[int] = []
    nums = entry['uniprot_nums']
    for i in range(1, len(nums)):
        if nums[i - 1] and nums[i]:
            dist = abs(nums[i] - nums[i - 1])
            distances.append(dist)
        else:
            distances.append(-1)
    return distances


def save_stitched(entries: List[Dict], output_file: str,
                  min_coverage: Optional[float] = None,
                  min_identity: Optional[float] = None,
                  max_evalue: Optional[float] = None,
                  min_alnscore: Optional[float] = None,
                  max_genomic_distance: Optional[int] = None,
                  use_aligned: bool = False) -> Tuple[int, List[Dict]]:
    if not entries:
        raise ValueError("No entries to save.")
    # Normalize percentages
    if min_coverage is not None and min_coverage > 1:
        min_coverage = min_coverage / 100
    if min_identity is not None and min_identity > 1:
        min_identity = min_identity / 100
    num_chains = len(entries[0]['sequences']) if entries else 0
    sequences_written = 0
    sequences_filtered = 0
    filtered_entries: List[Dict] = []
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    # Precompute per-chain aligned query lengths to stabilize block widths when use_aligned=True
    q_aln_len: List[int] = []
    if entries:
        for seq in entries[0]['sequences']:
            L = len(''.join(ch for ch in seq if not ch.islower())) if use_aligned else len(seq)
            q_aln_len.append(L)
    with open(output_file, 'w') as f:
        for entry in entries:
            if not entry.get('is_query', False):
                reasons: List[str] = []
                if min_coverage and any(c < min_coverage for c in entry['coverages']):
                    reasons.append('cov')
                if min_identity and any(i < min_identity for i in entry['identities']):
                    reasons.append('id')
                if max_evalue is not None and any(e is not None and e > max_evalue for e in entry['evalues']):
                    reasons.append('evalue')
                if min_alnscore is not None and any(a is not None and a < min_alnscore for a in entry['alnscores']):
                    reasons.append('aln')
                if max_genomic_distance is not None and entry['has_uniref']:
                    distances = _calc_genomic_distances(entry)
                    if num_chains == 2:
                        if distances and distances[0] != -1 and distances[0] > max_genomic_distance:
                            reasons.append('gdist')
                    else:
                        valid = [d for d in distances if d != -1]
                        if valid and all(d > max_genomic_distance for d in valid):
                            reasons.append('gdist')
                if reasons:
                    sequences_filtered += 1
                    continue
            # Build header and sequence
            if entry.get('is_query', False):
                # Use aligned lengths for header so downstream visualization can
                # accurately split concatenated positions into chains after
                # removing lowercase insertions.
                header = "query"
                if q_aln_len:
                    for L in q_aln_len:
                        header += f"_len{L}"
                else:
                    for seq in entry['sequences']:
                        L = len(''.join(ch for ch in seq if not ch.islower()))
                        header += f"_len{L}"
            else:
                parts: List[str] = []
                for i, uid in enumerate(entry['uids']):
                    parts.append(uid if uid else entry['headers'][i])
                header = '_'.join(parts)
                if entry['has_uniref'] and all(entry['uids']):
                    for dist in _calc_genomic_distances(entry):
                        if dist != -1:
                            header += f"_{dist}"
            if use_aligned:
                parts: List[str] = []
                for i, seq in enumerate(entry['sequences']):
                    s = ''.join(ch for ch in seq if not ch.islower())
                    # If a row is shorter than query (shouldn't happen often), right-pad with gaps
                    if i < len(q_aln_len) and len(s) < q_aln_len[i]:
                        s = s + ('-' * (q_aln_len[i] - len(s)))
                    parts.append(s)
                sequence = ''.join(parts)
            else:
                sequence = ''.join(entry['sequences'])
            sequence = sequence.replace('\x00', '')
            f.write(f">{header}\n{sequence}\n")
            sequences_written += 1
            filtered_entries.append(entry)
    return sequences_written, filtered_entries


def get_stats(entries: Optional[List[Dict]] = None) -> Dict:
    if not entries:
        return {}
    num_chains = len(entries[0]['sequences']) if entries else 0
    stats: Dict[str, object] = {
        'num_chains': num_chains,
        'num_entries': len(entries),
    }
    for i in range(num_chains):
        coverages = [e['coverages'][i] for e in entries[1:]]
        identities = [e['identities'][i] for e in entries[1:]]
        evalues = [e['evalues'][i] for e in entries[1:] if e['evalues'][i] is not None]
        alnscores = [e['alnscores'][i] for e in entries[1:] if e['alnscores'][i] is not None]
        chain_id = i + 101
        stats[f'chain_{chain_id}'] = {
            'query_length': len(entries[0]['sequences'][i]) if entries else 0,
            'avg_coverage': sum(coverages) / len(coverages) if coverages else 0,
            'avg_identity': sum(identities) / len(identities) if identities else 0,
            'avg_evalue': sum(evalues) / len(evalues) if evalues else 0,
            'avg_alnscore': sum(alnscores) / len(alnscores) if alnscores else 0,
            'min_evalue': min(evalues) if evalues else None,
            'max_alnscore': max(alnscores) if alnscores else None,
        }
    return stats
