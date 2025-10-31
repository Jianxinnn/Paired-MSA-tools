Paired MSA Tools (standalone)

Proposed GitHub name: paired-msa-tools

Overview
- Standalone Python CLI to fetch Multiple Sequence Alignments (MSA) for monomers and multimers from remote services and to inspect/compose paired multimer alignments.
- No dependency on the Protenix repository. Prefer the single entry `query_msa.py` in this folder. Advanced users can use the CLIs below.
  - `cli/fetch_msa.py` — remote retrieval from Protenix server or ColabFold server; outputs per-chain `non_pairing.a3m` and `pairing.a3m` ready for structure predictors.
  - `pair_msa.py` — optional inspector to export a single “concatenated paired multimer A3M” for visualization/debugging.

Highlights
- Supports single-chain and multi-chain input (colon-separated `A:B[:C...]` or multi-entry FASTA).
- Protenix backend (default host: https://protenix-server.com/api/msa): returns `0.a3m` plus `uniref_tax.m8`; we split into `pairing.a3m` and `non_pairing.a3m` using UniRef100+taxid mapping.
- ColabFold backend (default host: https://api.colabfold.com): returns `env`+`uniref` A3Ms for non-pairing and a `pair.a3m` container for pairing; we split these into per-chain files.

Install
```
pip install -r requirements.txt
```

Quick Start
0) One-file entry (recommended):
```
# Single job (mono or multimer)
python query_msa.py --seq AAAA:BBBB --server protenix --out_dir ./msa_job --auth_user <user> --auth_pass <pass>

# Batch (CSV/TSV with name,seq or name,fasta)
python query_msa.py --input jobs.csv --out_dir ./msa_batch --server colabfold --pairing greedy
```
1) Single sequence (Protenix server):
```
python cli/fetch_msa.py --seq AVEVL...K --out_dir ./msa_job --server protenix \
  --auth_user <username> --auth_pass <password> --email you@example.com
```

2) Multimer (use ':' to separate chains):
```
python cli/fetch_msa.py --seq AAAA:BBBBBBBB --out_dir ./msa_job --server protenix \
  --auth_user <username> --auth_pass <password>
```

3) From FASTA (multiple entries => multimer):
```
python cli/fetch_msa.py --fasta complex.fasta --out_dir ./msa_job --server colabfold
```

4) Build a single “paired multimer A3M” for inspection or third-party tools:
```
# by species (Protenix / local pairing logic)
python pair_msa.py msa_job/msa/0/pairing.a3m msa_job/msa/1/pairing.a3m --out msa_job/paired_concat.a3m --strategy species
# by index (ColabFold /ticket/pair row correspondence)
python pair_msa.py msa_job/msa/0/pairing.a3m msa_job/msa/1/pairing.a3m --out msa_job/paired_concat_idx.a3m --strategy index
```

Options
- `--server {protenix,colabfold}`: choose backend.
- `--host_url URL`: override server host.
- `--pairing {off,greedy,complete}`: request pairing alignment (remote), default `greedy` except `off` means non-pairing only.
- `--email EMAIL`: optional contact email passed to the server.
- `--auth_user/--auth_pass`: HTTP basic auth (Protenix server).

Output Layout
```
<out_dir>/
  msa/
    0/
      non_pairing.a3m
      pairing.a3m
    1/
      ...
  paired_concat.a3m   # (optional) produced by pair_msa.py
```

Single-Job Orchestrator (recommended for ColabFold)
```
# Automatically runs two-step for ColabFold (non-pairing per chain, then pairing for complex)
python cli/run_job.py --seq AAAA:BBBB --server colabfold --out_dir ./msa_ab --pairing greedy \
  --email you@example.com

# Protenix (one-step internally)
python cli/run_job.py --seq AVEVL...K --server protenix --out_dir ./msa_mono \
  --auth_user <username> --auth_pass <password>
```

Batch Processing (with rate limiting)
Prepare a CSV or TSV with columns:
```
name,seq
job1,AAAA:BBBB
job2,AVEVL...K
```
or
```
name,fasta
job3,/path/to/complex.fasta
```
Run the batch fetcher (safe defaults: one submit per 5 seconds, extra 2 seconds between jobs):
```
python cli/fetch_msa_batch.py --input jobs.csv --out_root ./msa_batch --server protenix \
  --auth_user <username> --auth_pass <password> --min_interval_s 5 --sleep_between_jobs 2
```

Single-Chain Retrieval
- Protenix server (recommended: provide BasicAuth and an email):
```
python cli/fetch_msa.py \
  --seq AVEVL...K \
  --server protenix \
  --out_dir ./msa_mono \
  --auth_user <username> --auth_pass <password> \
  --email you@example.com
```
Result: `./msa_mono/msa/0/{non_pairing.a3m, pairing.a3m}` (pairing derived locally using taxid mapping).

- ColabFold server (non-pairing):
```
python cli/fetch_msa.py --seq AVEVL...K --server colabfold --pairing off --out_dir ./msa_mono_cf
```
Result: `./msa_mono_cf/msa/0/non_pairing.a3m` and a minimal `pairing.a3m` stub.

Multi-Chain Retrieval (A:B)
- Protenix server (one command; pairing derived from tax mapping):
```
python cli/fetch_msa.py \
  --seq AAAA:BBBB \
  --server protenix \
  --out_dir ./msa_ab \
  --auth_user <username> --auth_pass <password>
```
Result: `./msa_ab/msa/0/` for A and `./msa_ab/msa/1/` for B, each containing `non_pairing.a3m` and `pairing.a3m`.

- ColabFold server (two-step recommended to mirror server semantics):
  1) Non-pairing per chain (env+uniref):
```
python cli/fetch_msa.py --seq AAAA --server colabfold --pairing off --out_dir ./msa_ab/A
python cli/fetch_msa.py --seq BBBB --server colabfold --pairing off --out_dir ./msa_ab/B
# Move or copy results to a common layout:
mkdir -p ./msa_ab/msa
cp -r ./msa_ab/A/msa/0 ./msa_ab/msa/0
cp -r ./msa_ab/B/msa/0 ./msa_ab/msa/1
```
  2) Pairing for the complex (server-side paired):
```
python cli/fetch_msa.py --seq AAAA:BBBB --server colabfold --pairing greedy --out_dir ./msa_ab
```
This writes `pairing.a3m` into `./msa_ab/msa/0` and `./msa_ab/msa/1`. You now have both `non_pairing.a3m` and `pairing.a3m` for each chain.

Optional: Export a Single “Paired Multimer A3M”
Most predictors expect per-chain directories; you don’t need a single merged file. For inspection/visualization you can export one:
```
# Species-based pairing (Protenix / local species logic)
python pair_msa.py ./msa_ab/msa/0/pairing.a3m ./msa_ab/msa/1/pairing.a3m \
  --out ./msa_ab/paired_concat.a3m --strategy species

# Index-based pairing (ColabFold pair semantics: row i across chains is paired)
python pair_msa.py ./msa_ab/msa/0/pairing.a3m ./msa_ab/msa/1/pairing.a3m \
  --out ./msa_ab/paired_concat_idx.a3m --strategy index
```

What “Pairing” Means
- Each chain has two inputs:
  - `non_pairing.a3m`: the chain’s own unpaired MSA (depth, background conservation).
  - `pairing.a3m`: the candidate set used to form cross-chain paired rows.
- In the model’s feature stage, cross-chain pairing is performed and merged into one assembly-level MSA:
  - For species-based pairing (Protenix/local): group by species ID (taxid), rank within each chain by similarity to the chain’s query, take K=min counts, and pair by rank to form K cross-chain rows per species. Overly large species (e.g., >600 rows) are skipped.
  - For index-based pairing (ColabFold): the server returns a `pair.a3m` container where row i of chain A’s block pairs with row i of chain B’s block; we split this into per-chain `pairing.a3m` and preserve row order.
- Final assembly-level MSA = [paired rows] stacked on top of [unpaired rows in block-diagonal form], with appropriate masks and cropping.

Notes & Tips
- Protenix backend requires BasicAuth; set `--auth_user/--auth_pass` and optionally `--email`.
- Public ColabFold has rate limits; run chains sequentially or with care.
- `paired_concat.a3m` is for inspection; for predictors, point each chain’s `precomputed_msa_dir` to its folder containing `non_pairing.a3m` and `pairing.a3m`.
- This project focuses on remote retrieval. For fully offline search (jackhmmer/colabfold_search with local DBs) you’d extend with local backends.
