Paired MSA Tools (standalone)

Proposed GitHub name: paired-msa-tools

Overview
- Standalone Python CLI to fetch Multiple Sequence Alignments (MSA) for monomers and multimers from remote services and to inspect/compose paired multimer alignments.
- No dependency on the Protenix repository. Prefer the single entry `query_msa.py` in this folder. Advanced users can use the CLIs below.
  - `cli/fetch_msa.py` — remote retrieval from Protenix server or ColabFold server; outputs per-chain `non_pairing.a3m` and `pairing.a3m` ready for structure predictors.
- `pair_msa.py` — optional inspector to export a single “concatenated paired multimer A3M” for visualization/debugging.
  - 支持 4 种配对策略：`species`（默认）、`index`、`rank`（非 species）、`stitch`（ColabFold 容器）。

Highlights
- Supports single-chain and multi-chain input (colon-separated `A:B[:C...]` or multi-entry FASTA).
- Protenix backend (default host: https://protenix-server.com/api/msa): returns `0.a3m` plus `uniref_tax.m8`; we split into `pairing.a3m` and `non_pairing.a3m` using UniRef100+taxid mapping.
- ColabFold backend (default host: https://api.colabfold.com): returns `env`+`uniref` A3Ms for non-pairing and a `pair.a3m` container for pairing;我们会拆分为每链文件；并可选择输出“Notebook 风格”的单一合并 A3M 与多种过滤。

Install
```
pip install -r requirements.txt
```

Input Formats
- Single job（二选一）：
  - `--seq` 直接传序列字符串，用冒号 `:` 连接多条链。示例：
    - 单链：`--seq AVEVL...K`
    - 多链二聚体：`--seq AAAA:BBBBBBBB`
    - 多链三聚体：`--seq AAAA:BBBB:CCCC`
    说明：程序会自动去掉空白并转为大写；链顺序即后续的链索引顺序（0,1,2,…）。
  - `--fasta` 传一个 FASTA 文件路径；文件中每个条目（每个 `>` header 开头的序列）视作一条链，按出现顺序确定链顺序。示例（双链 FASTA 文件 contents）：
    ```
    >chainA
    AAAA...
    >chainB
    BBBB...
    ```
    单链 FASTA 则只有一个条目。

- Batch（批量）
  - `--input` 读取 CSV/TSV，列头有两种形式（二选一）：
    - `name,seq`：其中 `seq` 使用与 `--seq` 相同的 `A:B[:C...]` 形式。
    - `name,fasta`：其中 `fasta` 是一个本地 FASTA 文件路径（可为绝对或相对路径）。
  - 默认分隔符是逗号，若为 TSV 则传 `--delimiter '\t'`。
  - 示例（CSV，按 `name,seq`）：
    ```
    name,seq
    job1,AAAA:BBBB
    job2,AVEVLK...
    ```
  - 示例（TSV，按 `name\tfasta`）：
    ```
    name	fasta
    complex1	/path/to/complex.fasta
    ```
  - 运行：
    ```
    # CSV 批量（Protenix）
    python query_msa.py --input jobs.csv --out_dir ./msa_batch --server protenix
    # TSV 批量（ColabFold），注意分隔符
    python query_msa.py --input jobs.tsv --delimiter '\t' --out_dir ./msa_batch --server colabfold
    ```

Notes
- 输入序列允许 `A–Z` 字母，程序会移除空白字符并转成大写；不需要也不建议在输入里包含 gap 或小写插入。
- 链顺序一旦确定，会体现在输出目录 `msa/<i>/` 的 `i` 上（0,1,...），也会体现在后续拼接的列顺序里。
- 运行完成后，可用 `pair_msa.py --auto_dir <作业目录>` 自动识别并导出合并 A3M（见下文 One-click autodetect）。

Quick Start
0) One-file entry (recommended):
```
# Single job (mono or multimer)
python query_msa.py --seq AAAA:BBBB --server protenix --out_dir ./msa_job \
  --auth_user <user> --auth_pass <pass>

# Batch (CSV/TSV with name,seq or name,fasta)
python query_msa.py --input jobs.csv --out_dir ./msa_batch --server colabfold
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

# non-species rank pairing (each chain sorted by similarity-to-query, zip by rank)
python pair_msa.py msa_job/msa/0/pairing.a3m msa_job/msa/1/pairing.a3m \
  --strategy rank --pad_missing --min_chains 2 --out msa_job/paired_concat_rank.a3m

# stitch (ColabFold pair.a3m + optional filters)
# Prefer raw pair.a3m; fall back to per-chain inputs when pair.a3m is unavailable
python pair_msa.py --strategy stitch \
  --pair_a3m ./msa_job/raw/pair.a3m \
  --out ./msa_job/paired_stitched.a3m

# With filters and Δgene (now default to get_msa-like values)
python pair_msa.py --strategy stitch \
  --pair_a3m ./msa_job/raw/pair.a3m \
  --min_coverage 0.75 --min_identity 0.15 --genomic_distance 1 \
  --out ./msa_job/paired_stitched_filt.a3m

# Fallback from per-chain pairing.a3m (index semantics + optional filters)
python pair_msa.py --strategy stitch \
  ./msa_job/msa/0/pairing.a3m ./msa_job/msa/1/pairing.a3m \
  --out ./msa_job/paired_stitched_from_index.a3m

Tip (common pitfall)
- 若使用 ColabFold 的合并容器 `pair.a3m`，请务必通过 `--pair_a3m` 传入；
  不要把 `pair.a3m` 当作位置参数给到 `a3m` 列表，否则会被当作“单链”输入，
  导致 Δgene 过滤被跳过，留下大量未过滤的条目（与 Notebook 行为不一致）。

One-click autodetect
- `--auto_dir <JOB_DIR>`：无需逐个传入文件，自动识别输入：
  - stitch：优先 `<JOB_DIR>/raw/pair.a3m`（ColabFold 配对容器）；若无，则回退到 `<JOB_DIR>/msa/<i>/pairing.a3m` 列表。
  - species/index/rank：使用 `<JOB_DIR>/msa/<i>/pairing.a3m`（按 i=0,1,... 排序）。
示例：
```
# Auto + stitch（推荐用于 ColabFold 配对任务）
python pair_msa.py --auto_dir ./msa_job --strategy stitch --out ./msa_job/paired_stitched.a3m

# Auto + species（或 index/rank）
python pair_msa.py --auto_dir ./msa_job --strategy species --out ./msa_job/paired_species.a3m
```
```

Options
- `--server {protenix,colabfold}`: choose backend.
- `--host_url URL`: override server host.
- `--pairing {off,greedy,complete}`: request pairing alignment (remote), default `complete`.
- `--email EMAIL`: optional contact email passed to the server.
- `--auth_user/--auth_pass`: HTTP basic auth (Protenix server).

Notebook-style stitched output and filters (new in query_msa.py)
- `--stitch_out PATH`: 额外写出一个“合并 A3M”（把多链在行上配对后串接为一条序列行），等价于 Notebook 的 `msa.a3m`。
- `--min_coverage FLOAT` (default: `0.75`): 覆盖度阈值（可用 0–1 或 0–100）。
- `--min_identity FLOAT` (default: `0.15`): 同一性阈值（0–1 或 0–100）。
- `--max_evalue FLOAT`: e-value 最大阈值。
- `--min_alnscore FLOAT`: 比对分数最小阈值。
- `--genomic_distance INT` (default: `1`): 仅对 ColabFold 多链+`pairing=complete` 生效，启用服务端 `paircomplete-pairfilterprox_Δgene`；同时在本地对 stitched 输出按 Δgene 做二次过滤。

Defaults aligned with get_msa.py
- 为了与 get_msa Notebook 的默认设置一致，本项目默认采用：`min_coverage=0.75`、`min_identity=0.15`、`genomic_distance=1`，且 ColabFold 多链默认 `--pairing complete`。
- 注意：若将 `min_coverage`/`min_identity` 写成百分数（如 75/15），程序会自动转换到 0–1 范围。

Notes
- 2025-11 Bugfix: 修正了 UniProt ID 提取在 `UniRef100_<ACC>_<taxid>` 头中的截断边界，避免把 `<taxid>` 混入 `<ACC>` 导致 Δgene 过滤失效。

pair_msa.py (stitch) options
- `--strategy stitch`: 启用基于 ColabFold pair 容器的行号配对语义，并输出单文件合并 A3M。
- `--pair_a3m PATH`: 指定原始 `pair.a3m`（优先）；未提供时可直接传每链 `pairing.a3m` 作为位置参数（将按索引对齐）。
- 过滤参数（默认 None 表示不过滤）：
  - `--min_coverage FLOAT`、`--min_identity FLOAT`、`--max_evalue FLOAT`、`--min_alnscore FLOAT`
  - `--genomic_distance INT`：在本地 stitched 输出上按“基因组邻近”进行二次过滤（需能解析到 UniRef/UniProt ID）。

non-species rank pairing（新）
- `--strategy rank`: 非 species 配对，每条链按“与 query 的相似度”降序排序，再按 rank 对齐。
- `--pad_missing`: 允许某些链缺失该 rank，缺失链将以全 gap 补齐；配合 `--min_chains` 限制至少出现的链数。
- `--min_chains`: 至少出现的链数（默认 2），低于该值的行会跳过。
- `--max_rows`: 限制输出的最大配对行数（可选）。

Which strategy to use
- `species`（Protenix 专长）：基于物种(TaxID)重新配对（本地构造配对关系）。
- `index`（最小拼接）：纯按行号对齐拼接，不做过滤。
- `stitch`（ColabFold 一致）：从 `pair.a3m` 还原“行号配对”并合并，可叠加过滤与 Δgene；无 `pair.a3m` 时等价于 `index` 语义但增强了过滤与统计。

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
  raw/                # 原始下载解包（用于 stitched 解析）

如果指定 `--stitch_out`，将额外生成：
```
<stitch_out>
```
一个 Notebook 风格的合并 A3M（多链序列行串接）。
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
python cli/fetch_msa.py --seq AAAA:BBBB --server colabfold --pairing complete --out_dir ./msa_ab
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

Notebook Parity（Notebook 流程在当前脚本中的映射）
- 获取/提交/轮询/下载：`utils/client.py` + `utils/retrieval.py::run_remote`（ColabFold 与 Protenix 通用，`--host_url` 可改 URL）。
- 解析与缝合（对应 Notebook 的 `_parse_paired_a3m` / `_parse_msa_lines` / `save_msa`）：`utils/stitch.py`
  - `parse_paired_a3m(<out_dir>/raw/pair.a3m)`：多链配对结果解析，保持“同位行配对”语义。
  - `parse_single_msas(<out_dir>/raw/)`：单链（env+uniref）合并并去重。
  - `save_stitched(..., min_coverage, min_identity, max_evalue, min_alnscore, max_genomic_distance)`：写出合并 A3M，并进行与 Notebook 相同的过滤与统计。
- `pair_msa.py --strategy stitch`：命令行等价封装，既可直接读 `raw/pair.a3m`，也可从每链 `pairing.a3m` 回退生成合并视图。
- 命令行入口：`query_msa.py --stitch_out` 会在常规每链落盘后，基于 `<out_dir>/raw/` 自动生成合并 A3M。

关于 Δgene（基因组邻近）
- ColabFold 多链 + `--pairing complete` 时：
  - 远端：通过 `--genomic_distance Δgene` 触发 API 模式 `paircomplete-pairfilterprox_Δgene`，服务端即按邻近约束构建候选配对。
  - 本地：若 `--stitch_out` 同时指定，会基于 UniRef/UniProt 头部解析出的 ID 计算链间“距离”，对合并行再做一次二次过滤。
- Protenix：`--genomic_distance` 会被忽略（当前服务端不支持该模式）；如需相似效果，请在 stitched 文件上另行筛选或告知我们整合策略需求。

使用示例（Stitched + 过滤）
- ColabFold（两条链、complete 配对、Δgene=1、基础过滤，并生成 `msa_merged.a3m`）：
```
python query_msa.py \
  --seq AAAA:BBBB \
  --server colabfold \
  --pairing complete \
  --genomic_distance 1 \
  --min_coverage 0.75 \
  --min_identity 0.15 \
  --stitch_out ./msa_job/msa_merged.a3m \
  --out_dir ./msa_job
```
- Protenix（两条链，生成合并 A3M；Δgene 忽略）：
```
python query_msa.py \
  --seq AAAA:BBBB \
  --server protenix \
  --auth_user <user> --auth_pass <pass> \
  --stitch_out ./msa_job/msa_merged.a3m \
  --out_dir ./msa_job
```

行为细节与差异
- 每链落盘始终可用（`<out_dir>/msa/<i>/{non_pairing.a3m, pairing.a3m}`），便于直接对接结构预测器。
- 指定 `--stitch_out` 后：
  - ColabFold 多链优先解析 `<out_dir>/raw/pair.a3m`；单链/非配对改解析 `<out_dir>/raw/{uniref.a3m,bfd.mgnify30.metaeuk30.smag30.a3m}`。
  - Protenix 没有 `pair.a3m`，使用每链 `pairing.a3m` 以“索引对齐”方式回退合并（行 i 跨链相接）。
- 过滤指标（coverage/identity/evalue/alnscore）依赖 mmseqs2 风格的头部 tab 字段；若缺失，对应阈值自动跳过。
- 统计信息会在终端打印“过滤前/过滤后”的链级均值与条目数量，便于快速评估阈值影响。

已知限制
- 未内置 Notebook 的本地缓存（可按需添加：基于序列+模式+Δgene 的 tar 包缓存目录）。
- Protenix 的合并当前使用“索引对齐”语义；如需基于“物种（taxid）”的行对齐与过滤（类似 `pair_msa.py --strategy species`），请提出需求，我们可在 `--stitch_out` 路径加入 `--stitch_strategy` 选项统一处理。

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
