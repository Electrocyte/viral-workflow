# Workflow Inventory — `viral-workflow` (ADVTIG untargeted viral screen)

- **Revision:** 8 (2026-10-02). It describes branch `workflow-audit` at commit `4af00a8`, and its new §0 summarises the current repository state. Earlier revisions:
  - revision 7 added §24, the current `run_coverage.py` interface decision;
  - revision 6 added §23, the B3/B4 code corrections;
  - revision 5 added §22, the repository packaging pass;
  - revision 4 (2026-10-01) added §21, the usability verification;
  - revision 3 incorporated the code-provenance, reference, Fourth-Pass and sequencing-input investigations.
- **Commits:**
  - `f1c0773`: the reconstructed files copied unchanged from SMART-CAMP: `pipeline/`, `q_control/`, `configs/`, `Pipfile`, `Pipfile.lock`, `CONFIGURATION-FILE-INPUTS.txt`.
  - `4af00a8`: the `run_coverage.py` interface correction, the stage-7 Fifth-Pass input correction, the protocol updates, the audit updates and the validation changes.
- **Companion files:** `WORKFLOW_INVENTORY.tsv` (per-file inventory, with per-file provenance columns since revision 5), `WORKFLOW_MISSING_DEPENDENCIES.tsv` (dependency completeness audit) and `USABILITY_BLOCKERS.tsv` (blockers, with status since revision 5).

**How to read this document.**
- **§0** gives the current state of the repository.
- **§1–§21** are the investigation and the 2026-10-01 usability audit. Their findings are kept as evidence. Where a later pass changed the repository, the statement is labelled **[Historical]** or **[Superseded]**, or carries a *Current:* note.
- **§22–§24** record what each remediation pass changed.

Four kinds of statement are kept apart:
- historical behaviour (what the code did in the past);
- earlier audit findings;
- current repository state;
- remaining remediation work.

**Method.** The investigation (§1–§21) was read-only: nothing in the workflow was executed, built, decompressed to disk or modified. The later passes (§22–§24) copied files from SMART-CAMP and made two narrow code corrections. Every validation in this document is static or interface-level only. **No biological analysis has been executed.** The sources examined:
- the repo copy at `/home/mangifera/viral-workflow`;
- the SMART-CAMP git repository at `/mnt/d/GitHub/SMART-CAMP` (working tree and full history; HEAD `65a94b3a`, branch `master`, clean for the ADVTIG files);
- `/home/mangifera/methylfilter` (working tree and history);
- the mounted Windows drives `/mnt/c,d,e,f`, using depth-bounded name searches;
- historical outputs under `/mnt/e/SequencingData/ADVTIG_primary` and `/mnt/e/SequencingData/ADVTIG`;
- reference files under `/mnt/e/SequencingData/ADVTIG_fa` and `/mnt/e/SequencingData/BlastN_libraries`, read as headers or streamed counts only.

**Notation.** **[Obs]** marks an observed fact. **[Inf]** marks an inference; its evidence is given alongside.

**Investigations that were not completed.** The tool-permission policy blocked these. They are listed so that absence of evidence is not read as evidence of absence.
1. **A forensic search for an uncommitted `run_coverage.py`:** git stashes, dangling objects, editor/Spyder autosaves, and shell history on the Windows drives.
2. **A dedicated Centrifuge/U-RVDB provenance search:** header-counting the surviving `U-RVDBvCurrent.fasta`, and so on. Parts of it were covered incidentally by the reference audit.
3. **Listing the contents of the two WSL backups:** `/mnt/f/ubuntu_backup.tar` (213 GB) and `/mnt/d/WSL_Backups/ubuntu.vhdx` (500 GB). These were judged too large and slow to scan in this pass.
4. *(Superseded.)* Committing and pushing these reports was blocked at revision 4. The work has since been committed and pushed to `origin/workflow-audit` as `f1c0773` and `4af00a8`.

---

## 0. Current repository state (branch `workflow-audit`, commit `4af00a8`)

**Verdict: NOT YET USABLE.** Several defects are now fixed: the missing packages, the misplaced module, the broken `run_coverage.py` parser and the stage 6→7 name mismatch. But the workflow still cannot be run by a new user. The environment is incomplete, paths are machine-specific, the reference resources are not prepared, two counter-screen panel sources are unresolved, and there is no current user-facing run documentation.

**[Current] What is now in the repository:**
- **Packages (`f1c0773`):** `pipeline/` holds 21 files (`__init__.py`, the full stage-1 import closure, and `map_genome_to_reads.py` with its helpers `bespoke_parse_map_reads.py` and `pre_trim_fq.py`). `q_control/` holds `adjust_q_scores.py` and `filter_q.py`. All are byte-identical copies from SMART-CAMP HEAD `65a94b3a`.
  - Every internal import of `mp_metagenomic_assessment_v4.py` (21/21) resolves inside this repository.
  - `pipeline/map_genome_to_reads.py` resolves its helpers (2/2).
  - **A SMART-CAMP checkout is no longer needed to satisfy these imports.** SMART-CAMP was the *provenance source* for the reconstructed files.
- **Sample sheets (`f1c0773`):** `configs/` holds the ADVTIG `viral_DNA_all9-2*` sheets and the Viral_FDA `viral_DNA_all16-2*` sheets. `CONFIGURATION-FILE-INPUTS.txt` (the existing column description) is at the repo root. All are copied unchanged.
  - The scripts and the run-book still refer to them by the old `/home/james/SMART-CAMP/configs/…` paths. This is part of the machine-specific path problem (B9).
- **Environment metadata (`f1c0773`):** `Pipfile` and `Pipfile.lock` are copied unchanged from SMART-CAMP as **historical** environment metadata.
  - They are **not** a complete current environment definition: they lack pysam and every command-line tool (B6, B7).
- **`run_coverage.py` interface (`4af00a8`, current workflow decision, §24):**
  - `-x` / `--extract-fastq` extracts FASTQ;
  - `--edit-distance VALUE` sets the edit-distance threshold, **default 0.15**, so the canonical reconstructed workflow uses 0.15;
  - there is **no `-e` option**.
  - Historically, `-e` meant FASTQ extraction and the threshold was fixed internally at 0.10. Commit `1c782e3e` then defined `-e` twice, and the script could not start (§6).
- **Stage 6 → 7 (`4af00a8`):** stage 6 (`count_the_screen.py --skip`) writes `MetaFilt-Fifth-Pass.csv`, and stage 7 (`collect_taxid_primary_counts.py`) now reads `MetaFilt-Fifth-Pass.csv`.
  - Historical `MetaFilt-Fourth-Pass*.csv` files were not renamed or modified. They remain historical artefacts.
  - Stage 7's final output name still contains Fourth-Pass (`Deduplicated-Read-Counts-MetaFilt-Fourth-Pass-FILTERED.csv`). That is a separate naming issue left for a later decision.
- **Run-book (`4af00a8`):** the current `run_coverage.py` commands in `ADVTIG-protocol.md` use `--edit-distance 0.15`, and the extraction commands use `--edit-distance 0.15 -x`. Historical commands are preserved as originally written and labelled HISTORICAL.

**[Current] Still present but misplaced:** the root-level `map_genome_to_reads.py` is a byte-identical duplicate of `pipeline/map_genome_to_reads.py`. It has not been removed. Its helper imports do not resolve from the repo root, so it should not be used (§22.6, D1).

**[Current] Blocker status** (detail in `USABILITY_BLOCKERS.tsv` and §21):

| ID | Blocker | Status |
|---|---|---|
| B1 | `pipeline/` and `q_control/` absent | **RESOLVED** (`f1c0773`) |
| B2 | `map_genome_to_reads.py` misplaced; its helpers absent | **RESOLVED** (`f1c0773`; root duplicate retained) |
| B3 | `run_coverage.py` duplicate `-e`, parser cannot be built | **RESOLVED** (`4af00a8`; new interface, §24) |
| B4 | Stage 6 writes Fifth-Pass, stage 7 read Fourth-Pass | **RESOLVED** for the stage 6→7 input dependency (`4af00a8`) |
| B5 | Sample-sheet format | **PARTIALLY RESOLVED**: sheets and the existing column description are present; no clean current user-facing specification yet |
| B6 | Python environment | **OPEN**: Pipfiles are historical and incomplete (no pysam) |
| B7 | Command-line tools and versions | **OPEN** |
| B8 | Current user-facing run-book, input layout, outputs | **OPEN** |
| B9 | Machine-specific hard-coded paths | **OPEN** |
| B10 | RVDB / Centrifuge reference preparation | **OPEN** |
| B11 | Mammalian counter-screen panels | **OPEN** |
| B12 | Viral and bacterial counter-screen panel source FASTAs | **OPEN** |
| B13 | Stage-1 `~` expansion, unchecked exit codes, sequencing-summary requirement | **OPEN** |

**[Current] Latest validation** (`audit/static_validation_results_cli.txt`; static and interface checks only, **not a biological execution test**):
- all 33 current Python files compile;
- the `run_coverage.py` parser constructs;
- `run_coverage.py --help` succeeds, using the validator's existing inert substitutes for pandas, pysam, numpy and Biopython, which are not installed (B6);
- the default edit-distance threshold is 0.15;
- `--edit-distance 0.15`, `-x` and `--edit-distance 0.15 -x` all parse correctly;
- all current supported documented commands parse (41 OK, 0 FAIL);
- historical commands are classified separately (26 still parse, 18 expected HFAIL);
- the only remaining FAIL is the `bash -n` syntax failure in `ADVTIG-untargeted/URVDB.md` (a Markdown table inside its code block). It is unrelated to these fixes and remains unresolved.

**Remaining remediation work:**
- B5: a clean current sample-sheet specification;
- B6/B7: a complete environment and pinned tool versions;
- B8: user-facing execution documentation;
- B9: configurable paths;
- B10–B12: reference preparation, the mammalian panels and the viral/bacterial panel sources;
- B13: the stage-1 robustness issues;
- decisions: the stage-7 output name and what to do with the root duplicate.

---

## 1. Executive summary

*The answers below are the investigation findings. Where the repository has since changed, a* Current: *note gives the present state (see §0).*

| Question | Answer |
|---|---|
| Is all workflow **source code** accounted for? | **Yes, in the repository.** *[Superseded finding at revision 4: `pipeline/` and `q_control/` were absent from the copy, `map_genome_to_reads.py` sat at the root, and no committed `run_coverage.py` revision accepted the documented command lines.]* **Current:** the packages and `pipeline/map_genome_to_reads.py` with its helpers are present (`f1c0773`), and imports resolve. `run_coverage.py` has a corrected interface (`4af00a8`, §24). The root `map_genome_to_reads.py` remains as a misplaced duplicate. |
| Is the workflow **operationally** complete? | **No.** These are absent and not regenerable from anything found: the Centrifuge index, the RVDB split collection and maps, the minimap2 panel indexes, the source FASTAs for the viral and bacterial "RefSeq" counter-screen panels, and **all raw sequencing reads**. **Current:** `Pipfile`/`Pipfile.lock` are now in the repo as historical metadata, but there is still no complete environment definition and no current user-facing run documentation. |
| Fourth Pass | **Generated by a script, now resolved.** `count_the_screen.py --skip` wrote `MetaFilt-Fourth-Pass.csv` from commit `dc1ef791` (2025-12-18) until `a68d1880` (2026-02-02), which renamed the output to `MetaFilt-Fifth-Pass.csv`. It is a merge that keeps every row. Every historical version is explained 100% by the code at its date. The spreadsheet copies were filtered and formatted only. |
| `glue.py` | Found at `methylfilter/Power/glue.py`. **It is not a workflow dependency.** The phrase "same as glue.py" is a help-text leftover. `deep_cov.py`'s metagenomic filter was adapted from methylfilter's `filter_short_hitlength_backend.py`. |
| 15% edit-distance behaviour | **[Historical]** No committed SMART-CAMP `run_coverage.py` ever accepted `-e 0.15`. The historical 15% outputs are nevertheless **plausibly explained by committed code**, through `deep_cov.py -e 0.15 -dd -x`, which calls `run_coverage` *functions* rather than its CLI (§6). An uncommitted file is possible but not required to explain them. **Current workflow decision (§24):** `run_coverage.py --edit-distance VALUE`, default **0.15**; extraction is `-x`. |
| Raw sequencing inputs | **Not on any searched drive.** ADVTIG: only derived outputs survive, and the logs record the exact raw paths on the unavailable `/mnt/usersData`. Viral_FDA 2024-02-16: no trace at all. The two unsearched WSL backups are the remaining candidates. |
| Exact historical reproduction | **Not currently possible.** |
| Functionally equivalent rebuild | **Possible, given decisions and inputs** (§20). |
| Current usability | **NOT YET USABLE** (§0, §21). B1–B4 are resolved and B5 is partially resolved; B6–B13 are open. |

---

## 2. Repository structure (current state at `4af00a8`)

```
viral-workflow/                          git: branch workflow-audit @ 4af00a8 (pushed); f1c0773 = copied files
├── README.md                            16 B   title only
├── ADVTIG-protocol.md                   26.3 KB run-book (SMART-CAMP HEAD + 4af00a8: current run_coverage lines use '--edit-distance 0.15' [-x]; historical lines labelled; §24)
├── mp_metagenomic_assessment_v4.py      33.5 KB stage 1 (= SMART-CAMP, md5 83e36075…)
├── map_genome_to_reads.py               39.3 KB ⚠ MISPLACED DUPLICATE of pipeline/map_genome_to_reads.py (not removed; do not use; §22.6)
├── pipeline/                            21 files (f1c0773): __init__.py + stage-1 closure + map_genome_to_reads.py and its 2 helpers
├── q_control/                           2 files (f1c0773); namespace package, no __init__.py, as in source
├── configs/                             10 sample sheets (f1c0773): ADVTIG viral_DNA_all9-2*, Viral_FDA viral_DNA_all16-2*
├── CONFIGURATION-FILE-INPUTS.txt        7.4 KB existing sample-sheet column description (f1c0773)
├── Pipfile, Pipfile.lock                historical SMART-CAMP environment metadata (f1c0773; unchanged; incomplete, no pysam/CLI tools)
├── ADVTIG-untargeted/                   9 files (= SMART-CAMP HEAD except run_coverage.py [interface, §24], collect_taxid_primary_counts.py [Fifth-Pass input, §23], URVDB.md [historical note]; 4af00a8)
├── audit/                               static checks and their results
├── WORKFLOW_INVENTORY.md / .tsv, WORKFLOW_MISSING_DEPENDENCIES.tsv, USABILITY_BLOCKERS.tsv   (these reports)
└── .git/
```

**[Historical: revision-4 finding, superseded by `f1c0773`.]** At revision 4 these were missing from the copy, though present in SMART-CAMP:
- `pipeline/` (`__init__.py` plus 12 imported modules and their second-level imports `dep_pre_chunker`, `updated_blastn_interpreter_v2`, `rank_BLAST_predictions_for_aa_v2`, `cat_nanostats`, `find_barcode_reads_for_QC`, `bespoke_parse_map_reads`, `pre_trim_fq`);
- `q_control/` (`adjust_q_scores.py`, `filter_q.py`; this is a namespace package with no `__init__.py`);
- `configs/` sample sheets;
- `Pipfile` and `Pipfile.lock`.

Before `f1c0773`, `mp_metagenomic_assessment_v4.py` failed with `ImportError: pipeline`, and `map_genome_to_reads.py` failed on `import bespoke_parse_map_reads, pre_trim_fq`. **Current:** all of these files are present, and their imports resolve (§22.4).

---

## 3. Reconstructed execution sequence (stage by stage)

Notation:
- `DIR` is the dataset root, e.g. `/mnt/usersData/ADVTIG_v3_untargeted/`.
- Sample names are `{date}_{NA}_{strain}_{CFU}_{batch}_{duration}`.
- Availability codes (current state at `4af00a8`): **P** present in this repo; **S** present only in SMART-CAMP; **E** exists elsewhere on this machine; **R** regenerable from available material; **G** generated by an earlier stage; **X** absent; **?** unknown.
- Stage commands are shown in their current form where one exists. Machine-specific paths are unchanged (B9).

| # | Command (protocol) | Required inputs | Reference / index | Config | Outputs → consumer | Availability |
|---|---|---|---|---|---|---|
| 0a | URVDB.md L45–107: clean `U-RVDBvCurrent.fasta` → `U-RVDB-clean-fix2.fasta`; build `fixed_URVDB_seqID2.map`; `centrifuge-build --conversion-table RVDB2_seqID.map --taxonomy-tree nodes.dmp --name-table names.dmp U-RVDB-clean-fix2.fasta mini-u-viral2 -p 8` | U-RVDBvCurrent.fasta; nucl_gb.accession2taxid → RVDB2_seqID.map; taxdump | — | — | `mini-u-viral2.{1,2,3[,4]}.cf` → stage 1; clean FASTA + seqID map → 0b and stages 2+ | Source FASTA **E** (candidate, unverified); accession2taxid **E** (2023-07 snapshot); taxdump **E** (2023-04/07); `.cf` **X** → **R** (not byte-identical) |
| 0b | `split_fa.py`, with the configuration at **commit `7f6e9993`** (2025-01-17), not the current configuration | U-RVDB-clean-fix2.fasta, fixed_URVDB_seqID2.map | — | — | `URVDB_split_fa/` (one FASTA per sequence) + `U-RVDB-clean-fix2-index.json` → stages 2–6 | **R** from 0a outputs; current `split_fa.py` is configured for PLSDB |
| 0c | ADVTIG-protocol.md L154–357: `cat` host genomes + viral targets → `all_mammalian*.fa`; `minimap2 -d`; header-derived TSVs | 5 host genomes, 11 viral FASTAs, NL63, lambda, `viral_input-sequences.fna`, `ss_bacteria_sequences.fna` | — | — | `.mmi` + `.tsv` panels → stage 6 | See §9: mammalian panels **R**; viral and bacterial RefSeq panels **X** (sources absent) |
| 1 | `mp_metagenomic_assessment_v4.py -d DIR -t 5 -ci "~/SequencingData/Centrifuge_libraries/viral/" -c -bl uviral25-2 -m -rs 1 -fd viral_DNA_all9-2.txt -skip` | Raw ONT reads at `{DIR}/{sample}/{RUN}/fastq_pass/*.fastq[.gz]` (exactly one run-dir level, because `**` is non-recursive) and `{DIR}/{sample}/{RUN}/*sequencing*summary*.txt` | Centrifuge `{ci}/mini-u-viral2.*.cf` | `.txt` sample list (P, `configs/`) | `trimmed/trimmed_{s}.fastq` and `centrifuge/{s}_mini-u-viral2_centrifuge_troubleshooting_report.tsv` → stages 2–3; NanoPlot output; `nanoplot_summary_data.csv` | Code **P** (imports resolve since `f1c0773`); reads **X**; index **X**; tools centrifuge 1.x, porechop and NanoPlot are not installed |
| 1b | `pipeline/map_genome_to_reads.py -d DIR -s …txt -t lambda -p 6` | `trimmed/trim*q` | `/home/james/SequencingData/Centrifuge_libraries/ecoli-lambda-phage.fasta` | `.txt` | `analysis/bespoke/…` → **no downstream consumer** (QC side branch) | Code **P** at `pipeline/map_genome_to_reads.py` with its helpers (`f1c0773`); root copy is a misplaced duplicate; lambda candidate **E** (3.56 kb, provenance uncertain) |
| 2 | `run_coverage.py -s …csv -d DIR -o DIR --database uviral25-2 --edit-distance 0.15` with modes `-x` (extract), `-r` (align), none (stats), `-l` (labels). *Historical form: `-e` = extract, ED fixed at 0.10.* | troubleshooting report, trimmed FASTQ | `URVDB_split_fa/`, index JSON, `fixed_URVDB_seqID2.map`, `RVDB2_seqID.map` | `.csv` sheet (P, `configs/`) | `fq_seqID_uv25_2/*.fastq`; `map-ont/{s}/{seqID}/*` (sam/bam/mpileup/coverage/editdist); `uviral25-2_mapping_stats{,_labelled}.csv` → stage 3 | **P.** Parser constructs and all current commands parse (`4af00a8`, §24). *[Historical: SMART-CAMP HEAD `1c782e3e` cannot start (duplicate `-e`); `aff9daba` worked with the old form at a fixed ED 0.10.]* |
| 3 | `deep_cov.py … -e 0.15 -dd`, then `-dd -x`, then `-dd -a` | stage-2 outputs + troubleshooting report | same DB files | `.csv` | `analysis/TaxID_SeqID_Sample.csv`, `MetaFilt-First-Pass…`, `map-ont/TaxID_uviral25-2_mapping_stats*`, `Third-Pass-FULL-…filtered.csv`, `MetaFilt-Third-Pass-…labelled.csv` → stage 4 | **P**; works as committed (since `1c782e3e`); its own CLI keeps `-e FLOAT` and `-x` |
| 4 | `filt_low_complexity.py -o DIR` | Third-Pass labelled CSV + primary BAMs | — | — | `MetaFilt-Third-Pass-…-complexity.csv` → stages 5–6 | **P** |
| 5 | `counter_screen.py -o DIR -s …csv -p 15 -c [--panel new_combined / refseq / bacterial_refseq] [-r] -e 0.15` | complexity CSV + primary BAMs | `all_mammalian_targets2.mmi` + `mammalian_ref_species2.tsv` (new_combined); `viral-input-sequences.mmi/.tsv` (refseq); `bacterial-input-sequences.mmi/.tsv`; the `split` default also needs `ADV5.mmi`, `all_mammalian.mmi`, lambda `.mmi` | `.csv` (needs `Spike_species`) | `counter_screen_*/…` PAFs/stats; `…-counter-screened-summary-{combined3, refseq, refseq-split, bacterial-refseq*}.csv` → stage 6 | Code **P**; indexes **X** (§9) |
| 6 | `count_the_screen.py -o DIR -s …csv -e 0.15`, then `--skip` | counter-screen outputs | viral/bacterial TSVs | `.csv` | `counter_screen_refseq/…-refseq-viral.csv`, `…bacterial.csv`; with `--skip`: **`MetaFilt-Fifth-Pass.csv`** (historically **`MetaFilt-Fourth-Pass.csv`**), `-VIRAL-ONLY`, `-BACTERIAL-ONLY` → stage 7 | **P** |
| 7 | `collect_taxid_primary_counts.py` ×3 | `map-ont/{s}/{seqID}/*taxID*.0*primary*p.bam`; **`MetaFilt-Fifth-Pass.csv`** (since `4af00a8`; historically `MetaFilt-Fourth-Pass.csv`) | — | 6 samples hard-coded | `sample_taxid_bams.json`, `sample_taxid_dedup_primary_counts.csv`, `Deduplicated-Read-Counts-MetaFilt-Fourth-Pass-FILTERED.csv` (final) | **P.** Its input now matches stage 6's output (`4af00a8`). It uses key columns only. ⚠ Its final output name still says Fourth-Pass; this naming issue awaits a later decision (§23.5). |

**[Obs] Facts about stage 1** (file:line citations refer to the SMART-CAMP originals, which are byte-identical to the copies now in `pipeline/` and the repo root):
- It runs NanoPlot (`dep_run_nanostat_analyses.py:60-72`), porechop (`mp_trim_reads_v2.py:309-312`), then Centrifuge.
- The Centrifuge call is `centrifuge -q -x {ci}/mini-u-viral2 … -p 25` (`dep_flat_sample_classifier_test.py:306,324-341`). There are 5 processes of 25 threads each, and each process loads its own copy of the index.
- With `-skip` and no `-hn`, there is no host removal, no BLAST and no KrakenUniq.
- The `-ci` value is double-quoted, so the literal `~` reaches Centrifuge. No code calls `expanduser`.
- No tool return codes are checked. The run ends with `pd.read_csv(nanoplot_summary_data.csv)` (`:721`), which fails if no sequencing summary exists.
- Trimming writes intermediates **into the raw-read tree**, which must therefore be writable.
- Only the `-bl uviral25-2` run produces the file names that downstream code consumes. The `uviral25` and `c-viral29` runs in the protocol are historical.

---

## 4. Code provenance

| Component | Status | Provenance evidence |
|---|---|---|
| `mp_metagenomic_assessment_v4.py` | **Found.** In the repo, = SMART-CAMP HEAD | md5 identical. *[Historical: its imports resolved only inside SMART-CAMP.]* **Current:** they resolve inside this repo (`f1c0773`). |
| 12 `pipeline.*` + 2 `q_control.*` modules, plus 5 second-level modules | **Found in SMART-CAMP; now copied into the repo** (`f1c0773`, byte-identical) | recursive import trace; all import only stdlib, pandas and numpy |
| `pipeline/map_genome_to_reads.py` | **Found.** Now at `pipeline/map_genome_to_reads.py`, with its helpers `bespoke_parse_map_reads.py` and `pre_trim_fq.py` (`f1c0773`) | md5 identical. *[Historical: it had been copied only to the repo root, without its helpers.]* The root copy remains as a duplicate. Its outputs are consumed by nothing downstream. |
| `ADVTIG-untargeted/*.py` | = SMART-CAMP HEAD, except `run_coverage.py` (interface) and `collect_taxid_primary_counts.py` (Fifth-Pass input), both changed in `4af00a8` | md5 identical for the other files; changes are documented in §23–§24 |
| `glue.py` | Found at `/home/mangifera/methylfilter/Power/glue.py`. **Not a dependency.** | Never present in SMART-CAMP history. No ADVTIG script imports it. The help text "(same as glue.py)" was copied on 2025-12-11 from methylfilter's `light_fast_aligner.py` (`df46570`). `glue.py` itself builds panel JSONs; its `-mhl/-mcfm/-msf` only name an output path. |
| `deep_cov.py` metagenomic filter | **Adapted** from methylfilter `Power/filter_short_hitlength_backend.py::summarize_troubleshooting` (`ca6fbd3`, 2025-08) | Same per-read collapse code and flag names. The "Middle" preset (mhl 50, mcfm 0.05, msf 0.75) equals methylfilter's defaults (`cat_align.py:276-277`). |
| Producer of `MetaFilt-Fourth-Pass.csv` | **Found in history:** `count_the_screen.py --skip`, from `dc1ef791` (2025-12-18) to `a68d1880` (2026-02-02) | `git log -S 'Fourth-Pass'`. HEAD writes the same table as `MetaFilt-Fifth-Pass.csv`. |
| `split_fa.py` configuration that produced `URVDB_split_fa/` + `U-RVDB-clean-fix2-index.json` | **Found in history:** `2e4de9c9` → `7f6e9993` (2025-01-17) | From `0ee9d9da` (2025-01-23) onward the file is configured for C-RVDB, U-RVDB v29 or PLSDB. |
| A SMART-CAMP `run_coverage.py` that accepts `-e 0.15 …` | **[Historical] Never committed** (§6). **Current:** this repo's `run_coverage.py` has a corrected interface (`--edit-distance`, `-x`; §24). | 228 revisions scanned |

---

## 5. Superseded conclusions from earlier revisions

| Earlier statement | Corrected finding |
|---|---|
| "`mp_metagenomic_assessment_v4.py` missing" | Present in the copy. At revision 4 its 14 imported modules (and 5 second-level modules) were in SMART-CAMP only. **Current:** all are in the repo (`f1c0773`). |
| "`pipeline/map_genome_to_reads.py` missing" | At revision 4 it was present only at the repo root, without its 2 sibling modules. **Current:** it is at `pipeline/` with its helpers (`f1c0773`); the root copy remains as a duplicate. |
| "`glue.py` missing" | Found. It is not a dependency, so this is documentation drift only. |
| "No producer of `MetaFilt-Fourth-Pass.csv`; probably manual" | **Wrong.** The producer was `count_the_screen.py --skip`, whose output was renamed to Fifth-Pass on 2026-02-02 (§7). |
| "U-RVDB source FASTA not found" | A candidate survives: `/mnt/e/SequencingData/BlastN_libraries/NCBI_RVDB/virus_RVDBs/U-RVDBvCurrent.fasta` (5.43 GB, mtime 2020-01-09), with a second copy in `/mnt/d/Dropbox/Data/ViralDB/`. Whether it is the snapshot actually used is unverified. |
| "NCBI taxonomy and accession2taxid not found" | Found: `/mnt/e/SequencingData/genbank/{taxonomy/names.dmp,nodes.dmp,nucl_gb.accession2taxid}` (2023-07-17) and `/mnt/e/SequencingData/BlastN_libraries/RVDBs/{names,nodes}.dmp` (2023-04-20). |
| "HCoV-NL63 FASTA missing" | NC_005831.2 can be extracted from `U-RVDBvCurrent.fasta` or `C-RVDBv21.0.fasta`; the header needs reformatting. |
| "Lambda FASTA missing" | A candidate exists: `/mnt/d/Dropbox/AA SMART/fourier_data/prior_runs/ecoli-lambda-phage.fasta`. It is 3,560 bp with header `>E.coli Lambda phage`. **[Inf]** It is the ONT DNA control-strand fragment, which matches the `DNA_CS` label. It is not the 48.5 kb NC_001416.1 genome. |
| "Historical 15% outputs require an uncommitted script" | **Not required** (§6). Plausibly produced by committed code through `deep_cov.py`. |

---

## 6. `run_coverage.py` history and the 15% edit-distance question

> **[Historical section.]** This describes SMART-CAMP history. The current interface in this repository is `-x`/`--extract-fastq` and `--edit-distance VALUE` (default 0.15), with no `-e` option. It is a workflow decision recorded in §24 (`4af00a8`), not a recovered historical flag name.

**[Obs]** The table below covers the 228 commits from `c6b3dd51` (2025-01-17) to `1c782e3e` (2026-02-12). The file was never renamed, and no other branch, tag or path contains it.

| Revisions | `-e` | ED threshold |
|---|---|---|
| `9142c431` → `aff9daba` (2025-01-20 → 2026-02-12 10:36) | one `-e`/`--extract-fastq`, `store_true` | hard-coded `target_threshold=0.10`. The threshold list `[…0.1, 0.15, 0.2, 0.5]` has existed since `16802532` (2025-11-07), so every `*-editdist.tsv` contains a 0.15 row. |
| **`aff9daba`** (2026-02-12 10:36, "update the edit distance to 15%") | unchanged | **Still 0.10.** It only changed the condition for regenerating the editdist TSV. This is the **last committed revision that starts.** |
| **`1c782e3e`** (2026-02-12 13:08, "update the edit distance to 15% – first test …") = HEAD | **adds** `parser.add_argument("-e", "--edit-distance-threshold", type=float, default=0.1)` beside the existing `-e/--extract-fastq` → `argparse.ArgumentError` | The threshold becomes a parameter. `needs_ed_rerun` compares the stored and requested threshold. **This is the commit that introduces configurable ED and the conflict.** The same commit rewrote the protocol to `-e 0.15 …`. |

**Other relevant findings:**
- **[Obs]** No revision accepts the protocol's `run_coverage.py … -e 0.15 [-e|-r|-l]` lines:
  - before `1c782e3e`, `-e 0.15` gives "unrecognized arguments: 0.15";
  - from `1c782e3e` onward, the parser cannot be constructed.

  The pre-`1c782e3e` protocol form (`-e`, `-r`, none, `-l`) works with `aff9daba`.
- **[Obs]** In the same commit `1c782e3e`, `deep_cov.py` resolved the identical clash correctly: it **renamed its extract flag from `-e` to `-x`** and added `-e` as a float. Its `-x` mode calls `run_coverage.extract_mapping_data(...)` with the threshold. `deep_cov` imports `run_coverage` functions only, so the broken `run_coverage` CLI is never built.
- **[Obs]** `a8951eed` (2026-02-10, "fix ED stringency") changed the ED ratio denominator from total reads to reads that carry an NM tag. That change affects results independently of the 0.10 → 0.15 change.
- **[Inf, moderate–high]** The historical 15% outputs can be explained entirely by committed code. These are `MetaFilt-Fifth-Pass-15percentED.csv` (2026-02-13) and `MetaFilt-Fifth-Pass.csv` (2026-03-09).
  - The first-pass `run_coverage` stages had already been run at ED 0.10.
  - The 15% re-extraction was then done by `deep_cov.py -e 0.15 -dd -x` and `-a`, which read the 0.15 row of the existing editdist TSVs.
  - After that came `counter_screen.py -e 0.15` and `count_the_screen.py -e 0.15`.
  - The `run_coverage.py -e 0.15` protocol lines were committed untested ("first test") and could never have run.
  - **Caveat:** the outputs were produced on a remote server, so matching commit times to file mtimes cannot prove which code ran.
- **[Obs]** The two 15% outputs also differ in counter-screen thresholds. `7e2073ab` (2026-02-10) relaxed Depth 1X to **≥ 0.35**, and `605181b8` (2026-03-09) restored **≥ 0.6**. `15percentED.csv` therefore used Depth ≥ 0.35, while the 03-09 Fifth Pass used ≥ 0.6.
- **Intended fix (revision-4 analysis).** The intent could be inferred with high confidence: a configurable float ED threshold, with extraction on a different flag. The author's exact flag name could not be, because it is not recorded anywhere, and the protocol line `-e 0.15 -e` cannot be interpreted literally under any parser. **Current:** resolved by an explicit workflow decision (§24): `-x` for extraction, following `deep_cov.py`, and `--edit-distance` (default 0.15) for the threshold.
- **Not searched (blocked):** stashes, dangling objects, editor autosaves and shell history for an uncommitted variant.

---

## 7. MetaFilt Fourth Pass: the data transformation

**[Obs] Producer.** `count_the_screen.py --skip` (`dc1ef791` → `a68d1880`) produces it by two left-merges on (Sample, SeqID, TaxID), with **no row filtering**:
1. the complexity table (`…-complexity.csv`) left-merged with the viral RefSeq counter-screen columns, producing `…-counter-screen-refseq.csv`;
2. that result left-merged with 53 `bcols` columns (panel hits, per-species `*_reads`/`*_pct`, `mammal_total_*`) from `…-counter-screened-summary-combined3.csv`.

Only rows that pass the `counter_screen.py` thresholds have non-empty panel columns.

**Historical versions.** All contain the same 6 ADVTIG samples, with no duplicate keys.

| File (location) | mtime | Rows | Cols | Rows with counter-screen data | Rule explaining which rows have counter-screen data (100% of rows) |
|---|---|---|---|---|---|
| `ADVTWG-MetaFilt-Fourth-Pass-2.csv` (Prior analyses) | 2025-12-22 | 1448 | 111 | 338 | lcf < 0.2, Depth1X ≥ 0.8, Ident ≥ 80, Mapped > 1, EDratio ≥ 0.7 (Dec-2025 code) |
| `MetaFilt-Fourth-Pass_v3.csv` (Prior analyses) | 2026-01-26 | 1513 | 112 | 629 | ≤ 0.2 / 0.7 / 70 / > 1 / 0.6 (`a2ae4441`) |
| `MetaFilt-Fourth-Pass_v4.csv` (ADVTIG_primary) | 2026-01-29 12:43 | 1513 | 112 | 650 | ≤ 0.2 / 0.6 / 60 / > 1 / 0.6 (`e4617a07`, 11:47) |
| `MetaFilt-Fifth-Pass.csv` (ADVTIG_primary) | 2026-03-09 | 1746 | 120 | 950 | the v4 rule, ED 0.15 |

**Changes between versions:**
- **Fourth-Pass-2 → v3:**
  - adds the column `HighComplexityGreater` (`17346697`);
  - adds 65 rows and removes none;
  - fills 64 rows whose ED values were previously NaN (the 01-23 "null ED rerun" commits);
  - The 65 new rows had empty complexity data while `low_complexity_frac` was 0.0, so 44 of them **wrongly passed** the counter-screen.
- **v3 → v4:**
  - same keys and columns;
  - complexity recomputed for those 65 rows (the BAM-glob fix in `e4617a07`);
  - counter-screen rows: 19 removed, 40 added (thresholds relaxed).
- **Lineage [Inf]:** v3 and v4 are **sibling runs**; v4 is not derived from the v3 file. Fifth-Pass is a later run of the same step, not computed from v4.

**Spreadsheet copies:**
- `v3_labelled.{csv,xlsx}` is a 72-column subset with identical values.
- `v4.xlsx` has the same values as `v4.csv`, sorted differently. Its yellow header row and colour scales are formatting only. Its AutoFilter (`Depth ≥ 0.6, Ident ≥ 60, Mapped ≥ 2, ED ≥ 0.6, lcf ≤ 0.2`) shows exactly the 650 rows the script rule selects.
- No hand-typed values or manual columns were found in any Fourth-Pass file.

**Consumer.** `collect_taxid_primary_counts.py` (stage 7) uses **key columns only**. *Current:* since `4af00a8` it reads `MetaFilt-Fifth-Pass.csv`, the file current stage 6 writes. The historical Fourth-Pass files described here were not renamed or modified. The observations below concern the historical Fourth-Pass runs:
- `Deduplicated-…_v3` and `_v4` are byte-identical: 594 (Sample, TaxID) pairs, equal to *all* Fourth-Pass keys, not only the 250 counter-screen-passing pairs.
- So "FILTERED" means "restricted to Third-Pass keys". It does not mean "passed the counter-screen".

**Most recent historical version:** **v4**. Its code state is `e4617a07`, and it is reproduced exactly inside `FifthPass_Charley3.xlsx` (2026-02-10; 1513 rows, ED 0.10). Whether v4 or the current Fifth Pass (ED 0.15, 1746 rows) is authoritative is a scientific decision for the user.

**Answers to the specific questions:**
- Generated manually? **No.**
- Omitted script? **No.** It came from a committed script, later renamed.
- Spreadsheet filtering? **Only as a view.** The filter reproduces the script rule and does not create the file.
- Reproducibly inferable? **Yes:** both by code (historical revisions) and empirically (100% rule match).

---

## 8. Centrifuge / U-RVDB provenance

**[Obs] Build recipe (URVDB.md L45–107; also `Seq commands.md:700-735`):**
1. `U-RVDBvCurrent.fasta` → awk keeps field 3 of the `|`-separated header → `U-RVDB-clean-fix2.fasta`.
2. `MAP=/mnt/usersData/RVDB2/RVDB2_seqID.map` (= `awk 'NR>1{print $2"\t"$3}' nucl_gb.accession2taxid`).
3. `centrifuge-build --conversion-table $MAP --taxonomy-tree nodes.dmp --name-table names.dmp U-RVDB-clean-fix2.fasta mini-u-viral2 -p 8`.

The `fixed_URVDB_seqID2.map` used later for labels comes from URVDB.md L72–91.

**[Obs] Surviving candidate material:**

| Item | Location | Notes |
|---|---|---|
| `U-RVDBvCurrent.fasta` | `/mnt/e/SequencingData/BlastN_libraries/NCBI_RVDB/virus_RVDBs/` and `/mnt/d/Dropbox/Data/ViralDB/` | 5,434,153,675 B, mtime 2020-01-09. It sits under the same `NCBI_RVDB/virus_RVDBs/` sub-path that URVDB.md uses. |
| `nucl_gb.accession2taxid` | `/mnt/e/SequencingData/genbank/` | 11.9 GB, 2023-07-17 |
| `names.dmp`, `nodes.dmp` | `/mnt/e/SequencingData/genbank/taxonomy/` (2023-07-17) and `/mnt/e/SequencingData/BlastN_libraries/RVDBs/` (2023-04-20) | — |
| `C-RVDBv21.0.fasta`, `C-RVDBv25.0.fasta` | `/mnt/e/SequencingData/BlastN_libraries/…` | — |
| `C-RVDBv29.0.fasta.gz` | `/mnt/e/SequencingData/ADVTIG_fa/` | — |
| `seqid2taxid.map` | `ADVTIG_fa/centrifuge-build/` | 8.08 M lines; **[Inf]** an unclustered RVDB of the ~v25 era (2023). It is *not* the `RVDB2_seqID.map`, whose format has no `acc\|` prefix. |

**Not found:** any `.cf` file for an RVDB index, `U-RVDB-clean-fix2.fasta`, `fixed_URVDB_seqID2.map`, `RVDB2_seqID.map`, `URVDB_split_fa/`, or the index JSON. The only `.cf` files on the machine are unrelated: `v_f_b4` and `fungus`.

**Which release?** **Unknown; ambiguous.**
- The key name `uviral25` suggests U-RVDB v25.0, which has 8,079,113 sequences per URVDB.md's own table.
- The recipe says "vCurrent", and URVDB.md tabulates "vCurrent (filtered)" as 2,957,858 sequences.
- That count is smaller than v18.0 (3,076,419), which is consistent with an older pre-v18 release such as the 2020-01 file on E:. **[Inf]**
- The build itself dates from about January 2025 (split_fa commits of 2025-01-17).
- **Checking whether the surviving file has 2,957,858 headers would settle this.** That cheap read-only check was not run in this pass (blocked investigation 2).

**Rebuild feasibility:**
- **Functionally equivalent:** yes, from any U-RVDB release + accession2taxid + taxdump.
- **Exact reproduction:** only if (i) the surviving vCurrent file is the snapshot that was used, (ii) the accession2taxid and taxdump snapshots match those used (unknown dates), and (iii) the Centrifuge version matches (unrecorded). Rebuilding from today's U-RVDB would give a **new, functionally similar database, not the historical one**.

---

## 9. Biological reference stack (status of each resource)

| Resource | Status | Rebuild class |
|---|---|---|
| Cat GCF_018350175.1, Human GCF_000001405.40 (full p14, 705 seqs incl. alt/patch), Monkey GCF_003339765.1, Pig GCF_000003025.6 | **Exist elsewhere, compressed only** (`/mnt/e/SequencingData/ADVTIG_fa/*.fa.gz`). gzip integrity OK; first headers match the expected assemblies. | Source available |
| Tamarin GCA_021498475.1 (*Saguinus midas*) | **Inside a zip only** (`Tamarin_GCA_021498475.1.zip`, tested OK). The internal path lacks the `GCA_…_ASM2149847v1/` top folder the protocol expects; `map_genome_to_reads` expects a third, flat path. | Needs extraction to a decided path |
| `Human_hg38_p14_Primary_Assembly.fasta` (map_genome_to_reads only) | Absent | Approximate: filter "Primary Assembly" headers |
| Viral targets: ADV5 (AY339865.1), EBV-B95-8 (V01555.2), FeLV (NC_001940.1), FeLV-KT (MT129531.1), PCV1 (NC_001792.2), RSV-A (JF920069.1), REOvirus 10 segments, MVM NC_001510.1, SMRV NC_001514.1, PERV AF038600.1, OC43 NC_006213.1 | **Exist elsewhere** (`/mnt/e/SequencingData/ADVTIG_fa/`). Some have CRLF line endings. | Source available |
| EBV/PCV1 doubled/split variants | Exist under `.2` filenames. `map_genome_to_reads` names them `.5`/`.6`, matching their header versions. | Rename or symlink |
| HCoV-NL63 NC_005831.2 | **Inside other files only** (U-RVDBvCurrent, C-RVDBv21) | Extract + reheader |
| Lambda / DNA CS (`ecoli-lambda-phage[.NC_001416.1].fasta`) | Candidate in Dropbox (3.56 kb). Identity to the original unverified. | Probably rebuildable |
| `viral_input-sequences.fna` (viral RefSeq panel source) | **Absent everywhere searched** | **Source missing** |
| `ss_bacteria_sequences.fna` (bacterial panel source) | **Absent** | **Source missing** |
| `fungal-seq-2022.fna`, SILVA 138.1, PA7, PLSDB | Absent or differently named | Not needed by the viral path |

**minimap2 indexes.** None present; all were absent in every search.

| Index | Intended source | Class |
|---|---|---|
| `ADV5.mmi` | `ADV5.fa` (inferred; no build command recorded) | Rebuildable from a clearly identified FASTA. The parameters are presumably defaults (unrecorded). |
| `all_mammalian.mmi` | `all_mammalian.fa` = 5 host genomes (cat recipe at protocol L154–160). The `minimap2 -d` command for it is not recorded. | Probably rebuildable; the FASTA recipe is explicit and the sources are present |
| `all_mammalian_targets.mmi` | + 11 viral targets incl. NL63 (L162–179); no `-d` command | Probably rebuildable (NL63 needs extraction). Only used by `old_combined`, which is unused. |
| **`all_mammalian_targets2.mmi`** | targets + AF038600.1 (again) + lambda NC_001416.1 (L181–183; `-d` recorded) | Probably rebuildable; **the lambda source version is uncertain** |
| `ecoli-lambda-phage.NC_001416.1.mmi` | lambda (L241) | Probably rebuildable (same caveat) |
| **`viral-input-sequences.mmi`** + `.tsv` | `viral_input-sequences.fna` (L256–271) | **Source data missing** |
| **`bacterial-input-sequences.mmi`** + `.tsv` | `ss_bacteria_sequences.fna` (L290–305) | **Source data missing** |
| `mammalian_ref_species2.tsv` | header-derived (L185–252) | Rebuildable, **but the protocol contains two conflicting recipes**: L231–236 appends `AF038600_1`, while L248–252 re-copies the base file and appends only `DNA_CS`. **[Obs]** Historical outputs contain columns for *both* `AF038600_1` and `DNA_CS`, so the historical file had both rows. |

Even with the same FASTAs, byte-identical `.mmi` files would also need the same minimap2 version and parameters, and neither is recorded. Having broadly similar genomes is **not** proof that the historical index can be reproduced exactly.

---

## 10. Sequencing inputs

| Dataset | Expected samples | Classification | Evidence |
|---|---|---|---|
| ADVTIG (legacy) and ADVTIG_v2/v3_untargeted (same 6 samples) | `20230119_DNA_advtig-Rep1B-E50_10CFU_41_36` (negative control); `20221222_DNA_advtig-Rep{2..6}B-E50_10CFU_41_36` (spiked with RSV, REO, EBV1, FeLV, PCV1, ADV5) | **Present only in derived form; raw data referenced on an unavailable filesystem** | NanoPlot 1.30.0 logs record the exact raw paths, e.g. `/mnt/usersData/ADVTIG/20230119_DNA_advtig-Rep1B-E50_10CFU_41_36/20230119_1550_X4_FAU43076_901473ea/sequencing_summary_FAU43076_901473ea_42144f8a.txt`. Reads per sample: 150k–1.99M. Surviving material: derived trees on `/mnt/e/SequencingData/ADVTIG` (7.2 GB: mpileups, NanoPlot, Recentrifuge from C-RVDBv25 runs, 16 small subset FASTQs) and Dropbox copies. No `trimmed/`, no `mini-u-viral2` Centrifuge reports, no BAMs. |
| Viral_FDA (`20240216_sDNA_WHO-colab2-…_10CFU_51_36` ×5) | from `viral_DNA_all16-2` | **Apparently lost or unknown** | No name hits on any drive. The Dropbox `Viral_FDA/` folder is a *different* run (20231019, batch 49, R10.4.1). The 16-2 sheet points all 5 samples to Identifier `S5B1` / `barcode144`, which is the code's "not barcoded" sentinel. |

**[Obs] Read format:**
- Single-end ONT 1D reads: R9.4.1, Guppy model `2021-05-17_dna_r9.4.1_minion_384_d37a2ab9`, GridION positions X1–X5, one flowcell per sample.
- Raw layout: `{DIR}/{sample}/{YYYYMMDD_HHMM_Xn_FLOWCELL_hash}/fastq_pass/`.
- Reads are untrimmed on input; porechop runs in stage 1.
- The sample sheets plus the NanoPlot logs are enough to reconstruct the expected ADVTIG input paths.

**[Inf] Remaining hope.** `/mnt/usersData` was most likely a directory on an old WSL root filesystem. `/mnt/f/ubuntu_backup.tar` (213 GB, a WSL root export) and `/mnt/d/WSL_Backups/ubuntu.vhdx` (500 GB) were not searched. They may contain the raw reads, the RVDB derivatives, the Centrifuge index, the panel FASTAs and `/home/james`.

---

## 11. Environment and software

| Item | Evidence | Gap |
|---|---|---|
| Python | `Pipfile` / `Pipfile.lock`, now copied into the repo root as **historical** metadata (`f1c0773`): python 3.10, pandas 1.4.3, numpy 1.23, biopython 1.85 | **Not a complete current environment.** pysam is missing from the lock, though `run_coverage` and `filt_low_complexity` import it, and no CLI tools are covered (B6, B7). SMART-CAMP's `requirements.txt` (2021, WinPython) is stale and was not copied. |
| Tools | `Seq commands.md` records downloads of minimap2 v2.17 and v2.24 and samtools 1.17. Centrifuge and porechop were installed by untagged git clone. NanoPlot 1.30.0 appears in logs. | **No record of which versions ran ADVTIG.** No conda environment, container or install document exists. |
| Hard-coded locations | ~40 paths across `/mnt/usersData`, `/home/james`, `~/` (unexpanded inside quotes) and `/mnt/e`; DB keys baked into filename globs; `github="/home/james/SMART-CAMP/"` | — |
| Resources | Stage 1: up to 5 Centrifuge processes × 25 threads, each holding the index in RAM. `run_coverage` uses `Pool(20)` and buffers whole SAM/mpileup output in memory. The split collection needs millions of inodes. Logged runtimes are 17–40 h per stage. | Not documented |

---

## 12. Historical outputs (evidence base)

**`/mnt/e/SequencingData/ADVTIG_primary`** (22 MB) holds:
- the Fourth-Pass versions (§7);
- `MetaFilt-Fifth-Pass*.csv/xlsx`, including `-15percentED` (02-13) and the current one (03-09);
- `_new_rows`, `VIRAL-ONLY`, `BACTERIAL-ONLY`;
- `Deduplicated-Read-Counts-*_{2025,v2,v3,v4}`;
- First/Second/Third-pass tables (some with an `ADVTWG-` prefix);
- reviewer workbooks (Charley3, Purple, edited): AutoFilters and hidden rows, but no added columns.

**Version suffixes.** `_v3`, `_v4` and `-15percentED` are manual renames made after rsync (`Seq commands.md:2991`).

**Other derived data.** `/mnt/e/SequencingData/ADVTIG` holds 2023 C-RVDBv25-era derived data.

---

## 13. Hard-coded paths and portability

The concerns in revision 2, §19, still stand. Additions:

- **Stage 1:**
  - the literal `~` in `-ci`;
  - `github`, `krakenuniq` and `hs-blastn` paths hard-coded to `/home/james` (unused on the `-skip` path);
  - non-recursive `**` globs fix the input layout at exactly one run-directory level;
  - intermediates are written into the raw tree;
  - `database_config.py` evaluates `/mnt/usersData` and `/home/james` paths at import time (prints only).
- **`map_genome_to_reads.py`:** about 60 hard-coded reference paths in `__main__`.

---

## 14. Documentation drift

These items from revision 2 still hold. The first two are now labelled HISTORICAL in the run-books (`4af00a8`):
- the `-a -c 2` commands (there is no `-c` option; it was removed in `aac5714b`);
- `deep_cov … -dd -e`;
- `Third-Pass-FULL-…TaxID_taxID-…` naming;
- `summarise_urvdb.py` misnamed;
- `split_fa.py` index name.

Revision-4 items, with their current status:
- **`-e 0.15` commands.** **[Resolved, `4af00a8`]** At revision 4, protocol `run_coverage.py … -e 0.15 …` could not run under any revision. The current protocol commands use `--edit-distance 0.15` [`-x`] and parse; each records its original `1c782e3e` text.
- **`collect_taxid_primary_counts.py` read `MetaFilt-Fourth-Pass.csv`.** **[Resolved, `4af00a8`]** It now reads `MetaFilt-Fifth-Pass.csv`, which stage 6 has written since 2026-02-02. Its output file name still says Fourth-Pass (open decision).
- **Relative paths.** `./pipeline/map_genome_to_reads.py` **now resolves** (`f1c0773`). `~/SMART-CAMP/mp_metagenomic_assessment_v4.py` **still points outside this repository** (B9, open).
- **Two recipes for `mammalian_ref_species2.tsv`.** They contradict each other (§9).
- **"(same as glue.py)"** refers to a methylfilter tool.

---

## 15. Unknown or ambiguous items

1. Which U-RVDB snapshot built `mini-u-viral2`. This is testable (§8).
2. The snapshot dates of the accession2taxid and taxdump files used to build the index.
3. The contents of `viral_input-sequences.fna` and `ss_bacteria_sequences.fna`.
4. The exact lambda / DNA CS sequence used.
5. The tool versions used (Centrifuge, porechop, minimap2, samtools).
6. Whether an uncommitted `run_coverage.py` was ever used. Not needed to explain the outputs; its search was blocked.
7. Viral_FDA 2024-02-16 raw data.
8. What the WSL backups contain.
9. *(Superseded.)* The author's intended option name for `run_coverage.py`'s extract mode is still unrecorded, but it no longer blocks anything: the current interface was set by a workflow decision (§24).

---

## 16. New User Readiness

> **Status.** Assessed at revision 4 (2026-10-01). Rows changed by `f1c0773` and `4af00a8` carry a *Current:* note. The verdict still stands.

The question: could a competent bioinformatician who has never seen this workflow determine each item below from the material available, without asking the author? General bioinformatics knowledge is assumed.

| Item | Determinable? | What is missing |
|---|---|---|
| Purpose of the workflow | **Partially.** Inferable from code and filenames (untargeted ONT adventitious-virus screen with host/panel counter-screen). | A README stating purpose, scope, spikes and controls |
| Stage order | **Partially.** The protocol lists commands, but mixed with historical runs, DB builds and unrelated experiments. | A clean run-book of current commands only |
| Required inputs | **No** without code reading | An input specification: raw ONT `fastq_pass` + `sequencing_summary` per sample, layout exactly `{DIR}/{sample}/{run}/fastq_pass` |
| Sample sheet format | **Partially.** `CONFIGURATION-FILE-INPUTS.txt` describes the columns. But two formats are needed (`.txt` for stage 1 and `map_genome`, `.csv` for stages 2–6), plus the `Spike_species` column and the hard-coded 6 samples in stage 7. *Current:* example sheets (`configs/`) and `CONFIGURATION-FILE-INPUTS.txt` are now in the repo (`f1c0773`). | A single, clean, current user-facing schema (B5, partial) |
| Where to put data | **No** | Directory-layout documentation; configurable roots |
| Which databases | **Partially.** The names are visible, but which of 4 Centrifuge DBs is current takes analysis. | A DB manifest: name, release, files, checksums |
| How to obtain or build them | **Partially.** Recipes are scattered in URVDB.md and the protocol, sometimes contradictory, and the `split_fa` configuration has to be recovered from git. | Build scripts; recorded RVDB release; panel FASTA sources |
| Software to install | **Partially** (from code) | An install guide |
| Software versions | **No** | Pinned environment incl. pysam and the CLI tools |
| Python environment | **Partially.** *Current:* the historical SMART-CAMP Pipfile pair is in the repo (`f1c0773`), but it is incomplete (no pysam). | A complete environment definition (B6) |
| Commands to run | **Partially.** *Current:* the current `run_coverage` commands now parse (`4af00a8`), but the protocol is not yet a clean current run-book. | A clean run-book of current commands (B8) |
| Historical versus reusable commands | **Partially.** *Current:* historical `run_coverage`/`deep_cov` commands are labelled HISTORICAL in the protocol and in `URVDB.md` (`4af00a8`), but other historical runs remain interleaved. | Full separation (B8) |
| Expected outputs per stage | **No** without code reading | An output table |
| Detecting stage success | **No.** Tool exit codes are unchecked; failures surface late with misleading errors. | Checks or documentation |
| Resuming after failure | **No.** It depends on undocumented cache files, `sys.exit`s, rerun-twice steps and SAM/BAM deletion. | Resume notes |
| Temporary versus final files | **No** | Output classification |
| CPU, memory, storage, runtime | **No.** Scattered `time` comments only. | A requirements statement |
| Comparing with historical results | **No.** The baselines are on E: with manual suffixes. | A baseline manifest naming the authoritative files, with checksums |

**Verdict.** A new user **still could not** run this workflow from the available material without the author.
- **No longer an obstacle:**
  - the `run_coverage` CLI problem (B3);
  - the stage 6→7 Fourth/Fifth-Pass mismatch (B4);
  - the hidden SMART-CAMP package dependencies, now packaged (B1/B2).
- **Still to rediscover, beyond the missing data:**
  - the `split_fa` configuration in git;
  - the input layout;
  - which of the protocol's remaining commands are current;
  - how to supply the environment, the references and the machine-specific paths.

---

## 17. Workflow Completeness

### Missing code
- **Current: none.** `pipeline/` (21 files), `q_control/` (2) and `pipeline/map_genome_to_reads.py` with its helpers are now in the repo (`f1c0773`). The root duplicate `map_genome_to_reads.py` remains (decision pending).
- **[Resolved, `4af00a8`]** Previously, no committed `run_coverage.py` combined a configurable ED threshold with an extract mode; the last working SMART-CAMP revision was `aff9daba` (ED fixed at 0.10). The current interface is `--edit-distance` (default 0.15) and `-x` (§24).
- **[Resolved, `4af00a8`]** The stage 6→7 name inconsistency: stage 7 now reads `MetaFilt-Fifth-Pass.csv`. Stage 7's output name still says Fourth-Pass (open decision).
- `glue.py` is not required, and the historical Fourth-Pass producer exists in git history.
- The `split_fa.py` configuration that built the URVDB split collection exists only in SMART-CAMP history (`7f6e9993`); see B10.

### Missing reference data
- **Missing sources:** `viral_input-sequences.fna`, `ss_bacteria_sequences.fna`.
- **Present but not yet verified as the historical source:** the exact U-RVDB snapshot (candidate present); the lambda DNA CS (candidate present).

### Missing generated indexes
These are all regenerable, in principle, from available sources:
- `mini-u-viral2.*.cf`, `U-RVDB-clean-fix2.fasta`, `RVDB2_seqID.map`, `fixed_URVDB_seqID2.map`, `URVDB_split_fa/`, `U-RVDB-clean-fix2-index.json`;
- `all_mammalian*.fa/.mmi`, `ADV5.mmi`, the lambda `.mmi`, `mammalian_ref_species2.tsv`.

These are **not** regenerable: `viral-input-sequences.{mmi,tsv}` and `bacterial-input-sequences.{mmi,tsv}`.

### Missing runtime inputs
- **ADVTIG raw reads:** on an unavailable filesystem; possibly in the WSL backups.
- **Viral_FDA 2024-02-16 raw reads:** no trace.
- **Sample sheets:** *Current:* present in `configs/` (`f1c0773`). The scripts still reference the old `/home/james/SMART-CAMP/configs/` paths (B9).

### Missing environment or software definitions
- A complete pinned environment: pysam plus the CLI tools (Centrifuge, porechop, NanoPlot, minimap2, samtools, ripgrep, GNU awk). The historical `Pipfile`/`Pipfile.lock` in the repo do not cover these.
- Tool versions.

### Missing documentation
- README, input specification, run-book, DB build guide, output guide, resource requirements, and a baseline definition (§16).

### Unknown provenance
- Which U-RVDB snapshot; accession2taxid and taxdump dates; the lambda sequence; the RefSeq panel contents; tool versions; the Tamarin and primary-assembly target paths. (The historical `run_coverage` extract-flag name is unrecorded but no longer relevant; see §24.)

### Potentially recoverable dependencies
- **Copy from SMART-CAMP:** packages, configs, Pipfile. **Done** (`f1c0773`).
- **Recover from git:** `split_fa.py` @ `7f6e9993`; `count_the_screen.py` Fourth-Pass naming @ `e4617a07`.
- **From E: and Dropbox:**
  - host genomes and viral targets (`ADVTIG_fa`);
  - `U-RVDBvCurrent.fasta`, `nucl_gb.accession2taxid`, taxdump;
  - NL63 (extract);
  - the lambda candidate.
- **Possibly in the unsearched WSL backups:** raw reads, the Centrifuge index, RVDB derivatives, the panel FASTAs, `/home/james`.

### Dependencies that appear genuinely lost
These are conditional on the WSL backups also lacking them:
- Viral_FDA 2024-02-16 raw reads;
- `viral_input-sequences.fna` and `ss_bacteria_sequences.fna`;
- the `mini-u-viral2` index as built;
- historical tool versions.

---

## 18. Recoverability by level

| Level | Status | Evidence |
|---|---|---|
| Source code recoverable | **Yes** (except an uncommitted CLI fix, which is not needed to explain the outputs) | §4, §6 |
| Workflow logic recoverable | **Yes** | §3, §7: every stage and transformation is traced, including Fourth Pass |
| Reference resources recoverable | **Partially** | Mammalian panels and the Centrifuge DB can be rebuilt from surviving sources (DB snapshot unverified). The viral and bacterial RefSeq panels cannot. |
| Historical execution environment recoverable | **No / Partially** | Python pins exist, minus pysam; CLI tool versions are unrecorded |
| Historical input data recoverable | **No** (pending the backup search) | §10 |
| Exact historical reproduction possible | **No** | Raw inputs, DB snapshot, tool versions, RefSeq panel contents; the code revision is inferable |
| Functionally equivalent new execution possible | **Yes, conditionally** | Needs: input reads (recovered or new), a rebuilt DB from a chosen RVDB release, substitute RefSeq panels, an environment, and configurable paths. (The `run_coverage` CLI has been corrected; `4af00a8`.) |

**Exact historical reproducibility: No.** Script revisions can be identified, and the parameters and thresholds are recoverable from git for each historical output. The Fourth-Pass step is fully reproducible. But the following cannot currently be recovered:
- the sample inputs;
- the database release (a candidate exists but is unverified);
- the panel reference sequences;
- the tool versions.

**Functional reproducibility: Achievable.** The scientific logic is fully recoverable and every stage can be rebuilt with documented contemporary resources. The required decisions are listed in §20.

---

## 19. Irreducible Unknowns

Searches performed:
- SMART-CAMP: working tree and all history;
- methylfilter: working tree and history;
- `/mnt/c/Users`, `/mnt/d`, `/mnt/e`, `/mnt/f` (depth-bounded) and `/home/mangifera`;
- historical outputs, NanoPlot logs, and `Seq commands.md`.

Because three searches were blocked or skipped (header block), **nothing below is certified irreducible**. Each item is irreducible *unless* the unsearched sources contain it.

| Unknown | Would be resolved by |
|---|---|
| Contents of the viral and bacterial "RefSeq" counter-screen panels (`viral_input-sequences.fna`, `ss_bacteria_sequences.fna`) | WSL backups (`/home/james/SequencingData/Centrifuge_libraries/…`) |
| Viral_FDA 2024-02-16 raw reads | WSL backups or the sequencing facility |
| ADVTIG raw reads (`/mnt/usersData/ADVTIG/…`) | WSL backups |
| Exact `mini-u-viral2` index and its source snapshot | Header count of the surviving vCurrent file (2,957,858?); WSL backups |
| Snapshot dates of the accession2taxid and taxdump files used | WSL backups (file mtimes); otherwise **irreducible** |
| Versions of Centrifuge, porechop, minimap2 and samtools actually used | WSL backups (conda envs or binaries under `/home/james`); otherwise **irreducible** |
| An exact uncommitted `run_coverage.py`, if one existed | Blocked forensic search; WSL backups. *Its intended behaviour is not an unknown* (§6). |

**Not** irreducible:
- the Fourth-Pass transformation: fully explained;
- the `split_fa` configuration: in git;
- the counter-screen thresholds per output: in git;
- the composition of `mammalian_ref_species2.tsv`: inferable from output columns.

---

## 20. Minimum recovery package (status at `4af00a8`)

**A. Exists; copy into the workflow**
- **[Done, `f1c0773`]** From SMART-CAMP: `pipeline/`, `q_control/`, the `configs/` sheets (ADVTIG 9-2 and Viral_FDA 16-2), `CONFIGURATION-FILE-INPUTS.txt`, and `Pipfile`/`Pipfile.lock` (historical, as a starting point). `pipeline/map_genome_to_reads.py` is in place; the root copy has not yet been removed.
- From E: `ADVTIG_fa` viral targets and host genomes; `U-RVDBvCurrent.fasta`; `nucl_gb.accession2taxid`; `names.dmp`/`nodes.dmp`.
- Historical baselines: `MetaFilt-Fourth-Pass_v4.csv`, `MetaFilt-Fifth-Pass-15percentED.csv`, `MetaFilt-Fifth-Pass.csv`, `Deduplicated-…_v4.csv`, the NanoPlot summaries.

**B. Regenerate deterministically (from A plus recorded recipes)**
- `U-RVDB-clean-fix2.fasta`, `RVDB2_seqID.map`, `fixed_URVDB_seqID2.map` (URVDB.md recipes).
- `URVDB_split_fa/` + index JSON (`split_fa.py` @ `7f6e9993`).
- `mini-u-viral2.*.cf` (recorded command). It is deterministic given the inputs and the Centrifuge version.
- `all_mammalian.fa`, `all_mammalian_targets.fa` (with NL63 extracted), `ADV5.mmi`, `all_mammalian*.mmi`.
- Decompressed or extracted host genomes.

**C. Reconstruct approximately**
- The lambda DNA CS reference and its `.mmi` (candidate file).
- `Human_hg38_p14_Primary_Assembly.fasta`.
- **[Done, `4af00a8`]** A `run_coverage.py` CLI with configurable ED: `-x` for extraction and `--edit-distance` (default 0.15) for the threshold (§24).

**D. Needs a human decision**
- Which U-RVDB release to use: the historical candidate (verify by header count) or a current release.
- The authoritative baseline: Fourth-Pass v4 (ED 0.10) versus Fifth Pass 15percentED (ED 0.15, Depth ≥ 0.35) versus Fifth Pass of 03-09 (ED 0.15, Depth ≥ 0.6).
- **[Decided, `4af00a8`]** `collect_taxid_primary_counts.py` reads the Fifth-Pass file. **Still open:** whether to rename its `…-Fourth-Pass-FILTERED.csv` output.
- Substitutes for the viral and bacterial RefSeq panels, and their sources and releases.
- **[Decided, §24]** The `run_coverage` flags (`-x`, `--edit-distance`, default 0.15). **Still open:** a canonical, clean protocol (B8).
- What to do with the root duplicate `map_genome_to_reads.py`.
- Whether Viral_FDA, the `map_genome_to_reads` lambda QC branch, and the P. aeruginosa / SILVA appendices are in scope.
- Tool versions to pin.
- Whether to search the WSL backups (213 GB tar / 500 GB vhdx) before declaring data lost.

**E. Apparently unavailable** (pending the WSL backup search)
- ADVTIG raw reads; Viral_FDA 2024-02-16 raw reads.
- The `viral_input-sequences.fna` and `ss_bacteria_sequences.fna` sources.
- The original `mini-u-viral2` index files.
- Historical tool versions.

---

## 21. Current Usability Verification

> **[Historical snapshot, 2026-10-01, commit `31c9ea1`.]** The checks, verdict text and blocker descriptions below record the state *before* `f1c0773` and `4af00a8`. The *Current (`4af00a8`)* columns give the present status. The overall verdict is unchanged: **NOT YET USABLE**. For the current summary see §0.

**Verified:** 2026-10-01, against the repository as committed at `31c9ea1` on branch `workflow-audit`, plus the files under `audit/` added in this pass. No workflow code was modified.

**How it was verified:** static checks only, using `audit/static_validation.py`. The full output is saved in `audit/static_validation_results.txt`. The script:
- byte-compiles every `.py` file into a temporary directory;
- resolves imports from the syntax tree, exactly as Python would when each script is run directly;
- rebuilds each script's argparse parser on its own (argparse only, with stand-in database keys for `mp_metagenomic_assessment_v4.py`);
- parses every command line documented in `ADVTIG-protocol.md` and `URVDB.md` against those parsers, after substituting the run-books' own `VAR="…"` assignments;
- syntax-checks the run-book code blocks with `bash -n`;
- greps for file names that one stage writes and another reads.

Nothing biological was executed: no Centrifuge, minimap2, samtools, database builds or decompression. No synthetic dry run was attempted, because every stage hard-codes absolute `/mnt/usersData` paths (see B9) and could not be pointed at a scratch directory without editing the code.

### Verdict: **NOT YET USABLE** (then and now)

*At 2026-10-01:* a new bioinformatician who clones this repository cannot run the intended workflow, even after supplying their own reads. Stage 1 and stage 1b fail on import, and stage 2 fails while building its argument parser. No environment is defined, no sample-sheet template exists, there is no current command sequence, and reference locations are hard-coded to a machine that no longer exists.

*Current (`4af00a8`):* the import failures and the `run_coverage.py` parser failure are fixed, and the sample sheets are present. Still missing: a complete environment, a clean user-facing command sequence and sample-sheet specification, configurable paths, and prepared references.

### Check results

| # | Check | Result at 2026-10-01 | Evidence | Current (`4af00a8`) |
|---|---|---|---|---|
| 1 | All Python files compile | **PASS** (10/10). There is a `SyntaxWarning` for an invalid escape `\{` at `counter_screen.py:812`; it is cosmetic. | results §1 | PASS (33/33) |
| 2 | `pipeline/` and `q_control/` modules present | **FAIL.** `mp_metagenomic_assessment_v4.py` imports 12 `pipeline.*` and 2 `q_control.*` modules, and none are in the repo (`pipeline/` and `q_control/` are ABSENT). They exist only in `/mnt/d/GitHub/SMART-CAMP`. | results §2, §7 | **PASS** (`f1c0773`): present; all internal imports resolve |
| 3 | `map_genome_to_reads.py` at the documented path, with its sibling imports resolving | **FAIL.** The protocol (L44–45) calls `./pipeline/map_genome_to_reads.py`, but the file is at the repo root. Its imports `bespoke_parse_map_reads` and `pre_trim_fq` are missing. | results §2, §4 | **PASS** (`f1c0773`): `pipeline/map_genome_to_reads.py` + helpers present and resolving; root duplicate remains |
| 4 | Sample sheets or example configs present | **FAIL.** There are no `.csv` or `.txt` sheets in the repo (§7 of the results: NONE). The documented sheets are `/home/james/SMART-CAMP/configs/viral_DNA_all9-2.{txt,csv}`, which live only in SMART-CAMP `configs/`. | results §7 | **PASS** for presence (`f1c0773`): `configs/` + `CONFIGURATION-FILE-INPUTS.txt`. No clean current spec yet (B5 partial) |
| 5 | A reproducible environment definition including pysam | **FAIL.** No `environment.yml`, `requirements.txt`, Pipfile, lock file or container exists in the repo. The required third-party Python packages are pandas, numpy, pysam and Biopython. The only candidate (SMART-CAMP `Pipfile.lock`) has no pysam. pandas is not installed in this machine's `python3` either. | results §2 | **FAIL**: historical Pipfile pair now in repo, but no pysam and no CLI tools |
| 6 | Command-line tools documented and installable at stated versions | **FAIL.** centrifuge, porechop, NanoPlot, minimap2, samtools, ripgrep and GNU awk are required, and none is listed with a version anywhere in the repo. Only `rg` is installed here. | `command -v` checks; §11 | FAIL |
| 7 | `run_coverage.py` builds its parser and the documented syntax is consistent | **FAIL.** `argparse.ArgumentError: argument -e/--edit-distance-threshold: conflicting option string: -e`. All 37 documented `run_coverage.py` invocations, active and commented, fail. The current protocol form `-e 0.15 -e` cannot be valid under any parser. | results §3–4 | **PASS** (`4af00a8`): parser constructs; current commands use `--edit-distance 0.15` [`-x`] and parse; no `-e` option |
| 8 | Other parsers and documented invocations | **PASS** for `mp_metagenomic_assessment_v4.py` (6 lines), `map_genome_to_reads.py` (2), `deep_cov.py` (active lines with `-e 0.15`), `filt_low_complexity.py`, `counter_screen.py` (12) and `count_the_screen.py` (4). Only the commented or historical `deep_cov.py … -dd -e` lines fail (protocol L126, URVDB.md L326). | results §4 | PASS; current commands 41 OK / 0 FAIL; historical lines classified separately |
| 9 | Fourth versus Fifth Pass naming | **FAIL.** `count_the_screen.py:564` writes `MetaFilt-Fifth-Pass.csv`, while `collect_taxid_primary_counts.py:91,179` reads `/mnt/usersData/ADVTIG_v3_untargeted/MetaFilt-Fourth-Pass.csv`. Run in the documented order, stage 7 fails with FileNotFoundError. | results §6 | **PASS** (`4af00a8`): both stages use `MetaFilt-Fifth-Pass.csv` |
| 10 | Each stage has one unambiguous current command | **FAIL.** The protocol interleaves three datasets, three Centrifuge databases, commented historical forms, DB builds and an unrelated P. aeruginosa appendix. Several stages appear more than once with different flags. `collect_taxid_primary_counts.py` must be run three times. `URVDB.md`'s code block also fails `bash -n` at L335 (a Markdown table inside the block), which is a sign that it is notes, not a script. | results §4–5 | FAIL (historical `run_coverage`/`deep_cov` lines now labelled; `URVDB.md` `bash -n` failure unchanged) |
| 11 | Required directory and input layout documented | **FAIL** in user-facing docs. The layout `{DIR}/{sample}/{run}/fastq_pass/*.fastq[.gz]` plus `{DIR}/{sample}/{run}/*sequencing*summary*.txt` (exactly one run-directory level) is documented only in this audit (§3, §10). There is no README or protocol text. | README is 16 bytes | FAIL |
| 12 | Database and reference build documented enough to reproduce | **PARTIAL.** Recipes exist in `URVDB.md` (L45–107) and the protocol (L154–357). But the `split_fa.py` settings that actually built `URVDB_split_fa/` exist only in git history (SMART-CAMP `7f6e9993`), not in the repo. The U-RVDB release is not stated. The protocol has two contradictory `mammalian_ref_species2.tsv` recipes. There is no build command for `ADV5.mmi` or `all_mammalian.mmi`. | §8–9 | PARTIAL |
| 13 | Reference assets supplied, generated or explicitly user-supplied | **FAIL.** The repo supplies **no** reference assets, and none are declared as user-supplied in any user document. The two RefSeq counter-screen panel sources (`viral_input-sequences.fna`, `ss_bacteria_sequences.fna`) have neither a build or download procedure nor a declaration as user-supplied. | §9 | FAIL |
| 14 | Paths configurable for a new user | **FAIL.** There are 43 hard-coded absolute paths outside comments: `run_coverage.py` 7, `deep_cov.py` 4, `counter_screen.py` 12, `count_the_screen.py` 2, `collect_taxid_primary_counts.py` 4 (plus 6 hard-coded sample names), `split_fa.py` 10, `mp_metagenomic_assessment_v4.py` 4. Database and panel locations are not command-line options. | `grep` count | FAIL (unchanged) |

### Blockers

| ID | Blocker (as found 2026-10-01) | What prevented execution | Class | Status at `4af00a8` |
|---|---|---|---|---|
| B1 | `pipeline/` and `q_control/` absent | Stage 1 raises `ImportError: No module named pipeline` | **Repository packaging** | **RESOLVED** (`f1c0773`) |
| B2 | `map_genome_to_reads.py` misplaced; `bespoke_parse_map_reads.py` and `pre_trim_fq.py` absent | Stage 1b: the documented path doesn't exist, and the script raises ImportError | **Repository packaging** | **RESOLVED** (`f1c0773`; root duplicate retained) |
| B3 | Duplicate `-e` in `run_coverage.py` (introduced in `1c782e3e`) | Stage 2 exits before parsing any arguments | **Code defect** | **RESOLVED** (`4af00a8`; `-x`, `--edit-distance`, default 0.15, no `-e`) |
| B4 | Stage 6 writes Fifth-Pass; stage 7 reads Fourth-Pass | Stage 7 FileNotFoundError (plus a project decision on which file is intended) | **Code defect** | **RESOLVED** for the stage 6→7 input (`4af00a8`); stage-7 output name still says Fourth-Pass |
| B5 | No sample-sheet template or schema in the repo | The user cannot build the two required formats: `.txt` for stage 1 and 1b, `.csv` with `date, NA, strain, concentration_CFU, batch, duration_h, Spike_species` for stages 2–6 | **Documentation** (+ packaging) | **PARTIALLY RESOLVED**: sheets + existing column description present; clean current spec pending |
| B6 | No environment definition (pysam, Biopython, pandas, numpy) | Python imports fail; nothing defines a reproducible install | **Environment** | OPEN (historical Pipfiles only) |
| B7 | CLI tools unspecified (centrifuge, porechop, NanoPlot, minimap2, samtools ≥1.10, rg, GNU awk) | Not installed; no versions chosen | **Environment** | OPEN |
| B8 | No current run-book, input-layout spec or output description | The user cannot tell which commands to run, in what order, or where to put reads | **Documentation** | OPEN |
| B9 | Hard-coded absolute paths (`/mnt/usersData`, `/home/james`) for the dataset root in stage 7, DB files, panels and the split collection | Stages cannot find user-supplied references without source edits | **Code defect** | OPEN |
| B10 | `mini-u-viral2` Centrifuge index, `URVDB_split_fa/`, index JSON, `fixed_URVDB_seqID2.map`, `RVDB2_seqID.map` | Not supplied. The recipe is partly in git only, and the RVDB release is unstated. These could be generated deterministically once a release is chosen and the split_fa configuration is restored. | **Missing reference resource** (needs a build script and a release decision) | OPEN |
| B11 | Mammalian panels: `all_mammalian_targets2.mmi`, `mammalian_ref_species2.tsv`, `ADV5.mmi`, `all_mammalian.mmi`, lambda DNA CS | Not supplied; the build recipe is partial or contradictory | **Missing reference resource** (deterministic build possible) | OPEN |
| B12 | `viral-input-sequences.*`, `bacterial-input-sequences.*` | Source FASTAs unknown or absent; no procedure | **Missing reference resource** (needs a project decision: substitute source or declare user-supplied) | OPEN |
| B13 | Literal `~` passed to Centrifuge (`-ci` in double quotes); unchecked tool exit codes; stage 1 needs a sequencing summary or it crashes at the end | Silent failure, then a misleading crash | **Code defect** (stage 1) / **Documentation** | OPEN |
| — | Raw ONT reads | Expected to be supplied by the user. **Not a blocker**, provided B8 documents the format and layout. | **Expected user-supplied input** | n/a |

### Minimum changes to reach USABLE WITH EXTERNAL INPUT DATA

Items marked **[Done]** were completed in `f1c0773` or `4af00a8`. The rest is the remaining remediation work.

1. **Packaging (B1, B2, B5):**
   - **[Done, `f1c0773`]** Copy SMART-CAMP `pipeline/` and `q_control/` into the repo root.
   - **[Done, `f1c0773`]** Place `map_genome_to_reads.py` in `pipeline/`, with its helpers. *Open:* remove or retire the root duplicate.
   - **[Done, `f1c0773`]** Add `configs/` with the ADVTIG and Viral_FDA `.txt`/`.csv` sample sheets. *Open:* a clean current sample-sheet specification (B5).
2. **Code (B3, B4, B9, B13):**
   - **[Done, `4af00a8`]** `run_coverage.py`: `-x` for extraction and `--edit-distance` (default 0.15) for the threshold. This supersedes the earlier suggestion of `-e <float>`; there is no `-e` option.
   - **[Done, `4af00a8`]** `collect_taxid_primary_counts.py` reads `MetaFilt-Fifth-Pass.csv`, the file stage 6 writes. *Open:* the name of its Fourth-Pass-labelled output.
   - Replace hard-coded `/mnt/usersData` and `/home/james` locations with command-line options or a single config file:
     - dataset root;
     - DB root (split collection, index JSON, seqID maps, accession→taxid map);
     - panel directory;
     - sample list for stage 7.
   - Expand `~`, or document an absolute `-ci` path.
3. **Environment (B6, B7):** add an `environment.yml` pinning python, pandas, numpy, pysam, biopython, centrifuge, minimap2, samtools (≥1.10), porechop, nanoplot and ripgrep. Choose versions explicitly, since the historical ones are unrecorded.
4. **References (B10–B12):**
   - Add a reference build script, or a `Makefile`, that turns user-downloaded inputs into every derived asset:
     - U-RVDB release → clean FASTA → seqID maps → split collection + JSON → `centrifuge-build`;
     - NCBI host assemblies → `all_mammalian*.fa/.mmi` → label TSVs, using a single resolved `mammalian_ref_species2.tsv` recipe that includes `AF038600_1` and `DNA_CS`.
   - Record the chosen accessions and releases in a manifest.
   - Decide and document the viral and bacterial RefSeq panel sources, or declare them explicitly user-supplied with a stated format.
5. **Documentation (B8):** write a README covering:
   - purpose;
   - the stage order with one current command per stage;
   - the input layout (`{DIR}/{sample}/{run}/fastq_pass/` + sequencing summary; ONT single-end);
   - the sample-sheet schema;
   - the reference prerequisites;
   - the outputs of each stage, which are final and which are temporary;
   - resume behaviour (caches, repeat runs);
   - resource expectations.

   Move historical commands and the P. aeruginosa appendix out of the run-book.

When items 1–5 are done, a new user supplying their own reads and running the reference build would meet the **USABLE WITH EXTERNAL INPUT DATA** criteria. Exact reproduction of the historical results stays out of reach for the reasons in §18–19.

---

## 22. Repository packaging pass (2026-10-02)

**Scope.** Copy only. Workflow files that the audit had located elsewhere were copied into their original SMART-CAMP locations in this repository. No code, path, sample sheet, environment file or reference was modified, and nothing was built. The copied files were later committed unchanged as `f1c0773`. The SMART-CAMP originals were left untouched; its working tree was still clean at HEAD `65a94b3a` after the copy. The WSL backups were not searched.

### 22.1 How the copy set was chosen

The set comes from the import graph, not from directory listings.
- An AST walk of the transitive local imports, run in SMART-CAMP from `mp_metagenomic_assessment_v4.py` and `pipeline/map_genome_to_reads.py`, gave 21 `pipeline/` files (including `__init__.py`) and 2 `q_control/` files. This matches the revision-4 audit: 12 + 2 direct imports, 5 second-level imports, and the 2 `map_genome_to_reads` siblings.
- **[Obs]** The pipeline modules' other file references were checked:
  - every hard-coded `/home/james/SMART-CAMP/configs/...` or seqid-CSV reference in the copied pipeline modules sits inside an `if __name__ == "__main__":` block (AST-verified), so it never runs when the workflow imports the module;
  - `{github}/all_seqids*.csv` (`dep_flat_sample_classifier_test.py` L108–110) is reached only on the BLAST branch, and the documented stage-1 commands pass `-skip`, which limits classification to Centrifuge.
- Sample sheets: every sheet referenced in `ADVTIG-protocol.md` or `URVDB.md` was copied, plus the TXT companions of the referenced 9-2-5/9-2-6 CSVs.
  - Viral_FDA `viral_DNA_all16-2.*` is included because the protocol runs it through the same scripts (L44, L102–140).
  - `viral_DNA_all2.txt` and `viral_DNA_all10.txt` appear only in commented lines for other projects (Viral_CHO, Viral_human), so they were not copied.
- The stage 2–7 scripts in `ADVTIG-untargeted/` were re-grepped and reference no other SMART-CAMP source.

### 22.2 What was copied (provenance)

All 36 files were copied with `cp -n --preserve=timestamps` from SMART-CAMP (`/mnt/d/GitHub/SMART-CAMP`, branch `master`, HEAD `65a94b3a`). Each one is git-tracked there, and its working-tree content equals the HEAD blob. "Source commit" below is the last commit that changed the file. Every copy was checked with `cmp` and is byte-identical. The same data is in `WORKFLOW_INVENTORY.tsv` (columns `Source Path`, `Source Repository`, `Source Commit`, `Byte Identical To Source`).

| Target path | Source path (under `/mnt/d/GitHub/SMART-CAMP/`) | Source commit | Byte-identical | sha256 (first 16) |
|---|---|---|---|---|
| `pipeline/__init__.py` | `pipeline/__init__.py` | (last changed d237e129) | yes |  `e3b0c44298fc1c14` |
| `pipeline/Convert_U2T_fastq.py` | `pipeline/Convert_U2T_fastq.py` | (last changed 57f47b92) | yes |  `54bde4fb8bdac730` |
| `pipeline/bespoke_parse_map_reads.py` | `pipeline/bespoke_parse_map_reads.py` | (last changed d1d06b9d) | yes |  `a0610da1a95a19f8` |
| `pipeline/cat_nanostats.py` | `pipeline/cat_nanostats.py` | (last changed 4fbec5e6) | yes |  `921f0b73d17a032a` |
| `pipeline/database_config.py` | `pipeline/database_config.py` | (last changed cdab6e94) | yes |  `d1667d68d89377bf` |
| `pipeline/dep_flat_sample_classifier_test.py` | `pipeline/dep_flat_sample_classifier_test.py` | (last changed 8f9dc1b3) | yes |  `3bcbd7b07d7ab099` |
| `pipeline/dep_pre_chunker.py` | `pipeline/dep_pre_chunker.py` | (last changed 83037469) | yes |  `6824ab295e5c0230` |
| `pipeline/dep_run_nanostat_analyses.py` | `pipeline/dep_run_nanostat_analyses.py` | (last changed b6cf93fc) | yes |  `fc5ed81b8ed930e9` |
| `pipeline/find_barcode_reads_for_QC.py` | `pipeline/find_barcode_reads_for_QC.py` | (last changed 3758dddd) | yes |  `293ba7531d4d18ec` |
| `pipeline/guppy_demux.py` | `pipeline/guppy_demux.py` | (last changed fa1d8406) | yes |  `e045081033e1ff2f` |
| `pipeline/map_genome_to_reads.py` | `pipeline/map_genome_to_reads.py` | (last changed b16899bb) | yes |  `f445a1a640b8eeb6` |
| `pipeline/mp_cent_interpreter_v3.py` | `pipeline/mp_cent_interpreter_v3.py` | (last changed 50fe201a) | yes |  `78efc48cf41bd7f8` |
| `pipeline/mp_comb_ce_BN_v3.py` | `pipeline/mp_comb_ce_BN_v3.py` | (last changed dd81b06f) | yes |  `3be89f5e53b3e868` |
| `pipeline/mp_demux_reads.py` | `pipeline/mp_demux_reads.py` | (last changed 619d90e1) | yes |  `393a50509d697fa7` |
| `pipeline/mp_host_remove.py` | `pipeline/mp_host_remove.py` | (last changed 0154a078) | yes |  `59c2266076119998` |
| `pipeline/mp_trim_reads_v2.py` | `pipeline/mp_trim_reads_v2.py` | (last changed 9ea2537c) | yes |  `191e275a61715e23` |
| `pipeline/pre_trim_fq.py` | `pipeline/pre_trim_fq.py` | (last changed 62418d33) | yes |  `5ebd84a64c57d648` |
| `pipeline/rank_BLAST_predictions_for_aa_v2.py` | `pipeline/rank_BLAST_predictions_for_aa_v2.py` | (last changed 62418d33) | yes |  `a05dc662432755a5` |
| `pipeline/rank_centrifuge_predictions_for_aa_v2.py` | `pipeline/rank_centrifuge_predictions_for_aa_v2.py` | (last changed 62418d33) | yes |  `4e22a336492086c0` |
| `pipeline/sample_name_sanity_check.py` | `pipeline/sample_name_sanity_check.py` | (last changed 62418d33) | yes |  `98671ad90f568bca` |
| `pipeline/updated_blastn_interpreter_v2.py` | `pipeline/updated_blastn_interpreter_v2.py` | (last changed 2400b8ec) | yes |  `72451a10f4a013ed` |
| `q_control/adjust_q_scores.py` | `q_control/adjust_q_scores.py` | (last changed 073a98fe) | yes |  `288562d88bcb19ce` |
| `q_control/filter_q.py` | `q_control/filter_q.py` | (last changed b6207bb1) | yes |  `1bfeb41b05e6fdfb` |
| `configs/viral_DNA_all9-2.txt` | `configs/viral_DNA_all9-2.txt` | (last changed 972e8db6) | yes |  `f8d792868e3dc449` |
| `configs/viral_DNA_all9-2.csv` | `configs/viral_DNA_all9-2.csv` | (last changed 0de43abe) | yes |  `41ad32a2dd3aed35` |
| `configs/viral_DNA_all9-2-1.csv` | `configs/viral_DNA_all9-2-1.csv` | (last changed a1012075) | yes |  `24676449b5f3a5ae` |
| `configs/viral_DNA_all9-2-4.csv` | `configs/viral_DNA_all9-2-4.csv` | (last changed b53f0bba) | yes |  `6507268bf8b92509` |
| `configs/viral_DNA_all9-2-5.csv` | `configs/viral_DNA_all9-2-5.csv` | (last changed 0e7689d5) | yes |  `36a97f3fbf7de1b4` |
| `configs/viral_DNA_all9-2-5.txt` | `configs/viral_DNA_all9-2-5.txt` | (last changed a2753ec5) | yes |  `f8d792868e3dc449` |
| `configs/viral_DNA_all9-2-6.csv` | `configs/viral_DNA_all9-2-6.csv` | (last changed 0e7689d5) | yes |  `84b4ca20fce20098` |
| `configs/viral_DNA_all9-2-6.txt` | `configs/viral_DNA_all9-2-6.txt` | (last changed a2753ec5) | yes |  `f8d792868e3dc449` |
| `configs/viral_DNA_all16-2.txt` | `configs/viral_DNA_all16-2.txt` | (last changed 8288c417) | yes |  `7d9391889b99ce8c` |
| `configs/viral_DNA_all16-2.csv` | `configs/viral_DNA_all16-2.csv` | (last changed d8649b03) | yes |  `ffd81db1e6bbd4be` |
| `Pipfile` | `Pipfile` | (last changed 2a305ed9) | yes |  `4ce2b104f1a5f513` |
| `Pipfile.lock` | `Pipfile.lock` | (last changed 2a305ed9) | yes |  `b8543dbf7fce4bd4` |
| `CONFIGURATION-FILE-INPUTS.txt` | `CONFIGURATION-FILE-INPUTS.txt` | (last changed e5360b0c) | yes |  `7501ed1bc6123a7b` |

The files that were already in the repo before this pass were re-checked the same way and are byte-identical to SMART-CAMP HEAD: `ADVTIG-protocol.md`, `mp_metagenomic_assessment_v4.py`, root `map_genome_to_reads.py` (= `pipeline/map_genome_to_reads.py`) and all 9 files in `ADVTIG-untargeted/`. Their provenance is also in the TSV.

### 22.3 Deliberately not copied

| Item | Reason |
|---|---|
| The other ~140 SMART-CAMP `pipeline/` modules, including `map_genome_to_reads2.py`, `mp_trim_reads.py`, `mp_cent_interpreter*.py` and `ML_tool/` | Not in the import closure |
| `q_control/basequal.py`, `bin_q_scores_length.py`, `fastq_parser.py` | Not in the import closure |
| `all_seqids.csv`, `all_seqids_incl_plasmid.csv`, `C_RVDB_seqids.csv` (8–13 MB) | Reference lookup tables. Only used on the BLAST branch (not reached with `-skip`) or in `__main__` blocks. |
| `glue.py` (methylfilter) | Not a workflow dependency (revision 4) |
| `split_fa.py` @ `7f6e9993`, `run_coverage.py` @ `aff9daba`, `count_the_screen.py` @ `e4617a07` | Historical revisions of files already in the repo. Restoring any of them would change behaviour; deferred. |
| `Seq commands.md`, SMART-CAMP `README.md` | Historical notes and parent-project documentation, not workflow support files |
| All FASTAs, indexes, databases, reads and historical outputs | Out of scope for this pass |

### 22.4 Verification

1. **The existing static checks, rerun unchanged** (`audit/static_validation.py` → `audit/static_validation_results_packaging.txt`):
   - §1: all 33 `.py` files compile.
   - §2: all 14 `pipeline.*`/`q_control.*` imports of `mp_metagenomic_assessment_v4.py` resolve, and `pipeline/map_genome_to_reads.py` resolves both of its helpers.
   - §2 still prints MISS lines, for two reasons:
     - the root duplicate `map_genome_to_reads.py`, whose MISSes are genuine;
     - five `from pipeline import …` lines inside `pipeline/*.py`. These are an artefact of the checker, which resolves each file against its own directory. At runtime they resolve against the entry script's directory (see check 2).
   - §4: the protocol's `./pipeline/map_genome_to_reads.py` invocations (L44–45) changed from `PATH MISSING` to `path OK`.
   - The 42 FAIL lines are otherwise identical to revision 4. They are the known `run_coverage.py` parser failure (B3), which is out of scope here and not counted against packaging.
   - §2 now reports `Pipfile`/`Pipfile.lock` present.
2. **The new runtime-faithful import check** (`audit/packaging_import_check.py` → `audit/packaging_import_check_results.txt`) uses `importlib.util.find_spec` with `sys.path` set as it is at runtime, and executes no workflow code. **OVERALL PASS:**
   - `mp_metagenomic_assessment_v4.py`: 21/21 local imports resolve, transitively, including every second-level `from pipeline import …`;
   - `pipeline/map_genome_to_reads.py`: 2/2;
   - `ADVTIG-untargeted/deep_cov.py`: 1/1;
   - the root `map_genome_to_reads.py` fails 0/2 (informational: misplaced duplicate).
   - The remaining non-stdlib imports are third-party (pandas, numpy, pysam, Bio). They are an environment matter (B6), not packaging.
3. **Byte identity:** 36/36 identical (`cmp`), and each source matches its git HEAD blob.

### 22.5 Blocker status after this pass

| ID | Status | Basis |
|---|---|---|
| B1 | **RESOLVED** | The packages are physically present, and every internal import resolves |
| B2 | **RESOLVED** | The canonical path exists, and its helpers are present and resolve. The root duplicate remains (§22.6). |
| B5 | **PARTIALLY RESOLVED** | Sheets and the schema file are present. The documented single schema is still open, and references still point at `/home/james/SMART-CAMP/configs` (B9). |
| B6 | OPEN | The Pipfile pair is present but unchanged and incomplete (no pysam, no CLI tools) |
| B3, B4, B7–B13 | OPEN at the time | Not packaging. B3 and B4 were later resolved (§23–§24) |

Verdict for §21 is unchanged: **NOT YET USABLE**. (For the current blocker status, see §0.)

### 22.6 Duplicated or misplaced files, for a later decision

- **D1** `./map_genome_to_reads.py` (repo root) is byte-identical to `pipeline/map_genome_to_reads.py`. Its location is wrong: the protocol calls `./pipeline/…`, and from the root its helper imports cannot resolve. It was retained as instructed.
- **D2** `configs/viral_DNA_all9-2-5.txt` and `configs/viral_DNA_all9-2-6.txt` are byte-identical to `configs/viral_DNA_all9-2.txt` (all six samples), despite their single-sample names. They were duplicated like this in the source and are preserved unchanged.
- **D3 (path, not file)** Several references still point outside this repository; none were changed (B9, deferred):
  - the run-book invokes `~/SMART-CAMP/mp_metagenomic_assessment_v4.py`;
  - sample sheets are referenced as `/home/james/SMART-CAMP/configs/…`;
  - `github = "/home/james/SMART-CAMP/"` (`mp_metagenomic_assessment_v4.py` L490).

### 22.7 Source or configuration dependencies still missing

- **Packaging:** none. Every source and configuration file that the audit found elsewhere, and that belongs in the repo, is now present.
- **Still outstanding**, all outside this pass:
  - a `run_coverage.py` revision that accepts the documented CLI (never committed; B3). *Later resolved* by the interface correction in `4af00a8` (§24).
  - the `split_fa.py` configuration for URVDB, which exists only as git revision `7f6e9993`;
  - a complete environment definition (B6/B7);
  - the source FASTAs for `viral-input-sequences` and `bacterial-input-sequences` (B12);
  - all reference data, indexes and raw reads (B10, B11, inputs).

---

## 23. B3 and B4 code corrections (2026-10-02)

**Scope.** Only B3 and B4 were corrected. No paths, configuration handling, sample sheets, references, environment, thresholds or scientific logic were changed, and the root `map_genome_to_reads.py` duplicate is still in place. The changes were committed in `4af00a8`, together with §24.

> **Interim B3 form, superseded by §24.** This pass first made `-e FLOAT` (default 0.1) the threshold. §24 then replaced it with `--edit-distance FLOAT` (default 0.15) and removed `-e` entirely. The `-e 0.15` interface rows below record that interim state only. The B4 content is current.

### 23.1 Source files modified

| File | Change | Lines |
|---|---|---|
| `ADVTIG-untargeted/run_coverage.py` | `parser.add_argument('-e', '--extract-fastq', …)` → `parser.add_argument('-x', '--extract-fastq', …)` | L1329 (1 line) |
| `ADVTIG-untargeted/collect_taxid_primary_counts.py` | `FILTER_CSV = ".../MetaFilt-Fourth-Pass.csv"` → `".../MetaFilt-Fifth-Pass.csv"` (directory unchanged) | L91, L179 (2 lines) |
| `ADVTIG-protocol.md` | *(Interim; now `--edit-distance 0.15 -x`, §24.)* Current commands L59 and L106: `-e 0.15 -e` → `-e 0.15 -x`. Each line has a trailing comment recording the original `-e 0.15 -e` text from `1c782e3e`. Historical commented lines L39, L119 and L126 got a trailing `[HISTORICAL …]` tag; their commands are unchanged. | 5 lines edited, none inserted, so line numbers cited elsewhere in this inventory stay valid |
| `ADVTIG-untargeted/URVDB.md` | A 3-line `# NOTE` at the end of the code block (L474–476) marks every `run_coverage.py`/`deep_cov.py` command in the log as historical, with the old flag meanings. No command was changed. | appended only |
| `audit/static_validation.py` | Validator extended; see §23.4 | — |

As a result, these four workflow files are no longer byte-identical to SMART-CAMP HEAD. `WORKFLOW_INVENTORY.tsv` records why.

### 23.2 Exact behavioural change

- **B3 (interim, superseded by §24):** the change is limited to the command-line flag name for FASTQ extraction, which is now `-x` (long form `--extract-fastq` and destination `extract_fastq` unchanged). `-e`/`--edit-distance-threshold` is now unambiguous: a float, with the default **unchanged at 0.1**.
  - The parser previously could not be built at all, so the script could not start. With the same flags it now does exactly what the `1c782e3e` code intended.
  - Unchanged: extraction (`finish_setup`), the alignment and edit-distance calculations (`extract_mapping_data`, `needs_ed_rerun`), and all thresholds.
- **B4:** stage 7 now reads the merged table that current stage 6 writes, `{out_dir}/MetaFilt-Fifth-Pass.csv` (`count_the_screen.py:564`).
  - The columns stage 7 reads (`Sample`, `TaxID`, `SeqID`) are present in the Fifth-Pass schema: the historical `/mnt/e/SequencingData/ADVTIG_primary/MetaFilt-Fifth-Pass.csv` has them at the same positions as `MetaFilt-Fourth-Pass_v4.csv` (cols 1, 2, 35).
  - Unchanged: filtering, joining and counting.
  - The **data content** of the stage-7 input does change, because current Fifth Pass holds more (Sample, TaxID) keys than historical Fourth Pass v4 (§7). This was always the consequence of running the current code in order, and it is why historical Fourth-Pass baselines are not directly comparable.

### 23.3 Provenance for the choices

- **B3: why `-x`.**
  - Commit `1c782e3e` (2026-02-12) added `-e/--edit-distance-threshold` (float) to `run_coverage.py` while leaving `-e/--extract-fastq`.
  - In the same commit, `deep_cov.py` resolved the identical clash by moving its extract flag to `-x` and keeping `-e` as the float (§6).
  - The protocol written in that commit uses `-e 0.15` on every stage, so `-e` = threshold was the intent. `-x` follows the existing `deep_cov.py` convention.
  - The 0.1 default was kept because it is the default in the `1c782e3e` code. Historical runs at 0.10 used the hard-coded value from before `1c782e3e`.
- **B4: why Fifth Pass.**
  - `a68d1880` (2026-02-02) renamed the stage-6 output from `MetaFilt-Fourth-Pass.csv` to `MetaFilt-Fifth-Pass.csv`, but `collect_taxid_primary_counts.py` was never updated to match (§7).
  - The user directed standardising on Fifth Pass for the current workflow.
  - No other executable code references either name: the only other `MetaFilt-*-Pass` strings are stage 6's `-VIRAL-ONLY`/`-BACTERIAL-ONLY` subset outputs, which are already Fifth-Pass.
  - Historical Fourth-Pass result files were not touched.

### 23.4 Validation evidence

`audit/static_validation_results_b3b4.txt` holds the output of `python3 audit/static_validation.py`.

**Changes to the validator** (classification only; no existing check removed or relaxed):
- §4 now tags each documented command:
  - **CURRENT** = an uncommented line in `ADVTIG-protocol.md`;
  - **HISTORICAL** = a commented line in `ADVTIG-protocol.md`, or any line in the `URVDB.md` log.
- Historical lines are still parsed against the current parsers. Their failures print as `HFAIL` and are counted separately.
- §6 gains an explicit stage 6 → 7 filename comparison.
- A new §8 holds targeted B3/B4 checks.

**Results:**

| Check | Result |
|---|---|
| All Python files compile (§1) | **PASS** 33/33 |
| `run_coverage.py` parser constructs (§3) | **PASS** (was `argparse.ArgumentError … conflicting option string: -e`) |
| Help text generates (§8) | **PASS** via `format_help()`. The real file's `run_coverage.py --help` also exits 0 through its own `__main__`. pandas, pysam, numpy and Biopython are not installed (B6), so inert stand-ins were injected for those imports only; `--help` exits before any analysis code. |
| `-e 0.15` → threshold 0.15, no extraction | **PASS** |
| `-x` → extraction, threshold default 0.1 | **PASS** (`--extract-fastq` also works) |
| `-e 0.15 -x` and `-x -e 0.15` coexist | **PASS** |
| `-e 0.15 -r`, `-e 0.15 -l`, `-e 0.15 -a -x` | **PASS** |
| Bare `-e` (old extract form) | rejected, **as expected** |
| Current documented commands (§4) | **41 OK, 0 FAIL**, including all 8 `run_coverage.py` lines (L59–62, L106–109) |
| Historical documented commands (§4) | 26 still parse; **18 HFAIL**, all explained (see list below) |
| Stage 6 output = stage 7 input (§6, §8) | **PASS**: both `MetaFilt-Fifth-Pass.csv` |
| Remaining `FAIL` lines in the whole report | **1**: the pre-existing `URVDB.md` `bash -n` failure (a Markdown table inside the code block, L335). It is a documentation issue, unrelated to B3/B4. |

**The 18 HFAIL lines:**
- 9 use the old bare `-e` = extract: protocol L39, L119 and L126 (`deep_cov`), and URVDB.md L251, L260, L311, L320, L326 (`deep_cov`) and L361.
- 6 use `-c 2` (model number), which was removed from `run_coverage.py` in `aac5714b` on 2025-03-27: protocol L40–42 and URVDB.md L312–314. L39 and L311 also use `-c 2` but are counted under the bare `-e` group above.
- 2 omit the `--database` argument, which later became required: URVDB.md L362–363.
- 1 is a commented `deep_cov.py` with no arguments: protocol L91.

None of these were made to parse. Supporting them would contradict the current interface (now `--edit-distance FLOAT` and `-x`, §24).

### 23.5 Issues found, not acted on

- **Stage-7 output name.** `collect_taxid_primary_counts.py` L182 still writes `Deduplicated-Read-Counts-MetaFilt-Fourth-Pass-FILTERED.csv`. After B4 its contents are filtered by Fifth Pass, but the name says Fourth Pass, and a rerun in `/mnt/usersData/ADVTIG_v3_untargeted/` would overwrite the historically named file. No later stage reads it, so it is not part of the stage-sequence link and was left unchanged. **Decision needed.**
- **Historical `-c` flag.** Protocol L39–42 (commented) and URVDB.md L311–314 use `-c 2`, which no longer exists. The revision-4 audit had attributed all `run_coverage` failures to the `-e` clash.
- **Cached outputs may skip the new threshold.** `run_coverage.py` reuses existing editdist TSVs unless `needs_ed_rerun` triggers. Runs into directories that already hold 0.10-era outputs depend on that logic, which was not changed or exercised here.

### 23.6 Blocker status after this pass

**B3 RESOLVED** (final interface in §24). **B4 RESOLVED.** B1 and B2 stay resolved (§22), and B5 is partially resolved. **B6–B13 are OPEN** and were not touched. The §21 verdict is still **NOT YET USABLE**.

---

## 24. Current workflow decision: `run_coverage.py` command line (2026-10-02, commit `4af00a8`)

**Decision.** In the reconstructed workflow, the **current** `run_coverage.py` interface is:

| Function | Current option | Default |
|---|---|---|
| FASTQ extraction | `-x` / `--extract-fastq` (flag) | off |
| Edit-distance threshold | `--edit-distance FLOAT` (long option only; internal name `edit_distance_threshold`) | **0.15** |
| `-e` | **does not exist** | — |

Omitting `--edit-distance` behaves exactly like `--edit-distance 0.15`. Any other float can still be supplied explicitly: the value flows unchanged into `main()` → `extract_mapping_data()` → `needs_ed_rerun()` / `get_editdist_metrics()`, as it did at `1c782e3e`. This supersedes the interim `-e FLOAT` (default 0.1) form recorded in §23. Only the CLI changed: the edit-distance calculation was neither investigated nor modified.

**How this differs from historical behaviour.**

| Revisions | `-e` meant | Edit-distance threshold |
|---|---|---|
| up to `aff9daba` (2026-02-12 10:36) | extract FASTQ (`store_true`) | hard-coded 0.10, not a CLI option |
| `1c782e3e` (2026-02-12 13:08), the SMART-CAMP HEAD | **both**: extract FASTQ *and* `-e/--edit-distance-threshold` (default 0.1), so argparse failed and the script could not start | intended 0.15 per the protocol's `-e 0.15` lines (never runnable) |
| **current (this repo)** | **nothing; no `-e`** | `--edit-distance`, default **0.15** |

Historical commands keep their original `-e` and `-c` syntax and are labelled HISTORICAL in the run-books (§23.1). The validator parses them against the current CLI and reports them as HFAIL.

**Files changed in this step:**
- `ADVTIG-untargeted/run_coverage.py` L1333: `parser.add_argument("-e", "--edit-distance-threshold", type=float, default=0.1, …)` → `parser.add_argument("--edit-distance", dest="edit_distance_threshold", type=float, default=0.15, …)`. L1329 keeps `-x/--extract-fastq` from §23.
- `ADVTIG-protocol.md`: all 8 current `run_coverage.py` lines (L59–62, L106–109) now use `--edit-distance 0.15`; the extraction lines L59 and L106 read `--edit-distance 0.15 -x`. Each line has a trailing note quoting its original `1c782e3e` text, for example `as committed in 1c782e3e this line read '-e 0.15 -r'`. No lines were inserted.
- `audit/static_validation.py` §8: updated to the new interface, and adds a scan of the current `run_coverage` commands.

**Validation.** `audit/static_validation_results_cli.txt` holds the output of `python3 audit/static_validation.py`.

| Check | Result |
|---|---|
| All Python files compile | PASS 33/33 |
| Parser constructs; `format_help()` lists `-x, --extract-fastq` and `--edit-distance`, with no `-e` | PASS |
| Real `run_coverage.py --help` (third-party imports stubbed, since B6 means they are not installed) | PASS, exit 0 |
| `--edit-distance 0.15` → 0.15 | PASS |
| `-x` → extraction on, threshold 0.15 | PASS |
| `--edit-distance 0.15 -x` and `-x --edit-distance 0.15` | PASS |
| `--edit-distance 0.15` with `-r`, with `-l`, and with `-a -x` | PASS |
| Default threshold is 0.15, and omitting the option equals `--edit-distance 0.15` | PASS |
| An explicit other value (`--edit-distance 0.10`) is accepted | PASS |
| `-e 0.15` and bare `-e` are rejected | PASS (as intended) |
| No current `run_coverage.py` command uses `-e`; all 8 pass `--edit-distance 0.15`; extraction lines end in `--edit-distance 0.15 -x` | PASS |
| §4 documented commands | CURRENT 41 OK / 0 FAIL; HISTORICAL 26 parse / 18 HFAIL (same classification as §23.4) |
| Remaining FAIL lines | 1, the pre-existing `URVDB.md` `bash -n` failure (unrelated) |

**Scope note.** The current `deep_cov.py`, `counter_screen.py` and `count_the_screen.py` commands still pass `-e 0.15`. That `-e` belongs to *their own* CLIs, where it is a valid float threshold (§4: OK). They are outside this `run_coverage.py`-only decision and were not changed.

