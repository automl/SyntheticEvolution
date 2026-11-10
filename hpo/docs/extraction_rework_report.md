# Report: DSSR step -> extraction step, RNA + protein

Date: 2026-11-10. Branch: `shs-analyzer`. Nothing is committed.

This work was first written against 72c8ba9 and then merged into HEAD (7ec3fd3),
which had meanwhile gained the persistent AF3 worker, per-RNA archives, scratch
directories, failure reports, retention, signals, environment checks and the
export rework. This report describes the merged state.

## Goal
Timeline note (28-10-26): for proteins, get C-alpha positions with Biopython, build a distance matrix, binary-mask it, and compare the resulting pair list with the generator's input. Rework `hpo/pipeline/dssr.py` into a type-aware **extraction step**, and make the dataset and config ready for several chain types (complexes later).

## Decisions
| Topic | Decision |
| --- | --- |
| Extractors | RNA -> DSSR, protein -> Biopython contacts. Both return `list[(i, j)]`, so `scoring.py` is unchanged. |
| Contact definition | CA atom, `distance <= 8.0 A`, `min_separation = 6`, all configurable. |
| Biopython env | Controller env, imported lazily (RNA-only runs never need it). |
| CSV | One row per chain, rows grouped by `id` (mirrors AF3 `sequences`). |
| Config | `generators:` block per type + `optimize:` key; only one generator is searched. |
| Trial parameters | Saved per chain type: `{"rna": {...}, "protein": {...}}` (only types present in the dataset). |
| Protein AF3 JSON | Fixed (it still emitted an `rna` chain). |

## What changed
**Extraction** ([extraction.py](../pipeline/extraction.py), renamed from `dssr.py` with `git mv`)
- Kept `parse_pairs` and HEAD's scratch-aware `run_dssr(..., scratch=...)` (log and DSSR's temporary directory in scratch, only `dssr.json` published). Removed the stray `DSSR` line that made `dssr.py` raise a `NameError` on import in 72c8ba9 and HEAD.
- New `contact_pairs()`: Biopython `MMCIFParser` -> CA coordinates -> numpy distance matrix -> boolean mask `np.triu(d <= cutoff, k=min_separation)` -> pair list. Strict validation (residue count, residue identity vs. sequence, chain, missing atom). Uses its own 3-to-1 table because Biopython's key case differs between versions.
- New `extract_pairs(chain_type, model, directory, sequence, settings, scratch=None)` with an `EXTRACTORS` registry. The published file is `dssr.json` (RNA) or `contacts.json` (protein, options and pairs); `RESULT_FILES` names both.

**Dataset and config** ([run.py](../pipeline/run.py))
- `create_dataset`: columns `id,chain_id,type,sequence,pairs,...`; `type` defaults to `rna` and `chain_id` to `A`, so old CSVs still load. Rows get `chains=[{id,type,sequence}]` plus `type`/`sequence` aliases. Per-type alphabet validation. More than one chain per `id` is rejected with "complexes ... not supported yet".
- `load_config`: new `generators`, `optimize`, `extraction` sections. Top-level `shs_python`, `fixed`, `search`, `dssr`, `dssr_timeout_seconds` are gone. A non-optimized generator's `search` is ignored with a warning. `validate_config_for_dataset` requires `python` / `extraction.dssr` only for chain types actually present. HEAD's `storage`, `failures`, `af3_buckets`, optional `evaluations` and `f1_soft` are unchanged.
- `fingerprint` hashes both generator dirs and the `extraction` settings; executable paths and timeouts are excluded (the DSSR path used to be included by mistake).
- `shs_generation_request(row, index, parameters, ...)` now takes the per-type parameters, selects the row's type, and passes only the fields the generators accept (`id, sequence, pairs, interactions, mutation_rates`). The trial and the exporter both use it, so the exported request still matches what the generator received.

**Trial** ([trial.py](../pipeline/trial.py)): HEAD's stage flow (scratch, archives, failure reports) with per-row generator and python by chain type. Sampled parameters go only to the optimized type; other types use their `fixed` values. Stage 3 calls `extraction.extract_pairs` with `scratch=` and publishes the archive again in a `finally` block, so `dssr.json` / `contacts.json` is kept even when the pairs are rejected. The progress stage is now called `extraction` (was `dssr`) in failure reports.

**Consumers of the new shapes**
- [environment.py](../pipeline/environment.py): checks every generator that has a `python` configured, DSSR if `extraction.dssr` is set, and `biopython` + `numpy` in the controller environment when a protein generator is configured.
- [export_results.py](../export_results.py): unchanged; it passes the saved per-type parameters to `shs_generation_request`.
- [trial_comparison.py](../../analysis/modules/trial_comparison.py): flattens per-type parameters to `rna.N` etc.
- [run_neps_hpo.py](../run_neps_hpo.py) reads `generators[optimize].search`; [run_standalone_trial.py](../run_standalone_trial.py) warning updated.

**Protein generator** ([json_generator.py](../../SHS-Generator-Protein/json_generator.py)): `build_input_json` now emits `{'protein': {id, sequence, unpairedMsa, pairedMsa: '', templates: []}}`.

**Configs and examples**: `helix.example.yaml` has both generator blocks, `optimize: rna`, the extraction section and HEAD's storage/failure/bucket settings. The protein search space only uses protein-valid parameters. `full-searchspace-example.yaml` is labelled RNA-only. Example CSVs migrated; added `protein_input.csv`. `biopython` added to `env/environment-controller.yml`. Docs updated (`pipeline_flow.mmd`, `storage_layout.md`, `hpo/README.md`, `analysis/INPUT_CONTRACT.md`).

## Tests
Run in a scratch uv venv (pytest, pyyaml, numpy, biopython, scipy, grakel). Result: **331 passed, 1 failed**. The failure, `test_submission.py::test_controller_script_starts_requested_module` (runs `controller.slurm`), fails in the same way on HEAD without these changes.

HEAD itself could not be tested as committed: the stray `DSSR` line in `dssr.py` broke the import and four test modules failed at collection.

Where the tests went: HEAD had reorganised the old files. The stash's changes were ported to their new homes: DSSR + protein/dispatch tests in `test_extraction.py`; config rules in `test_config.py`; chain-type CSV tests in `test_dataset.py`; fingerprint and request tests in `test_run.py`; the mixed RNA/protein trial in `test_trial.py`; per-type export in `test_export.py`; per-type comparison in `test_trial_comparison.py`; environment checks in `test_environment.py`; the protein AF3 JSON in `test_protein_json.py`.

## Not done / to know
- **Not run against real data**: no real AF3 CIF, DSSR binary, Slurm or NePS here. Protein extraction is tested on synthetic mmCIF only. Please run one protein trial on the cluster.
- Old workspaces fail the fingerprint check (new config shape): use a new `run_name`. Old `trial_result.json` files have flat parameters; `trial_comparison` still reads them, the exporter does not rebuild requests for them.
- `env/conda-lock.yml` is not regenerated; `biopython` is unpinned in `environment-controller.yml`.
- Task dirs and archives are still named `rna_<index>` for every chain type; the legacy CLI paths in the protein generator (`--rna-seq`, `attach_rna_msas`) still say RNA. The hpo pipeline only uses the request path.
- `analysis/modules/shs_best_config.py` is RNA-only: it loads `SHS-Generator-RNA/pair_map.py`, so exported protein tasks are not plotted by it.
- Generator parameter names are not validated at config load (that needs the generator's numpy env); a wrong key (e.g. `wobble_prob` for protein) fails when the generator runs.
- `hpo/docs/resulting_files.txt` no longer exists at HEAD, so there is nothing to update there.
