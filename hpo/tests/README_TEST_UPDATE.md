# Additional HPO tests

Copy the eight Python files from `hpo/tests/` into your repository's `hpo/tests/`.
Keep all existing tests and fixtures. No production files are included or replaced.
No conftest.py is required: fixtures are local to the test modules that use them.

From the SyntheticEvolution repository root, in your development environment
with pytest and the pipeline's Python dependencies installed:

```bash
python -m pytest -q hpo/tests
```

Run only these additions:

```bash
python -m pytest -q hpo/tests/test_controller.py hpo/tests/test_af3_execution.py hpo/tests/test_submission.py hpo/tests/test_neps_handoff.py hpo/tests/test_dssr_execution.py hpo/tests/test_evaluation_aggregation.py hpo/tests/test_input_validation.py hpo/tests/test_standalone_execution.py
```

## What is covered

| File | Coverage |
| --- | --- |
| test_controller.py | Fresh two-RNA trial; shared parameter and seed forwarding; request and AF3 manifest paths; GPU submission arguments; differing per-RNA losses and mean; cached completion; generation errors/timeouts; wait-timeout recovery without resubmission; uncertain submission and command timeout; DSSR partial failure; missing model/marker; invalid predicted pairs; changed experiment identity. |
| test_af3_execution.py | Select correct array task; AF3 interpreter, environment and CLI arguments; distinguish selected model from sample models; missing/empty/multiple models; inference failure; success-marker publication; scheduler terminal failures; incomplete accounting; wait timeout. |
| test_submission.py | Both controller submitters: resources, working directory, log paths, Python environment and shell-safe quoting of spaces, quotes and dollar signs. |
| test_neps_handoff.py | Search-space construction without changing config; two callback invocations with distinct parameters/directories; dataset, fingerprint and loss forwarding; failure propagation; version guard; duplicate-controller lock. |
| test_dssr_execution.py | DSSR arguments, isolated working directory, logging and timeout; parsed-output handoff; nonzero exit, timeout, malformed JSON and missing output. |
| test_evaluation_aggregation.py | Equal weight per RNA despite different pair counts; reject NaN and infinite losses; empty target with spurious predicted pairs. |
| test_input_validation.py | Invalid CSV IDs, sequences, pair JSON and pair indices; empty dataset; duplicate IDs across files; whitespace/case and duplicate/reversed pair normalization. |
| test_standalone_execution.py | Validation-only stops before workspace preparation/execution; normal standalone hands off the prepared dataset and prints the returned loss. |

All external execution in these additions is substituted. They do not invoke
Slurm, ws_find, the SHS generator, AlphaFold, DSSR or an actual NePS search.
Filesystem operations use pytest temporary directories. The controller-lock test
uses a real local advisory lock; like the controller, it requires Linux/WSL's
fcntl support. Other existing tests can still run the real local SHS generator.

The fresh-trial test deliberately uses two distinct RNA sequences and losses
0 and 1. This checks that both dataset rows are evaluated and contribute equally
to the reported 0.5 loss, rather than merely checking that execution finishes.
The AF3 model and generated JSON placeholders contain no scientific prediction;
the existing generator and DSSR parser tests remain responsible for their formats.

## Validation performed

56 new test cases passed in 0.25 seconds in a staged local package assembled
from the uploaded trial, AF3 worker, standalone and submission modules.

The uploaded run_neps_hpo.py still had the old positional call, so only the
validation copy was adjusted to the interface you said you fixed:

```python
run_pipeline(config, dataset, pipeline_directory / 'artifacts', fingerprint, parameters)
```

The new handoff test checks that contract and accepts an equivalent keyword call.

The current hpo/pipeline/dssr.py and scoring.py were not attached. Validation used
the earlier supplied implementations, adapting DSSR's import to the package
location. This is not verification of your full current repository. The new DSSR
execution tests substitute parse_pairs; your existing parser fixture tests still
check the real parsing behavior. No full existing-suite result is claimed.

The new test source was also parsed using Python 3.9's syntax rules. Execution
under your actual Python 3.9 environment still needs the command above.

These tests verify Python orchestration and failure handling. They do not verify
the AF3 Slurm shell script, real module loading, cluster filesystem permissions,
real NePS 0.16.0 integration, or scientific prediction quality. After the local
suite passes, use a small real cluster smoke run before spending the HPO budget.
