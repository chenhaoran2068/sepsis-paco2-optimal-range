# Validation Status

This public source update incorporates the accepted stacked-landmark analysis
and forward-only handling of time-varying covariates. Its executable path
passed a complete synthetic-demo run in a short-path Windows validation copy.

## Candidate checks

- Repository-candidate schema and structural validation: passed for the rebuilt
  candidate described by the release manifest.
- Prohibited-content, private-path, credential, and symlink scan: passed by the
  automated candidate validator and targeted review.
- Clean R 4.5.1 synthetic-demo run from a short Windows path: passed using
  `Rscript code/run_demo.R --run-id validation_20260911_08d` after restoring the
  recorded environment.
- Required-output contract and schema validation: passed. The run generated
  1,500 synthetic stays, 10,500 synthetic stay-day rows, and all 25 required
  outputs, with `qa/demo_qa_summary.csv` reporting `PASS`.
- Executable-input SHA-256 comparison between the reviewed candidate and the
  isolated validation copy: passed for the final rebuilt candidate inputs.

The first run stopped before analysis because dependencies had not yet been
restored, which is the documented prerequisite. A later run exposed a missing
synthetic `window_duration_min` field required by the updated completed-window
landmark code. The generator, data contract, and automated tests were updated,
and the complete run then passed. These failed runs are retained only as local
validation evidence and are not distributed.

## Output contract boundary

`expected/output-contract.csv` lists the 25 required files that form the
minimum stable public validation surface. The automated QA verifies that each
listed path exists and is non-empty. These files are the artifacts a reader may
use to confirm that the documented synthetic pipeline completed.

A successful demo run also creates additional synthetic support artifacts,
including XLSX and PDF alternatives, quality-control summaries, raw support
tables, simulated narratives, intermediate model objects, and run logs. Those
files are useful for local inspection but are not a stable public interface:
their presence, file names, or columns are not guaranteed across future
candidate builds. They are generated only at runtime, are ignored by version
control, and are not included in this candidate.

## Interpretation

Successful validation establishes that the synthetic demonstration can run
from a clean, short-path Windows copy after its documented environment restore.
It will not validate access to restricted study data, reproduce study
estimates, or authorize public release.

## Windows path note

The complete output tree includes long descriptive figure filenames. On Windows,
place a checkout in a short path before running the demo. The runner
intentionally stops before writing files when the planned output path would
exceed the conservative Windows path threshold.

## Future frozen-release gates

Before a future tagged GitHub Release or archive, confirm a version, release
date, archive DOI decision, and any preferred paper citation. No version tag,
GitHub Release, archive DOI, or paper DOI has been assigned.
