# Validation Status

This public source repository passed structural validation and a complete
synthetic-demo run before publication. The executable code sources were previously exercised in
an isolated short-path Windows R 4.5.1 mirror. Candidate-specific documentation,
licenses, and release metadata do not alter that executable analysis path.

## Completed checks for this repository

- Run the repository-candidate validator with no errors.
- Verify that no real data, generated outputs, logs, installed package library,
  private path, credential, or prohibited file type is included.
- Restore R 4.5.1 dependencies from `renv.lock` in an isolated short-path copy.
- Run `Rscript code/run_demo.R --run-id validation_20260908_04`.
- Confirm `qa/demo_qa_summary.csv` reports `PASS` and all 25 required outputs
  are present.
- Compare the SHA-256 hashes of every executable, configuration, test, and
  runtime documentation input used for the isolated run with this candidate.

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

Successful validation establishes that the synthetic demonstration can run from
a clean, short-path Windows copy after its documented environment restore. It
does not validate access to restricted study data, reproduce study estimates, or
authorize public release.

## Windows path note

The complete output tree includes long descriptive figure filenames. On Windows,
place a checkout in a short path before running the demo. The runner
intentionally stops before writing files when the planned output path would
exceed the conservative Windows path threshold.

## Future frozen-release gates

Before a future tagged GitHub Release or archive, confirm a version, release
date, archive DOI decision, and any preferred paper citation. No version tag,
GitHub Release, archive DOI, or paper DOI has been assigned.
