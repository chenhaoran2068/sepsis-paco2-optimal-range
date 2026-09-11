# Sepsis PaCO2 Optimal Range

**Status:** public source repository for a synthetic demonstration. No GitHub
Release, version tag, archive DOI, or paper DOI has been issued.

## Study identity

This repository accompanies the study, *Daily Time-Weighted PaCO2 and
28-Day Mortality in ICU Sepsis: A Multicohort Retrospective Study*. The study
examined the observational association between daily time-weighted arterial
carbon dioxide partial pressure (TWA-PaCO2) during ICU days 1 to 7 and 28-day
all-cause mortality in four ICU sepsis cohorts.

The synthetic demonstration intentionally uses three analytic labels: `MIMIC`,
`AmsterdamUMCdb`, and `Chinese cohort`. The synthetic `Chinese cohort` is a
single aggregate label used only to exercise a shared multi-source schema for
the two Chinese study components. It does not represent actual participant
allocation, sample size, source structure, or the study's real four-cohort
pooling.

"Optimal range" is a concise repository identity for the study's
model-estimated lower-risk range. It is not a recommended bedside treatment
target. The full study identity, research question, cohort descriptions, and
observational interpretation boundary are in
[protocol/study-summary.md](protocol/study-summary.md).

## Proposed GitHub topics

GitHub stores topics as a flat list. They are grouped here for readability.

### Clinical domain

- `sepsis`
- `critical-care`
- `paco2`

### Study design

- `observational-study`
- `retrospective-cohort`
- `multicohort-study`

### Analysis

- `longitudinal-analysis`
- `survival-analysis`

### Geographic coverage

- `north-america`
- `europe`
- `asia`

## What this candidate contains

- Parameterized R scripts that exercise the selected analytic workflow on newly
  generated synthetic data.
- A deterministic synthetic data generator and data contract.
- Configuration, a locked R environment, output contract, and automated
  validation.
- A map from public code stages to the corresponding manuscript display items
  and claim families.

The categorical analysis uses stacked landmark Cox models for overall
estimates. Eligible patient-landmark rows are stratified by ICU day and, for
pooled analyses, analytic cohort, with robust standard errors clustered by
patient. Separate day-specific landmark models are retained as supporting
analyses.

## What it does not contain

- MIMIC-IV, AmsterdamUMCdb, Chinese cohort, or any other participant-level
  data.
- Source-data extraction adapters, credentials, agreements, ethics material,
  internal logs, model objects, or internal paths.
- Manuscript tables, figures, supplements, PDF files, or their numerical
  results.

## Requirements

- R 4.5.1. The `renv.lock` file records package versions.
- Network access for the first environment restore.

## Run the synthetic demonstration

From the repository root, first restore the recorded R environment:

```powershell
Rscript code/install_dependencies.R
```

Then run the complete synthetic demonstration:

```powershell
Rscript code/run_demo.R --run-id demo_run_001
```

The run creates new synthetic inputs below `demo/input/` and writes only
synthetic outputs below `output/demo_runs/demo_run_001/`. It refuses to
overwrite an existing run identifier. The shortest verification command is the
second command after the environment has been restored.

## Expected validation evidence

The runner executes `tests/validate_demo_run.R`. A successful run produces:

- `output/demo_runs/demo_run_001/qa/demo_qa_summary.csv` with `status` equal to
  `PASS`.
- `output/demo_runs/demo_run_001/qa/demo_output_inventory.csv`.
- The outputs listed in [expected/output-contract.csv](expected/output-contract.csv).

All generated inputs and outputs are synthetic demonstrations. They verify code
paths and output schemas only. They are not a reproduction of the study's
participant data, numeric estimates, figures, or conclusions.

## Data access

No study data are distributed. Read [DATA_ACCESS.md](DATA_ACCESS.md) before
attempting any work with source datasets.

## License and citation status

The code and original repository documentation are licensed under the MIT
License in [LICENSE](LICENSE). Newly generated synthetic demo data are dedicated
under CC0-1.0 only as described in [LICENSE-DATA.md](LICENSE-DATA.md). The
distributed `renv/activate.R` autoloader is third-party MIT-licensed code; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Study text, figures, tables,
supplements, restricted source data, and other third-party materials are not
distributed and are not relicensed.

The canonical repository is
`https://github.com/chenhaoran2068/sepsis-paco2-optimal-range`. No version tag,
GitHub Release, archive DOI, or paper DOI exists yet. [CITATION.cff](CITATION.cff)
describes this unversioned source repository; cite a frozen tagged archive
release when one becomes available.

## Support boundary

This candidate supports only the synthetic demonstration. It does not provide an
authorized pipeline for restricted source data, nor does it replace
institutional, database-provider, ethics, privacy, or data-use requirements.
