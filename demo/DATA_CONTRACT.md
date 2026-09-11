# Synthetic Demonstration Data Contract

## Boundary

The demonstration generator creates new records using a fixed random seed and explicit illustrative distributions. It reads no source data, study result, model object, internal configuration, or locally stored cohort file. No distributional parameter in this generator is estimated from restricted study data.

## Generated inputs

`code/R/generate_demo_data.R` creates three Parquet files below `demo/input/pooled/`:

- `all_baseline_outcome.parquet`: 1,500 synthetic ICU stays.
- `all_day_long.parquet`: 10,500 synthetic stay-day observations for ICU days 1 to 7.
- `all_raw_paco2.parquet`: 10,500 synthetic PaCO2 observations used to exercise the TWA calculation path.

The Study has four cohorts. The synthetic demonstration intentionally uses three analytic labels: `MIMIC`, `AmsterdamUMCdb`, and `Chinese cohort`. The synthetic `Chinese cohort` label is an illustrative aggregate used only to exercise the shared multi-source schema for the two Chinese study components. It does not encode actual site membership, sample allocation, dates, variable availability, source structure, or any source-data pooling decision.

The synthetic cohort labels support a common multi-source schema. They do not represent source records, source identifiers, source dates, or source cohort counts.

Each synthetic stay-day row includes a 1,440-minute nominal ICU-day window and
an illustrative interval from the window end to the synthetic event or
censoring time. These fields exercise the completed-window and forward-looking
landmark code paths. They are generated independently of study records.

## Illustrative risk mechanism

The generator uses an explicitly illustrative U-shaped PaCO2 risk mechanism to exercise continuous, categorical, burden, subgroup, and pH-adjusted code paths. The demonstration must not be used to reproduce, estimate, validate, or infer the study's numerical results or clinical conclusions.

## Regeneration

The runner sets a fixed seed, refuses to overwrite an existing run identifier, and writes new synthetic data and outputs for each new identifier. The generated files are ignored by version control and are not included in this candidate.
