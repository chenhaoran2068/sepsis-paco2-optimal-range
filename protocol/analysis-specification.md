# Analysis Specification for the Synthetic Demonstration

The synthetic demonstration follows the public code sequence below:

1. Generate new synthetic baseline, stay-day, and PaCO2 inputs.
2. Derive person-period data for a 28-day follow-up framework.
3. Create descriptive baseline outputs.
4. Fit illustrative continuous association models using a piecewise exponential additive mixed modelling structure. Time-varying covariates are carried forward from previously observed values only.
5. Fit categorical stacked landmark Cox models for overall estimates, stratified by ICU day and, for pooled analyses, analytic cohort, with patient-clustered robust standard errors. Retain separate day-specific landmark models as supporting analyses.
6. Quantify cumulative high-risk exposure burden across completed ICU-day windows using the same stacked landmark structure.
7. Run illustrative subgroup analyses and a pH-adjusted sensitivity analysis. The pH comparison uses the same eligible person-period rows before and after adding daily pH.
8. Create synthetic display-item files and validate their schemas and presence.

The same executable order is used by `code/run_demo.R`. The analysis is configured for one thread and a fixed seed to make the synthetic demonstration deterministic on a given compatible R environment. Platform-specific binary package availability can still affect environment installation; the output contract validates structure rather than byte-identical graphics.

This is not a protocol for source-data analysis. The authorized internal workflow, source definitions, statistical decisions, and manuscript authority remain inside the governed Study and are not exported by this candidate.
