# Analysis Specification for the Synthetic Demonstration

The synthetic demonstration follows the public code sequence below:

1. Generate new synthetic baseline, stay-day, and PaCO2 inputs.
2. Derive person-period data for a 28-day follow-up framework.
3. Create descriptive baseline outputs.
4. Fit illustrative continuous association models using Poisson additive mixed modelling structure.
5. Run illustrative categorical landmark and cumulative burden analyses.
6. Run illustrative subgroup and pH-adjusted sensitivity analyses.
7. Create synthetic display-item files and validate their schemas and presence.

The same executable order is used by `code/run_demo.R`. The analysis is configured for one thread and a fixed seed to make the synthetic demonstration deterministic on a given compatible R environment. Platform-specific binary package availability can still affect environment installation; the output contract validates structure rather than byte-identical graphics.

This is not a protocol for source-data analysis. The authorized internal workflow, source definitions, statistical decisions, and manuscript authority remain inside the governed Study and are not exported by this candidate.
