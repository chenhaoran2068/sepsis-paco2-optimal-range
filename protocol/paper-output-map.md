# Code-to-Paper Output Map

This map documents conceptual correspondence only. Outputs created by this
package are synthetic demonstrations and do not reproduce manuscript estimates,
tables, figures, or claims.

| Public code stage | Required synthetic output families | Related manuscript display item or claim family |
|---|---|---|
| `code/R/00_functions/compute_twa_paco2.R` | No standalone display item | Daily TWA-PaCO2 exposure definition |
| `code/R/analysis/01_prepare_pamm_ped_data.R` | No standalone display item | 28-day follow-up modelling input and quality-control structure |
| `code/R/analysis/10_result1_population_baseline.R` | Table 1; Table S1 | Descriptive cohort characterization and missingness support |
| `code/R/analysis/11_result2_primary_pamm.R` | Figure 1; Figures S2-S4; Table S2 | Continuous association, TWA-PaCO2 distribution, adjustment hierarchy, day-specific curves, and model support |
| `code/R/analysis/13_result3_landmark_burden.R` | Tables 2-3; Tables S3-S6; Figure 2; Figure S5 | Landmark categories, high-risk burden, category sensitivity, descriptive event rates, secondary exposure definitions, and cohort/day support |
| `code/R/analysis/14_result4_subgroups_ph_sensitivity.R` | Figure 3; Figure S6; Tables S7-S8 | Subgroup analyses and pH-adjusted sensitivity |
| `code/R/analysis/20_build_latex_tables.R` | Synthetic LaTeX fragments derived from Tables 1-3 and Tables S1-S8 | Demonstration output formatting only; not manuscript-ready display files |
| `tests/validate_demo_run.R` | Presence and schema checks for the 25 required contract outputs | Public-package completion check only; not a numerical manuscript-results check |

## Output-contract boundary

`expected/output-contract.csv` enumerates the 25 files that form the minimum
stable public validation surface. The runner creates additional synthetic
support artifacts, including XLSX/PDF alternatives, quality-control files, raw
support tables, simulated narratives, intermediate objects, and logs. They aid
local inspection but are not part of the stable public contract and must not be
interpreted as study results.

## Selection-flow boundary

The source-cohort selection-flow material and final display graphic are not
distributed. The public Result 1 stage writes synthetic flow-data support only;
it intentionally does not create or promise a final Figure S1 graphic.

The manuscript PDF, final display files, and actual results are intentionally
not distributed. The public candidate cannot be used to independently verify
the manuscript's numerical findings without separate authorized access and a
future separately reviewed source-data route.
