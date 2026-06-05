## FVS-HI koa: stress test, CFI validation, and a survival fix

This branch adds a self-contained harness (`fvsOL/inst/extdata/koa_stress_test/`) that exercises the koa equations independent of the compiled variant, with a machine-precision base-R fidelity check against `HiGy.R`. Full write-up in `REPORT.md` and the dated `.docx`.

### What holds up
- `HiGy.R` is on the latest equation iteration (matches `koa_prediction_functions_FINAL.R` v2 for HT, HCB, dDBH, dHT, survival, and CFs).
- Components produce no NaN/Inf/negative values across 81,000 input combinations; increments clip cleanly; natural cohorts track Table 8; SDI stays within the self-thinning envelope.

### The blocker: survival (especially plantations)
- Applied per tree, the published survival cloglog collapses real stands. On the koa CFI/PSP network (first to last measurement) it drives basal area and density bias to about -100%. It is hypersensitive to crown ratio (0.80/yr at CR 0.5, ~0 at CR 0.7) and has no origin term, so dense high-crown plantation trees are wrongly predicted to die.
- We tried refitting. With only 280 deaths and collinear size/site/origin covariates, free GLMs overfit (separation; flipping origin flips survival) and minimal GLMs give absurd rates. The survival signal cannot support a free per-tree GLM.

### Proposed fix (tuned on the data, validated to 200 yr)
- `koa_survival_calibrated.R`: annual mortality = base(origin) + a linear self-thinning lift that starts at 0.65 relative density (SDI/SDImax), rises to 0.15/yr at 0.85 RD, and plateaus above. Base 0.5%/yr natural, 0.3%/yr plantation. BYI is deliberately NOT a direct mortality term: the data show its apparent effect is one cluster of high-BYI natural plots with no plantations, so it cannot be estimated for plantations; BYI instead shapes density through growth.
- Data diagnostics (`tune_survival.R`): mortality is lowest in small trees (not highest), so no small-tree term; origin and density are the only robust signals.
- 200-year stress test (`run_final_longterm.py`): no runaway (peak 70-88% of SDImax), no collapse, monotonic QMD across all origin x BYI. Natural QMD approaches ~90 cm, plantations plateau at the 60 cm cap; SDI tracks the self-thinning line.
- PAI/MAI (`run_mai_pai.py`): MAI peak rises with BYI (natural 4.3-8.2, plantation 5.4-10.3 m3/ha/yr); MAI culmination earlier on better sites and for plantations; PAI peaks before MAI culminates. Absolute culmination ages (9-15 yr natural) are on the young side and should be checked against observed koa yield.

### Also flagged
- Four code issues: first-cycle HCB unit-order bug in `customRun_fvsRunHi.R`; `Hi.GY()` calls `AcadianGYOneStand()`; tree-size cap units; koa-only parameterization vs the six-species GUI. See `discrepancies.md`.
- Add an ingrowth submodel (main residual gap in young-stand basal area).

Happy to walk through any of it.
