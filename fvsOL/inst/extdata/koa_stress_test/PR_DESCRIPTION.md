## FVS-HI koa: stress test, CFI validation, and a survival fix

This branch adds a self-contained harness (`fvsOL/inst/extdata/koa_stress_test/`) that exercises the koa equations independent of the compiled variant, with a machine-precision base-R fidelity check against `HiGy.R`. Full write-up in `REPORT.md` and the dated `.docx`.

### What holds up
- `HiGy.R` is on the latest equation iteration (matches `koa_prediction_functions_FINAL.R` v2 for HT, HCB, dDBH, dHT, survival, and CFs).
- Components produce no NaN/Inf/negative values across 81,000 input combinations; increments clip cleanly; natural cohorts track Table 8; SDI stays within the self-thinning envelope.

### The blocker: survival (especially plantations)
- Applied per tree, the published survival cloglog collapses real stands. On the koa CFI/PSP network (first to last measurement) it drives basal area and density bias to about -100%. It is hypersensitive to crown ratio (0.80/yr at CR 0.5, ~0 at CR 0.7) and has no origin term, so dense high-crown plantation trees are wrongly predicted to die.
- We tried refitting. With only 280 deaths and collinear size/site/origin covariates, free GLMs overfit (separation; flipping origin flips survival) and minimal GLMs give absurd rates. The survival signal cannot support a free per-tree GLM.

### Proposed fix
- `koa_survival_calibrated.R`: a calibrated mortality rate anchored to observed annual mortality by origin (natural 1.5%, plantation 0.5%), modulated by relative density and tree size. Stable, monotonic, origin-correct, and reproduces realistic 40-year retention (natural 53%, plantation 66%). Drop-in for `calc_mortality()`.

### Also flagged
- Four code issues: first-cycle HCB unit-order bug in `customRun_fvsRunHi.R`; `Hi.GY()` calls `AcadianGYOneStand()`; tree-size cap units; koa-only parameterization vs the six-species GUI. See `discrepancies.md`.
- Add an ingrowth submodel (main residual gap in young-stand basal area).

Happy to walk through any of it.
