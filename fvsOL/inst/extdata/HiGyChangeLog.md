# HiGy.R changelog 

### Version 0.4.2

Parameters and equation forms follow the koa v105final model.

#### Parameters
- `ht.pred.parm`, `ddbh.parm`, `dht.parm` and `mort.s3.parm` refit; `b3` of
  `ddbh.parm` and `dht.parm` fixed at 0.
- `ddbh()`: correction factor 1.57859; origin calibration 0.76421 (planted),
  0.23596 (natural).
- `dht()`: correction factor 2.39567; origin calibration 1.21046 (planted),
  0.17573 (natural).
- `mort_prob()` coefficients, `p.bar`, Garcia `beta`, `h.qmd.a` and `h.qmd.k`
  refit; natural stand mortality factor 2.43532.
- New: `sdi.max` (741.996871565), `rd.floor` (0.001130119195), `ingrowth.parm`.

#### Diameter and height increment
- `calc_plot_summary()` returns relative density `rd`: SDI on the 25.4 cm
  index diameter divided by `sdi.max`.
- BAL terms use the BAL share of plot BA times `rd`; `sqrt(ba*dbh)` and
  `sqrt(ba*ht)` replaced by `sqrt(rd*dbh)` and `sqrt(rd*ht)`.

#### Mortality
- Stage 3: `log(ba+1)` and `log(bal+1)` replaced by `log(rd)` (floored at
  `rd.floor`) and `log(bal/ba*rd+1)`.
- Stage 2: Garcia step from the previous year's `H_QMD` to the current
  year's, carried on the tree list as `h.qmd.prev`. Without it (first year)
  only the background floor applies.

#### BAL
- `calc_bal()`: plot BA times the expansion factor weighted share of the
  other trees with DBH at least as large (ties count as larger), as fitted.
- Tree BA and BAL recomputed after growth are now kept (the result was not
  assigned).

#### Crown ratio
- `make_ops(cr.update)`: 'engine' (default) predicts crown ratio from height
  to crown base each year, bounded 0.05 to 0.95; 'recession' keeps the 1%
  annual limit.
- `(pcr-cr)/cr<-0.01` was parsed as an assignment and stopped every
  `HiGYOneStand()` call; now `< -0.01`.

#### Ingrowth
- `pred_ingrowth()` and `calc_ingrowth()`, off by default
  (`make_ops(use.ingrowth)`). Recruits enter at 2.5 cm DBH, crown ratio 0.6,
  and are not returned to FVS.

#### Size limits, BYI, FVS hand-off
- `tree.size.cap` metric and by origin: DBH 90 cm natural, 69.7 cm planted;
  height 28.0416 m (92 ft). `max.dbh` 90 was converted as inches (228.6 cm).
  `make_fvs_calib()` gains `planted`.
- BYI cap of 600 removed; the model of record has none.
- `make_fvs_tree()`: `dg` returned inside bark (x 0.9, FVS-HI koa bark
  ratio). FVS adds `dg`/bark ratio to DBH, so diameter growth was 11% high.
- customRun_fvsRunHi.R: `h.qmd.prev` carried between cycles of a stand.

#### Checks
- Ten 60 year projections (natural and planted, BYI 142 to 418, one with
  ingrowth) against the v105final engine: largest relative difference
  1.6e-13 tree DBH, 7.7e-12 expansion factor, 6.0e-13 stand TPH.
- customRun_fvsRunHi.R with FVShi, 9 Hawaii stands, 50 years, 5 year
  cycles, NoTriple: no failures, NA or negative values. Tree DBH and height
  within 5.5e-7 and stand BA within 1.0e-6 (relative) of the engine
  projections of the same stands.

### Version 0.4.0

#### Height
- Relative diameter in the height equation changed from `dbh / qmd`
  (quadratic mean diameter) to `dbh / max.plot.dbh` (plot maximum
  diameter, capped at 1). `calc_plot_summary()` gains `max.plot.dbh`.
- All six `ht.pred.parm` coefficients (`a0`, `a1`, `b`, `c`, `g1`, `g2`)
  refit.
- Fixed a latent bug: the base/site parameter filter in `calc_ht()` was
  computed but never applied — the coefficient lookup read from the
  unfiltered table, so the site-quality row was silently never selected
  regardless of `byi`. Now applied correctly.

#### Height-to-crown-base
- Coefficients unchanged.
- Guards added to the equation: `ht/dbh` floored via
  `pmax(ht/pmax(dbh, 0.1), 0.5)`; `byi` itself floored at
  1 (`log(pmax(byi, 1)/100)`) instead of the ratio floored at 0.01.
- Same latent base/site filter bug as Height, fixed the same way.

#### Diameter increment
- All `ddbh.parm` coefficients refit.
- Stand density index and relative height dropped from the equation
  entirely. A competition term (`bal^2/log(dbh+5)`) and a size-competition
  term (`sqrt(ba*dbh)`) added; crown-ratio term simplified to
  `log(pmax(cr, 0.01))`.
- Duan correction factor: 1.026 → 1.36869.
- A multiplicative origin-calibration factor added: 0.40548 (natural) /
  1.43606 (planted).
- Output cap: 0-6 cm/yr → 0-4 cm/yr.
- A breast-height gate added: trees below 1.3716 m get zero diameter
  increment (previously ungated).

#### Height increment
- All `dht.parm` coefficients refit, mirroring the diameter-increment
  restructuring (competition term, simplified crown-ratio term, size term,
  origin-calibration factor).
- Duan correction factor: 0.863 → 1.030.
- Origin-calibration factor added: 0.51917 (natural) / 2.64739 (planted).
- Output cap: 0-6 m/yr → 0-2 m/yr.

#### Mortality
Replaces the single tree-level logistic survival equation with a
three-stage stand mortality model:

- **Stage 1 — plot probability.** A new stand-level logistic model
  (intercept, ln(SDI), planted offset) estimates the annual probability a
  stand experiences *any* mortality at all, and scales the Stage 2 rate
  by how that probability compares to the fitting population's average.
- **Stage 2 — whole-stand self-thinning.** Garcia (2009) self-thinning
  curve replaces the per-tree logistic survival equation. Given
  trees/ha and how much an allometric height-equivalent of QMD (`H_QMD`,
  not a literal top height) grew over the step, it returns the stand
  mortality fraction the self-thinning line implies, floored at a small
  origin-specific background rate (0.3%/yr natural, 0.6%/yr planted, even
  with no growth).
- **Origin-level factor.** Natural stands run this refit's mortality
  about 2.6x higher than planted stands, applied as a fixed multiplier
  after the Stage 1 gate.
- **Stage 3 — per-tree allocation.** The single stand-level death rate is
  spent across individual trees by a logistic weight scoring each tree's
  diameter, height relative to the plot's tallest tree, and local
  competition (plot BA, BA in larger trees) — smaller, more suppressed,
  more crowded trees get a larger share. Per-tree shares are capped at
  95%/yr and rescaled, so they still sum to the plot's target rate.
- **Result cap**: stand mortality rate capped at 95%/yr (was implicit,
  uncapped, in the 0.2.0 per-tree equation).
- Unlike every other equations, mortality does not use site index (`byi`) 
- `calc_plot_summary()` gains `max.plot.ht` (relative-height term) and
  `sdi` (stand density index, Reineke exponent 1.605.


### version 0.2.0
- updated equations- integration of biomass yield index (BYI) and planted indicator
- planted indicator is derived from FVS_Standinit.StdOrgCd (Stand Origin Code; also used by the FIAVBC keyword)
  Sourced from FVS `StdOrgCd` (Stand Origin Code) Event Monitor variable (requires
  SQLIN query and Event Monitor Compute)
    # Natural stand = 0 - established through natural regeneration
    # Plantation = 1 - established through planting

### version 0.1.0
- initial version 
- designed to work with FVS-HI and customRun_fvsRunHi.R
- contains equations for koa

