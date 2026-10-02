# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

R code for a Shiny-based interface to the USDA Forest Vegetation Simulator (FVS). The repo has no README, test suite, or linter config. It contains two R packages and two small standalone Shiny apps:

- `fvsOL/`: the main package, an R-Shiny UI for FVS run in "Online" or "Onlocal" configuration. Depends on `rFVS`.
- `rFVS/`: thin R wrappers (`fvs*` functions) around FVS compiled as a shared library (`dyn.load` in `fvsLoad.R`). Run-control functions such as `fvsRun` stop FVS at numbered stop points so R can read or modify tree and species attributes mid-simulation.
- `FVSPrjBldr/`, `FVSDataConvert/`: standalone Shiny apps (`ui.R` and `server.R`, no package). The first creates projects; the second converts input data.

## Build / install

Each package has a makefile that runs `devtools::document`, `build`, and `install`. A `*madeTag` file records a successful build. Install `rFVS` before `fvsOL`.

```
make -C rFVS            # document + build + install rFVS
make -C fvsOL           # also regenerates data/*.RData first
make -C fvsOL clean     # remove generated RData and tag
```

`fvsOL/makefile` also builds two generated data files that are loaded at runtime:
- `data/prms.RData` comes from `Rscript parms/mkpkeys.R`. It parses the `//start name` / `//end` sections of the `fvsOL/parms/*.kwd|*.prm` files into keyword and parameter-form definitions.
- `data/fvsOnlineHelpRender.RData` comes from `inst/extdata/mkhelp.R`, which renders `fvsOnlineHelp.html` and `databaseDescription.xlsx`.

Rebuild whenever you edit `parms/` or the help sources. Packages use roxygen (`Roxygen: list(markdown = TRUE)`), so edit the `#'` comments, not `NAMESPACE` or `man/`.

## Running

`fvsOL::fvsOL(prjDir=, runUUID=, fvsBin=)` is the entry point. It `setwd`s into the project directory and requires a directory of compiled FVS libraries (`FVSbin` in the project dir, or the `fvsBin` argument). It then sets globals via `<<-` and calls `shinyApp(FVSOnlineUI, FVSOnlineServer)`. Running FVS needs those platform-specific binaries, which are not in this repo.

`rFVS/tests/` holds only an example keyword/tree file pair (`iet01.key`, `.tre`), not automated tests.

## fvsOL architecture

- `R/server.R` (~9k lines) is a single `FVSOnlineServer` function holding nearly all reactive logic. `R/ui.R` defines `FVSOnlineUI`. Expect to edit both for any UI feature.
- A project is a directory containing `FVS_Data.db` (input SQLite, seeded from `inst/extdata/FVS_Data.db.default` if absent), a project database of saved runs, and `projectId.txt`.
- A run is the RefClass `fvsRun` (`mkfvsRun` in `server.R`), serialized into the project database via `storeFVSRun`. Its `uuid` identifies the run, and its components are `fvsCmp` keyword-component objects.
- Keyword components are built from the UI definitions in `parms/` (via `prms.RData`) using `mkInputElements.R` and `componentWins.R`. `writeKeyFile.R` turns the selected stands and components into an FVS keyword file.
- `fvsRunUtilities.R` has the database helpers, run execution, and output processing. `fvsOutUtilities.R` has the output and plot helpers. `svsTree.R` is the SVS 3D tree view (rgl).
- `externalCallable.R` is the exported `extn*` API (`extnMakeRun`, `extnAddStands`, `extnSimulateRun`, and others) for creating and running projects without the GUI. Keep it in sync with the internal functions it wraps.
- FVS runs in child processes (`parallel` cluster, with `rFVS` loaded in each). Custom run scripts live in `inst/extdata/customRun_*.R` and are registered in `inst/extdata/runScripts.R`. Only scripts whose file exists are offered in the UI. Never allow user-uploaded scripts in a client/server deployment.
- `inst/extdata/sqlQueries.R` and `sqlQueries_Metric.R` define the output queries (US and metric units), selected at runtime.
- The `AcadianGY`, `AdirondackGY`, and `HiGy` files in `inst/extdata` are growth-model extensions run alongside FVS (northeast US and Hawaii variants).

## Conventions

- Default code owners are listed in `.github/CODEOWNERS`.
- Package version strings are date-based (`YYYY.MM.DD`) in each `DESCRIPTION`.
