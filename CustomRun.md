# Adding FVS-HI Custom Run Options (HiGY)

## Overview

The `callCustom` mechanism in `server.R` (around line 3027) is used in the process of calling
R code for growth and mortality projections. Currently the options are Acadian and Adirondack,
defined in `fvsOL/inst/extdata/runScripts.R`. This document describes what is needed to add
the HiGY option.

---

## The System is Convention-Based

The server uses the script name (e.g., `fvsRunAcadian`) to automatically derive everything —
no hardcoded references to Acadian/Adirondack exist in `server.R` or `ui.R`. Adding HI
requires **only one code change** plus two new files.

---

## The Three-Part Mechanism in `server.R`

**1. UI rendering (`server.R:8680-8687`)** — `customRunOps()`:
```r
fn = paste0("customRun_", globals$fvsRun$runScript, ".R")  # customRun_fvsRunHi.R
source(fn, local=TRUE)
uiF = eval(parse(text=paste0(sub("fvsRun","ui", globals$fvsRun$runScript))))
# sub("fvsRun","ui","fvsRunHi") → "uiHi"  ← calls uiHi(globals$fvsRun)
output$uiCustomRunOps = renderUI(uiF(globals$fvsRun))
```

**2. Run-time execution (`server.R:4652-4666`)** — sources and calls:
```r
rsFn = paste0("customRun_", globals$fvsRun$runScript, ".R")  # customRun_fvsRunHi.R
clusterEvalQ(fvschild, source(rsFn))
# ...
cmd = paste0("clusterEvalQ(fvschild,", globals$fvsRun$runScript, "(runOps))")
# → fvsRunHi(runOps)
```

**3. Save/restore (`server.R:8967-8972`)** — generic, iterates `names(uiCustomRunOps)` automatically.

---

## Where `ui.R` Already Handles This

`ui.R:306` just embeds `customRunElements` which is built in `server.R:160-165`:
```r
customRunElements = list(
  selectInput("runScript", ..., choices=runScripts, ...),  # populated from runScripts.R
  uiOutput("uiCustomRunOps"))                              # dynamic panel
```

**No change needed in `ui.R`.**

---

## The One Required Code Change

**`fvsOL/inst/extdata/runScripts.R`** — add the HiGY entry:

```r
customRunScripts=list(
  "AcadianGY (Weiskittel et al.) normally run with FVSne"    = "fvsRunAcadian",
  "AdirondackGY (Weiskittel et al.) normally run with FVSne" = "fvsRunAdirondack",
  "HiGY normally run with FVShi"                             = "fvsRunHi"   # ADD THIS
)
```

---

## The Two New Files Required

**`fvsOL/inst/extdata/customRun_fvsRunHi.R`** — must contain exactly two things,
following the Acadian pattern (`customRun_fvsRunAcadian.R:6` and `:340`):

```r
# 1. The run function — name MUST be fvsRunHi
fvsRunHi <- function(runOps, logfile="Hi.log") {
  rFn = "HiGY.R"
  if (file.exists(rFn)) source(rFn) else {
    rFn = system.file("extdata", rFn, package="fvsOL")
    source(rFn)
  }
  # ... read runOps$uiHi* options, run FVS loop at stopPointCode=5 ...
}

# 2. The UI function — name MUST be uiHi (sub("fvsRun","ui","fvsRunHi") = "uiHi")
uiHi <- function(fvsRun) {
  # set defaults into fvsRun$uiCustomRunOps$uiHi* ...
  list(
    myRadioGroup("uiHiIngrowth", "Simulate ingrowth:", c("Yes","No"), ...),
    # ... etc.
  )
}
```

**`fvsOL/inst/extdata/HiGY.R`** — the Hawaii growth & yield model functions (sourced
inside `fvsRunHi`).

---

## Naming Constraints

| Piece | Required name | Why |
|---|---|---|
| Entry in `runScripts.R` | `"fvsRunHi"` | Used to derive all other names |
| Custom run file | `customRun_fvsRunHi.R` | `paste0("customRun_", runScript, ".R")` |
| Run function | `fvsRunHi` | Called as `fvsRunHi(runOps)` |
| UI function | `uiHi` | `sub("fvsRun","ui","fvsRunHi")` → `"uiHi"` |
| UI input IDs | `uiHi*` prefix (e.g. `uiHiIngrowth`) | Saved/restored via `names(uiCustomRunOps)` |
| GY file | `HiGY.R` | Sourced inside `fvsRunHi` |

---

## Naming Constraint Code References

### 1. Custom run filename → `customRun_fvsRunHi.R`

Two places, same formula:

**`server.R:4654`** (at run time):
```r
4654→  rsFn = paste0("customRun_",globals$fvsRun$runScript,".R")
```

**`server.R:8680`** (when UI option is selected):
```r
8680→  fn=paste0("customRun_",globals$fvsRun$runScript,".R")
```

---

### 2. UI function name → `uiHi`

**`server.R:8685`**:
```r
8685→  uiF = try(eval(parse(text=paste0(sub("fvsRun","ui",globals$fvsRun$runScript)))))
```
`sub("fvsRun", "ui", "fvsRunHi")` → `"uiHi"` — that string is then `eval`'d, so a function
named exactly `uiHi` must exist in the sourced file.

---

### 3. Run function call → `fvsRunHi(runOps)`

**`server.R:4680`**:
```r
4680→  cmd = paste0("clusterEvalQ(fvschild,",globals$fvsRun$runScript,"(runOps))")
```
`globals$fvsRun$runScript` is `"fvsRunHi"` → executes `fvsRunHi(runOps)`.

---

### 4. Script name `"fvsRunHi"` set when user picks from dropdown

**`server.R:8676`**:
```r
8676→  globals$fvsRun$runScript = input$runScript
```
`input$runScript` is the **value** (right-hand side) from `runScripts.R`. That value drives
all three derivations above.

---

### 5. Dropdown choices sourced from `runScripts.R`

**`server.R:154-163`**:
```r
154→  rsf <- "runScripts.R"
155→  if (file.exists(rsf)) source(rsf) else source(system.file("extdata", rsf, package="fvsOL"))
156→
157→  runScripts <- if (exists("customRunScripts") && length(customRunScripts))
158→                  append(x=customRunScripts,after=0,defaultRun) else defaultRun
159→
160→  customRunElements = list(
161→    selectInput("runScript",
162→                "Select run script (normally, use the default)",
163→                choices=runScripts,
```

The value side of `runScripts.R` (`"fvsRunAcadian"`, `"fvsRunAdirondack"`) is what becomes
`input$runScript` and `globals$fvsRun$runScript` — the root of all derived names.

---

### 6. UI input IDs saved/restored generically

**`server.R:8970-8971`** (`saveRun`):
```r
8970→  for (item in names(globals$fvsRun$uiCustomRunOps))
8971→    globals$fvsRun$uiCustomRunOps[[item]] = input[[item]]
```
`uiCustomRunOps` is populated by `uiHi()` when it sets defaults like
`fvsRun$uiCustomRunOps$uiHiIngrowth = "No"`. The key names (e.g. `uiHiIngrowth`) must match
the `inputId` strings used in the `myRadioGroup`/`myInlineTextInput` calls inside `uiHi()`.
