#!/usr/bin/env bash
# regenerate_all.sh -- rerun the entire koa FVS-HI analysis pipeline end to end,
# regenerate all CSVs and figures, rebuild the docx report, and log every step.
set -u
cd "$(dirname "$0")"
export R_LIBS_USER=/sessions/elegant-lucid-cannon/mnt/outputs/.Rlib
REPO=/sessions/elegant-lucid-cannon/mnt/outputs/fvshi_repo/fvsOL/inst/extdata/koa_stress_test
OUTDOC="/sessions/elegant-lucid-cannon/mnt/Koa/koa_fvs_stress_test/20260605_koa_FVS-HI_stress-test-report.docx"
LOG=results/processing-log.md
mkdir -p results
echo "# Processing log — koa FVS-HI regenerate ($(date))" > $LOG
run(){ echo ">>> $*"; echo "- \`$*\`" >> $LOG; "$@" >/tmp/step.out 2>&1; \
       if [ $? -eq 0 ]; then echo "  - OK" >> $LOG; else echo "  - FAILED" >> $LOG; tail -3 /tmp/step.out; fi; }

echo "=== R: data extraction + model fits ==="
run Rscript extract_ingrowth_psp.R
run Rscript fidelity_check.R
run Rscript fit_survival_candidates.R
run Rscript tune_survival.R
run Rscript check_origin_byi.R
run Rscript check_interaction.R
run Rscript fit_ingrowth.R
run Rscript fit_ingrowth_rd.R
run Rscript fit_ingrowth_byi.R

echo "=== Python: projections, stress tests, validation ==="
run python3 run_benchmarks.py
run python3 run_cfi.py
run python3 run_stress_comprehensive.py
run python3 run_final_longterm.py
run python3 run_mai_pai.py
run python3 run_survival_solution.py
run python3 run_origin_byi_final.py
run python3 run_alloc_test.py
run python3 run_ingrowth_longterm.py

echo "=== Figures ==="
for f in _figs.py _figs2.py _fig3.py _figsurv.py _fig_lt.py _figsurv2.py _fig_ing.py; do
  [ -f "$f" ] && run python3 "$f"
done

echo "=== Rebuild report ==="
run node build_report_docx.js "$OUTDOC" "$REPO/koa_survival_calibrated.R" "$REPO/koa_ingrowth.R"
echo ">>> validate"; python3 /sessions/elegant-lucid-cannon/mnt/.claude/skills/docx/scripts/office/validate.py "$OUTDOC" 2>&1 | tail -1
echo "Done. See $LOG"
