const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        ImageRun, AlignmentType, LevelFormat, HeadingLevel, BorderStyle,
        WidthType, ShadingType, PageNumber, Header, Footer } = require("docx");

const R = "/sessions/elegant-lucid-cannon/mnt/outputs/koa_harness/results";
const CW = 9360; // content width DXA (US Letter, 1in margins)
const border = { style: BorderStyle.SINGLE, size: 1, color: "B0B0B0" };
const borders = { top: border, bottom: border, left: border, right: border };
const HEAD = "1F4E79", ALT = "EAF1F8";

const P = (text, opts={}) => new Paragraph({ spacing:{after:120}, ...opts,
  children: (Array.isArray(text)?text:[new TextRun({text, ...(opts.run||{})})]) });
const H1 = t => new Paragraph({ heading: HeadingLevel.HEADING_1, children:[new TextRun(t)] });
const H2 = t => new Paragraph({ heading: HeadingLevel.HEADING_2, children:[new TextRun(t)] });
const bullet = t => new Paragraph({ numbering:{reference:"b",level:0}, spacing:{after:80},
  children:[new TextRun(t)] });

function table(headers, rows, widths){
  const cell = (txt, w, opts={}) => new TableCell({ borders, width:{size:w,type:WidthType.DXA},
    margins:{top:60,bottom:60,left:100,right:100},
    shading: opts.fill?{fill:opts.fill,type:ShadingType.CLEAR}:undefined,
    children:[new Paragraph({children:[new TextRun({text:String(txt),bold:opts.bold,
      color:opts.color, size:18})]})]});
  const head = new TableRow({ tableHeader:true, children: headers.map((h,i)=>
    cell(h,widths[i],{bold:true,fill:HEAD,color:"FFFFFF"})) });
  const body = rows.map((r,ri)=> new TableRow({ children: r.map((c,i)=>
    cell(c,widths[i],{fill: ri%2?ALT:undefined})) }));
  return new Table({ width:{size:CW,type:WidthType.DXA}, columnWidths:widths, rows:[head,...body] });
}
function img(file,w,h){ return new Paragraph({alignment:AlignmentType.CENTER, spacing:{before:80,after:80},
  children:[new ImageRun({type:"png",data:fs.readFileSync(`${R}/${file}`),
  transformation:{width:w,height:h}, altText:{title:file,description:file,name:file}})]}); }
const cap = t => new Paragraph({alignment:AlignmentType.CENTER, spacing:{after:160},
  children:[new TextRun({text:t,italics:true,size:16,color:"555555"})]});

const children = [];
// Title
children.push(new Paragraph({spacing:{after:60}, children:[new TextRun({
  text:"FVS-HI Koa Growth & Yield Model", bold:true, size:40, color:HEAD})]}));
children.push(new Paragraph({spacing:{after:40}, children:[new TextRun({
  text:"Comprehensive Stress Test, CFI Validation, and Survival Refinement", size:26, color:"333333"})]}));
children.push(P([new TextRun({text:"Prepared for: ",bold:true}), new TextRun("Ben Rice (Midgard Natural Resources), FVS-HI PR #30")]));
children.push(P([new TextRun({text:"Prepared by: ",bold:true}), new TextRun("Aaron Weiskittel, University of Maine  |  2026-06-05")]));
children.push(P([new TextRun({text:"Code basis: ",bold:true}), new TextRun("HiGy.R v0.2.0 (== koa_prediction_functions_FINAL.R v2, 2026-05-12); harness branch koa-stress-test")]));

// 1 Executive summary
children.push(H1("1. Executive summary"));
children.push(P("A comprehensive stress test of the current FVS-HI koa equations confirms the static and increment components are robust and track the manuscript Table 8 trajectories, but the survival component is not safe to ship unguarded. Applied per tree as FVS does, the published survival model collapses real stands."));
children.push(bullet("The harness reproduces HiGy.R to machine precision (max difference 1e-14 over 3,240 grid points), so these findings apply directly to the shipped code."));
children.push(bullet("HiGy.R is on the latest equation iteration: it matches koa_prediction_functions_FINAL.R (v2, post-pivot) for every component. The FVS_FINAL README and parameters.csv recommend a different (superseded) form and CFs and should be reconciled."));
children.push(bullet("Components produce no NaN, Inf, or negative values across 81,000 input combinations and clip cleanly (dDBH at 4 cm/yr, dHT at 2 m/yr)."));
children.push(bullet("Survival is the problem: 67% of the koa-plausible crown-ratio by relative-height domain yields annual survival below 0.90, and survival drops from 0.80 at CR 0.5 to ~0 at CR 0.7. Koa crowns are routinely 0.5 to 0.8."));
children.push(bullet("On the koa CFI/PSP network (first to last remeasurement), raw survival drives stand basal area and density bias to -99.9% (near-total mortality). A stabilized survival (clamp crown ratio and relative height to the calibration domain, floor annual survival at 0.90) recovers density bias to -9.3%."));
children.push(bullet("Recommendation: do not ship the raw survival cloglog. Ship the stabilizer as an interim or refit survival for individual-tree use. Keep the rest of the equations. Add an ingrowth submodel. Fix four code issues."));

// 2 Scope & method
children.push(H1("2. Scope and method"));
children.push(P("Two engines exercise the equations independent of the compiled FVS variant: a stand-level cohort projector (run with the operational guardrails of koa_projection.py, and with them stripped to expose raw behavior) and an individual-tree-list engine that mirrors HIGYOneStand for real plots. Survival is run in three modes: cohort (raw equation floored at 0.50/yr, the koa_projection.py crutch), raw (no floor, the FVS behavior), and stable (covariates clamped to the calibration domain and floored). Port fidelity to HiGy.R is verified in base R against the verbatim function bodies."));

// 3 Provenance
children.push(H1("3. Equation provenance (which version is current)"));
children.push(P("The latest finalized code, koa_prediction_functions_FINAL.R (v2, 2026-05-12 post-pivot), uses the NEW increment form, conditional Duan correction factors (1.026 and 1.030), and the published survival fit (intercept 18.133). HiGy.R v0.2.0 in PR #30 matches that file coefficient for coefficient across height, height-to-crown-base, diameter increment, height increment, survival, and the correction factors. The FVS_FINAL README TL;DR (old form, marginal factors) and koa_FVS_FINAL_parameters.csv (M1/S1) reflect a pre-pivot decision that was reversed and are stale. This contradiction is documented in discrepancies.md and should be resolved before release."));

// 4 Robustness
children.push(H1("4. Component robustness"));
children.push(P("Each component was evaluated over 81,000 input combinations spanning DBH 0.5 to 110 cm, height 1.4 to 35 m, basal area in larger trees 0 to 80, stand basal area 0.5 to 80, crown ratio 0.05 to 0.95, BYI 0 to 600, and both origins."));
children.push(table(["Component","NaN","Inf","Negative","Min","Max"],
  [["Height (m)","0","0","0","1.37","31.94"],
   ["Height to crown base (m)","0","0","0","0.08","30.40"],
   ["Diameter increment (cm/yr)","0","0","0","0.00","4.00"],
   ["Height increment (m/yr)","0","0","0","0.00","2.00"],
   ["Annual survival","0","0","0","0.00","1.00"]],
  [2600,900,900,1330,1815,1815]));
children.push(P("No numerical failures. Increment components clip at their documented ceilings. The survival range spanning the full 0 to 1 interval is the first signal of the instability quantified next."));

// 5 Survival instability
children.push(H1("5. Survival instability"));
children.push(P("The published survival cloglog has very large coefficients and is hypersensitive to crown ratio and BYI. At DBH 15, height 9, relative height 0.5, BYI 264:"));
children.push(table(["Crown ratio","0.4","0.5","0.6","0.7"],
  [["Annual survival","0.99","0.80","0.02","0.00"],
   ["30-year cumulative","0.82","0.001","~0","~0"]],
  [2760,1650,1650,1650,1650]));
children.push(img("fig_survival_heatmap.png", 430, 250));
children.push(cap("Figure 1. Annual survival across the crown-ratio by relative-height domain. Most of the koa-plausible space (67% of cells) falls below 0.90/yr; survival craters above crown ratio 0.55."));
children.push(P("The cohort projector only behaves because it fixes relative height at 0.5, clamps crown ratio at or above 0.20, and floors annual survival at 0.50. FVS applies the raw model with none of these crutches."));

// 6 Long horizon
children.push(H1("6. Long-horizon behavior"));
children.push(P("Natural-origin cohort QMD (cm) at age 40 (manuscript Table 8 expectation about 42 cm):"));
children.push(table(["BYI","Bounded cohort","Raw unbounded","Stable unbounded"],
  [["100","43.4","43.0","42.9"],["264","56.9","56.3","54.9"],["450","57.2","52.8","42.5"]],
  [2340,2340,2340,2340]));
children.push(P("For realistic natural starting densities the increment and height equations track Table 8 and do not run away. Runaway risk is confined to very small initial DBH and dense planted starts, the regime to confirm in true FVS runs."));
children.push(img("fig_cohort_trajectories.png", 470, 178));
children.push(cap("Figure 2. Natural cohort QMD and SDI trajectories (BYI 264). Bounded and raw increment paths are close; SDI stays within the envelope under self-thinning."));

// 7 Extrapolation
children.push(H1("7. Out-of-source extrapolation"));
children.push(P("Four off-distribution points show the survival contrast between raw and stabilized; increment predictions remain bounded."));
children.push(table(["Case","dDBH","dHT","Survival raw","Survival stable"],
  [["Kahikinui dry small","0.84","0.30","0.96","0.96"],
   ["Wet mature natural","0.34","0.06","0.00","0.90"],
   ["Old-growth low-BYI","0.20","0.04","0.00","0.90"],
   ["Dense plantation","2.34","0.55","0.00","0.90"]],
  [2760,1300,1300,2000,2000]));
children.push(P("Wherever crown ratio or BYI move off the cohort-average regime, raw survival collapses to zero while the stabilizer holds a defensible floor. (Illustrative dDBH/dHT shown; see extrapolation.csv.)"));

// 8 CFI validation
children.push(H1("8. CFI / PSP validation (first to last measurement)"));
children.push(P("Individual-tree projection of the koa CFI/PSP network (13 plots with usable remeasurements, about 9 to 14 year spans, BYI assumed 264) from first measurement to observed last:"));
children.push(table(["Variable","Survival","Mean bias","Bias %","RMSE"],
  [["QMD (cm)","raw","-3.2","-13.3%","5.4"],
   ["QMD (cm)","stable","-6.7","-27.5%","8.1"],
   ["Basal area (m2/ha)","raw","-9.0","-99.9%","9.7"],
   ["Basal area (m2/ha)","stable","-6.4","-71.0%","7.5"],
   ["Density (trees/ha)","raw","-186","-99.9%","199"],
   ["Density (trees/ha)","stable","-17","-9.3%","170"]],
  [2760,1500,1600,1750,1750]));
children.push(img("fig_survival_dialin.png", 470, 178));
children.push(cap("Figure 3. Left: survival versus crown ratio with the stabilizer clamp region. Right: CFI first-to-last bias, raw versus stabilized survival."));
children.push(P("Raw survival wipes out the stands. The stabilizer recovers realistic density (bias -9.3%). The residual basal-area shortfall is largely the missing ingrowth submodel: these are young building koa stands where observed basal area and density rise through recruitment that a no-ingrowth projection cannot reproduce."));

// 9 Dial-in
children.push(H1("9. Survival dial-in"));
children.push(P("Sweeping the stabilizer settings against CFI density bias shows the annual-survival floor is the dominant lever; the crown-ratio clamp ceiling has little effect once the floor is active."));
children.push(table(["Crown-ratio ceiling","Floor 0.85","Floor 0.90","Floor 0.95"],
  [["0.50","-35.1%","-9.2%","+29.1%"],["0.55","-35.2%","-9.3%","+29.1%"],["0.60","-35.2%","-9.3%","+29.1%"]],
  [2760,2200,2200,2200]));
children.push(P([new TextRun({text:"Recommended interim setting: ",bold:true}),
  new TextRun("clamp crown ratio to [0.20, 0.55] and relative height to [0.30, 0.70], floor annual survival near 0.90 to 0.92. This is an operational stopgap, not a refit; the proper fix is to refit survival for individual-tree application on the current AK.SURV records.")]));

// 10 Bakuzis
children.push(H1("10. Bakuzis / SDI envelope"));
children.push(P("Across both origins and three site classes, 100% of decade snapshots remain within 110% of the maximum stand density index (500) under self-thinning. The equations respect the self-thinning boundary; they do not push stands above the observed density envelope."));

// 11 Survival solution
children.push(H1("11. Survival solution for plantations and individual trees"));
children.push(P("The published survival model is not a usable per-tree solution, and a free refit does not fix it. With only 280 deaths and crown ratio, size, site, and origin strongly collinear, fitted GLMs overfit: coefficients are extreme (intercept 44, planted -24, BYI/1000 -87) and behave backwards out of sample, predicting near-total plantation mortality when origin is flipped at fixed size; minimal stable-coefficient GLMs instead predict absurd rates (about 34% per year). The survival signal cannot support a free per-tree GLM."));
children.push(H2("11.1 What the data actually show"));
children.push(P("Exposure-based (deaths per tree-year) diagnostics on the 6,489-record survival dataset reframe the problem:"));
children.push(bullet("Mortality is lowest in small trees (0 to 5 cm: 0.04%/yr; 5 to 10 cm: 0.56%/yr) and peaks in the 10 to 20 cm class (2.77%/yr). A small-tree mortality term is therefore wrong for koa."));
children.push(bullet("The apparent BYI effect is an artifact of one cluster: the top BYI tercile (above 408) shows 3.3%/yr versus about 0.15%/yr in the lower two, and that cluster is entirely natural stands. Plantations do not span high BYI, so a BYI mortality effect cannot be estimated for plantations at all. Forcing BYI into a GLM produces unstable 890-fold rate ratios."));
children.push(bullet("The robust, real signals are origin (plantations about 0.5%/yr versus natural) and density (mortality rises through the 0.55 to 0.8 relative-density self-thinning zone)."));
children.push(H2("11.2 Recommended model"));
children.push(P("Annual mortality is a low density-independent background that differs by origin, plus a density-dependent self-thinning lift that starts at 0.65 relative density (SDI/SDImax), increases linearly to a maximum at 0.85 RD, and plateaus above. BYI is deliberately not a direct mortality driver; it influences long-term density correctly through growth, by driving stands into self-thinning sooner. Drop-in code is koa_survival_calibrated.R; the constants are interpretable and meant to be refined as data accrue."));
children.push(P([new TextRun({text:"mortality = base(origin) + maxlift × clamp((SDI/SDImax − 0.65) / 0.20, 0, 1),",italics:true}),
  new TextRun("  with base 0.005/yr natural, 0.003/yr plantation, maxlift 0.15/yr reached at 0.85 RD, capped at 0.20/yr.")]));
children.push(H2("11.3 Final long-term stress test (200-year projections)"));
children.push(P("Cohort QMD (cm) / trees per hectare / percent of SDImax, final model, across origin and site:"));
children.push(table(["Origin / BYI","Age 25","Age 50","Age 100","Age 150","Age 200"],
  [["Natural BYI 100","27/307/69","40/157/67","59/83/65","72/59/65","83/46/63"],
   ["Natural BYI 264","34/217/70","50/110/67","72/60/65","87/44/65","90/34/53"],
   ["Natural BYI 450","38/178/70","56/91/67","80/51/65","90/38/59","90/29/46"],
   ["Plantation BYI 264","41/163/72","60/79/64","60/68/55","60/58/47","60/50/41"]],
  [2360,1400,1400,1400,1400,1400]));
children.push(img("fig_longterm_calibrated.png", 480, 148));
children.push(cap("Figure 4. 200-year behavior, final calibrated survival (ramp 0.65 to 0.85 RD). Left: SDI rises to 72-88% of SDImax then tracks the self-thinning line. Center: QMD approaches ~90 cm (natural) or the 60 cm cap (plantation). Right: stand trajectory in QMD by density space."));
children.push(P("Robustness checks across all six origin by BYI combinations pass: no runaway (peak 70 to 88% of SDImax), no collapse, and monotonic QMD over 200 years. Background annual mortality is 0.5% natural and 0.3% plantation, rising to about 15% per year at full self-thinning. Higher BYI reaches the self-thinning onset sooner and grows larger, so BYI shapes density through growth rather than a confounded mortality term. Plantations self-thin earlier from denser starts and plateau at the 60 cm cap. This is the long-term behavior the model was missing."));

// 11.4 PAI/MAI
children.push(H2("11.4 Volume increment: PAI and MAI culmination"));
children.push(P("Projecting from establishment (2 cm, age 1) to 200 years, periodic annual increment (PAI) and mean annual increment (MAI) of stem volume behave as expected, and MAI culmination (biological rotation age) responds correctly to site and origin:"));
children.push(table(["Origin / BYI","MAI culmination age","MAI peak (m3/ha/yr)","PAI peak age","Volume age 200"],
  [["Natural BYI 100","15","4.3","13","171"],
   ["Natural BYI 264","11","6.5","9","158"],
   ["Natural BYI 450","9","8.2","7","140"],
   ["Plantation BYI 100","10","5.4","8","110"],
   ["Plantation BYI 264","7","8.1","6","106"],
   ["Plantation BYI 450","7","10.3","5","104"]],
  [2360,1900,1850,1450,1800]));
children.push(img("fig_mai_pai.png", 470, 184));
children.push(cap("Figure 5. PAI (thin) and MAI (bold) for natural and plantation stands by BYI; dotted lines mark MAI culmination. PAI peaks before MAI culminates in every case."));
children.push(P("Trends are correct: MAI peak rises with BYI (natural 4.3 to 8.2, plantation 5.4 to 10.3 m3/ha/yr), culmination occurs earlier on better sites and earlier for plantations than natural stands, and PAI peaks before MAI culminates. The absolute culmination ages are on the young side for koa (9 to 15 years natural), reflecting rapid early increment and a constant form factor (0.40) in the simple volume proxy; these should be checked against observed koa yield data and the manuscript Table 8 before being used for rotation guidance. The qualitative site and origin ordering is robust."));

// 12 Code issues
children.push(H1("12. Open code issues"));
children.push(bullet("Unit-order bug in the first-cycle crown imputation in customRun_fvsRunHi.R: height is converted to feet before calc_hcb, which expects metres, biasing imputed crown ratios. Compute crown base while metric, convert to feet last."));
children.push(bullet("Hi.GY() calls AcadianGYOneStand() instead of HIGYOneStand(); latent but should be corrected or the wrapper removed."));
children.push(bullet("Tree-size cap units: tree.size.cap AK=(90,92) is multiplied by inch-to-cm and foot-to-metre; if 90 was cm the DBH cap is effectively disabled at 228 cm."));
children.push(bullet("Species scope: only koa is parameterized; ohia and sandalwood map to OT and are held static. Document or add surrogates."));

// 13 Recommendations
children.push(H1("13. Recommendations"));
children.push(new Paragraph({numbering:{reference:"n",level:0},children:[new TextRun("Replace the published survival cloglog with the calibrated rate (koa_survival_calibrated.R). Do not ship the raw cloglog; a free refit is not better.")]}));
children.push(new Paragraph({numbering:{reference:"n",level:0},children:[new TextRun("Keep the current height, height-to-crown-base, and increment equations with conditional correction factors (matches the FINAL code and tracks Table 8).")]}));
children.push(new Paragraph({numbering:{reference:"n",level:0},children:[new TextRun("Add an ingrowth submodel, or document the no-ingrowth limitation; it is the main residual gap in young-stand basal area and density.")]}));
children.push(new Paragraph({numbering:{reference:"n",level:0},children:[new TextRun("Fix the four code issues in section 11.")]}));
children.push(new Paragraph({numbering:{reference:"n",level:0},children:[new TextRun("Reconcile the FVS_FINAL README and parameters.csv with the actual FINAL code (see discrepancies.md).")]}));
children.push(new Paragraph({numbering:{reference:"n",level:0},children:[new TextRun("Validate in true FVS on the dense-planted, small-DBH regime that the equation-level harness cannot fully exercise.")]}));

children.push(H1("Appendix: reproduce"));
children.push(P("Branch koa-stress-test, folder fvsOL/inst/extdata/koa_stress_test: Rscript fidelity_check.R (port check); python3 run_stress_comprehensive.py (this report); python3 run_cfi.py (CFI validation). Python needs numpy, pandas, matplotlib; R needs no packages."));

const doc = new Document({
  styles: { default:{document:{run:{font:"Arial",size:22}}},
    paragraphStyles:[
      {id:"Heading1",name:"Heading 1",basedOn:"Normal",next:"Normal",quickFormat:true,
        run:{size:28,bold:true,color:HEAD,font:"Arial"},paragraph:{spacing:{before:260,after:140},outlineLevel:0}},
      {id:"Heading2",name:"Heading 2",basedOn:"Normal",next:"Normal",quickFormat:true,
        run:{size:24,bold:true,color:"2E5E8C",font:"Arial"},paragraph:{spacing:{before:180,after:100},outlineLevel:1}}]},
  numbering:{config:[
    {reference:"b",levels:[{level:0,format:LevelFormat.BULLET,text:"•",alignment:AlignmentType.LEFT,
      style:{paragraph:{indent:{left:540,hanging:280}}}}]},
    {reference:"n",levels:[{level:0,format:LevelFormat.DECIMAL,text:"%1.",alignment:AlignmentType.LEFT,
      style:{paragraph:{indent:{left:540,hanging:280}}}}]}]},
  sections:[{
    properties:{page:{size:{width:12240,height:15840},margin:{top:1440,right:1440,bottom:1440,left:1440}}},
    footers:{default:new Footer({children:[new Paragraph({alignment:AlignmentType.CENTER,
      children:[new TextRun({text:"FVS-HI koa stress test  |  ",size:16,color:"888888"}),
                new TextRun({children:["Page ",PageNumber.CURRENT],size:16,color:"888888"})]})]})},
    children }]
});
Packer.toBuffer(doc).then(b=>{ fs.writeFileSync(process.argv[2],b); console.log("wrote",process.argv[2]); });
