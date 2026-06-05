from koa_projector import project_cohort, SDI_MAX
from koa_survival_calibrated_py import make
print("maxlift x mort_max -> peak %SDImax for the stressors and natural")
print(f"{'maxlift':>7} {'mortmax':>7} {'nat450':>7} {'plt264':>7} {'plt450':>7}")
for ml in (0.12,0.15,0.18,0.22):
    for mm in (0.18,0.25):
        peaks={}
        for tag,(b,p) in {"nat450":(450,0),"plt264":(264,1),"plt450":(450,1)}.items():
            df=project_cohort(b,p,"A",bounded=True,surv_fn=make(maxlift=ml,mort_max=mm),max_age=200)
            peaks[tag]=df.SDI.max()/SDI_MAX*100
        print(f"{ml:7.2f} {mm:7.2f} {peaks['nat450']:7.0f} {peaks['plt264']:7.0f} {peaks['plt450']:7.0f}")
