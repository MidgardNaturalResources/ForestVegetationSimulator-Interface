"""
koa_projector.py
================
Two stepping engines for the koa component equations:

1. project_cohort(...)      -- stand-level even-aged cohort. With bounded=True it
   reproduces koa_projection.py's operational guardrails (self-thinning at
   60->55% SDImax, 0.65/0.35 dynamic-static HT blend, biological ceiling
   min(HT, 25+BYI/55), hard DBH max). With bounded=False it strips those
   guardrails to expose the *raw* increment behavior that FVS-HI/HiGy.R would
   show, since HiGYOneStand applies none of them (only per-step increment caps,
   FVS density mortality, and TreeSzCp size caps).

2. project_psp(...)         -- individual-tree list engine mirroring
   HiGYOneStand(): per-tree BA, BAL = cumsum(BA) over descending DBH, plot QMD /
   SDI / htmax, HCB->CR, annual dDBH/dHT, annual survival -> expf decay, 2.5%/yr
   crown-recession limit, size caps. Used to project real PSP tree lists and
   validate against observed remeasurements.

SDI (metric Reineke, used as a lineage-B covariate and for the envelope):
   SDI = TPH * (QMD/25)^1.6     (matches HiGy.R comment expf*(qmd/25)^1.6)
"""
import numpy as np, pandas as pd
from koa_equations import predict_HT, predict_HCB, bal_fraction, LINEAGES

SDI_MAX = 500.0            # upper boundary of FIA koa SDI (stress_test_summary)
FORM_FACTOR = 0.40         # V = BAPH * meanHT * FF
def sdi_of(tph, qmd): return tph * (np.maximum(qmd, 0.1) / 25.0) ** 1.6


# ---------------------------------------------------------------------------
def project_cohort(byi, planted, lineage, bounded=True, surv_mode="cohort", surv_fn=None,
                   init_dbh=None, init_tph=None, init_age=1, max_age=100,
                   ddbh_mult=1.0, dht_mult=1.0, mort_mult=1.0, sdi_max=SDI_MAX):
    # surv_mode: 'cohort'=raw eqn floored at 0.5/yr (koa_projection.py crutch),
    #            'raw'=raw eqn no floor (FVS behavior), 'stable'=clamped+floored
    L = LINEAGES[lineage]
    dbh_max = 60.0 if planted else 90.0
    if init_dbh is None: init_dbh = 6.0 if planted else 8.0
    if init_tph is None: init_tph = 1200.0 if planted else 500.0
    ages = list(range(init_age, max_age + 1)); n = len(ages)
    rec = {k: np.zeros(n) for k in ["QMD","HT","BAPH","TPH","CR","SDI","HCB","VOL"]}
    rec["age"] = np.array(ages)

    DBH, TPH = float(init_dbh), float(init_tph)
    QMD = DBH; BAPH = TPH * np.pi/4 * (DBH/100)**2
    HT = float(predict_HT(DBH, max(BAPH,0.1), QMD, byi)); CR = 0.65
    for step, age in enumerate(ages):
        SDI = sdi_of(TPH, QMD)
        BAL_avg = BAPH * bal_fraction(1.0)
        BAL_dom = BAPH * bal_fraction(1.5)
        HCB = float(predict_HCB(DBH, HT, BAL_avg, BAPH, byi))
        CR = max(0.20, (HT-HCB)/HT) if HT > 0.1 else 0.65
        rHT = 0.50
        rec["QMD"][step], rec["HT"][step], rec["BAPH"][step] = QMD, HT, BAPH
        rec["TPH"][step], rec["CR"][step], rec["SDI"][step] = TPH, CR, SDI
        rec["HCB"][step], rec["VOL"][step] = HCB, BAPH*HT*FORM_FACTOR
        if step == n-1: break

        if bounded and surv_fn is None:              # self-thinning guardrail (off when surv_fn drives mortality)
            pct = SDI / sdi_max
            if pct > 0.60:
                tr = 0.55/pct; tphn = TPH*tr
                DBH *= (TPH/max(tphn,1))**0.08; QMD = DBH; TPH = tphn
                BAPH = TPH*np.pi/4*(QMD/100)**2

        if surv_fn is not None:
            ps = float(np.atleast_1d(surv_fn(DBH, HT, CR, rHT, byi,
                                             bal=BAL_avg, baph=BAPH, planted=planted, sdi=SDI))[0])
            ps = np.clip(1.0 - mort_mult*(1.0-ps), 0.0, 1.0)
        elif surv_mode == "stable":
            ps = float(L.surv_annual_stable(DBH, HT, CR, rHT, byi))
            ps = np.clip(1.0 - mort_mult*(1.0-ps), 0.0, 1.0)
        elif surv_mode == "raw":
            ps = float(L.surv_annual(DBH, HT, CR, rHT, byi, sdi=SDI, planted=planted))
            ps = np.clip(1.0 - mort_mult*(1.0-ps), 0.0, 1.0)
        else:  # 'cohort' (koa_projection.py): floor at 0.50/yr
            ps = float(L.surv_annual(DBH, HT, CR, rHT, byi, sdi=SDI, planted=planted))
            ps = np.clip(1.0 - mort_mult*(1.0-ps), 0.50, 1.0)
        TPH *= ps
        if TPH < 5:
            for s2 in range(step+1, n):
                for k in rec:
                    if k != "age": rec[k][s2] = rec[k][step]
            break
        BAPH = TPH*np.pi/4*(QMD/100)**2

        dD = float(L.dDBH(DBH, BAPH, BAL_avg, CR, byi, planted, sdi=SDI, rht=rHT))*ddbh_mult
        dH = float(L.dHT(HT, BAPH, BAL_avg, CR, byi, planted, sdi=SDI, rht=rHT))*dht_mult
        DBH += dD; QMD = DBH; HT += dH
        if bounded:                                  # HT blend + ceiling + hard cap
            ht_static = float(predict_HT(DBH, BAPH, max(BAL_dom,0.1)*0+QMD, QMD, byi))
            HT = 0.65*HT + 0.35*ht_static
            HT = min(HT, 25 + byi/55)
            DBH = min(DBH, dbh_max); QMD = DBH
        BAPH = TPH*np.pi/4*(QMD/100)**2
    return pd.DataFrame(rec)


# ---------------------------------------------------------------------------
def project_psp(trees, byi, planted, lineage, n_years, surv_mode="raw", surv_fn=None,
                alloc_beta=3.0, ingrowth=False, ingrowth_byi_c=0.0, recruit_dbh=2.5,
                use_size_caps=True, dbh_max=None, ht_max=92*0.3048, return_traj=False):
    """trees: DataFrame with dbh(cm), ht(m), cr(0-1), expf(/ha). One plot.
    surv_mode: 'raw'|'stable'|'calib_alloc'. Numpy engine; return_traj gives the
    annual stand-summary trajectory."""
    from koa_survival_calibrated_py import surv_calibrated
    from koa_ingrowth import ingrowth_annual
    L = LINEAGES[lineage]
    surv = L.surv_annual_stable if surv_mode == "stable" else L.surv_annual
    if dbh_max is None: dbh_max = (60.0 if planted else 90.0)
    dbh = trees.dbh.to_numpy(float); ht = trees.ht.to_numpy(float)
    cr = trees.cr.to_numpy(float); expf = trees.expf.to_numpy(float)
    traj = []; recs = []
    for _ in range(int(n_years)):
        order = np.argsort(-dbh)
        dbh, ht, cr, expf = dbh[order], ht[order], cr[order], expf[order]
        ba = (dbh**2*0.00007854)*expf
        bal = np.cumsum(ba) - ba
        baph = ba.sum(); tph = expf.sum()
        qmd = np.sqrt(baph/(0.00007854*tph)) if tph > 0 else 0.0
        sdi = sdi_of(tph, qmd); htmax = ht.max() if len(ht) else 0.1
        if return_traj:
            traj.append((baph, tph, qmd, sdi))
        hcb = predict_HCB(dbh, ht, bal, baph, byi)
        cr = np.clip(1 - hcb/np.maximum(ht, 0.1), 0.05, 0.95)
        rht = ht/max(htmax, 0.1)
        below_bh = ht < 1.3716
        dD = np.where(below_bh, 0.0, L.dDBH(dbh, baph, bal, cr, byi, planted, sdi=sdi, rht=rht))
        dH = L.dHT(ht, baph, bal, cr, byi, planted, sdi=sdi, rht=rht)
        if surv_mode == "calib_alloc":
            m_stand = 1.0 - float(surv_calibrated(qmd, htmax, 0.5, 0.5, byi, baph=baph, planted=planted, sdi=sdi))
            w = np.exp(-alloc_beta*(dbh/max(qmd,0.1) - 1.0))
            wbar = np.sum(w*expf)/max(np.sum(expf), 1e-9)
            ps = 1.0 - np.clip(m_stand*w/max(wbar,1e-9), 0.0, 0.95)
        elif surv_fn is not None:
            ps = surv_fn(dbh, ht, cr, rht, byi, bal=bal, baph=baph, planted=planted, sdi=sdi)
        else:
            ps = surv(dbh, ht, cr, rht, byi, sdi=sdi, planted=planted)
        expf = np.maximum(expf*ps, 1e-5)
        nd = dbh + dD; nh = ht + dH
        if use_size_caps:
            nd = np.where(nd > dbh_max, dbh, nd); nh = np.where(nh > ht_max, ht, nh)
        dbh, ht = nd, nh
        if ingrowth:
            n_rec = ingrowth_annual(sdi=sdi, planted=planted, byi=byi, byi_c=ingrowth_byi_c)
            recs.append(n_rec)
            if n_rec > 1e-3:
                rh = float(predict_HT(recruit_dbh, baph, max(qmd,1.0), byi))
                dbh = np.append(dbh, recruit_dbh); ht = np.append(ht, rh)
                cr = np.append(cr, 0.6); expf = np.append(expf, n_rec)
        else:
            recs.append(0.0)
        # cap list size: bin to 0.5 cm DBH classes when large (keeps it fast)
        if len(dbh) > 300:
            keyb = np.round(dbh*2)/2.0
            uk = np.unique(keyb); ndbh=[]; nht=[]; ncr=[]; nexpf=[]
            for k in uk:
                m = keyb==k; e=expf[m].sum()
                ndbh.append(np.average(dbh[m],weights=expf[m])); nexpf.append(e)
                nht.append(np.average(ht[m],weights=expf[m])); ncr.append(np.average(cr[m],weights=expf[m]))
            dbh=np.array(ndbh); ht=np.array(nht); cr=np.array(ncr); expf=np.array(nexpf)
    baph = (dbh**2*0.00007854*expf).sum(); tph = expf.sum()
    qmd = np.sqrt(baph/(0.00007854*tph)) if tph > 0 else 0.0
    out = dict(QMD=qmd, BAPH=baph, TPH=tph, SDI=sdi_of(tph, qmd),
               HTmax=ht.max() if len(ht) else 0.0, VOL=baph*ht.mean()*FORM_FACTOR)
    if return_traj:
        tj = pd.DataFrame(traj, columns=["BAPH","TPH","QMD","SDI"])
        tj["ingrowth"] = recs[:len(tj)]
        out["traj"] = tj
    return out
