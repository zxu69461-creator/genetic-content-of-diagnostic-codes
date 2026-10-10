"""Reproduce the fixed revision analyses from small, frozen, local inputs."""
from pathlib import Path
import json,hashlib,resource,sys
import numpy as np,pandas as pd
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]; FROZEN_RESULTS=ROOT/'results'
OUT=ROOT/'recomputed' if (ROOT.parent/'SHA256.json').exists() else FROZEN_RESULTS
OUT.mkdir(exist_ok=True)
# The publication archive places these two frozen aggregate extracts beside other inputs.
EVIDENCE=ROOT/'inputs' if (ROOT/'inputs/g30_g4_fetch.parquet').exists() else ROOT.parent/'evidence'
L=pd.read_csv(ROOT/'protocol/frozen_locus_panel.tsv',sep='\t')
L['x']=np.where(L.ea==L.alt,L.beta,-L.beta);L=L.rename(columns={'se':'sx'})
EXCLUDED=[];OBS=[]
def align(d,cohort,label):
    rows=[]
    for r in L.itertuples():
        f=d[(d.chrom==r.chrom)&(d.pos==r.pos)]
        m=f[((f.ref==r.ref)&(f.alt==r.alt))|((f.ref==r.alt)&(f.alt==r.ref))]
        reason=None
        if len(m)!=1:reason='missing_or_ambiguous_alleles'
        else:
            z=m.iloc[0];beta=float(z.beta);se=float(z.se)
            if not np.isfinite(beta) or not np.isfinite(se) or se<=0:reason='invalid_effect_or_se'
        if reason:
            EXCLUDED.append({'cohort':cohort,'label':label,'rsid':r.rsid,'reason':reason});continue
        flip=(z.ref==r.alt)
        row={'cohort':cohort,'label':label,'rsid':r.rsid,'chrom':r.chrom,'pos':r.pos,'ref':r.ref,'alt':r.alt,
             'x':r.x,'sx':r.sx,'igas_p':r.p,'y':-beta if flip else beta,'sy':se,'flipped':bool(flip)}
        rows.append(row)
    OBS.extend(rows);return pd.DataFrame(rows)

datasets=[]
labels=['RHEUMA_SERONEG','RHEUMA_SEROPOS_WIDE','M13_ANKYLOSPON','RHEUMA_OTHER_WIDE']
for label in labels:
    p=ROOT/'inputs'/('ukb_targeted_'+label+'.tsv')
    if p.exists():
        raw=pd.read_csv(p,sep='\t',na_values=['NA'])
        d=raw.rename(columns={'CHR':'chrom','POS':'pos','REF':'ref','ALT':'alt','UKBB_beta':'beta','UKBB_sebeta':'se'})
        datasets.append(align(d,'UK Biobank',label))
    else:EXCLUDED.append({'cohort':'UK Biobank','label':label,'rsid':'all','reason':'source_unavailable'})
FG=pd.read_parquet(EVIDENCE/'g30_g4_fetch.parquet')
extra=ROOT/'inputs/finngen_seropos_panel.parquet'
if extra.exists():FG=pd.concat([FG,pd.read_parquet(extra)],ignore_index=True).drop_duplicates(['endpoint','chrom','pos','ref','alt'])
for label in labels:
    datasets.append(align(FG[FG.endpoint==label],'FinnGen R13',label))

def fit(d,R=None,eiv=False):
    x=d.x.to_numpy(); sx=d.sx.to_numpy();y=d.y.to_numpy();sy=d.sy.to_numpy()
    if R is None:R=np.eye(len(x))
    cov=sy[:,None]*R*sy[None,:];covx=sx[:,None]*R*sx[None,:]
    A=np.linalg.inv(cov);den=float(x@A@x)
    if eiv:den-=float(np.trace(A@covx))
    assert den>0
    theta=float(x@A@y/den)
    vy=float(x@A@cov@A@x/den**2);se=np.sqrt(vy)
    gx=(A@y-2*theta*A@x)/den
    sex=np.sqrt(max(0,float(gx@covx@gx)))
    resid=y-theta*x;q=float(resid@A@resid)
    return {'estimate':theta,'se':se,'lo':theta-1.96*se,'hi':theta+1.96*se,'p':2*stats.norm.sf(abs(theta/se)),
            'reference_error_se':sex,'unknown_overlap_se_bound':se+sex,'overlap_envelope_lo':theta-1.96*(se+sex),
            'overlap_envelope_hi':theta+1.96*(se+sex),'Q':q,'Q_df':len(x)-1,'Q_p':stats.chi2.sf(q,len(x)-1)}

results=[];influence=[];sims=[];contrasts=[]
ld={}
for group in ['EUR','FIN']:
    p=FROZEN_RESULTS/f'LD_{group}.tsv'
    if p.exists():ld[group]=pd.read_csv(p,sep='\t',index_col=0)
for d in datasets:
    if d.empty:continue
    cohort=d.cohort.iloc[0];label=d.label.iloc[0]
    for panel,select in [('21_no_chr6_20_40Mb',d.chrom!=6),('22_historical',np.ones(len(d),dtype=bool))]:
        a=d[select].reset_index(drop=True);n=len(a)
        if n<15:
            EXCLUDED.append({'cohort':cohort,'label':label,'rsid':panel,'reason':'fewer_than_15_loci'});continue
        for method,eiv in [('WLS',False),('EIV',True)]:
            res=fit(a,eiv=eiv)
            results.append({'cohort':cohort,'label':label,'panel':panel,'method':method,'n_loci':n,**res})
        primary=fit(a)
        for block in ['rsid','chrom']:
            vals=[]
            for b in a[block].unique():
                est=fit(a[a[block]!=b])['estimate'];vals.append(est)
                influence.append({'cohort':cohort,'label':label,'panel':panel,'block':block,'deleted':str(b),'estimate':est})
            v=np.array(vals);k=len(v);sj=np.sqrt((k-1)/k*np.sum((v-v.mean())**2));tc=stats.t.ppf(.975,k-1)
            results.append({'cohort':cohort,'label':label,'panel':panel,'method':'delete_'+block+'_t','n_loci':n,
                            'n_blocks':k,'estimate':primary['estimate'],'se':sj,'lo':primary['estimate']-tc*sj,
                            'hi':primary['estimate']+tc*sj,'p':2*stats.t.sf(abs(primary['estimate']/sj),k-1),
                            'leave_out_min':v.min(),'leave_out_max':v.max()})
        for group,mat in ld.items():
            z=a[a.rsid.isin(mat.index)].reset_index(drop=True)
            if len(z)<15:continue
            R=mat.loc[z.rsid,z.rsid].to_numpy()
            results.append({'cohort':cohort,'label':label,'panel':panel,'method':'GLS_'+group,'n_loci':len(z),**fit(z,R)})
            if len(z)!=n:results.append({'cohort':cohort,'label':label,'panel':panel,'method':'WLS_'+group+'_available','n_loci':len(z),**fit(z)})
            if group=='EUR':
                keep=[]
                for ix in np.argsort(z.igas_p.to_numpy(),kind='stable'):
                    if all(R[ix,j]**2<=.01 for j in keep):keep.append(ix)
                zz=z.iloc[sorted(keep)]
                results.append({'cohort':cohort,'label':label,'panel':panel,'method':'WLS_LD_pruned_EUR','n_loci':len(zz),**fit(zz)})
                if panel=='21_no_chr6_20_40Mb':
                    rng=np.random.default_rng(20261006);x=z.x.to_numpy();sy=z.sy.to_numpy()
                    C=sy[:,None]*R*sy[None,:];draw=rng.multivariate_normal(np.zeros(len(x)),C,size=10000)
                    for meth,A in [('WLS',np.diag(1/sy**2)),('GLS_EUR',np.linalg.inv(C))]:
                        den=x@A@x; se_nominal=1/np.sqrt(den);theta=draw@A@x/den
                        rate=float(np.mean(abs(theta/se_nominal)>stats.norm.ppf(.975)))
                        sims.append({'cohort':cohort,'label':label,'method':meth,'n_loci':len(z),'n_simulations':10000,
                                     'seed':20261006,'conditional_type1_error':rate,'monte_carlo_se':np.sqrt(rate*(1-rate)/10000),
                                     'scenario':'outcome null; observed IGAS effects fixed; EUR LD target errors'})
R=pd.DataFrame(results)
for cohort in R.cohort.unique():
    rows=R[(R.cohort==cohort)&(R.panel=='21_no_chr6_20_40Mb')&(R.method=='WLS')].set_index('label')
    if {'RHEUMA_SERONEG','RHEUMA_SEROPOS_WIDE'}.issubset(rows.index):
        a,b=rows.loc['RHEUMA_SERONEG'],rows.loc['RHEUMA_SEROPOS_WIDE'];diff=a.estimate-b.estimate;se=a.se+b.se
        contrasts.append({'cohort':cohort,'contrast':'M06.0 minus broad M05','difference':diff,'unknown_covariance_se_bound':se,
                          'envelope_lo':diff-1.96*se,'envelope_hi':diff+1.96*se})
for name,df in [('aligned_effects',pd.DataFrame(OBS)),('alignment_results',R),('influence',pd.DataFrame(influence)),
                ('conditional_null_simulation',pd.DataFrame(sims)),('serotype_contrast',pd.DataFrame(contrasts)),('exclusions',pd.DataFrame(EXCLUDED,columns=['cohort','label','rsid','reason']))]:
    df.to_csv(OUT/(name+'.tsv'),sep='\t',index=False,float_format='%.17g')

# Shared-reference HLA variance and worst-case first-order sampling-covariance envelope.
B=pd.read_parquet(EVIDENCE/'g28_atlas_b27.parquet');h=[]
for r in B.itertuples():
    if r.endpoint=='M13_ANKYLOSPON':continue
    f0=r.f_ref_b27;fA=r.f_as_b27;fL=r.f_case_b27;D=fA-f0;ratio=(fL-f0)/D
    grad=np.array([1,-ratio,ratio-1])/D
    sd=np.sqrt([fL*(1-fL)/(2*r.n_cases),fA*(1-fA)/(2*3875),f0*(1-f0)/(2*491870)])
    se=np.sqrt(np.sum((grad*sd)**2));bound=float(np.sum(abs(grad)*sd))
    h.append({'endpoint':r.endpoint,'n_cases':r.n_cases,'estimate':ratio,'historical_se':r.se,'shared_reference_se':se,
              'shared_reference_lo':ratio-1.96*se,'shared_reference_hi':ratio+1.96*se,
              'unknown_covariance_se_bound':bound,'covariance_envelope_lo':ratio-1.96*bound,'covariance_envelope_hi':ratio+1.96*bound})
H=pd.DataFrame(h);H.to_csv(OUT/'hla_overlap_sensitivity.tsv',sep='\t',index=False,float_format='%.17g')
decision=pd.read_parquet(ROOT/'inputs/prov_decision.parquet')
decision=decision[decision.gate.str.match(r'^G5:Q[1-4]$',na=False)][['gate','metric','value','threshold','passed','detail']]
assert len(decision)==4 and not decision.passed.any()
decision.to_csv(OUT/'G5_recorded_criteria.tsv',sep='\t',index=False,float_format='%.17g')
receipt={'revision_plan_sha256':hashlib.sha256((ROOT/'protocol/revision_analysis_plan.md').read_bytes()).hexdigest(),
         'n_aligned_rows':len(OBS),'n_exclusions':len(EXCLUDED),'n_fitted_rows':len(R),'G5_criteria_met':0,
         'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'rss_cap_mib':1024,
         'primary_result':R[(R.cohort=='UK Biobank')&(R.label=='RHEUMA_SERONEG')&(R.panel=='21_no_chr6_20_40Mb')&(R.method=='WLS')].to_dict('records')}
assert receipt['peak_rss_mib']<1024
(OUT/'analysis_receipt.json').write_text(json.dumps(receipt,indent=2,allow_nan=True))
print(R[(R.panel=='21_no_chr6_20_40Mb')&(R.method.isin(['WLS','GLS_EUR','delete_chrom_t']))].to_string(index=False))
print(H[H.endpoint=='RHEUMA_SERONEG'].to_string(index=False));print(json.dumps(receipt,indent=2))
