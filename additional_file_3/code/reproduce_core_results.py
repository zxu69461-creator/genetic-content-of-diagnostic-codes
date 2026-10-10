"""Offline numerical checks for the accompanying aggregate-data archive."""
from pathlib import Path
import hashlib,json,resource
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'inputs'
def read(t):return pd.read_parquet(D/(t+'.parquet'))

def verify_files():
 manifest=json.loads((ROOT/'SHA256.json').read_text())
 for name,expected in manifest.items():
  path=ROOT/name
  assert path.is_file(),f'Missing file: {name}'
  assert hashlib.sha256(path.read_bytes()).hexdigest()==expected,f'Checksum mismatch: {name}'
 return len(manifest)

checks={'verified_files':verify_files()}
M=read('ref_finngen_manifest').set_index('phenocode')
REF='RHEUMA_OTHER_WIDE';AS='M13_ANKYLOSPON'
a=read('g28_atlas_b27')
f0=a.f_ref_b27.to_numpy();fa=a.f_as_b27.to_numpy();f=a.f_case_b27.to_numpy()
n=a.n_cases.to_numpy();n0=float(M.loc[REF,'num_controls']);na=float(M.loc[AS,'num_cases'])
r=(f-f0)/(fa-f0)
se=np.sqrt((f*(1-f)/(2*n)+f0*(1-f0)/(2*n0)+r*r*fa*(1-fa)/(2*na))/(fa-f0)**2)
assert np.allclose(r,a.as_equivalent,rtol=1e-11,atol=1e-13)
assert np.allclose(se,a.se,rtol=1e-11,atol=1e-13)
z=r/se;p=2*norm.sf(abs(z))
counts={'endpoints':len(a),'positive_nominal':int(sum(z>3.29)),'negative_nominal':int(sum(z< -3.29)),
 'positive_bonferroni':int(sum((r>0)&(p*len(a)<.05))),'changed_at_exact_nominal_threshold':int(sum((z>3.29)!=(z>norm.isf(.0005))))}
assert counts=={'endpoints':1677,'positive_nominal':131,'negative_nominal':24,'positive_bonferroni':90,'changed_at_exact_nominal_threshold':0},counts
checks['atlas']=counts
j=a.index[a.endpoint=='RHEUMA_SERONEG'][0]
checks['HLA_seronegative_RA']={'estimate':float(r[j]),'se':float(se[j]),'lo':float(r[j]-1.96*se[j]),'hi':float(r[j]+1.96*se[j])}

E=read('g30_g4_fetch');L=read('g27_igas_loci');L=L[L.locus_set=='genome-wide']
gc=read('g0_commonref_lambda').set_index('endpoint').lambda_commonref
key=['chrom','pos','ref','alt']
ref=E[E.endpoint==REF][key+['af_alt_controls']].rename(columns={'af_alt_controls':'f0'})
x=E.merge(ref,on=key);v=x.f0*(1-x.f0);n=x.endpoint.map(M.num_cases)
g=x.endpoint.map(lambda e:max(1,round(gc.get(e,1.25),3)))
x['b']=(x.af_alt_cases-x.f0)/v
x['s']=np.sqrt(x.af_alt_cases*(1-x.af_alt_cases)/(2*n)+v/(2*M.loc[REF,'num_controls']))/v*np.sqrt(g)
loci=L.merge(E[key].drop_duplicates(),on=['chrom','pos'])
same=(loci.ea==loci.alt)&(loci.oa==loci.ref);flip=(loci.ea==loci.ref)&(loci.oa==loci.alt)
loci=loci[same|flip].copy();loci['b_as']=np.where(loci.ea==loci.alt,loci.beta,-loci.beta)
loci=loci.drop_duplicates(key).sort_values(['chrom','pos']).set_index(key)
B=x.pivot(index=key,columns='endpoint',values='b').reindex(loci.index)
S=x.pivot(index=key,columns='endpoint',values='s').reindex(loci.index)
bx=loci.b_as.to_numpy();sx=loci.se.to_numpy()
def slope(y,sy,m):
 w=1/sy[m]**2;den=np.sum(w*(bx[m]**2-sx[m]**2))
 if not np.isfinite(den) or den<=0:raise ValueError('Unstable errors-in-variables denominator')
 return np.sum(w*bx[m]*y[m])/den

yr=B[AS].to_numpy();sr=S[AS].to_numpy();lp=loci.reset_index()
expected=read('audit_mhc_sensitivity').set_index(['label','analysis'])
recalculated=[]
for label in sorted(expected.index.get_level_values(0).unique()):
 y=B[label].to_numpy();sy=S[label].to_numpy()
 for name,mask in [('original_22_loci',np.ones(len(loci),bool)),('exclude_chr6_20_40Mb',~((lp.chrom==6)&lp.pos.between(20000000,40000000)).to_numpy()),('exclude_all_chr6',lp.chrom.to_numpy()!=6)]:
  m=mask&np.isfinite(y)&np.isfinite(sy)&np.isfinite(yr)&(sy>0)
  estimate=slope(y,sy,m)/slope(yr,sr,m);reps=[]
  for k in np.flatnonzero(m):
   mm=m.copy();mm[k]=False;reps.append(slope(y,sy,mm)/slope(yr,sr,mm))
  reps=np.asarray(reps);err=np.sqrt((len(reps)-1)/len(reps)*sum((reps-reps.mean())**2))
  target=expected.loc[(label,name)]
  assert int(sum(m))==target.n_loci
  assert np.allclose([estimate,err],[target.estimate,target.se],rtol=1e-10,atol=1e-12),(label,name)
  recalculated.append({'label':label,'analysis':name,'n_loci':int(sum(m)),'estimate':float(estimate),'se':float(err),'lo':float(estimate-1.96*err),'hi':float(estimate+1.96*err)})
checks['AS_locus_results']=recalculated
rep=read('g10_replication');main=rep[rep.prediction!='B2-secondary']
assert len(main)==7 and int(main.met.sum())==3
checks['planned_external_criteria']={'met':3,'not_met':4,'secondary_polyarthrosis':'unstable'}
sim=read('g9_sim_grid');coverage_cols=[c for c in sim.columns if 'coverage' in c]
checks['simulation_coverage']={c:{'min':float(sim[c].min()),'max':float(sim[c].max())} for c in coverage_cols}
checks['peak_rss_mib_macos']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
checks['status']='PASS: selected numerical checks; not a full pipeline rerun'
print(json.dumps(checks,indent=2))
