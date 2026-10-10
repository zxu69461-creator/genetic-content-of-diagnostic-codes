"""Independent numerical and provenance checks for the revised targeted results."""
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from scipy import stats
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'recomputed' if (ROOT/'recomputed').exists() else ROOT/'results'
A=pd.read_csv(R/'aligned_effects.tsv',sep='\t',float_precision='round_trip')
F=pd.read_csv(R/'alignment_results.tsv',sep='\t',float_precision='round_trip')
checks=[]
for (cohort,label),d in A.groupby(['cohort','label']):
    d=d[d.chrom!=6];assert len(d)==21
    # Whitened least squares, independent of the analysis script's matrix inverse implementation.
    x=(d.x/d.sy).to_numpy()[:,None];y=(d.y/d.sy).to_numpy()
    b=np.linalg.lstsq(x,y,rcond=None)[0][0];se=1/np.linalg.norm(x)
    f=F[(F.cohort==cohort)&(F.label==label)&(F.panel=='21_no_chr6_20_40Mb')&(F.method=='WLS')].iloc[0]
    assert np.allclose([b,se],[f.estimate,f.se],atol=1e-12,rtol=1e-12)
    assert np.isclose(2*stats.norm.sf(abs(b/se)),f.p,rtol=1e-12,atol=1e-15)
    # Coefficient and its conditional error transform together under a common source rescaling.
    bscaled=np.linalg.lstsq(x*4,y,rcond=None)[0][0];assert np.isclose(bscaled,b/4)
    checks.append({'cohort':cohort,'label':label,'WLS_recomputed':'PASS','source_scale_invariance':'PASS'})
for g in ['EUR','FIN']:
    m=pd.read_csv(ROOT/'results'/f'LD_{g}.tsv',sep='\t',index_col=0).to_numpy()
    assert np.allclose(m,m.T) and np.allclose(np.diag(m),1) and np.linalg.eigvalsh(m).min()>0
source_path=ROOT/'inputs/IGAS_original_panel.tsv'
if not source_path.exists():source_path=ROOT/'inputs/IGAS_original_sumstats.tsv.gz'
if source_path.exists():
    raw=pd.read_csv(source_path,sep='\t',low_memory=False)
    frozen=pd.read_csv(ROOT/'protocol/frozen_locus_panel.tsv',sep='\t')
    raw=raw[raw.rsid.isin(frozen.rsid)]
    candidates=frozen.merge(raw,on='rsid',suffixes=('_frozen','_source'),validate='one_to_many')
    # The source has two rs1801274 records. Audit which record reproduces the already-frozen
    # harmonized values; do not change the frozen panel or select using target outcomes.
    match=np.isclose(abs(candidates.beta_frozen),abs(candidates.beta_source),atol=1e-14,rtol=1e-12)&np.isclose(candidates.se_frozen,candidates.se_source,atol=1e-14,rtol=1e-12)&np.isclose(candidates.p_frozen,candidates.p_source,rtol=1e-10,atol=0)
    candidates['matches_frozen_values']=match
    candidates.to_csv(R/'IGAS_source_candidate_records.tsv',sep='\t',index=False,float_format='%.17g')
    z=candidates[match].copy()
    assert len(z)==22
    assert z.rsid.nunique()==22
    assert np.allclose(abs(z.beta_frozen),abs(z.beta_source),atol=1e-14,rtol=1e-12)
    assert np.allclose(z.se_frozen,z.se_source,atol=1e-14,rtol=1e-12)
    assert np.allclose(z.p_frozen,z.p_source,rtol=1e-10,atol=0)
    assert np.allclose(z.beta_source,np.log(z.OR),atol=1e-12,rtol=1e-10)
    z.to_csv(R/'IGAS_source_numeric_audit.tsv',sep='\t',index=False,float_format='%.17g')
    source_check='All 22 frozen loci match original source effect magnitudes, standard errors and P values. There are 23 original records for these rsIDs, including an alternate rs1801274 record; all candidates are retained and the already-frozen matching record is identified. Source beta equals log(released OR); no clinical scale equivalence assumed.'
else:
    source_check='Original source not bundled; see frozen source audit and hashes.'
criteria=pd.read_csv(R/'G5_recorded_criteria.tsv',sep='\t');assert len(criteria)==4 and not criteria.passed.any()
assert F.shape[0]==112 and len(A)==176
receipt={'status':'PASS','independent_WLS_checks':checks,'IGAS_source_check':source_check,'G5_criteria_met':0,
         'effect_rows':len(A),'result_rows':len(F),'source_effect_unit':'released IGAS beta, not a patient fraction'}
if R.name=='recomputed':
    compared=[]
    for name in ['aligned_effects','alignment_results','influence','conditional_null_simulation','serotype_contrast','hla_overlap_sensitivity','G5_recorded_criteria','exclusions']:
        a=pd.read_csv(ROOT/'results'/(name+'.tsv'),sep='\t',float_precision='round_trip')
        b=pd.read_csv(R/(name+'.tsv'),sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-10,atol=1e-12)
        compared.append(name)
    receipt['recomputed_tables_matching_frozen']=compared
(R/'verification_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
