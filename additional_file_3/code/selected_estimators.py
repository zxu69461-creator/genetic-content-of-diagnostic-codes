"""Selected numerical estimators from the study code.
Coefficients are conditional genetic model summaries, not validated disease fractions.
These functions are supplied for method inspection; the compact reproduction script
checks selected headline results rather than rerunning every analysis.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import optimize, stats

def nearest_psd(A: np.ndarray, eps: float=1e-12) -> np.ndarray:
    A = (A + A.T) / 2
    w, V = np.linalg.eigh(A)
    return V * np.maximum(w, eps) @ V.T

def _moments(B, S, bL, sL, R, rL, w, Om=None, om=None):
    J, K = B.shape
    A = np.zeros((K, K))
    c = np.zeros(K)
    for j in range(J):
        s = S[j]
        Omj = Om[j] if Om is not None else R * np.outer(s, s)
        omj = om[j] if om is not None else rL * s * sL[j]
        A += w[j] * (np.outer(B[j], B[j]) - Omj)
        c += w[j] * (B[j] * bL[j] - omj)
    return (A, c)

def _solve(A, c, convex=True, free=None):
    K = len(c)
    free = np.zeros(K, bool) if free is None else np.asarray(free, bool)
    A = nearest_psd(A)
    if not convex:
        return np.linalg.solve(A, c)
    mix = ~free
    fun = lambda p: p @ A @ p - 2 * c @ p
    jac = lambda p: 2 * A @ p - 2 * c
    cons = [{'type': 'ineq', 'fun': lambda p: 1 - p[mix].sum(), 'jac': lambda p: -mix.astype(float)}]
    x0 = np.where(mix, 0.5 / max(mix.sum(), 1), 0.0)
    bounds = [(None, None) if f else (0, 1) for f in free]
    res = optimize.minimize(fun, x0, jac=jac, bounds=bounds, constraints=cons, method='SLSQP', options={'ftol': 1e-14, 'maxiter': 1000})
    return res.x

def deconvolve(B: np.ndarray, S: np.ndarray, bL: np.ndarray, sL: np.ndarray, R: np.ndarray | None=None, rL: np.ndarray | None=None, convex: bool=True, jackknife: bool=True, Om=None, om=None, free=None, blocks=None) -> dict:
    J, K = B.shape
    R = np.eye(K) if R is None else R
    rL = np.zeros(K) if rL is None else rL
    w = 1.0 / sL ** 2
    A, c = _moments(B, S, bL, sL, R, rL, w, Om, om)
    ev = np.linalg.eigvalsh((A + A.T) / 2)
    cond = float(ev.max() / max(ev.min(), 1e-300)) if ev.max() > 0 else np.inf
    free_m = np.zeros(K, bool) if free is None else np.asarray(free, bool)
    pi = _solve(A, c, convex, free_m)
    out = {'pi': pi, 'sum_pi': float(pi[~free_m].sum()), 'cond': cond, 'ill_conditioned': bool(cond > 1000000.0)}
    if jackknife:
        blk = np.arange(J) if blocks is None else np.asarray(blocks)
        ub = np.unique(blk)
        G = len(ub)
        reps = np.empty((G, K))
        for g, b in enumerate(ub):
            m = blk != b
            Aj, cj = _moments(B[m], S[m], bL[m], sL[m], R, rL, w[m], None if Om is None else Om[m], None if om is None else om[m])
            reps[g] = _solve(Aj, cj, convex, free_m)
        mean = reps.mean(0)
        se = np.sqrt((G - 1) / G * ((reps - mean) ** 2).sum(0))
        rs = reps[:, ~free_m].sum(1)
        ssum = np.sqrt((G - 1) / G * ((rs - rs.mean()) ** 2).sum())
        out.update(se=se, se_sum=float(ssum), reps=reps)
    fit = B @ pi
    if Om is not None:
        var = sL ** 2 + np.einsum('k,jkl,l->j', pi, Om, pi) - 2 * (om @ pi)
    else:
        var = sL ** 2 + np.einsum('jk,kl,jl->j', S * pi, R, S * pi) - 2 * (S * pi) @ rL * sL
    z = (bL - fit) / np.sqrt(np.maximum(var, 1e-30))
    out['resid_z'] = z
    out['Q'] = float((z ** 2).sum())
    out['Q_df'] = J - int((np.abs(pi) > 1e-06).sum())
    out['Q_p'] = float(stats.chi2.sf(out['Q'], max(out['Q_df'], 1)))
    return out

def component_heterogeneity(resid_z: np.ndarray, groups: pd.Series) -> pd.DataFrame:
    rows = []
    for g, idx in groups.groupby(groups).groups.items():
        zz = resid_z[np.asarray([groups.index.get_loc(i) for i in idx])]
        q = float((zz ** 2).sum())
        df = len(zz)
        rows.append({'sentinel_set': g, 'n': df, 'Q': q, 'Q_p': float(stats.chi2.sf(q, df)), 'I2': max(0.0, (q - df) / q) if q > 0 else 0.0})
    return pd.DataFrame(rows)

def null_z_corr(z1: np.ndarray, z2: np.ndarray, p_thresh: float=0.1) -> float:
    p1 = 2 * stats.norm.sf(np.abs(z1))
    p2 = 2 * stats.norm.sf(np.abs(z2))
    m = (p1 > p_thresh) & (p2 > p_thresh) & np.isfinite(z1) & np.isfinite(z2)
    return float(np.corrcoef(z1[m], z2[m])[0, 1])

def simulate(pi_true: np.ndarray, B_true: np.ndarray, S: np.ndarray, sL: np.ndarray, R: np.ndarray, rL: np.ndarray, delta_sd: float=0.0, rng=None):
    rng = np.random.default_rng(rng)
    J, K = B_true.shape
    Sig = np.zeros((J, K + 1, K + 1))
    for j in range(J):
        s = np.r_[S[j], sL[j]]
        C = np.eye(K + 1)
        C[:K, :K] = R
        C[:K, K] = rL
        C[K, :K] = rL
        Sig[j] = C * np.outer(s, s)
    E = np.stack([rng.multivariate_normal(np.zeros(K + 1), Sig[j]) for j in range(J)])
    B = B_true + E[:, :K]
    bL = B_true @ pi_true + rng.normal(0, delta_sd, J) + E[:, K]
    return (B, bL)

def winners_curse(beta: np.ndarray, se: np.ndarray, z_thresh: float=5.45) -> np.ndarray:
    out = np.array(beta, float).copy()
    z = np.asarray(beta) / np.asarray(se)
    for i, zi in enumerate(z):
        if not np.isfinite(zi) or abs(zi) <= z_thresh:
            continue

        def nll(mu):
            den = stats.norm.sf(z_thresh - mu) + stats.norm.cdf(-z_thresh - mu)
            return -(stats.norm.logpdf(zi - mu) - np.log(max(den, 1e-300)))
        r = optimize.minimize_scalar(nll, bounds=(0 if zi > 0 else zi, zi if zi > 0 else 0), method='bounded')
        out[i] = r.x * se[i]
    return out

def common_reference(eff: 'pd.DataFrame', n_cases: dict, ref_endpoint: str, n_ref: int, gc: dict | None=None, key=('chrom', 'pos', 'ref', 'alt')):
    import pandas as pd
    key = list(key)
    ref = eff[eff.endpoint == ref_endpoint][key + ['af_alt_controls']].rename(columns={'af_alt_controls': 'f_ref'})
    x = eff.merge(ref, on=key, how='inner')
    x = x[x.af_alt_cases.notna() & x.endpoint.isin(n_cases.keys())].copy()
    v = x.f_ref * (1 - x.f_ref)
    n = x.endpoint.map(n_cases).astype(float)
    g = x.endpoint.map(lambda e: (gc or {}).get(e, 1.0)).astype(float)
    x['beta'] = (x.af_alt_cases - x.f_ref) / v
    x['se'] = np.sqrt(x.af_alt_cases * (1 - x.af_alt_cases) / (2 * n) + v / (2 * n_ref)) / v * np.sqrt(g)
    return x[['endpoint'] + key + ['beta', 'se']]
