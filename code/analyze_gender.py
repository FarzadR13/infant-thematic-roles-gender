"""Does the gender of the agent modulate infants' pupil response to role switches?

Data: Papeo et al. (2024), Current Biology, Experiments 2 and 3 (final right-eye pupil data).
Run from the repository root:  python code/analyze_gender.py
Writes results/stats.csv and results/subject_scores.csv and prints a summary.
"""
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from common import ROOT, CLUSTER_RANGE, TIME_MS, load, surprise_table

RNG = np.random.default_rng(20261008)
N_BOOT, N_PERM, CLUSTER_P = 10_000, 5_000, 0.10


def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return (np.mean(a) - np.mean(b)) / sp * (1 - 3 / (4 * (na + nb) - 9))


def boot_ci(a, b):
    g = [hedges_g(RNG.choice(a, len(a)), RNG.choice(b, len(b))) for _ in range(N_BOOT)]
    return np.percentile(g, [2.5, 97.5])


def welch_t(A, B):
    """Welch t and two-sided p at every sample (A, B: infants x samples)."""
    t, p = stats.ttest_ind(A, B, axis=0, equal_var=False, nan_policy='omit')
    return np.asarray(t), np.asarray(p)


def max_cluster_mass(t, p):
    best, cur, sign = 0.0, 0.0, 0
    for ti, pi in zip(t, p):
        s = int(np.sign(ti)) if pi < CLUSTER_P else 0
        cur = cur + abs(ti) if (s != 0 and s == sign) else (abs(ti) if s != 0 else 0.0)
        sign = s
        best = max(best, cur)
    return best


def clusters(t, p):
    """List of (start, end, mass) of supra-threshold clusters."""
    out, start = [], None
    for i in range(len(t) + 1):
        s = int(np.sign(t[i])) if i < len(t) and p[i] < CLUSTER_P else 0
        if start is not None and (s == 0 or s != np.sign(t[start])):
            out.append((start, i - 1, float(np.sum(np.abs(t[start:i])))))
            start = None
        if s != 0 and start is None:
            start = i
    return out


def cluster_test(waves, is_boy, strata=None):
    """Between-group cluster-mass permutation test on surprise waves (infants x samples)."""
    W = waves[:, CLUSTER_RANGE]
    t, p = welch_t(W[is_boy], W[~is_boy])
    obs = clusters(t, p)
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        lab = is_boy.copy()
        if strata is None:
            lab = RNG.permutation(lab)
        else:
            for s in np.unique(strata):
                m = strata == s
                lab[m] = RNG.permutation(lab[m])
        null[k] = max_cluster_mass(*welch_t(W[lab], W[~lab]))
    t0 = TIME_MS[CLUSTER_RANGE]
    res = [(t0[a], t0[b], mass, (np.sum(null >= mass) + 1) / (N_PERM + 1)) for a, b, mass in obs]
    return sorted(res, key=lambda r: -r[2])


def describe_clusters(res):
    if not res:
        return 'no supra-threshold samples'
    a, b, mass, p = res[0]
    return f'largest cluster {a:.0f}-{b:.0f} ms, p = {p:.3f}'


def main():
    subjects, info, X = load()
    S = surprise_table(subjects, info, X)
    S.drop(columns='wave').to_csv(ROOT / 'results' / 'subject_scores.csv', index=False)
    rows = []

    for (exp, cond), d in S.groupby(['experiment', 'condition'], sort=False):
        t0, p0 = stats.ttest_1samp(d.surprise, 0)
        boy, girl = d.surprise[d.agent_gender == 'boy'].to_numpy(), d.surprise[d.agent_gender == 'girl'].to_numpy()
        t, p = stats.ttest_ind(boy, girl, equal_var=False)
        df = (boy.var(ddof=1) / len(boy) + girl.var(ddof=1) / len(girl)) ** 2 / (
            (boy.var(ddof=1) / len(boy)) ** 2 / (len(boy) - 1) + (girl.var(ddof=1) / len(girl)) ** 2 / (len(girl) - 1))
        lo, hi = boot_ci(boy, girl)
        W = np.vstack(d.wave.to_numpy())
        cl = cluster_test(W, (d.agent_gender == 'boy').to_numpy())
        tb, pb = stats.ttest_1samp(boy, 0)
        tg, pg = stats.ttest_1samp(girl, 0)
        rows.append(dict(data=f'{exp} {cond}', n_boy=len(boy), n_girl=len(girl),
                         all_mean=d.surprise.mean(), all_t=t0, all_df=len(d) - 1, all_p_one_tailed=p0 / 2 if t0 > 0 else 1 - p0 / 2,
                         boy_mean=boy.mean(), boy_sd=boy.std(ddof=1), boy_t=tb, boy_p=pb,
                         girl_mean=girl.mean(), girl_sd=girl.std(ddof=1), girl_t=tg, girl_p=pg,
                         welch_t=t, welch_df=df, welch_p=p, g=hedges_g(boy, girl), g_lo=lo, g_hi=hi,
                         cluster=describe_clusters(cl)))

    # pooled: Exp2 face-to-face + Exp3; each score divided by the SD of its experiment
    P = S[(S.condition == 'face-to-face') | (S.experiment == 'exp3')].copy()
    P['z'] = P.surprise / P.groupby('experiment').surprise.transform('std')
    P['boy'] = (P.agent_gender == 'boy').astype(int)
    fit = smf.ols('z ~ boy + C(experiment)', data=P).fit(cov_type='HC3')
    ci = fit.conf_int().loc['boy']
    inter = smf.ols('z ~ boy * C(experiment)', data=P).fit(cov_type='HC3')
    W = np.vstack(P.wave.to_numpy())
    cl = cluster_test(W, P.boy.to_numpy().astype(bool), strata=P.experiment.to_numpy())
    rows.append(dict(data='pooled exp2 face-to-face + exp3', n_boy=int(P.boy.sum()), n_girl=int((1 - P.boy).sum()),
                     welch_t=fit.tvalues['boy'], welch_df=fit.df_resid, welch_p=fit.pvalues['boy'],
                     g=fit.params['boy'], g_lo=ci[0], g_hi=ci[1],
                     interaction_p=inter.pvalues['boy:C(experiment)[T.exp3]'], cluster=describe_clusters(cl)))

    R = pd.DataFrame(rows)
    R.to_csv(ROOT / 'results' / 'stats.csv', index=False)
    pd.set_option('display.width', 200, 'display.max_columns', 40)
    for r in R.itertuples():
        print(f'\n{r.data}  (boy-agent n = {r.n_boy}, girl-agent n = {r.n_girl})')
        if not r.data.startswith('pooled'):
            print(f'  all infants: surprise = {r.all_mean:.3f} mm, t({r.all_df}) = {r.all_t:.2f}, one-tailed p = {r.all_p_one_tailed:.3f}')
            print(f'  boy-agent : {r.boy_mean:.3f} mm (SD {r.boy_sd:.3f}), vs 0: t = {r.boy_t:.2f}, p = {r.boy_p:.3f}')
            print(f'  girl-agent: {r.girl_mean:.3f} mm (SD {r.girl_sd:.3f}), vs 0: t = {r.girl_t:.2f}, p = {r.girl_p:.3f}')
            print(f'  boy vs girl: Welch t({r.welch_df:.1f}) = {r.welch_t:.2f}, p = {r.welch_p:.3f}, '
                  f"Hedges' g = {r.g:.2f} [{r.g_lo:.2f}, {r.g_hi:.2f}]")
        else:
            print(f'  boy - girl = {r.g:.2f} SD [{r.g_lo:.2f}, {r.g_hi:.2f}], t({r.welch_df:.0f}) = {r.welch_t:.2f}, '
                  f'p = {r.welch_p:.3f}; gender x experiment p = {r.interaction_p:.3f}')
        print(f'  time course: {r.cluster}')


if __name__ == '__main__':
    main()
