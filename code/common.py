"""Shared loading and preprocessing for the gender x thematic-role analyses."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FS = 60                      # Tobii T60XL sampling rate (Hz)
N_SAMPLES = 181              # 0-3000 ms
BASELINE = slice(0, 6)       # first 100 ms (samples 1-6)
WINDOW = slice(11, 118)      # 183-1950 ms (samples 12-118), the window of the original paper
CLUSTER_RANGE = slice(6, 176)  # samples 7-176, range tested in the original permutation script
TIME_MS = np.arange(N_SAMPLES) * 1000 / FS

# surprise contrast per data set: (role-violating trials, reference trials)
CONTRAST = {'exp2': ('deviant', 'standard'), 'exp3': ('role_switch', 'same_role')}


def load():
    """Return subjects table, trial table and baseline-corrected traces (trials x samples)."""
    subjects = pd.read_csv(ROOT / 'data' / 'subjects.csv', dtype={'subject_id': str})
    trials = pd.read_csv(ROOT / 'data' / 'pupil_trials.csv', dtype={'subject_id': str})
    tcols = [c for c in trials.columns if c.startswith('t') and c[1:].isdigit()]
    X = trials[tcols].to_numpy(float)
    X = X - np.nanmean(X[:, BASELINE], axis=1, keepdims=True)
    info = trials[['experiment', 'condition', 'subject_id', 'trial_type', 'trial_n', 'agent_in_trial']].copy()
    return subjects, info, X


def infant_waves(info, X):
    """Mean baseline-corrected trace per infant and trial type: dict[(exp, cond, id, type)] -> trace."""
    waves = {}
    for key, idx in info.groupby(['experiment', 'condition', 'subject_id', 'trial_type']).groups.items():
        waves[key] = np.nanmean(X[np.asarray(idx)], axis=0)
    return waves


def surprise_table(subjects, info, X):
    """One row per infant: window means per trial type, surprise score and surprise wave."""
    info = info.assign(window=np.nanmean(X[:, WINDOW], axis=1))
    w = info.groupby(['experiment', 'condition', 'subject_id', 'trial_type'])['window'].mean().unstack()
    waves = infant_waves(info, X)
    rows = []
    for (exp, cond, sid), r in w.iterrows():
        a, b = CONTRAST[exp]
        rows.append(dict(experiment=exp, condition=cond, subject_id=sid, surprise=r[a] - r[b],
                         wave=waves[(exp, cond, sid, a)] - waves[(exp, cond, sid, b)]))
    out = pd.DataFrame(rows).merge(subjects, on=['experiment', 'condition', 'subject_id'])
    return out
