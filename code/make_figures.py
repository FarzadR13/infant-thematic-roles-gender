"""Figures for the gender x thematic-role analyses. Run from the repository root:
    python code/make_figures.py
Saves PNG (300 dpi) and SVG files in figures/.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT, TIME_MS, WINDOW, load, infant_waves, surprise_table

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'svg.fonttype': 'none'})
BLUE, RED, TEAL, ORANGE = '#3166b8', '#cf4d54', '#17646d', '#c0661f'
OUT = ROOT / 'figures'


def band(ax, W, color, label, lw=1.6):
    W = W[:, :180]  # samples 1-180 (0-2983 ms); the last sample is often missing
    m = np.nanmean(W, axis=0)
    se = np.nanstd(W, axis=0, ddof=1) / np.sqrt(np.sum(~np.isnan(W), axis=0))
    t = TIME_MS[:180]
    ax.fill_between(t, m - se, m + se, color=color, alpha=0.18, lw=0)
    ax.plot(t, m, color=color, lw=lw, label=label)


def frame(ax, zero=False):
    ax.axvspan(TIME_MS[WINDOW.start], TIME_MS[WINDOW.stop - 1], color='#ececec', zorder=0, lw=0)
    if zero:
        ax.axhline(0, color='#777777', lw=0.8)
    ax.set_xlim(0, 3000)
    ax.set_xticks(range(0, 3001, 500))
    ax.set_xlabel('Time from image onset (ms)')


def save(fig, name):
    fig.savefig(OUT / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT / f'{name}.svg', bbox_inches='tight')
    plt.close(fig)


def group_waves(waves, subjects, exp, cond, gender, ttype):
    ids = subjects[(subjects.experiment == exp) & (subjects.condition == cond) & (subjects.agent_gender == gender)].subject_id
    return np.vstack([waves[(exp, cond, i, ttype)] for i in ids]), len(ids)


def main():
    subjects, info, X = load()
    waves = infant_waves(info, X)
    S = surprise_table(subjects, info, X)

    # Figure 1: Experiment 2 time courses
    fig, axes = plt.subplots(2, 2, figsize=(6.6, 4.8), sharey=True, sharex=True)
    for r, cond in enumerate(['face-to-face', 'back-to-back']):
        for c, g in enumerate(['boy', 'girl']):
            ax = axes[r, c]
            Ws, n = group_waves(waves, subjects, 'exp2', cond, g, 'standard')
            Wd, _ = group_waves(waves, subjects, 'exp2', cond, g, 'deviant')
            band(ax, Ws, BLUE, 'Standard'); band(ax, Wd, RED, 'Deviant')
            frame(ax)
            ax.set_title(f'{cond.capitalize()}, {g} agent (n = {n})')
            if c == 0:
                ax.set_ylabel('Pupil change (mm)')
            if r == 0:
                ax.set_xlabel('')
    axes[0, 0].legend(loc='lower left')
    fig.tight_layout()
    save(fig, 'fig1_exp2_time_courses')

    # Figure 2: Experiment 2 difference waves
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), sharey=True)
    for ax, cond in zip(axes, ['face-to-face', 'back-to-back']):
        for g, col in (('boy', TEAL), ('girl', ORANGE)):
            d = S[(S.experiment == 'exp2') & (S.condition == cond) & (S.agent_gender == g)]
            band(ax, np.vstack(d.wave), col, f'{g.capitalize()} agent (n = {len(d)})')
        frame(ax, zero=True)
        ax.set_title(cond.capitalize())
        ax.legend(loc='upper left', fontsize=8)
    axes[0].set_ylabel('Deviant − standard (mm)')
    fig.tight_layout()
    save(fig, 'fig2_exp2_difference_waves')

    # Figure 3: Experiment 3 time courses
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6), sharey=True)
    for ax, g in zip(axes, ['boy', 'girl']):
        Wsr, n = group_waves(waves, subjects, 'exp3', 'face-to-back', g, 'same_role')
        Wrs, _ = group_waves(waves, subjects, 'exp3', 'face-to-back', g, 'role_switch')
        band(ax, Wsr, BLUE, 'Same role'); band(ax, Wrs, RED, 'Role switch')
        frame(ax)
        ax.set_title(f'{g.capitalize()} agent (n = {n})')
    axes[0].set_ylabel('Pupil change (mm)')
    axes[0].legend(loc='upper left')
    fig.tight_layout()
    save(fig, 'fig3_exp3_time_courses')

    # Figure 4: Experiment 3 and pooled difference waves
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
    sets = [(S[S.experiment == 'exp3'], 'Experiment 3', 'Role switch − same role (mm)'),
            (S[(S.condition == 'face-to-face') | (S.experiment == 'exp3')], 'Pooled: Exp. 2 face-to-face + Exp. 3',
             'Surprise response (mm)')]
    for ax, (d, title, ylab) in zip(axes, sets):
        for g, col in (('boy', TEAL), ('girl', ORANGE)):
            dd = d[d.agent_gender == g]
            band(ax, np.vstack(dd.wave), col, f'{g.capitalize()} agent (n = {len(dd)})')
        frame(ax, zero=True)
        ax.set_title(title)
        ax.set_ylabel(ylab)
        ax.legend(loc='upper left', fontsize=8)
    fig.tight_layout()
    save(fig, 'fig4_exp3_pooled_difference_waves')


if __name__ == '__main__':
    main()
