# Agent gender and thematic roles in infants: pupillometry data and analysis code

This repository contains the pupillometry data of Experiments 2 and 3 of:

> Papeo, L., Vettori, S., Serraille, E., Odin, C., Rostami, F., & Hochmann, J.-R. (2024). Abstract thematic roles in infants' representation of social events. *Current Biology, 34*, 1–7. https://doi.org/10.1016/j.cub.2024.07.081

It also contains the code for a follow-up analysis. That analysis asks whether 7-month-old infants' pupil response to a reversal of agent and patient roles depends on whether the agent is a boy or a girl. Every infant and every trial is labelled with the gender of the character in the agent role.

## Repository structure

```
data/
  subjects.csv          one row per infant (80 infants)
  pupil_trials.csv      one row per trial, with the pupil trace of that trial
code/
  common.py             loading and preprocessing
  analyze_gender.py     statistics  -> results/
  make_figures.py       figures     -> figures/
results/                output of analyze_gender.py
figures/                output of make_figures.py
requirements.txt
```

## Data

The infants and trials are the final samples of the published analyses:

- Experiment 2: 24 infants in the face-to-face condition and 24 in the back-to-back condition.
- Experiment 3: 32 infants.

The pupil values are the final right-eye data of the original study (60 Hz, Tobii T60XL), with gaps already interpolated.

### Experimental design

Each infant saw one boy–girl dyad; one child wore blue and the other orange. The gender of the agent, the agent's colour and the agent's side on the screen were counterbalanced across infants. About 5 of 6 trials were *standard* trials, in which one child was in an agent-like posture and the other in a patient-like posture.

| Experiment | Trial type | What is shown |
|---|---|---|
| 2 | `standard` | Agent-like and patient-like postures. The children face each other (face-to-face condition) or turn away from each other (back-to-back condition). |
| 2 | `deviant` | Same children with the two postures reversed. |
| 3 | `standard` | As in Experiment 2, face-to-face. |
| 3 | `same_role` | Both children in a neutral posture. The standard agent faces the back of the other child. |
| 3 | `role_switch` | Both children in a neutral posture. The standard patient faces the back of the other child. |

### `data/subjects.csv`

| Column | Description |
|---|---|
| `experiment` | `exp2` or `exp3` |
| `condition` | `face-to-face` or `back-to-back` (Exp. 2); `face-to-back` (Exp. 3) |
| `subject_id` | Infant number used in the original file names. Numbers repeat across experiments and conditions, so identify an infant by `experiment` + `condition` + `subject_id`. |
| `agent_gender` | `boy` or `girl`: the child in the agent role in the standard trials |
| `agent_colour` | `blue` or `orange`: clothes of that child. The other child wears the other colour. |
| `agent_side` | `left` or `right`: side of that child on the screen in the standard trials |
| `n_standard`, `n_deviant`, `n_role_switch`, `n_same_role` | Number of trials of each type in the data |

### `data/pupil_trials.csv`

| Column | Description |
|---|---|
| `experiment`, `condition`, `subject_id` | As above |
| `trial_type` | See the table above |
| `trial_n` | Order of the trial among the trials of that type for that infant |
| `agent_in_trial` | `boy` or `girl`: the child in the agent role on this trial. In Exp. 2 this is the child in the agent-like posture. In Exp. 3 deviants it is the child facing the other. |
| `t0000` … `t3000` | Pupil diameter (mm) from image onset to 3000 ms, one column per sample at 60 Hz. The number in the column name is the time in ms. Values are raw diameters, not baseline-corrected. Empty cells are missing samples. |

### Stimulus codes

In the original stimulus logs, the images are named `SI_<item>_<fa|fc>_m<side><posture><colour>.jpg`. The code describes the man (boy):
- `side`: `l` = left, `r` = right
- `posture`: `a` = agent-like, `p` = patient-like
- `colour`: `b` = blue, `o` = orange

For example, `mlab` means the boy is on the left, in the agent-like posture, wearing blue. The girl is on the other side, in the other posture, wearing the other colour. The Exp. 3 deviants are named `SI_25_pp_m<side><colour>_<l|r>.jpg`, where the last letter gives the direction the boy faces. The `agent_gender` and `agent_in_trial` columns were derived from these codes.

## Analysis

1. **Baseline:** for each trial, subtract the mean of the first 100 ms (samples 1–6).
2. **Trial mean:** average the pupil change over 183–1950 ms, the window of the original study.
3. **Infant mean:** average the trials of each type for each infant.
4. **Surprise score:** deviant − standard (Exp. 2), or role-switch − same-role (Exp. 3).
5. **Group comparison:** compare boy-agent and girl-agent infants with Welch's t-test and Hedges' g (bootstrap 95% CI), for each experiment and condition.
6. **Time-course test:** a cluster-based permutation test on the difference waves, over 100–2917 ms with 5,000 permutations.
7. **Pooled test (Exp. 2 face-to-face + Exp. 3):** divide each score by the SD of its experiment, then fit a regression with agent gender and experiment (HC3 standard errors).

### Main results

| Data set | Boy agent, mean (mm) | Girl agent, mean (mm) | Welch's t (df), p | Hedges' g [95% CI] |
|---|---|---|---|---|
| Exp. 2, face-to-face | 0.028 (n = 9) | 0.049 (n = 15) | −0.57 (21.9), .58 | −0.21 [−0.98, 0.52] |
| Exp. 2, back-to-back | 0.036 (n = 14) | −0.029 (n = 10) | 0.83 (15.3), .42 | 0.35 [−0.51, 1.24] |
| Exp. 3 | 0.060 (n = 14) | 0.103 (n = 18) | −0.55 (28.1), .59 | −0.19 [−0.80, 0.61] |
| Pooled (Exp. 2 face-to-face + Exp. 3) | n = 23 | n = 33 | −0.76 (53), .45 | −0.21 SD [−0.73, 0.32] |

No cluster-based test showed a difference between the groups. The same code also reproduces the published effects:
- Exp. 2 face-to-face, deviant > standard: t(23) = 2.10, one-tailed p = .024
- Exp. 3, role-switch > same-role: t(31) = 2.19, one-tailed p = .018

### Running the code

```
pip install -r requirements.txt
python code/analyze_gender.py
python code/make_figures.py
```

`data/*.csv` can also be read in MATLAB with `readtable`.

## Contact

Farzad Rostami, Institut des Sciences Cognitives Marc Jeannerod, CNRS UMR 5229, Lyon.
