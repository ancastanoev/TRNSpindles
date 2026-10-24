# MatCont Continuation Analysis

This directory contains the **numerical continuation analysis** of the
cortico-thalamo-reticular Wilson–Cowan model, performed in
[MatCont](https://sourceforge.net/projects/matcont/) **7.4**. It is the
applied-mathematics companion to the Julia simulation package in the repository
root: where the Julia code *simulates* the dynamics, this component *traces the
bifurcation structure* — the equilibrium branches, limit-cycle envelopes, Hopf
boundaries, and the torus — that organise that dynamics.

The analysis is split into the two model variants studied in the thesis:

| Folder | Model | State dimension | Thesis chapter |
|--------|-------|-----------------|----------------|
| `WC_TRN/` | Wilson–Cowan CTRC, **no delay** | 4 (`P, I, T, R`) | §4.2 |
| `WC_TRN_delay/` | + distributed delay (first-order relaxation filter) | 8 (adds delayed copies) | §4.3 |

---

## 1. At a glance — figure map

Every figure in the bifurcation chapters is produced by one script here. Each
script lives in its own folder together with the continuation data it plots.

| Thesis figure | Script | Folder | What it shows |
|---------------|--------|--------|---------------|
| **4.1** | `wPP_bifurcation_multiple_wTP.m` | `WC_TRN/wPP_bifurcation_multiple_wTP/` | Equilibria + limit-cycle envelopes vs cortical self-coupling `wPP`, for three drives `wTP` |
| **4.2** | `wTP_bifurcation_multiple_wPP.m` | `WC_TRN/wTP_bifurcation_multiple_wPP/` | Equilibria + cycles vs corticothalamic drive `wTP`, for three `wPP` |
| **4.3** | `wTR_bifurcation_multiple_wPP.m` | `WC_TRN/wTR_bifurcation_multiple_wPP/` | Equilibria + cycles vs reticular inhibition `wTR`, for three `wPP` |
| **4.4** | `invariant_torus.m` | `WC_TRN/invariant_torus/` | Neimark–Sacker invariant torus in the `(P, T, R)` phase space |
| **4.5** | `hopf_boundary_multiple_rho.m` | `WC_TRN_delay/hopf_boundary_multiple_rho/` | Hopf stability boundary in `(wTR, wRT)`, traced for several mean delays `rho` |
| **4.6** | `wTR_delay_bifurcation_multiple_rho.m` | `WC_TRN_delay/wTR_delay_bifurcation_multiple_rho/` | Delayed `wTR` bifurcation with limit cycles, four panels across `rho` |
| **4.7** | `spindle_proximity_panel.m` | `WC_TRN_delay/spindle_proximity_panel/` | Spindle-proximity diagrams in `(rho, R)` for three connectivity settings |
| **4.8** | `spindle_cycles_TR.m` | `WC_TRN_delay/spindle_cycles_TR/` | TRN limit cycles in `(T, R)`, coloured by physiological proximity |

---

## 2. Why the continuation data (`.mat`) is included

Each figure folder ships the `.mat` continuation curves it plots. **These files
are not regenerable build artifacts — they are the primary numerical results of
the analysis**, and they are included deliberately.

A continuation curve is the output of an *interactive, expert-guided* numerical
continuation run inside the MatCont GUI: selecting an initial equilibrium or
cycle, choosing the continuation parameters and direction, detecting bifurcation
points, and steering the continuer through folds and branch switches. This is a
manual scientific procedure — there is no single "run this and the curve
appears" command, and reproducing it requires the MatCont toolbox and the same
guided workflow. The curves therefore represent the substantive numerical work
of the bifurcation analysis, in the same way the figures in the thesis body do.

The curves can be reconstructed from scratch in MatCont — the model system
files in each `system/` folder define the equations MatCont continues — but
doing so means repeating the guided continuation in the GUI (selecting the
initial point, setting parameters and direction, and steering through the
bifurcations), rather than running a single script. The shipped `.mat` files
capture the result of that process so the figures can be regenerated directly.

Including them makes the analysis **immediately inspectable and the figures
regenerable on the spot**: a reader can load any `.mat`, examine the traced
branch, and reproduce the corresponding figure by running its script — without
re-deriving the continuation. The data is compact (≈57 MB total across both
systems; the largest single curve is ~5 MB), so the entire bifurcation analysis
ships with the code it belongs to.

> This is the deliberate counterpart to the Julia side, where the large
> simulation caches are *excluded* precisely because they **are** regenerable
> from a single script. Here the opposite holds: the curves are the result, so
> they stay.

---

## 3. Directory layout

```text
matcont/
├── WC_TRN/                                  no-delay model (4 state variables)
│   ├── system/
│   │   ├── Sigmo.m                          shifted sigmoid activation
│   │   ├── TRN_WC_system.m                  ODE right-hand side
│   │   └── TRN_WC_solver.m                  driver: integrate + plot time series
│   ├── wPP_bifurcation_multiple_wTP/        Fig 4.1  (script + .mat curves)
│   ├── wTP_bifurcation_multiple_wPP/        Fig 4.2  (script + .mat curves)
│   ├── wTR_bifurcation_multiple_wPP/        Fig 4.3  (script + .mat curves)
│   └── invariant_torus/                     Fig 4.4  (script + torus .mat)
│
└── WC_TRN_delay/                            delay model (8 state variables)
    ├── system/
    │   ├── Sigmo.m                          shifted sigmoid activation
    │   ├── TRN_WC_system_delays.m           delayed ODE right-hand side
    │   └── TRN_WC_solver_delays.m           driver: integrate + plot time series
    ├── hopf_boundary_multiple_rho/          Fig 4.5  (script + .mat curves)
    ├── wTR_delay_bifurcation_multiple_rho/  Fig 4.6  (script + H1_/EP1_ .mat)
    ├── spindle_proximity_panel/             Fig 4.7  (script + rho_LC_ .mat)
    └── spindle_cycles_TR/                   Fig 4.8  (script + cycle .mat)
```

Each figure folder is self-contained: the script and the `.mat` files it loads
live together, so MATLAB resolves every `load` from the script's own folder with
no path setup.

---

## 4. Reproducing a figure

Each figure regenerates by running its script from inside its own folder:

```matlab
cd matcont/WC_TRN_delay/hopf_boundary_multiple_rho
hopf_boundary_multiple_rho
```

The figure scripts plot pre-computed continuation curves and require only base
MATLAB — **MatCont does not need to be on the path** to reproduce the figures.

The two time-series drivers (`system/TRN_WC_solver.m` and
`system/TRN_WC_solver_delays.m`) integrate the model directly with `ode45` and
also need only base MATLAB.

> **If `ode45` errors on `odeget` / `MinStep`:** MatCont ships its own
> `odeget.m` / `odeset.m`, and if MatCont is on the MATLAB path they shadow the
> built-in versions and break `ode45`. Clear the path before running:
>
> ```matlab
> restoredefaultpath
> rehash toolboxcache
> ```
>
> Re-add MatCont only when re-running a continuation in the GUI.

---

## 5. Parameter notation (code symbol → thesis symbol)

The model files use the thesis notation. The sigmoid response is
`Sigmo(alpha, theta, beta, u)` for amplitude `alpha`, threshold `theta`, slope
`beta`, and net input `u`. The mapping to the thesis symbols is:

| Group | Code variable | Thesis symbol | Role |
|-------|---------------|---------------|------|
| Sigmoid `S_E` (P, T) | `alpha_E`, `theta_E`, `beta_E` | `α_E`, `θ_E`, `β_E` | Excitatory response |
| Sigmoid `S_I` (I) | `alpha_I`, `theta_I`, `beta_I` | `α_I`, `θ_I`, `β_I` | Interneuron response |
| Sigmoid `S_R` (R) | `alpha_R`, `theta_R`, `beta_R` | `α_R`, `θ_R`, `β_R` | Reticular response |
| Coupling | `wPP`, `wPI`, `wPT` | `w_PP`, `w_PI`, `w_PT` | Inputs to `P` |
| Coupling | `wIP`, `wII`, `wIT` | `w_IP`, `w_II`, `w_IT` | Inputs to `I` |
| Coupling | `wTP`, `wTR` | `w_TP`, `w_TR` | Inputs to `T` |
| Coupling | `wRP`, `wRT`, `wRR` | `w_RP`, `w_RT`, `w_RR` | Inputs to `R` |
| Drive | `B` | `B` | Common background drive |
| Drive | `IC`, `IT`, `IR` | `I_C`, `I_T`, `I_R` | Pyramidal / relay / reticular offsets |
| Saturation | `r_P`, `r_I`, `r_T`, `r_R` | `r_P`, `r_I`, `r_T`, `r_R` | Saturation coefficients |
| Timing | `tau` | `τ` | Time constant |
| Timing | `rho` | `ρ` | Mean coupling delay (delay model only) |

Each weight `w_XY` is the directed connection **from** population `Y` **to**
population `X`; the sign of the interaction is carried in the sigmoid argument
(inhibitory terms enter negatively). Baseline values are listed in the thesis
parameter tables (Chapter 3).

---

## 6. Environment

- **MATLAB** R2025a (any recent release should work; the scripts use only base
  MATLAB and `ode45`).
- **MatCont 7.4** — required only to *recompute* the continuation curves, not to
  reproduce the figures from the shipped `.mat` data.