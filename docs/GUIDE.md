# TRNSpindles — Full Project Build Guide

A professional Julia implementation of the corticothalamic Wilson–Cowan model
with distributed delays and gap-junction-coupled TRN sectors for studying
NREM sleep spindles.

---

## 1. Scientific context

| File | Model variant | Paper reference |
|------|--------------|----------------|
| `WilsonCowan.jl` | 4-node WC CTRC (no delay) | — |
| `DelayedWC.jl` | + weak Gamma delay (linear-chain trick) | Macdonald (1978) |
| `RinzelGolomb.jl` | Conductance-based TRN neurons | Golomb & Rinzel (1994); trnlib.py |
| `Bridge.jl` | Open-loop WC→RG drive bridge (novel core) | this work |

The gap-junction topology and cluster-detection logic are ported directly from
[alexfikl/2025-trn-gap-junctions-experiments](https://github.com/alexfikl/2025-trn-gap-junctions-experiments)
into idiomatic Julia.

---

## 2. Requirements

| Software | Version |
|----------|---------|
| Julia | ≥ 1.10 |
| Git   | any    |

No C compiler is required (unlike the Python reference, which uses JiTCODE);
DifferentialEquations.jl handles all numerics natively in Julia.

---

## 3. Installation

```bash
# Clone / enter the project directory
cd TRNSpindles

# Activate and instantiate the environment
julia --project -e "using Pkg; Pkg.instantiate()"
```

This downloads and precompiles all registered dependencies (~2–5 minutes on
first run).

### Optional: pre-compile system image (for faster startup)

```bash
julia --project -e "
using PackageCompiler
create_sysimage([:DifferentialEquations, :OrdinaryDiffEq, :CairoMakie],
                sysimage_path=\"TRNSpindles.so\",
                precompile_execution_file=\"scripts/wc_baseline.jl\")
"
# Then run scripts with:
julia --project --sysimage TRNSpindles.so scripts/wc_baseline.jl
```

---

## 4. Running the scripts

All scripts read configuration from environment variables with sensible defaults.

### 4.1 Baseline 4-node WC simulation

```bash
julia --project scripts/wc_baseline.jl

# Override key parameters:
WC_wPP=1.1 WC_wTR=0.8 WC_T=1000 julia --project scripts/wc_baseline.jl
```

**Outputs:**
- `results/data/wc_baseline.jld2`
- `results/figures/wc_timeseries.png`
- `results/figures/wc_phase_plane.png`

---

### 4.2 Delay parameter sweep (ρ sweep)

```bash
julia --project scripts/delay_sweep.jl
```

**Outputs:**
- `results/data/delay_sweep.jld2`
- `results/figures/delay_spindle_conditions.png`
- `results/figures/delay_freq_profile.png`

**Regime labelling:**

| ρ range | Frequency | Regime label |
|---------|-----------|--------------|
| small   | 7–15 Hz   | N2_SPINDLE   |
| medium  | 1–7 Hz    | TRANSITION   |
| large   | 0.1–1 Hz  | N3_SLOW_WAVE |
| very large | < 0.1 Hz | INFRASLOW |

---

### 4.3 WC-driven Rinzel–Golomb cluster sweep

The production experiment drives the Rinzel–Golomb TRN network with the
afferent spindle drive computed from a Wilson–Cowan operating point (Bridge.jl)
and sweeps the gap-junction plane (g_syn, g_el) to characterise synchronisation
clustering.

```bash
# Production driven sweep (11 g_syn × 16 g_el × 11 seeds = 1936 simulations)
julia --project scripts/sweep_rg_driven.jl
```

This writes a timestamped publication file under `results/data/`
(`sweep_driven_publication_*.jld2`), which the figure generators
`scripts/figures/cluster_count_field.jl`, `cluster_decomposition.jl`, `multistability_maps.jl`
and `synchrony_chi_field.jl` consume (see the top-level README.md for the exact
figure commands).

For a fast pre-flight that checks cluster detection without running the full
sweep:

```bash
julia --project scripts/cluster_check.jl
```

---

## 5. Running the tests

```bash
julia --project test/runtests.jl
```

Or via the package manager:
```julia
using Pkg; Pkg.test()
```

Test suite covers:
- Parameter struct construction
- Sigmoid evaluation and derivatives
- ODE dimension and finiteness
- Jacobian shape
- Graph utilities (Laplacian, adjacency sets)
- Analysis functions (frequency, spindle score, regime classification)
- Cluster detection (Pfeuty χ, fixed-point detection)
- Short integration tests for all model variants
- Bridge open-loop drive (g_drive=0 recovers the undriven RHS exactly)
- Verification suite: linear-chain ≈ Gamma convolution, delayed(ga=0) ≈ base WC,
  analytic vs finite-difference Jacobian, limit-cycle period gate

---

## 6. Project structure

```
TRNSpindles/
├── Project.toml            # package metadata + pinned dependencies
├── Manifest.toml           # exact resolved environment (kept for reproducibility)
├── LICENSE                 # MIT, with third-party attribution
├── README.md               # project front door
├── Makefile                # common entry points
├── src/
│   ├── TRNSpindles.jl      # main module
│   ├── Parameters.jl       # typed parameter structs
│   ├── Sigmoids.jl         # WC sigmoid + derivative
│   ├── WilsonCowan.jl      # 4-node WC ODE + Jacobian
│   ├── DelayedWC.jl        # Gamma-delay via the linear-chain trick
│   ├── RinzelGolomb.jl     # conductance-based TRN neuron
│   ├── GraphUtils.jl       # gap-junction adjacency / Laplacian
│   ├── Bridge.jl           # open-loop WC→RG drive bridge (novel core)
│   ├── Simulate.jl         # unified ODE solve interface
│   ├── Analysis.jl         # frequency / spindle scoring / regime labelling
│   ├── ClusterDetect.jl    # Pfeuty chi index + cluster grouping
│   └── Plotting.jl         # CairoMakie figure helpers
├── scripts/
│   ├── cluster_check.jl # pre-flight cluster sanity check
│   ├── wc_baseline.jl   # baseline WC simulation
│   ├── delay_sweep.jl   # delay (rho) sweep
│   ├── sweep_rg_driven.jl  # production driven sweep (publication data)
│   └── figures/            # thesis figure generators (content-named)
├── test/
│   ├── runtests.jl         # entry point + unit/integration/verification tests
│   ├── test_regression.jl  # WC-core regression vs golden fixture
│   ├── test_kernel_bridge.jl  # Gamma-kernel known-answer + Bridge coverage
│   ├── test_bridge_units.jl   # Bridge read-out unit tests
│   └── golden/             # pinned reference fixtures
├── results/
│   ├── figures/
│   └── data/
└── docs/
    └── GUIDE.md            ← you are here
```

---

## 7. Key design decisions

### Why Julia (not Python/MATLAB)?

- Native type stability eliminates per-call heap allocation in the ODE right-hand
  side, and the solver runs as compiled native code rather than through an interpreter.
- StaticArrays give an allocation-free solve for the 4-node WC system.
- DifferentialEquations.jl provides automatic solver selection, dense output, events.
- No external C compiler is required (unlike JiTCODE in the reference Python code).

### Linear-chain trick vs. true delay DDE

The paper uses weak Gamma-distributed delays (`m=2`, shape parameter).  The
linear-chain trick converts this exactly into `m` auxiliary ODEs, so the full
system remains an ODE (not a DDE).  This avoids history-function storage and
allows standard adaptive Runge-Kutta solvers (Tsit5, Rodas5P).

### Immutable parameter structs

All parameter types are `Base.@kwdef struct` with concrete `{F<:AbstractFloat}`
type parameter.  This allows:
- Zero-overhead field access in hot loops.
- Easy `Float32` simulations for large sweeps.
- Clean differentiation (Zygote/ForwardDiff compatible).

### Gap-junction Laplacian

The electrical coupling term `−g_gap (L R)ᵢ = −g_gap ∑ⱼ Aᵢⱼ (Rᵢ − Rⱼ)` is
computed as a dense matrix-vector product `L * R`.  For large N (>50), replace
with `make_laplacian_sparse()` and use `SparseArrays` for a ~N-fold speedup.

---

## 8. Extending the project

### Add a new coupling topology

```julia
# In GraphUtils.jl, define:
function make_adjacency_ring(N; dtype=Float64)
    A = zeros(dtype, N, N)
    for i in 1:N
        j = mod1(i+1, N)
        A[i,j] = A[j,i] = 1
    end
    return A, [N], [0]
end
```

### Run a 2D bifurcation sweep

```julia
using TRNSpindles

results = Dict{Tuple{Float64,Float64}, SimResult}()
for wPP in 0.8:0.05:1.3, wTR in 0.5:0.05:1.0
    p = WCParams{Float64}(wPP=wPP, wTR=wTR)
    results[(wPP, wTR)] = solve_wc(p; tspan=(0.0, 300.0))
end
```

### Sensitivity analysis (ForwardDiff)

```julia
using ForwardDiff

function loss(params::AbstractVector{T}) where T
    p  = WCParams{T}(wPP=params[1], wTR=params[2])
    res = solve_wc(p; tspan=(T(0), T(200.0)))
    n_skip = size(res.u, 2) ÷ 2
    return -maximum(res.u[4, n_skip:end])  # maximise TRN amplitude
end

∇ = ForwardDiff.gradient(loss, [1.0, 0.8])
```

---

## 9. Biological parameter reference

| Parameter | Symbol | Value | Reference |
|-----------|--------|-------|-----------|
| PYR → PYR | wPP | 0.8–1.2 | Brunel (1999); Harris & Shepherd (2015) |
| TRN → TR  | wTR | 0.6–1.0 | Schellenberger Costa et al. (2016) |
| TR → TRN  | wRT | 0.2–0.4 | Pinault (2004); Schellenberger Costa et al. |
| TRN → TRN | wRR | 0.1–0.3 | Deleuze & Huguenard (2006) |
| Spindle target T | Tmax | 4–6 Hz | Mak-McCully et al. (2017); Pinault (2003) |
| Spindle target R | Rmax | 10–15 Hz | Halassa et al. (2011); Steriade et al. |
| Infraslow period | — | ~50 s  | Lázár et al. (2019); Lecci et al. (2017) |


---

## 10. Citation

If you use this code, please cite:

```bibtex
@article{Radulescu2022,
  author    = {Rădulescu, Anca and Kaslik, Eva and Fikl, Alexandru and
               Anderson, Michael and Norwood, Alex},
  title     = {Gap Junctions and Synchronization Clusters in the Thalamic Reticular Nuclei},
  journal   = {arXiv},
  year      = {2022},
  eprint    = {2209.00384},
}
```
