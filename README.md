# TRNSpindles

**A computational study of sleep-spindle dynamics in the cortico-thalamo-reticular system**, accompanying a BSc Computer Science thesis at the West University of Timișoara.

The project links two levels of description of the same neural system:

- a **mean-field Wilson–Cowan** model of the cortico-thalamo-reticular loop (with an optional distributed delay), studied both by direct simulation and by numerical bifurcation analysis; and
- a **detailed, gap-junction-coupled network** of conductance-based reticular (TRN) neurons, in which synchronisation clusters emerge.

The **novel contribution** is the bridge between them: an open-loop WC→RG drive (`Bridge.jl`) that takes the afferent spindle rhythm from a Wilson–Cowan operating point and uses it to drive the detailed TRN network, then characterises how gap-junction coupling shapes synchronisation clustering under that drive. The coupling is strictly one-way — the network read-out is never fed back into the mean-field equations — which makes the population rhythm and the single-neuron clustering directly comparable.

---

## Repository at a glance

The repository has **two components**, each self-documented:

| Component | What it is | Language / tool | Guide |
| --- | --- | --- | --- |
| **Julia package** (repo root) | The simulation engine: Wilson–Cowan, distributed delay, Rinzel–Golomb network, and the WC→RG bridge, plus analysis and the driven parameter sweep | Julia ≥ 1.10 | [`docs/GUIDE.md`](docs/GUIDE.md) |
| **`matcont/`** | The bifurcation analysis: equilibrium/limit-cycle continuation, Hopf boundaries, and the invariant torus for both model variants | MATLAB + MatCont 7.4 | [`matcont/README.md`](matcont/README.md) |

The Julia package is the computer-science core of the thesis; the MatCont component is the applied-mathematics companion that traces the bifurcation structure organising the simulated dynamics.

---

## Quick start (5 minutes)

```bash
# 1. enter the project and install the pinned environment
cd TRNSpindles
julia --project -e 'using Pkg; Pkg.instantiate()'

# 2. confirm it loads
julia --project -e 'using TRNSpindles; println("OK")'

# 3. run the baseline Wilson–Cowan simulation (writes a figure)
julia --project scripts/wc_baseline.jl

# 4. run the full test suite
julia --project test/runtests.jl
```

`instantiate` installs and precompiles every dependency from the **pinned `Manifest.toml`**, so the environment is reproduced exactly (a few minutes on first run). A `Makefile` wraps the common entry points: `make install`, `make test`, `make baseline`, `make delay`.

---

## What the package does

The Julia source (`src/`) is organised as a layered package:

- **Model layer** — `WilsonCowan.jl` (4-population mean-field model + analytic Jacobian), `DelayedWC.jl` (distributed delay realised as a finite ODE via the linear-chain construction, so standard adaptive solvers apply), and `RinzelGolomb.jl` (conductance-based TRN neuron with pre-allocated, thread-safe network scratch).
- **Bridge** — `Bridge.jl`, the novel core: extracts the spindle drive from a WC limit cycle, calibrates it, drives the TRN network, and sweeps the gap-junction plane to characterise clustering.
- **Analysis** — `Analysis.jl` (frequency, spindle-proximity scoring, regime classification) and `ClusterDetect.jl` (Pfeuty χ synchrony index and correlation-based cluster grouping).
- **Infrastructure** — `Simulate.jl` (a unified solve interface), `GraphUtils.jl` (gap-junction adjacency and Laplacian, dense and sparse), `Parameters.jl` (typed, immutable parameter structs), and `Plotting.jl`.

### Engineering highlights

- **Reproducible by construction** — a pinned `Manifest.toml` fixes every dependency version; the test suite includes a regression test that pins core output against a committed reference.
- **A 155-case automated test suite**, including formal **verification tests**: the linear-chain delay reduction is checked against a direct Gamma convolution, the delayed model with zero feedback is checked against the base model, and the analytic Jacobian is checked against finite differences.
- **Performance-conscious numerics** — `StaticArrays` give an allocation-free right-hand side for the mean-field model; the Rinzel–Golomb network uses in-place RHS evaluation with pre-allocated buffers and a sparse adjacency matrix; the parameter sweep is **multithreaded** with an independent, reproducible RNG per grid cell.
- **Native Julia throughout** — no external C compiler is required; published SciML benchmarks report order-of-magnitude speed-ups for the Julia ODE stack over the SciPy equivalent.

---

## Reproducing the thesis figures

### Julia figures (synchronisation analysis)

The driven parameter sweep (`scripts/sweep_rg_driven.jl`) characterises clustering across an `11 × 16` grid of `(g_syn, g_el)` with 11 seeds (1936 simulations) and writes a timestamped publication file under `results/data/`. The figure generators read that file as a command-line argument:

```bash
DATA=results/data/sweep_driven_publication_.jld2

julia --project scripts/figures/cluster_count_field.jl    "$DATA"
julia --project scripts/figures/cluster_decomposition.jl  "$DATA"
julia --project scripts/figures/synchrony_chi_field.jl    "$DATA"
julia --project scripts/figures/multistability_maps.jl    "$DATA"

# voltage traces and the read-out overlay generate their own data through the bridge:
julia --project scripts/figures/voltage_traces.jl
julia --project scripts/figures/readout_vs_cycle.jl
```

To regenerate the sweep from scratch (expensive — 1936 simulations, several hours):

```bash
julia --project scripts/sweep_rg_driven.jl
```

Generated simulation data and output figures are not tracked (they are regenerable from the scripts above); the figures shown in the thesis body are the canonical reference.

### MatCont figures (bifurcation analysis)

The bifurcation figures (thesis §4.1–4.8) are produced from the continuation curves in `matcont/`, which ships both the figure scripts and the continuation data. See [`matcont/README.md`](matcont/README.md) for the full figure map and reproduction steps.

---

## Project structure

```text
TRNSpindles/
├── Project.toml            package metadata + pinned dependencies
├── Manifest.toml           exact resolved environment (reproducibility)
├── LICENSE                 MIT, with third-party attribution (see below)
├── Makefile                common entry points
├── README.md               this file
├── src/                    the Julia package
│   ├── TRNSpindles.jl         main module
│   ├── Parameters.jl          typed parameter structs
│   ├── Sigmoids.jl            WC sigmoid + derivative
│   ├── WilsonCowan.jl         4-population WC model + Jacobian
│   ├── DelayedWC.jl           distributed delay via the linear-chain construction
│   ├── RinzelGolomb.jl        conductance-based TRN neuron network
│   ├── GraphUtils.jl          gap-junction adjacency / Laplacian
│   ├── Bridge.jl              open-loop WC→RG drive bridge (novel core)
│   ├── Simulate.jl            unified ODE solve interface
│   ├── Analysis.jl            frequency / spindle scoring / regime labelling
│   ├── ClusterDetect.jl       Pfeuty χ index + cluster grouping
│   └── Plotting.jl            CairoMakie figure helpers
├── scripts/
│   ├── cluster_check.jl       pre-flight cluster sanity check
│   ├── wc_baseline.jl         baseline WC simulation
│   ├── delay_sweep.jl         delay (rho) sweep
│   ├── sweep_rg_driven.jl     production driven sweep
│   └── figures/               thesis figure generators (content-named)
├── test/                   automated test + verification suite
│   ├── runtests.jl
│   ├── test_regression.jl
│   ├── test_kernel_bridge.jl
│   ├── test_bridge_units.jl
│   └── golden/                pinned reference fixtures
├── docs/
│   └── GUIDE.md            extended build / usage / extension guide
└── matcont/                MATLAB + MatCont bifurcation analysis (see its README)
```

---

## Provenance and license

TRNSpindles is released under the MIT License (see [`LICENSE`](LICENSE)).

`RinzelGolomb.jl`, `ClusterDetect.jl` and `GraphUtils.jl` are Julia ports of, or mirror the logic in, the MIT-licensed project [alexfikl/2025-trn-gap-junctions-experiments](https://github.com/alexfikl/2025-trn-gap-junctions-experiments) (© 2025 Alexandru Fikl). That project's copyright and MIT notice are reproduced in `LICENSE`, and per-file provenance is noted in each ported source file. All other code — the Wilson–Cowan model, the distributed-delay extension, the WC→RG bridge, the driven sweep, and the analysis — is original to this work.