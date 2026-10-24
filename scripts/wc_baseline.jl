"""
wc_baseline.jl
─────────────────
Baseline 4-node Wilson–Cowan CTRC simulation (no delays).

Reproduces representative time-series and phase-plane behaviour.
Saves:
  results/figures/wc_timeseries.png
  results/figures/wc_phase_plane.png
  results/data/wc_baseline.jld2

Usage
─────
    julia --project scripts/wc_baseline.jl

Optional environment overrides
───────────────────────────────
    WC_wPP=1.1  WC_wTR=0.8  WC_T=500  julia --project scripts/wc_baseline.jl
"""

using TRNSpindles
using CairoMakie, JLD2, Printf

# resolve runtime overrides
wPP = parse(Float64, get(ENV, "WC_wPP", "1.0"))
wRT = parse(Float64, get(ENV, "WC_wTR", "0.8"))
wTR = parse(Float64, get(ENV, "WC_wRT", "0.4"))
wRR = parse(Float64, get(ENV, "WC_wRR", "0.2"))
T_end = parse(Float64, get(ENV, "WC_T", "500.0"))

@info "Running baseline WC" wPP wRT wTR wRR T_end

# build parameters
p = WCParams{Float64}(; wPP=wPP, wRT=wRT, wTR=wTR, wRR=wRR)

# initial conditions (near literature rest values)
u0 = [2.0, 1.5, 1.5, 5.0]   # [P, I, T, R]

# simulate
res = solve_wc(p; u0=u0, tspan=(0.0, T_end), reltol=1e-9, abstol=1e-11)

# frequency analysis (TRN channel, index 4)
dt = res.t[2] - res.t[1]
stats = analyze_timeseries(res, 4, dt)
@info "TRN frequency analysis" stats.freq_fft stats.regime stats.amplitude

# spindle proximity of final cycle
n_skip = round(Int, 0.5 * length(res.t))
Tmax = maximum(res.u[3, n_skip:end])
Rmax = maximum(res.u[4, n_skip:end])
score = spindle_proximity_score(Tmax, Rmax)
@info @sprintf("Spindle proximity score: %.3f", score)

# save data
mkpath("results/data")
jldsave(
    "results/data/wc_baseline.jld2";
    t=res.t,
    u=res.u,
    p_wPP=wPP,
    p_wTR=wRT,
    freq_fft=stats.freq_fft,
    score=score,
)

# figures
mkpath("results/figures")

plt_ts = plot_timeseries(
    res;
    title=@sprintf(
        "WC baseline  wPP=%.2f  wRT=%.2f  freq≈%.2f Hz", wPP, wRT, stats.freq_fft
    )
)
save("results/figures/wc_timeseries.png", plt_ts)

# Phase-plane (T, R)
n_skip = round(Int, 0.5 * length(res.t))
T_sig = res.u[3, n_skip:end]
R_sig = res.u[4, n_skip:end]
fig_pp = Figure(; size=(700, 500))
ax_pp = Axis(
    fig_pp[1, 1]; xlabel="Relay T [a.u.]", ylabel="TRN R [a.u.]", title="Phase plane (T, R)"
)
lines!(ax_pp, T_sig, R_sig; linewidth=1.5, color=:blue)
save("results/figures/wc_phase_plane.png", fig_pp)

@info "Done." files=[
    "results/data/wc_baseline.jld2",
    "results/figures/wc_timeseries.png",
    "results/figures/wc_phase_plane.png",
]
