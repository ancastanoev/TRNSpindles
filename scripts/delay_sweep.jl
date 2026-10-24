"""
delay_sweep.jl
─────────────────
Sweep the Gamma mean-delay ρ for the distributed-delay WC system.

For each ρ value the steady-state TRN oscillation is characterised
(Rmin, Rmax, dominant frequency, spindle-proximity score) and saved.

Outputs
───────
  results/data/delay_sweep.jld2
  results/figures/delay_spindle_conditions.png
  results/figures/delay_freq_profile.png

Usage
─────
    julia --project scripts/delay_sweep.jl
"""

using TRNSpindles
using CairoMakie, JLD2, Printf, ProgressMeter

# grid
ρ_vals = vcat(
    collect(0.02:0.02:0.2),    # fine resolution in spindle range
    collect(0.3:0.1:1.0),
    collect(1.5:0.5:5.0),
    collect(6.0:1.0:12.0),
)
n_rho = length(ρ_vals)

# fixed coupling parameters
wc_p = WCParams{Float64}(; wPP=1.1, wRT=0.5, wTR=0.3, wRR=0.2)
m = 2         # Gamma shape
ga = 0.0       # no active feedback (pure delay effect)
T_end = 500.0

# preallocate results
freq_fft_arr = zeros(n_rho)
Rmin_arr = zeros(n_rho)
Rmax_arr = zeros(n_rho)
score_arr = zeros(n_rho)
regime_arr = Vector{SpindleRegime}(undef, n_rho)

@info "Delay sweep" n_rho=n_rho

prog = Progress(n_rho; desc="Delay sweep: ")
for (k, ρ) in enumerate(ρ_vals)
    dp = WCDelayParams{Float64}(; wc=wc_p, m=m, ρ=ρ, ga=ga)
    res = solve_delayed_wc(
        dp;
        u0_4=[2.0, 1.5, 1.5, 5.0],
        tspan=(0.0, T_end),
        init_chain=:equilibrium,
        reltol=1e-9,
        abstol=1e-11,
    )

    # Compute dt from the saved time grid
    dt = length(res.t) > 1 ? (res.t[end] - res.t[1]) / (length(res.t) - 1) : 0.01
    n_skip = round(Int, 0.5 * length(res.t))
    R_ss = res.u[4, n_skip:end]   # TRN steady-state

    Rmin_arr[k] = minimum(R_ss)
    Rmax_arr[k] = maximum(R_ss)
    freq_fft_arr[k] = dominant_frequency(R_ss, dt)
    score_arr[k] = spindle_proximity_score(maximum(res.u[3, n_skip:end]), Rmax_arr[k])
    regime_arr[k] = classify_regime(freq_fft_arr[k])

    next!(prog)
end

# save
mkpath("results/data")
jldsave(
    "results/data/delay_sweep.jld2";
    rho_vals=ρ_vals,
    freq_fft=freq_fft_arr,
    Rmin=Rmin_arr,
    Rmax=Rmax_arr,
    score=score_arr,
    regime=string.(regime_arr),
)

# print regime table
@info "Regime summary:"
for (ρ, f, reg) in zip(ρ_vals, freq_fft_arr, regime_arr)
    @printf("  ρ = %5.2f  f = %7.4f Hz  → %s\n", ρ, f, reg)
end

# figure: joint spindle-conditions (ρ, R)
mkpath("results/figures")
plt_sc = plot_spindle_conditions(
    ρ_vals,
    Rmin_arr,
    Rmax_arr;
    title="Joint Spindle Conditions in the (ρ, R) Plane",
    size=(800, 500),
)
save("results/figures/delay_spindle_conditions.png", plt_sc)

# figure: frequency vs delay
fig_freq = Figure(; size=(800, 500))
ax_freq = Axis(
    fig_freq[1, 1];
    xlabel="Mean delay ρ",
    ylabel="Dominant frequency (Hz)",
    title="TRN oscillation frequency vs delay",
    yscale=log10,
)
lines!(ax_freq, ρ_vals, freq_fft_arr; linewidth=2, label="FFT frequency")
scatter!(ax_freq, ρ_vals, freq_fft_arr; markersize=4)
hlines!(ax_freq, [7.0, 15.0]; linestyle=:dash, color=:blue, label="Spindle band (7–15 Hz)")
hlines!(ax_freq, [1.0]; linestyle=:dash, color=:green, label="1 Hz boundary")
hlines!(ax_freq, [0.1]; linestyle=:dash, color=:red, label="0.1 Hz boundary")
axislegend(ax_freq; position=:rt)
save("results/figures/delay_freq_profile.png", fig_freq)

@info "Done." files=[
    "results/data/delay_sweep.jld2",
    "results/figures/delay_spindle_conditions.png",
    "results/figures/delay_freq_profile.png",
]
