"""
    Bridge.jl

Open-loop WC→RG bridge.

Computes the afferent spindle drive from a Wilson–Cowan operating point,
delivers it to a Rinzel–Golomb TRN population, calibrates the coupling
strength, and sweeps (g_syn, g_el) to characterise gap-junction clustering
under that drive.

Coupling is strictly one-way (open-loop): R_eff is a read-out; it is
never fed back into the WC equations.
"""

"""
    WCDriveResult{F}

Output of `extract_wc_drive`.  Holds the uniformly-sampled, normalised
afferent drive D̂(t) and all metadata the driven RG population needs.

Fields
──────
- `t`       uniform WC time grid (model units)                         (Nt,)
- `D`       raw drive  D = wRP·P + wRT·T + B + IR                     (Nt,)
- `D̂`       normalised  D̂ = D / αE, clamped to [0,1]                  (Nt,)
- `R`       R_WC(t) — TRN population rate from WC model                (Nt,)
- `dt_wc`   uniform step of the WC grid (enables O(1) interpolation)
- `P_WC`    WC limit-cycle period (model time units) — FFT estimate
- `P_WC_pp` WC limit-cycle period — peak-to-peak estimate
- `κ`       P_phys / P_WC  (maps RG time → WC time)
- `α_E`     normalisation constant (= p.αE from WCParams)
- `mean_R`  ⟨R_WC⟩ time-mean
"""
struct WCDriveResult{F<:AbstractFloat}
    t::Vector{F}
    D::Vector{F}
    D̂::Vector{F}
    R::Vector{F}
    dt_wc::F
    P_WC::F
    P_WC_pp::F
    κ::F
    α_E::F
    mean_R::F
    P_t::Vector{F}   # P(t) — cortical pyramidal  (for drive decomposition)
    T_t::Vector{F}   # T(t) — thalamic relay       (for drive decomposition)
end

"""
    extract_wc_drive(p; u0, tspan, P_phys, transient_frac, dt_save) → WCDriveResult

Run the 4-node WC model and extract the afferent TRN drive time series.

The drive is the total excitatory input to the R equation (from WilsonCowan.jl):
    D(t) = wRP·P(t) + wRT·T(t) + B + IR
The self-inhibition term −wRR·R is excluded (it is internal to R).

D̂(t) = D(t) / αE is the normalised drive, clamped to [0, 1].

κ = P_phys / P_WC maps RG time to WC time.  P_phys is caller-supplied
(default 3.0 s); P_WC is measured from the simulation via FFT and
independently via peak-to-peak.  A mismatch > 20 % triggers a warning.

Parameters
──────────
- `P_phys`         physical spindle period (s); sets the time-scale κ
- `dt_save`        uniform WC output step (default 0.01 model units);
                   controls D̂ lookup resolution
- `transient_frac` fraction of tspan discarded as transient (default 0.5)
"""
function extract_wc_drive(
    p::WCParams{F};
    u0::AbstractVector{F}=F[1.0, 5.0, 2.0, 5.0],
    tspan::Tuple{F,F}=(F(0), F(2000)),
    P_phys::F=F(3.0),
    transient_frac::F=F(0.5),
    dt_save::F=F(0.01),
) where {F<:AbstractFloat}
    res = solve_wc(p; u0=u0, tspan=tspan, saveat=dt_save)

    n_skip = round(Int, transient_frac * length(res.t))
    t = res.t[(n_skip + 1):end]
    Pv = res.u[1, (n_skip + 1):end]
    Tv = res.u[3, (n_skip + 1):end]
    Rv = res.u[4, (n_skip + 1):end]

    # amplitude guard — catches saturated fixed points before any period measurement
    let R_amp = maximum(Rv) - minimum(Rv)
        R_amp > one(F) || error(
            "extract_wc_drive: R_WC amplitude $(round(R_amp; digits=4)) ≤ 1.0 — " *
            "WC is at a fixed point, not a limit cycle. " *
            "Use an oscillating IC, e.g. u0 = F[1.0, 5.0, 2.0, 5.0].",
        )
    end

    # drive = excitatory part of uR (WilsonCowan.jl line: uR = wRP*P + wRT*T - wRR*R + B + IR)
    D = p.wRP .* Pv .+ p.wRT .* Tv .+ p.B .+ p.IR
    α_E = p.αE
    D̂ = clamp.(D ./ α_E, zero(F), one(F))
    dt_wc = t[2] - t[1]

    # P_WC via FFT
    f_fft = dominant_frequency(Rv, dt_wc)
    P_WC = f_fft > eps(F) ? one(F) / f_fft : F(Inf)

    # P_WC via peak-to-peak
    peaks = extract_peaks(Rv)
    P_WC_pp = length(peaks) >= 2 ? mean_cycle_period(t, peaks) : F(Inf)

    # consistency check
    if isfinite(P_WC) && isfinite(P_WC_pp)
        rel = abs(P_WC - P_WC_pp) / P_WC
        if rel > 0.2
            @warn "extract_wc_drive: FFT and peak-to-peak periods disagree " *
                "by $(round(100 * rel; digits=1)) % — operating point may " *
                "not have a clean limit cycle. Check WC parameters."
        end
    end

    !isfinite(P_WC) && error(
        "extract_wc_drive: WC model produced no oscillation " *
        "(dominant_frequency returned 0).  Use an oscillating " *
        "WCParams operating point.",
    )

    κ = P_phys / P_WC
    mean_R = mean(Rv)

    return WCDriveResult{F}(t, D, D̂, Rv, dt_wc, P_WC, P_WC_pp, κ, α_E, mean_R, Pv, Tv)
end

"""
    DrivenRGParams{F, IM}

Packages a `RinzelGolombNetwork` with the WC drive signal and coupling
strength g_drive.

The extra current added to each neuron i is (from the spec):
    ΔdVᵢ/dt = −g_drive · (Vᵢ − E_exc) · D̂(t_RG / κ) / C

with E_exc = 0 (excitatory reversal potential) so the term simplifies to:
    ΔdVᵢ/dt = −g_drive · Vᵢ · D̂(t_RG / κ) / C

Setting g_drive = 0 recovers the undriven `rinzel_golomb_rhs!` exactly.
"""
struct DrivenRGParams{F<:AbstractFloat,IM<:AbstractMatrix{F}}
    net::RinzelGolombNetwork{F,IM}
    drive::WCDriveResult{F}
    g_drive::F
    E_exc::F      # = 0  (excitatory reversal potential, mV)
end

function DrivenRGParams(
    net::RinzelGolombNetwork{F,IM}, drive::WCDriveResult{F}; g_drive::F=zero(F)
) where {F,IM}
    return DrivenRGParams{F,IM}(net, drive, g_drive, zero(F))
end

"""
    _interp_D̂(drive, t_rg) → F

Look up D̂ at RG time t_rg by mapping to WC time t_wc = t_rg / κ with
periodic wrapping over the stored window.  O(1) on the uniform WC grid.
"""
@inline function _interp_D̂(drive::WCDriveResult{F}, t_rg::Real) where {F<:AbstractFloat}
    t_wc_raw = F(t_rg) / drive.κ
    t0 = drive.t[1]
    t_dur = drive.t[end] - t0
    t_dur < eps(F) && return zero(F)
    t_wc = t0 + mod(t_wc_raw, t_dur)
    idx_f = (t_wc - t0) / drive.dt_wc + one(F)
    idx_lo = max(1, min(length(drive.D̂) - 1, floor(Int, idx_f)))
    α = idx_f - F(idx_lo)
    return (one(F) - α) * drive.D̂[idx_lo] + α * drive.D̂[idx_lo + 1]
end

"""Same as `_interp_D̂` but for the R_WC channel (used in Stage C residual)."""
@inline function _interp_R_WC(drive::WCDriveResult{F}, t_rg::Real) where {F<:AbstractFloat}
    t_wc_raw = F(t_rg) / drive.κ
    t0 = drive.t[1]
    t_dur = drive.t[end] - t0
    t_dur < eps(F) && return zero(F)
    t_wc = t0 + mod(t_wc_raw, t_dur)
    idx_f = (t_wc - t0) / drive.dt_wc + one(F)
    idx_lo = max(1, min(length(drive.R) - 1, floor(Int, idx_f)))
    α = idx_f - F(idx_lo)
    return (one(F) - α) * drive.R[idx_lo] + α * drive.R[idx_lo + 1]
end

"""
    driven_rg_rhs!(dx, x, dp::DrivenRGParams, t)

In-place ODE for the WC-driven Rinzel–Golomb population.

Calls the unmodified `rinzel_golomb_rhs!` first (Stage B wrapping, not forking),
then adds the drive current to the N voltage equations.

State layout: x = [V₁…Vₙ | h₁…hₙ | s₁…sₙ]
"""
function driven_rg_rhs!(
    dx::AbstractVector{F}, x::AbstractVector{F}, dp::DrivenRGParams{F}, t
) where {F<:AbstractFloat}
    # 1. undriven base RHS (modifies dx in-place)
    rinzel_golomb_rhs!(dx, x, dp.net, t)

    # 2. add drive current; short-circuit if g_drive == 0 for exact equivalence
    iszero(dp.g_drive) && return nothing

    N = dp.net.N
    D̂_t = _interp_D̂(dp.drive, t)
    C = dp.net.param.C
    g_drv = dp.g_drive
    E_exc = dp.E_exc    # = 0

    @inbounds for i in 1:N
        dx[i] -= g_drv * (x[i] - E_exc) * D̂_t / C
    end
    return nothing
end

"""
    solve_driven_rg(dp; u0, tspan, rng, kwargs...) → SimResult

Simulate the WC-driven Rinzel–Golomb population.
Wraps `_run_ode` from Simulate.jl.
"""
function solve_driven_rg(
    dp::DrivenRGParams{F};
    u0::Union{Nothing,AbstractVector{F}}=nothing,
    tspan::Tuple{F,F}=(F(0), F(3000)),
    rng::AbstractRNG=Random.GLOBAL_RNG,
    alg=Tsit5(),
    reltol::Real=1e-6,
    abstol::Real=1e-8,
    saveat=1.0,
    maxiters::Int=Int(1e8),
) where {F<:AbstractFloat}
    u0_actual = isnothing(u0) ? make_rg_ic(dp.net; rng=rng) : F.(u0)
    f! = (du, u, params, t) -> driven_rg_rhs!(du, u, params, t)
    return _run_ode(
        f!,
        u0_actual,
        tspan,
        dp;
        alg=alg,
        reltol=reltol,
        abstol=abstol,
        saveat=saveat,
        maxiters=maxiters,
        model_type="DrivenRG-N=$(dp.net.N)",
    )
end

"""
    _mean_sbar(res, N) → (mean_s̄, modulation_depth)

Compute the time-mean and peak-to-trough of the population-mean synaptic
variable s̄(t) from a DrivenRG simulation.
State layout: s = x[2N+1 : 3N].
"""
function _mean_sbar(res::SimResult{F}, N::Int) where {F<:AbstractFloat}
    s_mat = res.u[(2N + 1):3N, :]
    s_bar = vec(mean(s_mat; dims=1))
    mean_s = mean(s_bar)
    mod_dep = maximum(s_bar) - minimum(s_bar)
    return mean_s, mod_dep
end

"""
    _rg_with_params(net, g_syn, g_el) → RinzelGolombNetwork

Return a copy of `net` with `g_syn` and `g_el` overridden.  All other
parameters (including the adjacency matrix and gating) are preserved.
"""
function _rg_with_params(
    net::RinzelGolombNetwork{F}, g_syn::F, g_el::F
) where {F<:AbstractFloat}
    p_new = RinzelGolombParams{F}(;
        C=net.param.C,
        g_Ca=net.param.g_Ca,
        g_L=net.param.g_L,
        g_syn=g_syn,
        g_el=g_el,
        V_Ca=net.param.V_Ca,
        V_L=net.param.V_L,
        V_syn=net.param.V_syn,
        kr=net.param.kr,
        kf=net.param.kf,
        phi=net.param.phi,
        V_threshold=net.param.V_threshold,
    )
    return RinzelGolombNetwork(p_new, net.gating, net.A)
end

# Calibration & sweep driver

"""
    find_g_drive(net, drive; ...) → F

Coarse geometric sweep over `g_drive` ∈ g_range to find a value where:
  • mean(s̄) ∈ [mean_s_lo, mean_s_hi]   (not silent / saturated)
  • peak-to-trough(s̄) > ε_mod           (drive modulates firing)

Prints a line for each candidate.  Errors if no suitable value is found.
"""
function find_g_drive(
    net::RinzelGolombNetwork{F},
    drive::WCDriveResult{F};
    g_range::Tuple{F,F}=(F(1e-3), F(5.0)),
    n_coarse::Int=20,
    ε_mod::F=F(0.02),
    mean_s_lo::F=F(0.05),
    mean_s_hi::F=F(0.90),
    n_cycles::Int=6,
    rng::AbstractRNG=Random.MersenneTwister(0),
    kwargs...,
) where {F<:AbstractFloat}
    T_drive = drive.P_WC * drive.κ
    tspan = (F(0), F(n_cycles * T_drive))

    g_vals = exp.(range(log(g_range[1]), log(g_range[2]); length=n_coarse))
    found = nothing

    for g in g_vals
        dp = DrivenRGParams(net, drive; g_drive=F(g))
        u0 = make_rg_ic(net; rng=rng)
        res = solve_driven_rg(dp; u0=u0, tspan=tspan, kwargs...)
        μ_s, mod_s = _mean_sbar(res, net.N)

        if μ_s < mean_s_lo
            @printf("  g_drive=%.5f → SILENT      mean s̄=%.4f  mod=%.4f\n", g, μ_s, mod_s)
            continue
        elseif μ_s > mean_s_hi
            @printf("  g_drive=%.5f → SATURATED   mean s̄=%.4f  mod=%.4f\n", g, μ_s, mod_s)
            break
        elseif mod_s < ε_mod
            @printf(
                "  g_drive=%.5f → UNMODULATED mean s̄=%.4f  mod=%.4f < ε_mod=%.4f\n",
                g,
                μ_s,
                mod_s,
                ε_mod
            )
            continue
        else
            @printf("  g_drive=%.5f ✓  mean s̄=%.4f  peak-to-trough=%.4f\n", g, μ_s, mod_s)
            found = F(g)
            break
        end
    end

    isnothing(found) && error(
        "find_g_drive: no valid g_drive found in range $g_range. " *
        "Check that the WC drive oscillates and the RG parameters are compatible.",
    )
    return found
end

"""
    CalibrationResult{F}

Output of the Stage C calibration.

Fields
──────
- `g_drive`   calibrated drive coupling
- `g_el_star` synchronized anchor (minimum g_el where χ ≥ χ_threshold for all seeds)
- `a`         FROZEN amplitude scale ⟨R_WC⟩ / ⟨s̄⟩_sync — computed once, never refit
- `residual`  self-consistency ‖R_eff − R_WC‖_rms at the anchor (diagnostic only)
"""
struct CalibrationResult{F<:AbstractFloat}
    g_drive::F
    g_el_star::F
    a::F
    residual::F
end

"""
    calibrate(net_base, drive; ...) → CalibrationResult

Full Stage C calibration:
  C.1  find_g_drive
  C.2  sweep g_el → find synchronized anchor g_el*
  C.2  compute a = ⟨R_WC⟩ / ⟨s̄⟩_synchronized  (FROZEN — never refit per grid cell)
  C.3  report self-consistency residual ‖R_eff − R_WC‖_rms (diagnostic, not minimised)
"""
function calibrate(
    net_base::RinzelGolombNetwork{F},
    drive::WCDriveResult{F};
    g_drive_range::Tuple{F,F}=(F(1e-3), F(5.0)),
    anchor_gel_range::Tuple{F,F}=(F(0.0), F(1.0)),
    n_gel_steps::Int=20,
    n_seeds::Int=5,
    χ_threshold::F=F(0.95),
    n_cycles::Int=8,
    rng_seed::Int=0,
    kwargs...,
) where {F<:AbstractFloat}
    T_drive = drive.P_WC * drive.κ
    tspan = (F(0), F(n_cycles * T_drive))

    println("\n── Stage C.1: searching g_drive ──────────────────────────────")
    g_drive = find_g_drive(
        net_base,
        drive;
        g_range=g_drive_range,
        n_cycles=n_cycles,
        rng=Random.MersenneTwister(rng_seed),
        kwargs...,
    )
    @printf("  → g_drive = %.5f\n\n", g_drive)

    println("── Stage C.2: sweeping g_el for synchronized anchor ──────────")
    g_el_vals = range(anchor_gel_range[1], anchor_gel_range[2]; length=n_gel_steps)
    g_el_star = F(anchor_gel_range[2])  # fallback: use maximum

    for g_el in g_el_vals
        net_ij = _rg_with_params(net_base, net_base.param.g_syn, F(g_el))
        dp_ij = DrivenRGParams(net_ij, drive; g_drive=g_drive)

        chis = F[]
        for seed in 1:n_seeds
            u0 = make_rg_ic(net_ij; rng=Random.MersenneTwister(seed))
            res = solve_driven_rg(dp_ij; u0=u0, tspan=tspan, kwargs...)
            n_skip = round(Int, 0.5 * size(res.u, 2))
            V_win = res.u[1:net_ij.N, (n_skip + 1):end]
            push!(chis, pfeuty_chi(V_win))
        end

        min_χ = minimum(chis)
        @printf(
            "  g_el=%.3f  min_χ=%.4f  (seeds: %s)\n",
            g_el,
            min_χ,
            join([@sprintf("%.3f", x) for x in chis], ", ")
        )

        if min_χ >= χ_threshold
            g_el_star = F(g_el)
            @printf(
                "  → anchor g_el* = %.4f  (min χ = %.4f ≥ threshold %.2f)\n\n",
                g_el_star,
                min_χ,
                χ_threshold
            )
            break
        end
    end

    println("── Stage C.2: computing amplitude scale a at anchor ──────────")
    net_anc = _rg_with_params(net_base, net_base.param.g_syn, g_el_star)
    dp_anc = DrivenRGParams(net_anc, drive; g_drive=g_drive)

    sbar_vals = F[]
    for seed in 1:n_seeds
        u0 = make_rg_ic(net_anc; rng=Random.MersenneTwister(seed))
        res = solve_driven_rg(dp_anc; u0=u0, tspan=tspan, kwargs...)
        μ_s, _ = _mean_sbar(res, net_anc.N)
        push!(sbar_vals, μ_s)
    end
    mean_sbar = mean(sbar_vals)
    a = drive.mean_R / mean_sbar
    @printf("  ⟨R_WC⟩ = %.4f,  ⟨s̄⟩_sync = %.4f,  a = %.4f\n\n", drive.mean_R, mean_sbar, a)

    println("── Stage C.3: self-consistency residual ──────────────────────")
    u0_sc = make_rg_ic(net_anc; rng=Random.MersenneTwister(0))
    res_sc = solve_driven_rg(dp_anc; u0=u0_sc, tspan=tspan, kwargs...)
    N = net_anc.N
    n_skip = round(Int, 0.5 * size(res_sc.u, 2))

    s_bar_t = vec(mean(res_sc.u[(2N + 1):3N, (n_skip + 1):end]; dims=1))
    R_eff_t = a .* s_bar_t
    t_window = res_sc.t[(n_skip + 1):end]
    R_wc_aligned = [_interp_R_WC(drive, t) for t in t_window]

    residual = sqrt(mean((R_eff_t .- R_wc_aligned) .^ 2))
    @printf("  ‖R_eff − R_WC‖_rms = %.4f  (diagnostic only; loop is open)\n\n", residual)

    return CalibrationResult{F}(g_drive, g_el_star, a, residual)
end

"""
    compute_R_eff(res, N, a) → Vector{F}

R_eff(t) = a · s̄(t)  where  s̄(t) = (1/N) Σᵢ sᵢ(t).

`a` MUST be the frozen `CalibrationResult.a` — do not pass a refit value.
"""
function compute_R_eff(res::SimResult{F}, N::Int, a::F) where {F<:AbstractFloat}
    s_mat = res.u[(2N + 1):3N, :]
    s_bar = vec(mean(s_mat; dims=1))
    return a .* s_bar
end

"""
    sweep_driven(net_base, drive, calib, g_syns, g_els, A; ...) → Matrix{Vector{ClusterResult}}

Sweep the (g_syn, g_el) plane for the WC-driven RG population.
Reuses `ClusterResult` / `detect_clusters` from ClusterDetect.jl.

The amplitude scale `calib.a` is FROZEN — the same value is applied at every
grid cell.  Refitting `a` to match R_WC at each cell is a bug; the function
asserts this cannot happen.

Returns the same `Matrix{Vector{ClusterResult}}` layout as `sweep_gsyn_gel`,
so `cluster_count_fraction` and plotting helpers work unchanged.
"""
function sweep_driven(
    net_base::RinzelGolombNetwork{F},
    drive::WCDriveResult{F},
    calib::CalibrationResult{F},
    g_syns::AbstractVector{F},
    g_els::AbstractVector{F},
    A::AbstractMatrix{F};
    n_runs::Int=5,
    n_cycles::Int=8,
    show_progress::Bool=true,
    rng::AbstractRNG=Random.MersenneTwister(0),
    kwargs...,
) where {F<:AbstractFloat}

    # a is frozen — assert it is a plain scalar (immutable CalibrationResult)
    a = calib.a
    @assert a > zero(F) "calib.a must be positive"
    # Matching R_WC everywhere by refitting a per grid cell would be a bug:
    # the purpose of the sweep is to vary (g_syn, g_el) under a fixed scale.

    T_drive = drive.P_WC * drive.κ
    tspan = (F(0), F(n_cycles * T_drive))

    ng = length(g_syns)
    ne = length(g_els)
    results = Matrix{Vector{ClusterResult{F}}}(undef, ng, ne)

    base_seed = rand(rng, UInt64)
    prog = show_progress ? Progress(ng * ne; desc="Driven sweep: ") : nothing

    Threads.@threads for k in 1:(ng * ne)
        i = (k - 1) % ng + 1
        j = (k - 1) ÷ ng + 1

        net_ij = _rg_with_params(net_base, g_syns[i], g_els[j])
        # re-attach the chosen adjacency
        net_ij = RinzelGolombNetwork(net_ij.param, net_ij.gating, A)
        dp_ij = DrivenRGParams(net_ij, drive; g_drive=calib.g_drive)

        local_rng = Random.MersenneTwister(
            base_seed ⊻ (UInt64(i) * 0x9e3779b97f4a7c15) ⊻ UInt64(j)
        )

        run_results = ClusterResult{F}[]
        for _ in 1:n_runs
            u0 = make_rg_ic(net_ij; rng=local_rng)
            res = solve_driven_rg(dp_ij; u0=u0, tspan=tspan, kwargs...)
            cr = detect_clusters(res, net_ij)
            push!(run_results, cr)
        end
        results[i, j] = run_results
        isnothing(prog) || next!(prog)
    end

    return results
end
