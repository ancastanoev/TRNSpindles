"""
    Simulate.jl

Unified simulation interface for all model variants.

Wraps DifferentialEquations.jl with sensible defaults and provides
the `SimResult` container for downstream analysis.
"""

"""
    SimResult{F}

Output of any solve function.  Stores the time grid and the solution matrix.

Fields
──────
- `t`  — time grid, shape (Nt,)
- `u`  — solution array, shape (n_state × Nt)
- `model_type` — human-readable string identifying the model variant
- `params` — the parameter struct used (Any, for flexibility)
"""
struct SimResult{F<:AbstractFloat}
    t::Vector{F}
    u::Matrix{F}
    model_type::String
    params::Any
end

# generic solver helper

function _run_ode(
    f!,
    u0,
    tspan,
    params;
    alg=Tsit5(),
    reltol=1e-8,
    abstol=1e-10,
    saveat=nothing,
    maxiters=Int(1e7),
    model_type="unknown",
)
    F = eltype(u0)
    prob = ODEProblem(f!, u0, tspan, params)
    kw = Dict{Symbol,Any}(:reltol => F(reltol), :abstol => F(abstol), :maxiters => maxiters)
    if !isnothing(saveat)
        kw[:saveat] = saveat
    end
    sol = solve(prob, alg; kw...)

    # Convert to dense matrix (n_state × Nt)
    u_mat = hcat(sol.u...)

    return SimResult{F}(sol.t, u_mat, model_type, params)
end

# 4-node Wilson-Cowan (no delay)

"""
    solve_wc(p; u0, tspan, kwargs...) → SimResult

Simulate the 4-node WC CTRC system.

Arguments
─────────
- `p`     — `WCParams`
- `u0`    — initial state [P₀, I₀, T₀, R₀]; default near zero
- `tspan` — (t_start, t_end)

Keyword arguments passed to DifferentialEquations:
- `alg`       default `Tsit5()`
- `reltol`    default `1e-8`
- `abstol`    default `1e-10`
- `saveat`    time points to save (nothing = adaptive)
"""
function solve_wc(
    p::WCParams{F};
    u0::AbstractVector{F}=F[2.0, 1.0, 1.5, 5.0],
    tspan::Tuple{F,F}=(F(0), F(500)),
    alg=Tsit5(),
    reltol::Real=1e-8,
    abstol::Real=1e-10,
    saveat=nothing,
    maxiters::Int=Int(1e7),
) where {F<:AbstractFloat}
    f! = (du, u, params, t) -> wc_rhs!(du, u, params, t)
    return _run_ode(
        f!,
        F.(u0),
        tspan,
        p;
        alg=alg,
        reltol=reltol,
        abstol=abstol,
        saveat=saveat,
        maxiters=maxiters,
        model_type="WilsonCowan-4node",
    )
end

# Delayed Wilson-Cowan

"""
    solve_delayed_wc(p; u0, tspan, init_chain, kwargs...) → SimResult

Simulate the Gamma-distributed delayed WC system.

`u0` is the 4-node IC [P₀, I₀, T₀, R₀]; chain variables are
appended automatically (zero or equilibrium-warm-start).
"""
function solve_delayed_wc(
    p::WCDelayParams{F};
    u0_4::AbstractVector{F}=F[2.0, 1.0, 1.5, 5.0],
    tspan::Tuple{F,F}=(F(0), F(500)),
    init_chain::Symbol=:zero,
    alg=Tsit5(),
    reltol::Real=1e-8,
    abstol::Real=1e-10,
    saveat=nothing,
    maxiters::Int=Int(1e7),
) where {F<:AbstractFloat}
    u0 = make_delay_ic(p, F.(u0_4); init_chain=init_chain)
    f! = (du, u, params, t) -> delayed_wc_rhs!(du, u, params, t)
    return _run_ode(
        f!,
        u0,
        tspan,
        p;
        alg=alg,
        reltol=reltol,
        abstol=abstol,
        saveat=saveat,
        maxiters=maxiters,
        model_type="DelayedWC-Gamma",
    )
end

# Rinzel–Golomb network

"""
    solve_rinzel_golomb(net; u0, tspan, kwargs...) → SimResult

Simulate N coupled Rinzel–Golomb TRN neurons.

- `net` — `RinzelGolombNetwork`
- `u0`  — (3N,) initial condition; if `nothing`, random IC is generated.
"""
function solve_rinzel_golomb(
    net::RinzelGolombNetwork{F};
    u0::Union{Nothing,AbstractVector{F}}=nothing,
    tspan::Tuple{F,F}=(F(0), F(2000)),
    rng::AbstractRNG=Random.GLOBAL_RNG,
    alg=Tsit5(),
    reltol::Real=1e-6,
    abstol::Real=1e-8,
    saveat=nothing,
    maxiters::Int=Int(1e8),
) where {F<:AbstractFloat}
    u0_actual = isnothing(u0) ? make_rg_ic(net; rng=rng) : F.(u0)
    f! = (du, u, params, t) -> rinzel_golomb_rhs!(du, u, params, t)
    return _run_ode(
        f!,
        u0_actual,
        tspan,
        net;
        alg=alg,
        reltol=reltol,
        abstol=abstol,
        saveat=saveat,
        maxiters=maxiters,
        model_type="RinzelGolomb-N=$(net.N)",
    )
end

# Convenience: parameter sweeps

"""
    parameter_sweep_wc(p_base, key, values; kwargs...) → Vector{SimResult}

Run `solve_wc` for each value in `values`, varying the single parameter `key`.
Returns a vector of `SimResult`s.

Example
───────
```julia
results = parameter_sweep_wc(p, :wRT, 0.6:0.1:1.0)
```
"""
function parameter_sweep_wc(
    p_base::WCParams{F}, key::Symbol, values::AbstractVector; kwargs...
) where {F<:AbstractFloat}
    map(values) do v
        p_new = _setprop(p_base, key, F(v))
        return solve_wc(p_new; kwargs...)
    end
end

# Helper: reconstruct a WCParams with one field changed (since it's immutable)
function _setprop(p::WCParams{F}, key::Symbol, val::F) where {F}
    kwargs = Dict{Symbol,Any}()
    for fname in fieldnames(WCParams)
        kwargs[fname] = (fname == key) ? val : getfield(p, fname)
    end
    return WCParams{F}(; kwargs...)
end
