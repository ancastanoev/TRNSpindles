"""
    RinzelGolomb.jl

Conductance-based Rinzel–Golomb TRN neuron model.

This is a direct Julia port of the `model_source` function from
`trnlib.py` (alexfikl/2025-trn-gap-junctions-experiments), which
implements the Golomb–Rinzel (1994) TRN neuron with:
  • Low-threshold T-type Ca²⁺ current (I_Ca)
  • Leak current (I_L)
  • Chemical synaptic current (I_syn, all-to-all, inhibitory)
  • Electrical (gap-junction) current (I_el, via adjacency matrix A)

State vector for N neurons:
  x = [V₁…Vₙ, h₁…hₙ, s₁…sₙ]   (length 3N)

Current equations
─────────────────
  I_Ca,i = g_Ca · m_∞(Vᵢ)³ · hᵢ · (Vᵢ − V_Ca)
  I_L,i  = g_L  · (Vᵢ − V_L)
  I_syn,i = g_syn / (N−1) · (Vᵢ − V_syn) · ∑ⱼ≠ᵢ sⱼ
  I_el,i  = g_el / Mᵢ · ∑ⱼ Aᵢⱼ (Vᵢ − Vⱼ)   [gap junction]

  C · dVᵢ/dt = −(I_Ca,i + I_L,i + I_syn,i + I_el,i)
  dhᵢ/dt     = k_h(Vᵢ) · (h_∞(Vᵢ) − hᵢ)
  dsᵢ/dt     = k_f · s_∞(Vᵢ) · (1 − sᵢ) − k_r · sᵢ

where
  k_h(V) = phi · k_∞.exp(V) / h_∞(V)   (closing rate of h-channel)
  X.exp(V) = exp(−(V − θ_X) / σ_X)     (raw exponential without +1 in denominator)

References
──────────
Golomb & Rinzel (1994), J. Neurophysiol. 71(1):272-280.
Destexhe, Contreras, Steriade, Sejnowski, Huguenard (1994).
trnlib.py, model_source() — alexfikl/2025-trn-gap-junctions-experiments.
"""

# ── low-level gating helpers ─────────────────────────────────────────────────

"""Raw exponential: exp(−(V − θ) / σ)"""
@inline gate_exp(V, θ, σ) = exp(-(V - θ) / σ)

"""Sigmoid gating value: 1 / (1 + exp(−(V−θ)/σ))"""
@inline gate_inf(V, θ, σ) = 1.0 / (1.0 + gate_exp(V, θ, σ))

"""
    gating_m_inf(V, g) → Float

m_∞(V) using parameters from `g::RinzelGolombGating`.
"""
gating_m_inf(V, g::RinzelGolombGating) = gate_inf(V, g.m_inf.θ, g.m_inf.σ)

"""h_∞(V)"""
gating_h_inf(V, g::RinzelGolombGating) = gate_inf(V, g.h_inf.θ, g.h_inf.σ)

"""s_∞(V)"""
gating_s_inf(V, g::RinzelGolombGating) = gate_inf(V, g.s_inf.θ, g.s_inf.σ)

"""
    gating_k_h(V, p, g) → Float

h-channel closing rate: k_h(V) = phi · k_∞.exp(V) / h_∞(V).
"""
function gating_k_h(V, p::RinzelGolombParams, g::RinzelGolombGating)
    kexp = gate_exp(V, g.k_inf.θ, g.k_inf.σ)
    h = gating_h_inf(V, g)
    return p.phi * kexp / h
end

# ── in-place ODE for N coupled neurons ──────────────────────────────────────

"""
    RinzelGolombNetwork

Pre-computed derived quantities for an N-neuron Rinzel–Golomb network.
Constructed once and reused during integration.
"""
struct RinzelGolombNetwork{F<:AbstractFloat,IM<:AbstractMatrix{F}}
    N::Int
    param::RinzelGolombParams{F}
    gating::RinzelGolombGating{F}
    A::IM              # adjacency matrix     (N × N)
    M::Vector{F}       # degree vector        (N,)  = sum(A, dims=2)
    # Pre-allocated scratch — eliminates per-RHS-call allocations.
    # THREAD-SAFETY: each simulation must own its own RinzelGolombNetwork instance.
    # Never share one network object across concurrent ODE solves.
    _m_inf::Vector{F}   # gating m_inf(V)   (N,)
    _h_inf::Vector{F}   # gating h_inf(V)   (N,)
    _s_inf::Vector{F}   # gating s_inf(V)   (N,)
    _kh::Vector{F}   # gating k_h(V)     (N,)
    _AV::Vector{F}   # gap-junction A*V  (N,)
end

function RinzelGolombNetwork(
    param::RinzelGolombParams{F}, gating::RinzelGolombGating{F}, A::AbstractMatrix{F}
) where {F<:AbstractFloat}
    N = size(A, 1)
    @assert size(A) == (N, N) "A must be square"
    M = vec(sum(A; dims=2))
    return RinzelGolombNetwork{F,typeof(A)}(
        N,
        param,
        gating,
        A,
        M,
        Vector{F}(undef, N),
        Vector{F}(undef, N),
        Vector{F}(undef, N),
        Vector{F}(undef, N),
        Vector{F}(undef, N),
    )
end

"""Default all-to-all network with N neurons."""
function RinzelGolombNetwork(
    N::Int; g_syn::AbstractFloat=0.2, g_el::AbstractFloat=0.0, dtype::Type{F}=Float64
) where {F<:AbstractFloat}
    A = ones(F, N, N)
    fill!(view(A, diagind(A)), zero(F))
    p = RinzelGolombParams{F}(; g_syn=F(g_syn), g_el=F(g_el))
    g = RinzelGolombGating{F}()
    return RinzelGolombNetwork(p, g, A)
end

"""
    rinzel_golomb_rhs_ref!(dx, x, net, t)

Reference (allocating) ODE RHS for N coupled Rinzel–Golomb TRN neurons.
Preserved for correctness comparison against `rinzel_golomb_rhs!`.

State layout:  x = [V₁…Vₙ | h₁…hₙ | s₁…sₙ]
"""
function rinzel_golomb_rhs_ref!(
    dx::AbstractVector{F}, x::AbstractVector{F}, net::RinzelGolombNetwork{F}, t
) where {F<:AbstractFloat}
    (; N, param, gating, A, M) = net
    (; C, g_Ca, g_L, g_syn, g_el, V_Ca, V_L, V_syn, kr, kf) = param

    V = view(x, 1:N)
    h = view(x, (N + 1):2N)
    s = view(x, (2N + 1):3N)

    dV = view(dx, 1:N)
    dh = view(dx, (N + 1):2N)
    ds = view(dx, (2N + 1):3N)

    m_inf = [gating_m_inf(V[i], gating) for i in 1:N]
    h_inf = [gating_h_inf(V[i], gating) for i in 1:N]
    s_inf = [gating_s_inf(V[i], gating) for i in 1:N]
    kh = [gating_k_h(V[i], param, gating) for i in 1:N]

    s_total = sum(s)
    AV = A * V

    @inbounds for i in 1:N
        Vi, hi, si = V[i], h[i], s[i]

        I_Ca = g_Ca * m_inf[i]^3 * hi * (Vi - V_Ca)
        I_L = g_L * (Vi - V_L)
        I_syn = (N > 1) ? g_syn / (N - 1) * (Vi - V_syn) * (s_total - si) : zero(F)

        if g_el < 1e-14
            I_el = zero(F)
        elseif M[i] > 0
            LVi = M[i] * Vi - AV[i]
            I_el = g_el / M[i] * LVi
        else
            I_el = zero(F)
        end

        dV[i] = -(I_Ca + I_L + I_syn + I_el) / C
        dh[i] = kh[i] * (h_inf[i] - hi)
        ds[i] = kf * s_inf[i] * (1 - si) - kr * si
    end

    return nothing
end

"""
    rinzel_golomb_rhs!(dx, x, net, t)

Zero-allocation in-place ODE RHS for N coupled Rinzel–Golomb TRN neurons.
Uses pre-allocated scratch buffers stored in `net` (_m_inf, _h_inf, _s_inf, _kh, _AV).

THREAD-SAFETY: `net` must not be shared across concurrent ODE solves.
Each simulation must own its own `RinzelGolombNetwork` instance.

State layout:  x = [V₁…Vₙ | h₁…hₙ | s₁…sₙ]
"""
function rinzel_golomb_rhs!(
    dx::AbstractVector{F}, x::AbstractVector{F}, net::RinzelGolombNetwork{F}, t
) where {F<:AbstractFloat}
    (; N, param, gating, A, M, _m_inf, _h_inf, _s_inf, _kh, _AV) = net
    (; C, g_Ca, g_L, g_syn, g_el, V_Ca, V_L, V_syn, kr, kf) = param

    V = view(x, 1:N)
    h = view(x, (N + 1):2N)
    s = view(x, (2N + 1):3N)

    dV = view(dx, 1:N)
    dh = view(dx, (N + 1):2N)
    ds = view(dx, (2N + 1):3N)

    @inbounds for i in 1:N
        _m_inf[i] = gating_m_inf(V[i], gating)
        _h_inf[i] = gating_h_inf(V[i], gating)
        _s_inf[i] = gating_s_inf(V[i], gating)
        _kh[i] = gating_k_h(V[i], param, gating)
    end

    s_total = sum(s)

    if g_el >= 1e-14
        mul!(_AV, A, V)
    end

    @inbounds for i in 1:N
        Vi, hi, si = V[i], h[i], s[i]

        I_Ca = g_Ca * _m_inf[i]^3 * hi * (Vi - V_Ca)
        I_L = g_L * (Vi - V_L)
        I_syn = (N > 1) ? g_syn / (N - 1) * (Vi - V_syn) * (s_total - si) : zero(F)

        if g_el < 1e-14
            I_el = zero(F)
        elseif M[i] > 0
            LVi = M[i] * Vi - _AV[i]
            I_el = g_el / M[i] * LVi
        else
            I_el = zero(F)
        end

        dV[i] = -(I_Ca + I_L + I_syn + I_el) / C
        dh[i] = _kh[i] * (_h_inf[i] - hi)
        ds[i] = kf * _s_inf[i] * (1 - si) - kr * si
    end

    return nothing
end

"""
    make_rg_ic(net; rng) → Vector{F}

Random initial condition for a Rinzel–Golomb network.
  V ~ Uniform(-60, 0) mV
  h ~ Uniform(0, 1)
  s ~ Uniform(0, 1)
"""
function make_rg_ic(
    net::RinzelGolombNetwork{F}; rng::AbstractRNG=Random.GLOBAL_RNG
) where {F<:AbstractFloat}
    N = net.N
    V0 = rand(rng, F, N) .* 60 .- 60     # [-60, 0]
    h0 = rand(rng, F, N)                   # [0, 1]
    s0 = rand(rng, F, N)                   # [0, 1]
    return vcat(V0, h0, s0)
end
