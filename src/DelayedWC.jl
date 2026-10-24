"""
    DelayedWC.jl

4-node Wilson–Cowan system with **weak Gamma-distributed delays** on the
TRN feedback, implemented via the **linear-chain trick**.

Mathematical formulation
────────────────────────
The Gamma kernel  h_{m,θ}(s) = s^{m-1} exp(-s/θ) / (θ^m Γ(m))  has mean ρ = m·θ.

Via the linear chain trick (integer shape m), the convolution

    A(t) = ∫₀^∞ h_{m,θ}(s) R(t-s) ds

is replaced by m auxiliary ODEs:

    θ · ȧ₁ = −a₁ + R(t)
    θ · ȧⱼ = −aⱼ + aⱼ₋₁     j = 2…m
    A(t)  = aₘ(t)

Full state vector (4 + m components):  u = [P, I, T, R, a₁, a₂, …, aₘ]

The slow negative feedback then enters the TRN equation as  −gₐ · A(t).

References
──────────
Macdonald (1978) — linear-chain trick for integral-delay equations.
Linear-chain construction (Macdonald 1978).
"""

# State-vector indexing helpers
@inline _base_idx() = 1:4
@inline _aux_idx(m::Int) = 5:(4 + m)

"""
    delayed_wc_rhs!(du, u, p, t)

In-place ODE for the delay-augmented 4-node WC system.
`p` must be a `WCDelayParams`.

State layout:
  u[1:4]      → [P, I, T, R]
  u[5:4+m]    → auxiliary chain [a₁, …, aₘ]
"""
function delayed_wc_rhs!(
    du::AbstractVector{T}, u::AbstractVector{T}, p::WCDelayParams{T}, t
) where {T<:AbstractFloat}
    wcp = p.wc
    m = p.m
    θ = gamma_scale(p)
    ga = p.ga

    P, I, Tv, R = u[1], u[2], u[3], u[4]

    # slow delayed feedback signal A(t) = aₘ
    A = u[4 + m]

    (;
        wPP,
        wPI,
        wII,
        wIP,
        wPT,
        wIT,
        wTP,
        wTR,
        wRP,
        wRT,
        wRR,
        αE,
        θE,
        βE,
        αI1,
        θI1,
        βI1,
        αI2,
        θI2,
        βI2,
        τ,
        rP,
        rI,
        rT,
        rR,
        B,
        IP,
        IT,
        IR,
        IC,
    ) = wcp

    uP = wPP * P - wPI * I + wPT * Tv + B + IP + IC
    uI = wIP * P - wII * I + wIT * Tv + B
    uT = wTP * P - wTR * R + B + IT
    uR = wRP * P + wRT * Tv - wRR * R + B + IR - ga * A  # ← slow feedback

    SP = wc_sigmoid(uP, αE, θE, βE)
    SI = wc_sigmoid(uI, αI1, θI1, βI1)
    ST = wc_sigmoid(uT, αE, θE, βE)
    SR = wc_sigmoid(uR, αI2, θI2, βI2)

    du[1] = (-P + (1 - rP * P) * SP) / τ
    du[2] = (-I + (1 - rI * I) * SI) / τ
    du[3] = (-Tv + (1 - rT * Tv) * ST) / τ
    du[4] = (-R + (1 - rR * R) * SR) / τ

    # linear chain auxiliaries
    # θ · ȧ₁ = −a₁ + R
    du[5] = (-u[5] + R) / θ

    # θ · ȧⱼ = −aⱼ + aⱼ₋₁    j = 2…m
    for j in 2:m
        du[4 + j] = (-u[4 + j] + u[4 + j - 1]) / θ
    end

    return nothing
end

"""
    make_delay_ic(p, u0_wc; init_chain=:zero) → Vector{T}

Construct the full initial-condition vector for the delayed system.

- `u0_wc`  — 4-element initial condition for [P, I, T, R].
- `init_chain` — `:zero` (cold start) or `:equilibrium` (chain initialised at R₀).
"""
function make_delay_ic(
    p::WCDelayParams{T}, u0_wc::AbstractVector{T}; init_chain::Symbol=:zero
) where {T<:AbstractFloat}
    m = p.m
    u0 = zeros(T, 4 + m)
    u0[1:4] .= u0_wc

    if init_chain === :equilibrium
        # warm start: all auxiliary variables equal to R₀
        R0 = u0_wc[4]
        u0[5:end] .= R0
    end
    # else: zero (cold start)

    return u0
end

"""
    delayed_wc_ndim(p) → Int

Return the dimension of the state space for the delay-augmented system.
"""
delayed_wc_ndim(p::WCDelayParams) = 4 + p.m
