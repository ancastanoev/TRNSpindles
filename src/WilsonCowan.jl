"""
    WilsonCowan.jl

Right-hand side of the 4-node Wilson–Cowan CTRC system (no delays).

State vector layout:  u = [P, I, T, R]

Equations:

    τ Ṗ = −P + (1 − rP P) · S_E( wPP·P − wPI·I + wPT·T + B + IP + IC )
    τ İ = −I + (1 − rI I) · S_I( wIP·P − wII·I + wIT·T + B )
    τ Ṫ = −T + (1 − rT T) · S_E( wTP·P − wTR·R + B + IT )
    τ Ṙ = −R + (1 − rR R) · S_R( wRP·P + wRT·T − wRR·R + B + IR )

Population interpretation
  P — cortical pyramidal      (excitatory)
  I — cortical interneurons   (inhibitory)
  T — thalamic relay          (excitatory)
  R — thalamic reticular (TRN)(inhibitory)

Note on sign convention:
  Inhibitory connections appear with a MINUS sign inside the sigmoid argument;
  excitatory connections appear with a PLUS sign.
"""

"""
    wc_rhs!(du, u, p, t)

In-place ODE right-hand side for the 4-node Wilson–Cowan system.
`p` must be a `WCParams`.
"""
function wc_rhs!(
    du::AbstractVector{T}, u::AbstractVector{T}, p::WCParams{T}, t
) where {T<:AbstractFloat}
    P, I, Tv, R = u[1], u[2], u[3], u[4]
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
    ) = p

    uP = wPP * P - wPI * I + wPT * Tv + B + IP + IC
    uI = wIP * P - wII * I + wIT * Tv + B
    uT = wTP * P - wTR * R + B + IT
    uR = wRP * P + wRT * Tv - wRR * R + B + IR

    SP = wc_sigmoid(uP, αE, θE, βE)
    SI = wc_sigmoid(uI, αI1, θI1, βI1)
    ST = wc_sigmoid(uT, αE, θE, βE)
    SR = wc_sigmoid(uR, αI2, θI2, βI2)

    du[1] = (-P + (1 - rP * P) * SP) / τ
    du[2] = (-I + (1 - rI * I) * SI) / τ
    du[3] = (-Tv + (1 - rT * Tv) * ST) / τ
    du[4] = (-R + (1 - rR * R) * SR) / τ

    return nothing
end

"""
    wc_rhs_static(u, p, t) → SVector{4}

Out-of-place (StaticArrays) version.  ~10× faster for single-trajectory work.
"""
function wc_rhs_static(u::SVector{4,T}, p::WCParams{T}, t) where {T<:AbstractFloat}
    P, I, Tv, R = u[1], u[2], u[3], u[4]
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
    ) = p

    uP = wPP * P - wPI * I + wPT * Tv + B + IP + IC
    uI = wIP * P - wII * I + wIT * Tv + B
    uT = wTP * P - wTR * R + B + IT
    uR = wRP * P + wRT * Tv - wRR * R + B + IR

    SP = wc_sigmoid(uP, αE, θE, βE)
    SI = wc_sigmoid(uI, αI1, θI1, βI1)
    ST = wc_sigmoid(uT, αE, θE, βE)
    SR = wc_sigmoid(uR, αI2, θI2, βI2)

    return SVector{4,T}(
        (-P + (1 - rP * P) * SP) / τ,
        (-I + (1 - rI * I) * SI) / τ,
        (-Tv + (1 - rT * Tv) * ST) / τ,
        (-R + (1 - rR * R) * SR) / τ,
    )
end

"""
    wc_jacobian(u, p) → Matrix{T}

Analytic 4×4 Jacobian of the WC ODE at state `u` with parameters `p`.
Useful for equilibrium classification and Hopf-bifurcation detection.
"""
function wc_jacobian(u::AbstractVector{T}, p::WCParams{T}) where {T<:AbstractFloat}
    P, I, Tv, R = u[1], u[2], u[3], u[4]
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
    ) = p

    uP = wPP*P - wPI*I + wPT*Tv + p.B + p.IP + p.IC
    uI = wIP*P - wII*I + wIT*Tv + p.B
    uT = wTP*P - wTR*R + p.B + p.IT
    uR = wRP*P + wRT*Tv - wRR*R + p.B + p.IR

    SP = wc_sigmoid(uP, αE, θE, βE)
    SI = wc_sigmoid(uI, αI1, θI1, βI1)
    ST = wc_sigmoid(uT, αE, θE, βE)
    SR = wc_sigmoid(uR, αI2, θI2, βI2)
    dSP = sigmoid_derivative(uP, αE, θE, βE)
    dSI = sigmoid_derivative(uI, αI1, θI1, βI1)
    dST = sigmoid_derivative(uT, αE, θE, βE)
    dSR = sigmoid_derivative(uR, αI2, θI2, βI2)

    # effective gains including refractory factor
    gP = (1 - rP * P) * dSP
    gI = (1 - rI * I) * dSI
    gT = (1 - rT * Tv) * dST
    gR = (1 - rR * R) * dSR

    J = zeros(T, 4, 4)

    # ∂Ṗ/∂(P,I,T,R)
    J[1, 1] = (-1 + gP * wPP - rP * SP) / τ
    J[1, 2] = (-gP * wPI) / τ
    J[1, 3] = (gP * wPT) / τ
    J[1, 4] = zero(T)

    # ∂İ/∂(P,I,T,R)
    J[2, 1] = (gI * wIP) / τ
    J[2, 2] = (-1 - gI * wII - rI * SI) / τ
    J[2, 3] = (gI * wIT) / τ
    J[2, 4] = zero(T)

    # ∂Ṫ/∂(P,I,T,R)
    J[3, 1] = (gT * wTP) / τ
    J[3, 2] = zero(T)
    J[3, 3] = (-1 - rT * ST) / τ
    J[3, 4] = (-gT * wTR) / τ

    # ∂Ṙ/∂(P,I,T,R)
    J[4, 1] = (gR * wRP) / τ
    J[4, 2] = zero(T)
    J[4, 3] = (gR * wRT) / τ
    J[4, 4] = (-1 - gR * wRR - rR * SR) / τ

    return J
end

"""
    wc_equilibrium_stability(u_star, p) → (eigenvalues, is_stable)

Return the eigenvalues of the Jacobian at equilibrium `u_star` and whether
all real parts are negative (stable equilibrium).
"""
function wc_equilibrium_stability(
    u_star::AbstractVector{T}, p::WCParams{T}
) where {T<:AbstractFloat}
    J = wc_jacobian(u_star, p)
    λs = eigvals(J)
    return λs, all(real(λ) < 0 for λ in λs)
end
