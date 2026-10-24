"""
    Parameters.jl

Typed, immutable parameter structs for every model variant.
All structs use `@kwdef` so callers can override only what they need.

Biological reference ranges
───────────────────────────
Coupling weights.
Sigmoid parameters: Wilson & Cowan (1972); Schellenberger Costa et al. (2016).
Rinzel–Golomb parameters: Golomb & Rinzel (1994); trnlib.py (alexfikl).
"""

# Wilson–Cowan 4-node CTRC model (no delay)

"""
    WCParams

All parameters for the 4-node Wilson–Cowan CTRC system.

Populations
-----------
P — cortical pyramidal cells
I — cortical interneurons
T — thalamic relay cells
R — thalamic reticular nucleus (TRN)

Coupling convention
-------------------
wXY  = weight of the directed connection FROM node Y TO node X.
All weights wXY ≥ 0.  Sign is handled inside the ODE (I inhibits, P/T excite).

CALIBRATION: matched to the MATLAB/MatCont delayed-system definition (see matcont/WC_TRN_delay/)
─────────────────────────────────────────────────────────────────────
The previous defaults did NOT reproduce the paper: the network saturated
(P≈14, I≈20, T≈11, R≈17), gain→0, no oscillations.  The mapping below was
verified by reproducing the MatCont continuation files exactly (equilibrium
curve and a stable limit cycle).  MATLAB symbol → field mapping:

    ce1 → wPP      ci1 → wPI      cet1 → wPT   (T→P)
    ce2 → wIP      ci2 → wII      cet2 → wIT
    ce3 → wTP      (P→T)
    cir → wTR      (R→T inhibition; bifurcation parameter)
    cet3→ wRT      (T→R excitation; bifurcation parameter)
    crr → wRR      ce2 → wRP      (P→R collateral)
    PPT → B        (common background, was 0.0, now 0.6; IP folded into B → IP=0)

Default coupling here is set to a VERIFIED OSCILLATING operating point
(wPP=1.0, wRT=0.2, wTR=1.0): a stable limit cycle with R∈[3.7,12.9],
T∈[0.1,1.5], f≈0.34 Hz, matching wTR_cycles_wPP_1_0_wTP_1_2_wRT_0_2.mat.
To recover the RESTING (non-oscillating) baseline used in the solver script,
set wTR=0.8 (cir) — the equilibrium is then P≈2.7, I≈5.5, T≈1.5, R≈6.8.
"""
Base.@kwdef struct WCParams{F<:AbstractFloat}
    # coupling weights
    wPP::F = 1.0        # ce1   PYR → PYR  (recurrent excitation)
    wIP::F = 0.7        # ce2   PYR → IN   (E→I drive)
    wPI::F = 0.4        # ci1   IN  → PYR  (I→E feedback)
    wII::F = 0.1        # ci2   IN  → IN   (self-inhibition)
    wPT::F = 0.4        # cet1  TR  → PYR  (thalamocortical drive)
    wTP::F = 1.2        # ce3   PYR → TR   (corticothalamic drive)
    wIT::F = 0.7        # cet2  TR  → IN   (thalamic drive on interneurons)
    wRT::F = 0.2        # cet3  TR  → TRN  (relay→reticular excitation)   [bif param]
    wTR::F = 1.0        # cir   TRN → TR   (reticular→relay inhibition)   [bif param]
    wRR::F = 0.2        # crr   TRN → TRN  (intra-reticular self-inhibition)
    wRP::F = 0.7        # ce2   PYR → TRN  (corticothalamic collateral)

    # sigmoidal parameters (αX, θX, βX)
    # Excitatory (PYR and TR share the same class)
    αE::F = 50.0       # Ae1 = Ae2
    θE::F = 4.0        # te1 = te2
    βE::F = 1.2        # be
    # Inhibitory class 1 (IN)
    αI1::F = 50.0       # Ai1
    θI1::F = 3.8        # ti
    βI1::F = 2.5        # bi
    # Inhibitory class 2 (TRN — higher maximal output)
    αI2::F = 75.0       # Ai2  (MatCont continuation runs use 75)
    θI2::F = 3.8        # ti
    βI2::F = 2.5        # bi

    # time constant
    τ::F = 1.0

    # refractory/saturation factors
    rP::F = 0.05       # r1
    rI::F = 0.03       # r2
    rT::F = 0.03       # r3
    rR::F = 0.03       # r4  (MatCont P0 = 0.03; solver script uses 0.05)

    # external inputs
    B::F = 0.6        # PPT  common background drive
    IP::F = 0.0        # folded into B (PPT)
    IT::F = 2.8        # IT   additional input to TR
    IR::F = 1.5        # IR   additional input to TRN
    IC::F = 0.1        # IC   small PYR bias term
end

"""Return the default WC parameter set (verified oscillating point; see WCParams docstring)."""
default_wc_params(::Type{F}=Float64) where {F} = WCParams{F}()

# Distributed-delay Wilson–Cowan (linear-chain trick)

"""
    WCDelayParams

Extends WCParams with a weak Gamma-distributed delay on the TRN→TR feedback loop.

The Gamma kernel with integer shape m and scale θ has mean delay ρ = m·θ.
Via the linear-chain trick this is implemented as m auxiliary ODEs:

    θ · ȧ₁ = −a₁ + R̄(t)
    θ · ȧⱼ = −aⱼ + aⱼ₋₁    (j = 2…m)
    A(t)  = aₘ(t)

The delayed feedback enters the TRN equation as:  −gₐ · A(t)
"""
Base.@kwdef struct WCDelayParams{F<:AbstractFloat}
    wc::WCParams{F} = WCParams{F}()

    # delay kernel (Gamma distribution)
    m::Int = 2        # integer shape parameter
    ρ::F = 1.0      # mean delay  ρ = m · θ  (θ = ρ/m)

    # feedback gain
    ga::F = 0.0      # strength of slow negative feedback on TRN
end

"""Theta (scale) of the Gamma kernel: θ = ρ / m."""
gamma_scale(p::WCDelayParams) = p.ρ / p.m

"""Mean delay."""
gamma_mean_delay(p::WCDelayParams) = p.ρ

default_delay_params(::Type{F}=Float64) where {F} = WCDelayParams{F}()

# Rinzel–Golomb conductance-based TRN neuron

"""
    GatingFunction

Sigmoid gating function used by the Rinzel–Golomb model:

    S(V; θ, σ) = 1 / (1 + exp(−(V − θ) / σ))

A negative σ produces a decreasing sigmoid (used for h_∞ and k_∞).
"""
Base.@kwdef struct GatingFunction{F<:AbstractFloat}
    θ::F   # half-activation voltage  (mV)
    σ::F   # slope factor              (mV); negative → decreasing
end

(g::GatingFunction)(V::Number) = 1.0 / (1.0 + exp(-(V - g.θ) / g.σ))
(g::GatingFunction)(V::AbstractArray) = @. 1.0 / (1.0 + exp(-(V - g.θ) / g.σ))

"""
    RinzelGolombParams

Dimensional parameters for the Golomb–Rinzel TRN neuron model as used in
trnlib.py (alexfikl/2025-trn-gap-junctions-experiments).

State variables per neuron: (V, h, s)
  V — membrane potential       (mV)
  h — inactivation of I_Ca     (dimensionless ∈ [0,1])
  s — synaptic gate             (dimensionless ∈ [0,1])

Currents
  I_Ca = g_Ca · m_∞(V)³ · h · (V − V_Ca)
  I_L  = g_L  · (V − V_L)
  I_syn = g_syn / (N−1) · (V − V_syn) · Σⱼ≠ᵢ sⱼ
  I_el  = g_el / Mᵢ · Σⱼ Aᵢⱼ (Vᵢ − Vⱼ)   [gap junctions]

  C · dV/dt = −(I_Ca + I_L + I_syn + I_el)
  dh/dt     = k_h(V) · (h_∞(V) − h)
  ds/dt     = k_f · s_∞(V) · (1 − s) − k_r · s
"""
Base.@kwdef struct RinzelGolombParams{F<:AbstractFloat}
    # membrane & conductances
    C::F = 1.0e-3      # membrane capacitance       (F/cm²)
    g_Ca::F = 0.5         # max Ca²⁺ conductance       (mS/cm²)
    g_L::F = 0.05        # leak conductance            (mS/cm²)
    g_syn::F = 0.2         # synaptic conductance        (mS/cm²)
    g_el::F = 0.0         # gap-junction conductance    (mS/cm²)

    # reversal potentials
    V_Ca::F = 120.0       # Ca²⁺ reversal potential    (mV)
    V_L::F = -60.0       # leak reversal potential     (mV)
    V_syn::F = -80.0       # synaptic reversal potential (mV)

    # channel kinetics
    kr::F = 58.0        # channel closing rate        (1/s)
    kf::F = 1.0e3       # channel opening rate        (1/s)
    phi::F = 2.0e3       # h-kinetics scaling factor   (1/s)

    # voltage threshold for "firing" classification
    V_threshold::F = -40.0  # (mV)
end

"""
    RinzelGolombGating

Gating-function parameters for a single Rinzel–Golomb neuron.
"""
Base.@kwdef struct RinzelGolombGating{F<:AbstractFloat}
    m_inf::GatingFunction{F} = GatingFunction{F}(θ=-65.0, σ=7.8)
    h_inf::GatingFunction{F} = GatingFunction{F}(θ=-81.0, σ=-11.0)
    s_inf::GatingFunction{F} = GatingFunction{F}(θ=-45.0, σ=2.0)
    k_inf::GatingFunction{F} = GatingFunction{F}(θ=-162.3, σ=17.8)
end

function default_rinzel_golomb(::Type{F}=Float64) where {F}
    return (RinzelGolombParams{F}(), RinzelGolombGating{F}())
end
