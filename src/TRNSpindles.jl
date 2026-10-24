"""
    TRNSpindles

Sleep-related corticothalamic rhythm simulation package.

Models
──────
1. 4-node Wilson–Cowan CTRC (no delay)
2. Gamma-distributed delayed WC (linear-chain trick)
3. Rinzel–Golomb conductance-based TRN neurons
"""
module TRNSpindles

# ── stdlib (always available, no Project.toml entry needed) ─────────────────
using LinearAlgebra
using SparseArrays
using Printf
using Logging
using Random
using Dates
using Statistics
using LaTeXStrings

# ── registered packages ──────────────────────────────────────────────────────
using DifferentialEquations
using OrdinaryDiffEq
using StaticArrays
using FFTW
using DSP
using StatsBase
using Clustering
using JLD2
using ProgressMeter

# ── Wilson–Cowan core ──────────────────────────────────────────────────────
include("Parameters.jl")
include("Sigmoids.jl")
include("WilsonCowan.jl")
include("DelayedWC.jl")
# ── Rinzel–Golomb TRN network ────────────────────────────────────────────────
include("RinzelGolomb.jl")
include("GraphUtils.jl")

# ── simulation, analysis, plotting ───────────────────────────────────────────
include("Simulate.jl")
include("Analysis.jl")
include("ClusterDetect.jl")
include("Plotting.jl")

# ── WC→RG bridge ─────────────────────────────────────────────────────────────
include("Bridge.jl")

# ── public API ────────────────────────────────────────────────────────────────
export
    # Parameters
    WCParams,
    WCDelayParams,
    RinzelGolombParams,
    default_wc_params,
    default_delay_params,
    default_rinzel_golomb,
    gamma_scale,
    gamma_mean_delay,

    # Sigmoids
    wc_sigmoid,
    sigmoid_derivative,

    # Wilson-Cowan
    wc_rhs!,
    wc_rhs_static,
    wc_jacobian,
    wc_equilibrium_stability,

    # Delayed WC
    delayed_wc_rhs!,
    make_delay_ic,
    delayed_wc_ndim,

    # Rinzel-Golomb
    RinzelGolombNetwork,
    rinzel_golomb_rhs!,
    rinzel_golomb_rhs_ref!,
    gating_m_inf,
    gating_h_inf,
    gating_s_inf,
    gating_k_h,
    make_rg_ic,

    # Graph utilities
    make_adjacency_all,
    make_adjacency_groups,
    make_laplacian,
    make_laplacian_sparse,
    make_partition_from_fractions,
    make_adjacency_set1,
    make_adjacency_set2,
    make_adjacency_set3,
    make_adjacency_set4,
    make_adjacency_random_clusters,

    # Simulation
    SimResult,
    solve_wc,
    solve_delayed_wc,
    solve_rinzel_golomb,

    # Analysis
    dominant_frequency,
    spindle_proximity_score,
    joint_spindle_score,
    classify_regime,
    analyze_timeseries,
    extract_peaks,
    mean_cycle_period,
    autocorrelation_freq,
    SpindleRegime,
    N2_SPINDLE,
    TRANSITION,
    N3_SLOW_WAVE,
    INFRASLOW,
    UNKNOWN,

    # Clustering
    ClusterResult,
    detect_clusters,
    pfeuty_chi,
    is_fixed_point,
    find_clusters,
    sweep_gsyn_gel,
    cluster_count_fraction,

    # Plotting
    plot_timeseries,
    plot_phase_plane,
    plot_spindle_conditions,
    plot_bifurcation_overlay,
    plot_cluster_grid,

    # Bridge (WC→RG open-loop drive)
    WCDriveResult,
    extract_wc_drive,
    DrivenRGParams,
    driven_rg_rhs!,
    solve_driven_rg,
    find_g_drive,
    CalibrationResult,
    calibrate,
    compute_R_eff,
    sweep_driven

end
