"""
    Plotting.jl

Publication-quality figure helpers for TRN spindle analysis.
All functions return CairoMakie Figure objects.

- `plot_timeseries`         — time series for all 4 WC populations
- `plot_phase_plane`        — (T, R) phase-plane limit cycles coloured by spindle proximity
- `plot_spindle_conditions` — (ρ, R) joint-spindle-condition diagram
- `plot_bifurcation_overlay`— overlay bifurcation boundaries on parameter plane
- `plot_cluster_grid`       — (g_syn, g_el) cluster-count heatmap
"""

using CairoMakie
using ColorSchemes
using LaTeXStrings

const SPINDLE_CMAP = :plasma
const ESCAPE_COLOR = :yellow
const SPINDLE_COLOR = :blue

function _time_label()
    return "Time (model units)"
end

function _spindle_cbar_label()
    return "Physiological Proximity (Spindle Health)"
end

"""
    plot_timeseries(res; pop_names, transient_frac, kwargs...) → Figure

Plot all WC population time series.  The first four rows of `res.u` are
assumed to be [P, I, T, R].
"""
function plot_timeseries(
    res::SimResult{F};
    pop_names::Vector{String}=["PYR (P)", "IN (I)", "Relay (T)", "TRN (R)"],
    transient_frac::F=F(0.5),
    title::String="",
    figure_size::Tuple{Int,Int}=(900, 600),
    kwargs...,
) where {F<:AbstractFloat}
    n_skip = round(Int, transient_frac * Base.size(res.u, 2))
    t = res.t[(n_skip + 1):end]
    U = res.u[1:min(4, Base.size(res.u, 1)), (n_skip + 1):end]
    npop = Base.size(U, 1)

    fig = Figure(; size=figure_size)
    axes = Axis[]
    for i in 1:npop
        ax = Axis(fig[i, 1]; ylabel="Activity")
        (i == 1 && !isempty(title)) && (ax.title = title)
        lines!(ax, t, U[i, :]; label=pop_names[i], linewidth=1.5)
        push!(axes, ax)
    end
    axes[npop].xlabel = _time_label()
    linkxaxes!(axes...)
    return fig
end

"""
    plot_phase_plane(results; T_target, R_target, title) → Figure

Overlay limit cycles in the (T, R) plane, coloured by spindle proximity score.
`results` is a vector of `SimResult`.


"""
function plot_phase_plane(
    results::Vector{SimResult{F}};
    T_idx::Int=3,
    R_idx::Int=4,
    T_target::Tuple{F,F}=(F(4), F(6)),
    R_target::Tuple{F,F}=(F(10), F(15)),
    transient_frac::F=F(0.5),
    title::String="Wilson–Cowan Limit Cycles in the (T,R) Plane",
    figure_size::Tuple{Int,Int}=(700, 500),
    kwargs...,
) where {F<:AbstractFloat}
    fig = Figure(; size=figure_size)
    ax = Axis(
        fig[1, 1];
        xlabel="Relay Population Rate T [Hz]",
        ylabel="TRN Population Rate R [Hz]",
        title=title,
    )

    for res in results
        n_skip = round(Int, transient_frac * Base.size(res.u, 2))
        T_sig = res.u[T_idx, (n_skip + 1):end]
        R_sig = res.u[R_idx, (n_skip + 1):end]
        Tmax = maximum(T_sig)
        Rmax = maximum(R_sig)
        score = spindle_proximity_score(Tmax, Rmax; T_target=T_target, R_target=R_target)
        col = get(colorschemes[SPINDLE_CMAP], score)
        lines!(ax, T_sig, R_sig; color=col, linewidth=0.8, alpha=0.85)
    end

    Colorbar(fig[1, 2]; colormap=SPINDLE_CMAP, limits=(0, 1), label=_spindle_cbar_label())
    return fig
end

"""
    plot_spindle_conditions(ρ_vals, Rmin_vals, Rmax_vals;
                            ρ_target, R_target, Rmin_target,
                            title) → Figure

Draw the (ρ, R) joint-spindle-condition diagram.

For each delay value ρ, draws a vertical bar from Rmin to Rmax, coloured
by the joint spindle score.  Cycles satisfying ALL conditions are circled.
"""
function plot_spindle_conditions(
    ρ_vals::AbstractVector{F},
    Rmin_vals::AbstractVector{F},
    Rmax_vals::AbstractVector{F};
    ρ_target::Tuple{F,F}=(F(0.3), F(0.6)),
    R_target::Tuple{F,F}=(F(10), F(15)),
    Rmin_target::Tuple{F,F}=(F(0), F(4)),
    title::String="Wilson–Cowan Dynamics: Joint Spindle Conditions in (ρ, R)",
    size::Tuple{Int,Int}=(800, 500),
    kwargs...,
) where {F<:AbstractFloat}
    fig = Figure(; size=size)
    ax = Axis(fig[1, 1]; xlabel="ρ", ylabel="TRN Rate R [Hz]", title=title)

    satisfies_all = Bool[]
    for (ρ, Rmin, Rmax) in zip(ρ_vals, Rmin_vals, Rmax_vals)
        score = joint_spindle_score(
            ρ, Rmin, Rmax; ρ_target=ρ_target, R_target=R_target, Rmin_target=Rmin_target
        )
        col = get(colorschemes[SPINDLE_CMAP], score)
        lines!(ax, [ρ, ρ], [Rmin, Rmax]; color=col, linewidth=2)
        scatter!(ax, [ρ], [Rmin]; color=col, markersize=4)
        scatter!(ax, [ρ], [Rmax]; color=col, markersize=4)
        sat = (
            ρ_target[1] <= ρ <= ρ_target[2] &&
            R_target[1] <= Rmax <= R_target[2] &&
            Rmin_target[1] <= Rmin <= Rmin_target[2]
        )
        push!(satisfies_all, sat)
    end

    for (i, sat) in enumerate(satisfies_all)
        sat || continue
        scatter!(
            ax,
            [ρ_vals[i]],
            [Rmax_vals[i]];
            color=:black,
            markersize=7,
            strokewidth=2,
            strokecolor=:black,
        )
        scatter!(
            ax,
            [ρ_vals[i]],
            [Rmin_vals[i]];
            color=:black,
            markersize=7,
            strokewidth=2,
            strokecolor=:black,
        )
    end

    Colorbar(
        fig[1, 2];
        colormap=SPINDLE_CMAP,
        limits=(0, 1),
        label="Joint Physiological Proximity Score",
    )
    return fig
end

"""
    plot_cluster_grid(g_syns, g_els, mean_counts; kwargs...) → Figure

Heatmap of mean cluster count over the (g_syn, g_el) parameter plane.

"""
function plot_cluster_grid(
    g_syns::AbstractVector{F},
    g_els::AbstractVector{F},
    mean_counts::AbstractMatrix{F};
    title::String="Mean Number of Synchronisation Clusters",
    size::Tuple{Int,Int}=(700, 550),
    kwargs...,
) where {F<:AbstractFloat}
    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1];
        xlabel=L"g_{\mathrm{syn}}\;(\mathrm{mS/cm}^2)",
        ylabel=L"g_{\mathrm{el}}\;(\mathrm{mS/cm}^2)",
        title=title,
    )
    hm = heatmap!(ax, g_syns, g_els, mean_counts'; colormap=:turbo, colorrange=(0, 4))
    Colorbar(fig[1, 2], hm)
    return fig
end

"""
    plot_bifurcation_overlay(param_vals, eq_branches, cycle_ranges;
                             param_label, pop_label) → Figure

Generic 1D bifurcation diagram overlaying equilibrium branches and cycle
amplitude bands.

Arguments
─────────
- `param_vals`    — bifurcation parameter values
- `eq_branches`   — list of (values, stability) pairs for equilibrium branches
- `cycle_ranges`  — list of (param_vals, min_vals, max_vals, color) for limit cycles
"""
function plot_bifurcation_overlay(
    param_vals::AbstractVector{F},
    eq_branches::Vector{<:Tuple},
    cycle_ranges::Vector{<:Tuple};
    param_label::String="Connectivity parameter",
    pop_label::String="Component of equilibrium",
    title::String="",
    size::Tuple{Int,Int}=(700, 450),
    kwargs...,
) where {F<:AbstractFloat}
    fig = Figure(; size=size)
    ax = Axis(fig[1, 1]; xlabel=param_label, ylabel=pop_label, title=title)

    for (vals, is_stable, col) in eq_branches
        ls = is_stable ? :solid : :dash
        lines!(ax, param_vals, vals; color=col, linestyle=ls, linewidth=2)
    end

    for (p_vals, mins, maxs, col) in cycle_ranges
        band!(ax, p_vals, mins, maxs; color=(col, 0.25))
        lines!(ax, p_vals, maxs; color=col, linewidth=1)
    end

    return fig
end
