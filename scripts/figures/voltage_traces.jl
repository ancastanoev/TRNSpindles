# voltage_traces.jl — representative TRN voltage traces per dynamical regime.
# Exports ONE four-panel, title-less figure:
# A = synchronized / 1-cluster
# B = 2-cluster
# C = 3-cluster
# D = multistable
#
# The LaTeX caption should be supplied outside the image.

using TRNSpindles,
    CairoMakie, Colors, SparseArrays, Random, Statistics, Printf, LaTeXStrings

latexnum(vs) = [latexstring(rstrip(rstrip(@sprintf("%.3f", v), '0'), '.')) for v in vs]

set_theme!(;
    fonts=(; regular="TeX Gyre Heros", bold="TeX Gyre Heros Bold"),
    Axis=(;
        xtickformat=latexnum,
        ytickformat=latexnum,
        spinewidth=1.0,
        xgridvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        xlabelsize=15,
        ylabelsize=15,
        xticklabelsize=12,
        yticklabelsize=12,
    ),
)

const PALETTE = [
    colorant"#D55E00",
    colorant"#0072B2",
    colorant"#009E73",
    colorant"#E69F00",
    colorant"#CC79A7",
]

const LW = 1.5

const G_DRIVE = 0.001
const P_PHYS = 3.0
const N_CYCLES = 8

const N = 1000
const NCLUST = 50
const ADJ_SEED = 42

# x-axis crop width.
# Use 1.5 for "one-something"; use 2.0 for "two-something".
const T_VIEW = 2.0

drive = extract_wc_drive(
    default_wc_params();
    u0=Float64[1, 5, 2, 5],
    tspan=(0.0, 2000.0),
    P_phys=P_PHYS,
    dt_save=0.01,
)

T_DRIVE = drive.P_WC * drive.κ

A, _, _ = make_adjacency_random_clusters(MersenneTwister(ADJ_SEED), N, NCLUST; gap=0)

A = sparse(A)
Pb, Gb = default_rinzel_golomb(Float64)

function run_cell(g_syn, g_el, seed)
    p = RinzelGolombParams{Float64}(;
        g_syn=g_syn,
        g_el=g_el,
        C=Pb.C,
        g_Ca=Pb.g_Ca,
        g_L=Pb.g_L,
        V_Ca=Pb.V_Ca,
        V_L=Pb.V_L,
        V_syn=Pb.V_syn,
        kr=Pb.kr,
        kf=Pb.kf,
        phi=Pb.phi,
    )

    net = RinzelGolombNetwork(p, Gb, A)
    dp = DrivenRGParams(net, drive; g_drive=G_DRIVE)
    u0 = make_rg_ic(net; rng=MersenneTwister(seed))

    res = solve_driven_rg(
        dp; u0=u0, tspan=(0.0, N_CYCLES * T_DRIVE), saveat=0.001, reltol=1e-6, abstol=1e-8
    )

    return res, net
end

function pick_spread(V::AbstractMatrix, k::Int)
    n = size(V, 1)
    k <= 1 && return [1]

    cand = unique(round.(Int, range(1, n; length=min(n, 80))))

    D = 1 .- cor(permutedims(V[cand, :]))
    D[isnan.(D)] .= 0.0

    chosen = [1]

    while length(chosen) < min(k, length(cand))
        s = [minimum(D[i, chosen]) for i in 1:length(cand)]
        s[chosen] .= -Inf
        push!(chosen, argmax(s))
    end

    return cand[chosen]
end

# Panel order:
# A = 1-cluster / synchronized
# B = 2-cluster
# C = 3-cluster
# D = multistable
panels = [
    (label="(a)", xfrac=1.0, gs=0.60, ge=0.020, seed=101, k=2),
    (label="(b)", xfrac=1.0, gs=0.30, ge=0.000, seed=101, k=2),
    (label="(c)", xfrac=0.55, gs=0.40, ge=0.000, seed=101, k=3),
    (label="(d)", xfrac=1.0, gs=0.45, ge=0.005, seed=202, k=4),
]

# Run all simulations first, then use common x/y limits.
sims = map(panels) do P
    @printf("running panel %s …\n", P.label)
    flush(stdout)

    res, net = run_cell(P.gs, P.ge, P.seed)

    t0 = max(res.t[1], res.t[end] - T_VIEW)
    keep = res.t .>= t0

    tw = res.t[keep] .- res.t[keep][1]
    Vk = res.u[1:net.N, keep]

    return (; tw, Vk, reps=pick_spread(Vk, P.k))
end

ymin = minimum(minimum(s.Vk[s.reps, :]) for s in sims)
ymax = maximum(maximum(s.Vk[s.reps, :]) for s in sims)

pad = 0.04 * (ymax - ymin)
ymin -= pad
ymax += pad

xmax = maximum(maximum(s.tw) for s in sims)

mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))

fig = Figure(; size=(900, 620), figure_padding=(16, 18, 12, 14), backgroundcolor=:white)

axes = Axis[]

for (idx, (P, s)) in enumerate(zip(panels, sims))
    row = idx <= 2 ? 1 : 2
    col = isodd(idx) ? 1 : 2

    ax = Axis(fig[row, col]; xlabel=L"\text{time (s)}", ylabel=L"V\ (\text{mV})")

    push!(axes, ax)

    for (j, ni) in enumerate(s.reps)
        alpha = if P.k <= 2
            1.0
        elseif P.k == 3
            0.92
        else
            0.85
        end

        lines!(ax, s.tw, s.Vk[ni, :]; color=(PALETTE[j], alpha), linewidth=LW)
    end

    xlims!(ax, 0, P.xfrac * xmax)
    ylims!(ax, ymin, ymax)

    text!(
        ax,
        0.035,
        0.94;
        text=L"\textbf{%$(P.label)}",
        space=:relative,
        align=(:left, :top),
        fontsize=15,
        color=:black,
    )
end

# Keep panels compact and aligned.
colgap!(fig.layout, 24)
rowgap!(fig.layout, 20)

save(
    joinpath(
        @__DIR__, "..", "..", "results", "figures", "fig6_voltage_traces_four_panel.pdf"
    ),
    fig,
)
save(
    joinpath(
        @__DIR__, "..", "..", "results", "figures", "fig6_voltage_traces_four_panel.png"
    ),
    fig;
    px_per_unit=6,
)

println("saved results/figures/fig6_voltage_traces_four_panel.{pdf,png}")
println("done — four-panel figure exported.")
