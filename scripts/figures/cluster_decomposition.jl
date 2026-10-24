# cluster_decomposition.jl — per-category seed-fraction decomposition.

using CairoMakie, JLD2, Colors, Statistics, Printf

infile = isempty(ARGS) ? error("Pass the publication JLD2 path.") : ARGS[1]
d = load(infile)
g_syn = Float64.(d["g_syn_grid"]);
g_el = Float64.(d["g_el_grid"])
N_GS, N_GE = length(g_syn), length(g_el)

fc = d["frac_cat"]   # N_GS × N_GE × 5 — [quiescent, sync, 2, 3, ≥4]

set_theme!(;
    fonts=(; regular="TeX Gyre Heros", bold="TeX Gyre Heros Bold"),
    Axis=(;
        spinewidth=0.9,
        xgridvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        xlabelsize=14,
        ylabelsize=14,
        xticklabelsize=11,
        yticklabelsize=11,
        titlesize=15,
        titlefont=:bold,
        titlealign=:left,
    ),
)

frac_cmap = cgrad(:dense)

fig = Figure(; size=(1160, 760), figure_padding=16, backgroundcolor=:white)

panel_grid = fig[1, 1] = GridLayout()
cbar_grid = fig[1, 2] = GridLayout()

xtick_args = (1:2:N_GS, [@sprintf("%.2f", g_syn[k]) for k in 1:2:N_GS])
ytick_args = (1:3:N_GE, [@sprintf("%.3f", g_el[k]) for k in 1:3:N_GE])

panel_defs = [
    ("(a)  Quiescent", fc[:, :, 1], 1, 1:2, false, true),
    ("(b)  Synchronized", fc[:, :, 2], 1, 3:4, false, false),
    ("(c)  Two clusters", fc[:, :, 3], 1, 5:6, true, false),
    ("(d)  Three clusters", fc[:, :, 4], 2, 2:3, true, true),
    ("(e)  Four or more", fc[:, :, 5], 2, 4:5, true, false),
]

hm_ref = nothing
for (title_str, mat, r, cspan, show_x, show_y) in panel_defs
    global hm_ref
    ax = Axis(
        panel_grid[r, cspan];
        title=title_str,
        xlabel=show_x ? L"g_\mathrm{syn}\ \mathrm{(mS/cm^2)}" : "",
        ylabel=show_y ? L"g_\mathrm{el}\ \mathrm{(mS/cm^2)}" : "",
        xticklabelsvisible=show_x,
        yticklabelsvisible=show_y,
        xticksvisible=show_x,
        yticksvisible=show_y,
    )
    ax.xticks = xtick_args;
    ax.yticks = ytick_args;
    ax.xticklabelrotation = π/4

    hm = heatmap!(
        ax,
        1:N_GS,
        1:N_GE,
        mat;
        colormap=frac_cmap,
        colorrange=(0.0, 1.0),
        interpolate=false,
    )
    ax.limits = (0.5, N_GS+0.5, 0.5, N_GE+0.5);
    ax.aspect = 1.0
    isnothing(hm_ref) && (hm_ref = hm)
end

Colorbar(
    cbar_grid[1, 1],
    hm_ref;
    vertical=true,
    width=16,
    label="fraction of seeds",
    labelsize=13,
    ticks=([0.0, 0.5, 1.0], ["0", "0.5", "1"]),
    ticklabelsize=11,
    spinewidth=0.9,
    height=Relative(0.92),
    valign=:center,
)

colgap!(panel_grid, 8);
rowgap!(panel_grid, 16)
colgap!(fig.layout, 14)
colsize!(fig.layout, 2, Fixed(74))

mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))
save(joinpath(@__DIR__, "..", "..", "results", "figures", "fig4_fraction_panels.pdf"), fig)
save(
    joinpath(@__DIR__, "..", "..", "results", "figures", "fig4_fraction_panels.png"),
    fig;
    px_per_unit=4,
)
println("Saved results/figures/fig4_fraction_panels.{pdf,png}")
