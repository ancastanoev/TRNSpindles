# synchrony_chi_field.jl — mean synchrony χ heatmap over (g_syn, g_el).

using CairoMakie, JLD2, Statistics, Printf, LaTeXStrings

infile = isempty(ARGS) ? error("Pass the publication JLD2 path.") : ARGS[1]
d = load(infile)
g_syn = Float64.(d["g_syn_grid"]);
g_el = Float64.(d["g_el_grid"])
N_GS, N_GE = length(g_syn), length(g_el)
chi = Float64.(d["mean_chi"])

set_theme!(; fonts=(; regular="TeX Gyre Heros", bold="TeX Gyre Heros Bold"))
chi_cmap = cgrad(:magma; rev=true)

fig = Figure(; size=(740, 580), figure_padding=14, backgroundcolor=:white)
ax = Axis(
    fig[1, 1];
    xlabel=L"\text{Chemical coupling }\; g_{\mathrm{syn}}\;\,(\text{mS/cm}^2)",
    ylabel=L"\text{Electrical coupling }\; g_{\mathrm{el}}\;\,(\text{mS/cm}^2)",
    xlabelsize=14,
    ylabelsize=14,
    xticklabelsize=11,
    yticklabelsize=11,
    xgridvisible=false,
    ygridvisible=false,
    topspinevisible=false,
    rightspinevisible=false,
    spinewidth=0.9,
)

ax.xticks = (1:N_GS, [latexstring(@sprintf("%.2f", v)) for v in g_syn]);
ax.xticklabelrotation = π/4
ystep = N_GE > 14 ? 2 : 1
ax.yticks = (
    collect(1:ystep:N_GE), [latexstring(@sprintf("%.3f", g_el[k])) for k in 1:ystep:N_GE]
)

hm = heatmap!(
    ax, 1:N_GS, 1:N_GE, chi; colormap=chi_cmap, colorrange=(0, 1), interpolate=false
)

ax.limits = (0.5, N_GS + 0.5, 0.5, N_GE + 0.5)
ax.aspect = 1.3

Colorbar(
    fig[1, 2],
    hm;
    label=L"\chi\;\;(1 = \text{synchronized})",
    labelsize=13,
    width=13,
    ticklabelsize=11,
    spinewidth=0.9,
    ticks=(0:0.25:1, [latexstring(@sprintf("%.2f", t)) for t in 0:0.25:1]),
)
colgap!(fig.layout, 10)

mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))
save(joinpath(@__DIR__, "..", "..", "results", "figures", "fig7_chi_field.pdf"), fig)
save(
    joinpath(@__DIR__, "..", "..", "results", "figures", "fig7_chi_field.png"),
    fig;
    px_per_unit=6,
)
println("Saved results/figures/fig7_chi_field.{pdf,png}")
