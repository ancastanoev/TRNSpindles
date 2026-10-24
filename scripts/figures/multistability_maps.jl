# multistability_maps.jl — seed-entropy and coexisting-regime count maps.

using CairoMakie, JLD2, Colors, Statistics, Printf, LaTeXStrings

latexnum(vs) = [latexstring(rstrip(rstrip(@sprintf("%.3f", v), '0'), '.')) for v in vs]

infile = isempty(ARGS) ? error("Pass the publication JLD2 path.") : ARGS[1]
d = load(infile)
g_syn = Float64.(d["g_syn_grid"]);
g_el = Float64.(d["g_el_grid"])
N_GS, N_GE = length(g_syn), length(g_el)
H = Float64.(d["cat_entropy"])   # Shannon bits per cell
nd = Int.(d["n_distinct"])        # distinct categories per cell

set_theme!(;
    fonts=(; regular="TeX Gyre Heros", bold="TeX Gyre Heros Bold"),
    Axis=(;
        spinewidth=1.0,
        xgridvisible=false,
        ygridvisible=false,
        topspinevisible=false,
        rightspinevisible=false,
        xlabelsize=14,
        ylabelsize=14,
        xticklabelsize=10.5,
        yticklabelsize=10.5,
        titlesize=15,
        titlefont=:bold,
        titlealign=:left,
    ),
)

H_cmap = cgrad([
    colorant"#0b0a1a",
    colorant"#3b0f70",
    colorant"#8c2981",
    colorant"#de4968",
    colorant"#fe9f6d",
    colorant"#fcfdbf",
])
nd_cmap = cgrad(
    [colorant"#0b0a1a", colorant"#2a9d8f", colorant"#9be7d8"], 3; categorical=true
)

function setticks!(ax)
    return (
        ax.xticks=(1:2:N_GS, [latexstring(@sprintf("%.2f", g_syn[k])) for k in 1:2:N_GS]);
        ax.xticklabelrotation=π/4;
        ax.yticks=(1:2:N_GE, [latexstring(@sprintf("%.3f", g_el[k])) for k in 1:2:N_GE])
    )
end

fig = Figure(; size=(1120, 540), figure_padding=16, backgroundcolor=:white)

ax1 = Axis(
    fig[1, 1];
    xlabel=L"g_{\mathrm{syn}}\ (\text{mS/cm}^2)",
    ylabel=L"g_{\mathrm{el}}\ (\text{mS/cm}^2)",
    title=L"\textbf{(a)}\;\; \text{Seed entropy }\, H\ (\text{bits})",
)
setticks!(ax1)
hm_H = heatmap!(
    ax1,
    1:N_GS,
    1:N_GE,
    H;
    colormap=H_cmap,
    colorrange=(0, max(1e-9, maximum(H))),
    interpolate=false,
)
ax1.limits = (0.5, N_GS+0.5, 0.5, N_GE+0.5);
ax1.aspect = 1.0
Colorbar(
    fig[1, 2],
    hm_H;
    label=L"H\ (\text{bits})",
    labelsize=13,
    width=13,
    tickformat=latexnum,
    ticklabelsize=11,
    spinewidth=1.0,
)

ax2 = Axis(
    fig[1, 3];
    xlabel=L"g_{\mathrm{syn}}\ (\text{mS/cm}^2)",
    ylabel="",
    title=L"\textbf{(b)}\;\; \text{Coexisting regimes across seeds}",
)
setticks!(ax2)
hm_nd = heatmap!(
    ax2, 1:N_GS, 1:N_GE, nd; colormap=nd_cmap, colorrange=(0.5, 3.5), interpolate=false
)
ax2.limits = (0.5, N_GS+0.5, 0.5, N_GE+0.5);
ax2.aspect = 1.0
Colorbar(
    fig[1, 4],
    hm_nd;
    label=L"\text{coexisting regimes}",
    labelsize=13,
    width=13,
    ticks=(1:3, [L"\text{monostable}", L"\text{bistable}", L"\text{tristable}"]),
    ticklabelsize=11,
    spinewidth=1.0,
)

linkaxes!(ax1, ax2)
hideydecorations!(ax2; grid=false)

colgap!(fig.layout, 10)

nmulti = count(>(1), nd)
@printf(
    "multistable cells (>1 regime): %d / %d  (%.1f%%)\n",
    nmulti,
    length(nd),
    100nmulti/length(nd)
)

mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))
save(joinpath(@__DIR__, "..", "..", "results", "figures", "fig5_multistability.pdf"), fig)
save(
    joinpath(@__DIR__, "..", "..", "results", "figures", "fig5_multistability.png"),
    fig;
    px_per_unit=5,
)
println("Saved results/figures/fig5_multistability.{pdf,png}")
