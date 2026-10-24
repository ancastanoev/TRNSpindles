# cluster_count_field.jl — mean-cluster-count heatmap over (g_syn, g_el).

using CairoMakie, JLD2, Colors, Statistics, Printf, LaTeXStrings

infile = isempty(ARGS) ? error("Pass the publication JLD2 path.") : ARGS[1]
d = load(infile)
g_syn = Float64.(d["g_syn_grid"]);
g_el = Float64.(d["g_el_grid"])
N_GS, N_GE = length(g_syn), length(g_el)

sg = Float64.(d["sim_g_syn"]);
se = Float64.(d["sim_g_el"])
snc = Int.(d["sim_n_cl"]);
sfp = Bool.(d["sim_is_fp"])

meanfield = fill(NaN, N_GS, N_GE);
fp_frac = zeros(N_GS, N_GE);
maxobs = 1.0
for i in 1:N_GS, j in 1:N_GE
    sel = findall(
        k -> isapprox(sg[k], g_syn[i]; atol=1e-9) && isapprox(se[k], g_el[j]; atol=1e-9),
        eachindex(sg),
    )
    isempty(sel) && continue
    fps = sfp[sel];
    fp_frac[i, j] = mean(fps)
    active = snc[sel][.!fps]
    if !isempty(active)
        ;
        m = mean(active);
        meanfield[i, j] = m;
        global maxobs = max(maxobs, m);
    end
end
quiescent_mask = fp_frac .>= 0.5
@printf("Mean-field range: 1.0 – %.2f clusters\n", maxobs)

set_theme!(; fonts=(; regular="TeX Gyre Heros", bold="TeX Gyre Heros Bold"))
field_cmap = cgrad([
    colorant"#1f6fb2",
    colorant"#4bb6c4",
    colorant"#9ed6a8",
    colorant"#f6d23b",
    colorant"#f49a1a",
    colorant"#d2480a",
    colorant"#7a0403",
])
FP_COLOR = colorant"#34394a"

fig = Figure(; size=(690, 930), figure_padding=(60, 60, 75, 45), backgroundcolor=:white)
ax = Axis(
    fig[1, 1];
    xlabel=L"\text{Chemical coupling }\; g_{\mathrm{syn}}\;\,(\text{mS/cm}^2)",
    ylabel=L"\text{Electrical coupling }\; g_{\mathrm{el}}\;\,(\text{mS/cm}^2)",
    xlabelsize=16,
    ylabelsize=16,
    xticklabelsize=11.5,
    yticklabelsize=11.5,
    xgridvisible=false,
    ygridvisible=false,
    topspinevisible=false,
    rightspinevisible=false,
    spinewidth=1.0,
)
ax.alignmode = Outside()   # axis cell = content area; labels protrude into figure_padding
ax.xticks=(1:1:N_GS, [@sprintf("%.2f", v) for v in g_syn]);
ax.xticklabelrotation=π/4
ystep = N_GE>14 ? 2 : 1
ax.yticks=(collect(1:ystep:N_GE), [@sprintf("%.3f", g_el[k]) for k in 1:ystep:N_GE])

fpmat = fill(NaN, N_GS, N_GE);
for i in 1:N_GS, j in 1:N_GE
    ;
    quiescent_mask[i, j] && (fpmat[i, j]=1.0);
end
heatmap!(
    ax,
    1:N_GS,
    1:N_GE,
    fpmat;
    colormap=cgrad([FP_COLOR, FP_COLOR]),
    colorrange=(0, 1),
    interpolate=false,
)
fd = copy(meanfield);
for i in 1:N_GS, j in 1:N_GE
    ;
    quiescent_mask[i, j] && (fd[i, j]=NaN);
end
hm = heatmap!(
    ax,
    1:N_GS,
    1:N_GE,
    fd;
    colormap=field_cmap,
    colorrange=(1.0, ceil(maxobs)),
    interpolate=false,
)
ax.limits=(0.5, N_GS+0.5, 0.5, N_GE+0.5);
ax.aspect = DataAspect()

Colorbar(
    fig[1, 2],
    hm;
    label=L"\text{Mean cluster count (uncapped)}",
    labelsize=14.5,
    ticklabelsize=12,
    width=16,
    spinewidth=1.0,
)
colsize!(fig.layout, 1, Fixed(500))
rowsize!(fig.layout, 1, Aspect(1, N_GE/N_GS))   # colorbar height = heatmap content height
colgap!(fig.layout, 14)
resize_to_layout!(fig)

mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))
save(joinpath(@__DIR__, "..", "..", "results", "figures", "fig3_cluster_field.pdf"), fig)
save(
    joinpath(@__DIR__, "..", "..", "results", "figures", "fig3_cluster_field.png"),
    fig;
    px_per_unit=5,
)
println("Saved results/figures/fig3_cluster_field.{pdf,png}")
