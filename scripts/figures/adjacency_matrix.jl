# adjacency_matrix.jl — block-diagonal gap-junction adjacency matrix (publication).
#
# Reproduces the EXACT electrical-coupling topology used by the (g_syn, g_el)
# sweep in scripts/sweep_rg_driven.jl:
#       N_NEURONS = 1000, N_ADJ_CLUSTERS = 50, gap = 0, ADJ_SEED = 42
# i.e. 50 all-to-all gap-junction clusters of 20 cells, no between-cluster edges.
#
# All text is set in LaTeX (Computer Modern) via CairoMakie's MathTeXEngine,
# matching the thesis body font. Output is vector PDF.
#
# Usage (from the TRNSpindles project root):
#       julia --project scripts/figures/adjacency_matrix.jl
#
using TRNSpindles, CairoMakie, Random, Printf

const N_NEURONS = 1000
const N_ADJ_CLUSTERS = 50
const ADJ_SEED = 42
const NZ = 80                 # zoom: first 4 clusters (4 × 20)
const NAVY = RGBf(0.094, 0.212, 0.408)

A, ms, gaps = make_adjacency_random_clusters(
    MersenneTwister(ADJ_SEED), N_NEURONS, N_ADJ_CLUSTERS; gap=0
)

@printf(
    "N=%d  clusters=%d  block sizes: min=%d max=%d  total=%d  edges=%d\n",
    N_NEURONS,
    length(ms),
    minimum(ms),
    maximum(ms),
    sum(ms),
    Int(sum(A) ÷ 2)
)

cmap = cgrad([:white, NAVY]; categorical=true)

# Thesis-style theme: Computer Modern serif, thin frames, LaTeX everywhere.
set_theme!(
    Theme(;
        fontsize=11,
        fonts=(; regular="CMU Serif", bold="CMU Serif Bold"),  # falls back to LM Roman
        Axis=(
            xgridvisible=false,
            ygridvisible=false,
            spinewidth=0.8,
            xtickwidth=0.8,
            ytickwidth=0.8,
            xticksize=3,
            yticksize=3,
            titlealign=:left,
            titlegap=6,
        ),
    ),
)

fig = Figure(; size=(640, 320))

# (a) full matrix
ax1 = Axis(
    fig[1, 1];
    title=L"\mathbf{(a)}\;\; \mathrm{full\ adjacency}\ A",
    xlabel=L"\mathrm{neuron\ index}\ j",
    ylabel=L"\mathrm{neuron\ index}\ i",
    aspect=1,
    yreversed=true,
    xticks=([1, 250, 500, 750, 1000], string.([1, 250, 500, 750, 1000])),
    yticks=([1, 250, 500, 750, 1000], string.([1, 250, 500, 750, 1000])),
)
heatmap!(
    ax1, 1:N_NEURONS, 1:N_NEURONS, A; colormap=cmap, colorrange=(0, 1), interpolate=false
)
xlims!(ax1, 0.5, N_NEURONS + 0.5);
ylims!(ax1, N_NEURONS + 0.5, 0.5)

# (b) zoom on first four clusters
ax2 = Axis(
    fig[1, 2];
    title=L"\mathbf{(b)}\;\; \mathrm{first\ four\ clusters}",
    xlabel=L"\mathrm{neuron\ index}\ j",
    ylabel=L"\mathrm{neuron\ index}\ i",
    aspect=1,
    yreversed=true,
    xticks=([1, 20, 40, 60, 80], string.([1, 20, 40, 60, 80])),
    yticks=([1, 20, 40, 60, 80], string.([1, 20, 40, 60, 80])),
)
heatmap!(
    ax2, 1:NZ, 1:NZ, A[1:NZ, 1:NZ]; colormap=cmap, colorrange=(0, 1), interpolate=false
)
xlims!(ax2, 0.5, NZ + 0.5);
ylims!(ax2, NZ + 0.5, 0.5)

colgap!(fig.layout, 18)

mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))
out = joinpath(@__DIR__, "..", "..", "results", "figures", "fig8_adjacency.pdf")
save(out, fig; pt_per_unit=1)
@printf("saved %s\n", out)

set_theme!()   # restore defaults for any later scripts
fig
