# sweep_rg_driven.jl — WC-driven RG (g_syn × g_el) parameter-plane sweep.
#
# Production script — same code for dry run and full publication sweep.
# Set DRY_RUN = true for the 60-sim coarse grid (~5-10 min at 16 threads).
# Set DRY_RUN = false for the 1936-sim publication grid (~3 h at 16 threads).
#
# Frozen bridge parameters (Gate 4): a=15.149, g_drive=0.001, one-way open-loop.
#
# Launch:  julia --threads 16 --project scripts/sweep_rg_driven.jl

using TRNSpindles, SparseArrays, Random, Printf, Statistics, StatsBase
using JLD2, ProgressMeter, Dates

# CONFIGURATION  (only block that differs between dry-run and full sweep)

const DRY_RUN = false   # set false for full publication sweep

if DRY_RUN
    const G_SYN_GRID = Float64[0.20, 0.35, 0.50, 0.65]
    const G_EL_GRID = Float64[0.00, 0.01, 0.03, 0.06, 0.12]
    const SEEDS = Int[101, 202, 303]
    const RUN_TAG = "dryrun"
else
    # Full publication sweep  (11 g_syn × 16 g_el × 11 seeds = 1936 sims)
    # g_el: dense in [0,0.02] where transition + multistability live,
    #        medium in (0.02,0.06], sparse sentinels beyond.
    const G_SYN_GRID = collect(Float64, 0.20:0.05:0.70)
    const G_EL_GRID = vcat(
        collect(Float64, 0.0000:0.0025:0.0200),  # 9 pts
        collect(Float64, 0.0300:0.0100:0.0600),  # 4 pts
        Float64[0.09, 0.12, 0.15],
    )                   # 3 pts
    const SEEDS = Int[101, 202, 303, 404, 505, 606, 707, 808, 909, 1001, 1102]
    const RUN_TAG = "publication"
end

# FIXED PARAMETERS (do not change)

const G_DRIVE_LOCKED = 0.001
const P_PHYS = 3.0
const N_CYCLES = 8
const N_NEURONS = 1000
const N_ADJ_CLUSTERS = 50
const ADJ_SEED = 42
const SAVEAT = 0.01
const RELTOL = 1e-6
const ABSTOL = 1e-8

const N_GS = length(G_SYN_GRID)
const N_GE = length(G_EL_GRID)
const N_SS = length(SEEDS)
const N_CATS = 5   # bins: 0=quiescent 1=sync 2=2-clust 3=3-clust 4=≥4-clust

# build shared read-only infrastructure

println("═"^70)
println("WC-driven RG sweep  [$(DRY_RUN ? "DRY RUN (coarse)" : "PUBLICATION")]")
@printf(
    "  Grid    : %d g_syn × %d g_el × %d seeds = %d sims\n",
    N_GS,
    N_GE,
    N_SS,
    N_GS*N_GE*N_SS
)
@printf("  Threads : %d\n", Threads.nthreads())
println("═"^70);
println()

print("WC drive … ");
flush(stdout)
const DRIVE = extract_wc_drive(
    default_wc_params();
    u0=Float64[1, 5, 2, 5],
    tspan=(0.0, 2000.0),
    P_phys=P_PHYS,
    dt_save=0.01,
)
const T_DRIVE = DRIVE.P_WC * DRIVE.κ
@printf("T_drive = %.4f s\n\n", T_DRIVE)

print("Sparse adjacency … ");
flush(stdout)
_A_tmp, _, _ = make_adjacency_random_clusters(
    MersenneTwister(ADJ_SEED), N_NEURONS, N_ADJ_CLUSTERS; gap=0
)
const A_SPARSE = sparse(_A_tmp);
_A_tmp = nothing
const P_BASE, G_BASE = default_rinzel_golomb(Float64)
@printf(
    "nnz=%d  (%.2f%% fill)\n\n",
    nnz(A_SPARSE),
    100.0*nnz(A_SPARSE)/(N_NEURONS*(N_NEURONS-1))
)

# cell / sim lists

struct SweepCell
    gi::Int
    ei::Int
    si::Int
    g_syn::Float64
    g_el::Float64
    seed::Int
end

# Flat list: iteration order gi→ei→si (gi varies fastest in Julia comprehension)
const CELL_LIST = vec([
    SweepCell(gi, ei, si, G_SYN_GRID[gi], G_EL_GRID[ei], SEEDS[si]) for
    gi in eachindex(G_SYN_GRID), ei in eachindex(G_EL_GRID), si in eachindex(SEEDS)
])

const N_SIMS = length(CELL_LIST)
const N_CELLS = N_GS * N_GE   # unique (g_syn, g_el) pairs

# Linear cell index: ci = (gi-1)*N_GE + ei  (1-based, row=gi, col=ei)
cell_idx(gi, ei) = (gi - 1)*N_GE + ei

# Precompute: for each cell, which sim indices belong to it
const CELL_SIM_MAP = [Int[] for _ in 1:N_CELLS]
for k in eachindex(CELL_LIST)
    c = CELL_LIST[k]
    push!(CELL_SIM_MAP[cell_idx(c.gi, c.ei)], k)
end

# core simulation function (pure — all state is local to this call)

function sim_one(g_syn::Float64, g_el::Float64, seed::Int)
    p = RinzelGolombParams{Float64}(;
        g_syn=g_syn,
        g_el=g_el,
        C=P_BASE.C,
        g_Ca=P_BASE.g_Ca,
        g_L=P_BASE.g_L,
        V_Ca=P_BASE.V_Ca,
        V_L=P_BASE.V_L,
        V_syn=P_BASE.V_syn,
        kr=P_BASE.kr,
        kf=P_BASE.kf,
        phi=P_BASE.phi,
    )
    net = RinzelGolombNetwork(p, G_BASE, A_SPARSE)
    dp = DrivenRGParams(net, DRIVE; g_drive=G_DRIVE_LOCKED)
    u0 = make_rg_ic(net; rng=MersenneTwister(seed))
    res = solve_driven_rg(
        dp;
        u0=u0,
        tspan=(0.0, N_CYCLES*T_DRIVE),
        saveat=SAVEAT,
        reltol=RELTOL,
        abstol=ABSTOL,
    )
    cr = detect_clusters(res, net)
    return (chi=cr.pfeuty_chi, n_cl=cr.n_clusters, is_fp=cr.is_fixed_point)
end

# Map n_clusters → category bin (0=quiescent, 1=sync, 2=2-clust, 3=3-clust, 4=≥4-clust)
cat_bin(n_cl::Int)::Int = n_cl <= 0 ? 0 : min(n_cl, 4)

const CAT_LABELS = ["quiescent", "sync", "2-clust", "3-clust", "≥4-clust"]
cat_bin_label(b::Int)::String = CAT_LABELS[b + 1]

# Shannon entropy over a 5-bin probability vector
function bin_entropy(frac::AbstractVector{Float64})::Float64
    h = 0.0
    for p in frac
        p > 0.0 && (h -= p * log2(p))
    end
    return h
end

# JIT warmup
print("JIT warmup … ");
flush(stdout)
sim_one(G_SYN_GRID[1], G_EL_GRID[1], SEEDS[1])
println("done.\n")

# parallel sweep

struct SimRecord
    chi::Float64
    n_cl::Int
    is_fp::Bool
    wall::Float64
end
sim_records = Vector{SimRecord}(undef, N_SIMS)

cell_seeds_done = zeros(Int, N_CELLS)   # protected by lock_prog

prog = Progress(N_SIMS; desc="Sweep: ", dt=0.5, barlen=40, showspeed=true)
lock_prog = ReentrantLock()

t_start = time()

Threads.@threads for k in eachindex(CELL_LIST)
    c = CELL_LIST[k]
    ci = cell_idx(c.gi, c.ei)
    wall = @elapsed val = sim_one(c.g_syn, c.g_el, c.seed)

    lock(lock_prog) do
        sim_records[k] = SimRecord(val.chi, val.n_cl, val.is_fp, wall)
        cell_seeds_done[ci] += 1

        if cell_seeds_done[ci] == N_SS
            idx = CELL_SIM_MAP[ci]
            ncls = [sim_records[j].n_cl for j in idx]
            bins = [cat_bin(nc) for nc in ncls]
            counts = [count(==(b), bins) for b in 0:4]  # per-bin seed counts
            dom_b = argmax(counts) - 1                  # dominant bin (0-based)
            fracs = counts ./ N_SS
            h = bin_entropy(fracs)
            n_dist = count(>(0), fracs)
            mean_t = mean(sim_records[j].wall for j in idx)

            # compact distribution string: only non-zero bins
            dist_str = join(
                ["$(cat_bin_label(b)[1])$(counts[b+1])" for b in 0:4 if counts[b + 1] > 0],
                " ",
            )
            multi_flag = n_dist > 1 ? @sprintf("  H=%.2f n=%d ←MULTI", h, n_dist) : ""

            @printf(
                "  cell %3d/%-3d  (g_syn=%.2f, g_el=%.4f)  dom=%-9s  [%s]%s  %.1fs\n",
                ci,
                N_CELLS,
                c.g_syn,
                c.g_el,
                cat_bin_label(dom_b),
                dist_str,
                multi_flag,
                mean_t
            )
            flush(stdout)
        end
        ProgressMeter.next!(prog)
    end
end

t_total = time() - t_start
println()
@printf("Sweep done in %.1f s  (%.1f s/sim mean)\n\n", t_total, t_total/N_SIMS)

# aggregate (n_gs × n_ge matrices over seeds)
#
# frac_cat[gi, ei, b+1]  — fraction of seeds in bin b  (b = 0:4)
# cat_entropy[gi, ei]    — Shannon entropy of the bin distribution (bits)
# n_distinct[gi, ei]     — number of bins with at least one seed
# dom_cat_bin[gi, ei]    — bin of the dominant (most-frequent) category
# mean_chi / std_chi     — Pfeuty χ statistics across seeds

frac_cat = zeros(Float64, N_GS, N_GE, N_CATS)
cat_entropy = zeros(Float64, N_GS, N_GE)
n_distinct = zeros(Int, N_GS, N_GE)
dom_cat_bin = zeros(Int, N_GS, N_GE)
mean_chi = zeros(Float64, N_GS, N_GE)
std_chi = zeros(Float64, N_GS, N_GE)
mean_wall = zeros(Float64, N_GS, N_GE)

for gi in eachindex(G_SYN_GRID), ei in eachindex(G_EL_GRID)
    ci = cell_idx(gi, ei)
    idx = CELL_SIM_MAP[ci]
    isempty(idx) && continue

    ncls = [sim_records[k].n_cl for k in idx]
    chis = [sim_records[k].chi for k in idx]
    walls = [sim_records[k].wall for k in idx]

    bins = [cat_bin(nc) for nc in ncls]
    for b in 0:4
        frac_cat[gi, ei, b + 1] = count(==(b), bins) / N_SS
    end

    fracs = view(frac_cat, gi, ei, :)
    cat_entropy[gi, ei] = bin_entropy(fracs)
    n_distinct[gi, ei] = count(>(0.0), fracs)
    dom_cat_bin[gi, ei] = argmax(fracs) - 1   # back to 0-based bin

    mean_chi[gi, ei] = mean(chis)
    std_chi[gi, ei] = length(chis) > 1 ? std(chis) : 0.0
    mean_wall[gi, ei] = mean(walls)
end

# summary table (sorted: most multistable cells first)
println("─"^88)
println("CELL SUMMARY  (per-seed category distribution, $(N_SS) seeds)")
println("─"^88)
@printf(
    "  %-6s  %-7s  %-9s  %-5s  %-5s  %s\n",
    "g_syn",
    "g_el",
    "dom_cat",
    "H",
    "n_d",
    "  Q    S    2    3   ≥4  (counts/$(N_SS))"
)
println("  " * "─"^80)

rows = [
    (cat_entropy[gi, ei], gi, ei) for
    gi in eachindex(G_SYN_GRID), ei in eachindex(G_EL_GRID)
]
rows = vec(rows)
sort!(rows; rev=true)   # most multistable first

for (h, gi, ei) in rows
    dom = dom_cat_bin[gi, ei]
    cnts = [round(Int, frac_cat[gi, ei, b + 1]*N_SS) for b in 0:4]
    nd = n_distinct[gi, ei]
    flag = nd > 1 ? " ←" : ""
    @printf(
        "  %-6.2f  %-7.4f  %-9s  %4.2f  %3d   %3d  %3d  %3d  %3d  %3d%s\n",
        G_SYN_GRID[gi],
        G_EL_GRID[ei],
        cat_bin_label(dom),
        h,
        nd,
        cnts[1],
        cnts[2],
        cnts[3],
        cnts[4],
        cnts[5],
        flag
    )
end

println()
n_multi = sum(n_distinct .> 1)
@printf(
    "Multistable cells (n_distinct > 1): %d / %d  (%.0f%%)\n",
    n_multi,
    N_CELLS,
    100.0*n_multi/N_CELLS
)
n_mono = sum(n_distinct .== 1)
@printf(
    "Monostable cells (n_distinct = 1):  %d / %d  (%.0f%%)\n\n",
    n_mono,
    N_CELLS,
    100.0*n_mono/N_CELLS
)

# ETA projection (for full sweep if this was a dry run)
all_walls = [sim_records[k].wall for k in eachindex(CELL_LIST)]
med_wall = median(all_walls)
seq_sum = sum(all_walls)
obs_speedup = seq_sum / t_total

n_full = 11 * 16 * 11   # full publication sweep: 11 g_syn × 16 g_el × 11 seeds = 1936
eta_seq_h = med_wall * n_full / 3600
eta_perf_min = med_wall * n_full / Threads.nthreads() / 60
eta_obs_min = med_wall * n_full / obs_speedup / 60

@printf(
    "Per-sim wall-time:  mean=%.1f s  median=%.1f s  min=%.1f s  max=%.1f s\n",
    mean(all_walls),
    med_wall,
    minimum(all_walls),
    maximum(all_walls)
)
@printf(
    "Thread speedup (this run): %.1fx  (sum_seq=%.0f s / wall=%.0f s)\n",
    obs_speedup,
    seq_sum,
    t_total
)
println()
@printf("Full-sweep ETA (%d sims @ median %.1f s/sim):\n", n_full, med_wall)
@printf("  Sequential              : %.1f h\n", eta_seq_h)
@printf(
    "  %2d threads (ideal %2dx)  : %.0f min\n",
    Threads.nthreads(),
    Threads.nthreads(),
    eta_perf_min
)
@printf(
    "  %2d threads (%.1fx obs.)  : %.0f min  (%.1f h)  ← realistic\n\n",
    Threads.nthreads(),
    obs_speedup,
    eta_obs_min,
    eta_obs_min/60
)

# save
mkpath("results/data")
ts = Dates.format(now(), "yyyy-mm-dd_HH-MM")
outfile = "results/data/sweep_driven_$(RUN_TAG)_$(ts).jld2"

jldsave(
    outfile;
    # grid
    g_syn_grid=G_SYN_GRID,
    g_el_grid=G_EL_GRID,
    seeds=SEEDS,
    # per-sim raw (flat, same order as CELL_LIST)
    sim_g_syn=[CELL_LIST[k].g_syn for k in eachindex(CELL_LIST)],
    sim_g_el=[CELL_LIST[k].g_el for k in eachindex(CELL_LIST)],
    sim_seed=[CELL_LIST[k].seed for k in eachindex(CELL_LIST)],
    sim_chi=[sim_records[k].chi for k in eachindex(CELL_LIST)],
    sim_n_cl=[sim_records[k].n_cl for k in eachindex(CELL_LIST)],
    sim_cat_bin=[cat_bin(sim_records[k].n_cl) for k in eachindex(CELL_LIST)],
    sim_is_fp=[sim_records[k].is_fp for k in eachindex(CELL_LIST)],
    sim_wall_s=[sim_records[k].wall for k in eachindex(CELL_LIST)],
    # aggregated matrices  (N_GS × N_GE, seeds collapsed)
    frac_cat=frac_cat,      # N_GS × N_GE × 5  (bins 0:4)
    cat_entropy=cat_entropy,   # N_GS × N_GE  Shannon entropy (bits)
    n_distinct=n_distinct,    # N_GS × N_GE  number of distinct bins
    dom_cat_bin=dom_cat_bin,   # N_GS × N_GE  dominant bin (0:4)
    mean_chi=mean_chi,
    std_chi=std_chi,
    mean_wall=mean_wall,
    # run metadata
    N_neurons=N_NEURONS,
    N_adj_clust=N_ADJ_CLUSTERS,
    g_drive=G_DRIVE_LOCKED,
    reltol=RELTOL,
    abstol=ABSTOL,
    n_cycles=N_CYCLES,
    saveat=SAVEAT,
    T_drive=T_DRIVE,
    dry_run=DRY_RUN,
    run_tag=RUN_TAG,
    julia_threads=Threads.nthreads(),
    total_wall_s=t_total,
    timestamp=string(now()),
)

println("Saved → $outfile")
println("Plot:   julia --project scripts/plot_rg_driven.jl \"$outfile\"")
