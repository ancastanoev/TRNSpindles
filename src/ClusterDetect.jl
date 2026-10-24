"""
    ClusterDetect.jl

Synchronisation-cluster detection for the Rinzel–Golomb and extended models.

Mirrors the cluster-detection logic of trnlib.py
(alexfikl/2025-trn-gap-junctions-experiments):

1.  **Pfeuty χ** — population-synchrony measure (Pfeuty 2007).
2.  **Fixed-point detection** — is the solution essentially stationary?
3.  **Cluster grouping** — partition neurons into groups with coherent dynamics.

The number of synchronisation clusters is the key output used to
produce the synchronisation-cluster analysis.
"""

# Pfeuty synchrony measure

"""
    pfeuty_chi(V) → Float ∈ [0, 1]

Global synchrony index of Pfeuty et al. (2007).

    χ = Var(⟨V⟩_i) / ⟨Var(Vᵢ)⟩_i

where the outer angle brackets are time averages and inner are neuron averages.
χ = 1 → perfect synchrony;  χ = 0 → fully incoherent.

`V` is a (N × Nt) matrix of voltages.
"""
function pfeuty_chi(V::AbstractMatrix{F})::F where {F<:AbstractFloat}
    N, Nt = size(V)
    N == 1 && return one(F)

    V_mean = vec(mean(V; dims=1))  # population mean at each time (Nt,)
    var_mean = var(V_mean)        # variance of population mean
    mean_vars = mean([var(V[i, :]) for i in 1:N])  # mean of individual variances

    mean_vars < eps(F) && return one(F)  # all neurons identical → synchronised
    return var_mean / mean_vars
end

# Fixed-point detection

"""
    is_fixed_point(V; eps=1e-4) → Bool

True if the voltage matrix `V` (N × Nt) appears to be converged to a
fixed point, i.e. all trajectories are essentially constant in time.
"""
function is_fixed_point(
    V::AbstractMatrix{F}; eps_tol::F=F(1e-4)
)::Bool where {F<:AbstractFloat}
    for i in 1:size(V, 1)
        row = view(V, i, :)
        rng = maximum(row) - minimum(row)
        mx = maximum(abs, row)
        mx < eps(F) && continue
        rng / mx > eps_tol && return false
    end
    return true
end

# Cluster detection via correlation

"""
    find_clusters(V; corr_thresh=0.95) → Vector{Vector{Int}}

Partition neurons into synchronisation clusters based on pairwise
Pearson correlation of their voltage time series.

Two neurons i, j belong to the same cluster if  |corr(Vᵢ, Vⱼ)| ≥ corr_thresh.

Returns a list of index vectors, one per cluster.  Isolated neurons form
singleton clusters.
"""
function find_clusters(
    V::AbstractMatrix{F}; corr_thresh::F=F(0.95)
)::Vector{Vector{Int}} where {F<:AbstractFloat}
    N = size(V, 1)

    # Build correlation matrix
    C = cor(V')   # cor(X) treats columns as observations; V' has neurons as cols

    # Union-Find / greedy grouping
    assigned = fill(-1, N)
    clusters = Vector{Vector{Int}}()
    cluster_id = 0

    for i in 1:N
        assigned[i] != -1 && continue
        cluster_id += 1
        grp = [i]
        assigned[i] = cluster_id
        for j in (i + 1):N
            assigned[j] != -1 && continue
            if abs(C[i, j]) >= corr_thresh
                push!(grp, j)
                assigned[j] = cluster_id
            end
        end
        push!(clusters, grp)
    end
    return clusters
end

# ClusterResult container

"""
    ClusterResult

Summary of the clustering analysis for a single simulation.

Fields
──────
- `n_clusters`      — number of distinct synchronisation clusters (0 = fixed point)
- `is_fixed_point`  — true if the solution converged to a fixed point
- `above_threshold` — true if the voltage is above the firing threshold
- `cluster_indices` — list of neuron-index vectors, one per cluster
- `pfeuty_chi`      — global synchrony index ∈ [0, 1]
- `g_syn`           — synaptic conductance (for sweep bookkeeping)
- `g_el`            — gap-junction conductance
"""
struct ClusterResult{F<:AbstractFloat}
    n_clusters::Int
    is_fixed_point::Bool
    above_threshold::Bool
    cluster_indices::Vector{Vector{Int}}
    pfeuty_chi::F
    g_syn::F
    g_el::F
end

"""
    detect_clusters(res, net; window_frac=0.125, corr_thresh=0.95,
                    V_threshold=-40.0) → ClusterResult

Detect synchronisation clusters in a Rinzel–Golomb simulation.

- `res`   — `SimResult` from `solve_rinzel_golomb`
- `net`   — `RinzelGolombNetwork`
- `window_frac` — fraction of total time used as analysis window (last portion)
"""
function detect_clusters(
    res::SimResult{F},
    net::RinzelGolombNetwork{F};
    window_frac::F=F(0.4),
    corr_thresh::F=F(0.95),
    V_threshold::F=F(-40.0),
) where {F<:AbstractFloat}
    N = net.N
    Nt = size(res.u, 2)
    n_window = max(2, round(Int, window_frac * Nt))

    V = res.u[1:N, (end - n_window + 1):end]   # voltage only, last window

    fp = is_fixed_point(V)
    χ = pfeuty_chi(V)

    if fp
        return ClusterResult{F}(
            0, true, false, [collect(1:N)], χ, net.param.g_syn, net.param.g_el
        )
    end

    clusters = find_clusters(V; corr_thresh=corr_thresh)
    nc = length(clusters)
    above = maximum(V) >= V_threshold

    return ClusterResult{F}(nc, false, above, clusters, χ, net.param.g_syn, net.param.g_el)
end

# Parameter-plane sweep (g_syn × g_el)

"""
    sweep_gsyn_gel(N_neurons, g_syns, g_els, A;
                   tspan, n_runs, rng, kwargs...) → Matrix{ClusterResult}

Sweep the (g_syn, g_el) parameter plane for a fixed gap-junction topology `A`.

For each grid point, `n_runs` simulations with independent random ICs are run,
and the most-frequent cluster count is returned.

This produces the synchronisation-cluster analysis
(alexfikl/2025-trn-gap-junctions-experiments).
"""
function sweep_gsyn_gel(
    N_neurons::Int,
    g_syns::AbstractVector{F},
    g_els::AbstractVector{F},
    A::AbstractMatrix{F};
    tspan::Tuple{F,F}=(F(0), F(2000)),
    n_runs::Int=10,
    rng::AbstractRNG=Random.default_rng(42),
    show_progress::Bool=true,
    kwargs...,
) where {F<:AbstractFloat}
    ng = length(g_syns)
    ne = length(g_els)
    results = Matrix{Vector{ClusterResult{F}}}(undef, ng, ne)

    p_base, g_base = default_rinzel_golomb(F)

    # Derive a base seed so each grid point gets an independent, reproducible RNG.
    base_seed = rand(rng, UInt64)

    prog = show_progress ? Progress(ng * ne; desc="(g_syn × g_el) sweep: ") : nothing

    Threads.@threads for k in 1:(ng * ne)
        i = (k - 1) % ng + 1
        j = (k - 1) ÷ ng + 1
        g_syn_k = g_syns[i]
        g_el_k = g_els[j]

        p_ij = RinzelGolombParams{F}(;
            g_syn=g_syn_k,
            g_el=g_el_k,
            C=p_base.C,
            g_Ca=p_base.g_Ca,
            g_L=p_base.g_L,
            V_Ca=p_base.V_Ca,
            V_L=p_base.V_L,
            V_syn=p_base.V_syn,
            kr=p_base.kr,
            kf=p_base.kf,
            phi=p_base.phi,
        )
        net_ij = RinzelGolombNetwork(p_ij, g_base, A)

        # Each grid point gets its own RNG — no shared state between threads.
        local_rng = Random.MersenneTwister(
            base_seed ⊻ (UInt64(i) * 0x9e3779b97f4a7c15) ⊻ UInt64(j)
        )

        run_results = ClusterResult{F}[]
        for _ in 1:n_runs
            res_ij = solve_rinzel_golomb(net_ij; tspan=tspan, rng=local_rng, kwargs...)
            cr = detect_clusters(res_ij, net_ij)
            push!(run_results, cr)
        end
        results[i, j] = run_results
        isnothing(prog) || next!(prog)
    end
    return results
end

"""
    cluster_count_fraction(results) → (mean_count_matrix, fraction_matrix)

Post-process the output of `sweep_gsyn_gel`:
- `mean_count_matrix[i,j]` — average cluster count over runs at grid point (i,j)
- `fraction_fixed_point[i,j]` — fraction of runs that converged to fixed point
"""
function cluster_count_fraction(results::Matrix{Vector{ClusterResult{F}}}) where {F}
    mean_counts = map(runs -> F(mean([cr.n_clusters for cr in runs])), results)
    frac_fp = map(
        runs -> F(count(cr.is_fixed_point for cr in runs)) / length(runs), results
    )
    return mean_counts, frac_fp
end
