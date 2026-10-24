"""
    GraphUtils.jl

Adjacency-matrix and graph-Laplacian utilities for the TRN gap-junction topology.

All functions mirror the adjacency construction in
`trnlib.py` (alexfikl/2025-trn-gap-junctions-experiments).

Graph Laplacian convention
──────────────────────────
L = D − A    where  Dᵢᵢ = Σⱼ Aᵢⱼ  (degree matrix)

The gap-junction current for sector i is then:

    I_gap,i = −g_gap · (L R)ᵢ = −g_gap · Σⱼ Aᵢⱼ (Rᵢ − Rⱼ)
"""

"""
    make_laplacian(A) → Matrix

Combinatorial graph Laplacian L = D − A of an adjacency matrix `A`.
"""
function make_laplacian(A::AbstractMatrix{F}) where {F<:AbstractFloat}
    D = Diagonal(vec(sum(A; dims=2)))
    return Matrix(D) - A
end

"""
    make_laplacian_sparse(A) → SparseMatrixCSC

Sparse version — preferred for large N.
"""
function make_laplacian_sparse(A::AbstractMatrix{F}) where {F<:AbstractFloat}
    D = spdiagm(0 => vec(sum(A; dims=2)))
    return D - sparse(A)
end

# Standard adjacency constructors

"""
    make_adjacency_all(N; dtype=Float64) → Matrix

All-to-all gap-junction adjacency matrix for N sectors (no self-coupling).
"""
function make_adjacency_all(N::Int; dtype::Type{F}=Float64) where {F<:AbstractFloat}
    A = ones(F, N, N)
    fill!(view(A, diagind(A)), zero(F))
    return A
end

"""
    make_partition_from_fractions(N_total, ps) → Vector{Int}

Partition N_total neurons into groups with sizes proportional to fractions `ps`.
Rounding is handled by distributing remainder to the groups with largest fractional parts.

Mirrors `make_partition_from_fractions` in trnlib.py.
"""
function make_partition_from_fractions(
    N_total::Int, ps::AbstractVector{F}
) where {F<:AbstractFloat}
    @assert abs(sum(ps) - 1) < 1e-10 "Fractions must sum to 1, got $(sum(ps))"
    weights = ps .* N_total
    ms = floor.(Int, weights)
    leftover = N_total - sum(ms)
    fracs = weights .- ms
    order = sortperm(fracs; rev=true)
    for k in 1:leftover
        ms[order[k]] += 1
    end
    return ms
end

"""
    make_adjacency_groups(ms, gaps; dtype=Float64) → (A, ms, gaps)

Block-diagonal adjacency for gap-junction clusters of sizes `ms` separated
by `gaps` isolated neurons.  Within each block neurons are all-to-all connected.

Returns `(A, ms, gaps)` for downstream bookkeeping.
"""
function make_adjacency_groups(
    ms::AbstractVector{Int}, gaps::AbstractVector{Int}; dtype::Type{F}=Float64
) where {F<:AbstractFloat}
    @assert length(ms) == length(gaps)
    N = sum(ms) + sum(gaps)
    A = zeros(F, N, N)
    i = 1
    for (m, g) in zip(ms, gaps)
        if m > 1
            A[i:(i + m - 1), i:(i + m - 1)] .= one(F)
        end
        i += m + g
    end
    fill!(view(A, diagind(A)), zero(F))
    return A, ms, gaps
end

# Preset topologies (mirrors trnlib.py sets 1–4)

"""
    make_adjacency_set1(N) → (A, ms, gaps)

All-to-all — set 1 from alexfikl/2025-trn-gap-junctions-experiments.
"""
function make_adjacency_set1(N::Int; dtype::Type{F}=Float64) where {F<:AbstractFloat}
    ms = [N]
    gaps = [0]
    return make_adjacency_groups(ms, gaps; dtype=F)
end

"""
    make_adjacency_set2(N) → (A, ms, gaps)

Two groups (90 % / 10 %) with 2-neuron gaps — set 2.
"""
function make_adjacency_set2(N::Int; dtype::Type{F}=Float64) where {F<:AbstractFloat}
    ps = [0.9, 0.1]
    gaps = [2, 2]
    m = N - sum(gaps)
    ms = make_partition_from_fractions(m, F.(ps))
    return make_adjacency_groups(ms, gaps; dtype=F)
end

"""
    make_adjacency_set3(N) → (A, ms, gaps)

Three groups (27 % / 40 % / 33 %) — set 3.
"""
function make_adjacency_set3(N::Int; dtype::Type{F}=Float64) where {F<:AbstractFloat}
    ps = [0.27, 0.40, 0.33]
    gap0 = max(1, N ÷ 15)
    gaps = fill(gap0, 3)
    m = N - sum(gaps)
    ms = make_partition_from_fractions(m, F.(ps))
    return make_adjacency_groups(ms, gaps; dtype=F)
end

"""
    make_adjacency_set4(N) → (A, ms, gaps)

Four groups (29 % / 25 % / 33 % / 13 %) — set 4.
"""
function make_adjacency_set4(N::Int; dtype::Type{F}=Float64) where {F<:AbstractFloat}
    ps = [0.29, 0.25, 0.33, 0.13]
    gap0 = max(1, N ÷ 15)
    gaps = fill(gap0, 4)
    m = N - sum(gaps)
    ms = make_partition_from_fractions(m, F.(ps))
    return make_adjacency_groups(ms, gaps; dtype=F)
end

"""
    make_adjacency_random_clusters(rng, N, n_clusters; gap=1) → (A, ms, gaps)

Generate a random cluster structure with `n_clusters` groups of roughly equal size.
"""
function make_adjacency_random_clusters(
    rng::AbstractRNG, N::Int, n_clusters::Int; gap::Int=1, dtype::Type{F}=Float64
) where {F<:AbstractFloat}
    gaps = fill(gap, n_clusters)
    m = N - sum(gaps)
    ps = ones(F, n_clusters) ./ n_clusters
    ms = make_partition_from_fractions(m, ps)
    return make_adjacency_groups(ms, gaps; dtype=F)
end
