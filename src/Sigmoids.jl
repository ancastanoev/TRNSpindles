"""
    Sigmoids.jl

Wilson–Cowan sigmoidal transfer functions and their derivatives.

The standard WC sigmoid used throughout:

    S_{α,θ,β}(u) = α / (1 + exp(−β(u − θ))) − α / (1 + exp(β θ))

The subtracted term ensures S(0) = 0.

Reference: Wilson & Cowan (1972), Eq. 5; Marreiros et al. (2008).
"""

"""
    wc_sigmoid(u, α, θ, β) → Float

Evaluate the normalised Wilson–Cowan sigmoid at input `u`.
"""
@inline function wc_sigmoid(u, α::T, θ::T, β::T) where {T<:AbstractFloat}
    α / (1 + exp(-β * (u - θ))) - α / (1 + exp(β * θ))
end

"""
    wc_sigmoid(u, p, node)

Convenience dispatch.  `node` ∈ {:P, :I, :T, :R}.
"""
function wc_sigmoid(u::T, p::WCParams{T}, node::Symbol)::T where {T<:AbstractFloat}
    if node === :P || node === :T
        return wc_sigmoid(u, p.αE, p.θE, p.βE)
    elseif node === :I
        return wc_sigmoid(u, p.αI1, p.θI1, p.βI1)
    elseif node === :R
        return wc_sigmoid(u, p.αI2, p.θI2, p.βI2)
    else
        throw(ArgumentError("Unknown node: $node.  Must be :P, :I, :T, or :R."))
    end
end

"""
    sigmoid_derivative(u, α, θ, β) → Float

Derivative  dS/du  of the WC sigmoid.
"""
@inline function sigmoid_derivative(u::T, α::T, θ::T, β::T)::T where {T<:AbstractFloat}
    e = exp(-β * (u - θ))
    α * β * e / (1 + e)^2
end

function sigmoid_derivative(u::T, p::WCParams{T}, node::Symbol)::T where {T<:AbstractFloat}
    if node === :P || node === :T
        return sigmoid_derivative(u, p.αE, p.θE, p.βE)
    elseif node === :I
        return sigmoid_derivative(u, p.αI1, p.θI1, p.βI1)
    elseif node === :R
        return sigmoid_derivative(u, p.αI2, p.θI2, p.βI2)
    else
        throw(ArgumentError("Unknown node: $node."))
    end
end
