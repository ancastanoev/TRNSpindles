# test/test_kernel_bridge.jl — Step 7: kernel known-answer + Bridge coverage
using Test, TRNSpindles

@testset "Gamma kernel integrates to 1 (known-answer)" begin
    # h_{m,θ}(s) = s^{m-1} e^{-s/θ} / (θ^m Γ(m));  integer m ⇒ Γ(m) = (m-1)!
    kernel(s, m, θ) = s^(m - 1) * exp(-s / θ) / (θ^m * factorial(m - 1))
    for m in (1, 2, 3, 5), θ in (0.5, 1.0, 2.0)
        S = 30.0 * m * θ          # integrate to 30 means (tail negligible)
        ds = (m * θ) / 4000        # fine step relative to the mean
        grid = 0.0:ds:S
        I =
            ds * (
                sum(kernel.(grid, m, θ)) -
                0.5 * (kernel(grid[1], m, θ) + kernel(grid[end], m, θ))
            )
        @test I ≈ 1.0 atol=1e-3    # trapezoid; mean ρ = m·θ as documented
    end
end

@testset "Bridge: extract_wc_drive coverage (network-free)" begin
    d = extract_wc_drive(default_wc_params())
    @test d isa WCDriveResult
    @test all(isfinite, d.R)
    @test all(isfinite, d.D)
    @test all(isfinite, d.D̂)
    # default params sit in the limit-cycle regime the f-gate pins → R oscillates,
    # so the R_amp>1 guard is satisfied (this is WHY extract returns rather than aborting)
    @test (maximum(d.R) - minimum(d.R)) > 1.0
    @test d.P_WC > 0
    @test d.P_WC_pp > 0
    @test d.P_WC ≈ d.P_WC_pp rtol=0.02   # two period estimators agree
    @test isfinite(d.κ)
end
