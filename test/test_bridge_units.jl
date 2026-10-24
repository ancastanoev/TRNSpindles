# test/test_bridge_units.jl — Step 7: pure unit tests for Bridge readout
using Test, TRNSpindles
const _T = TRNSpindles

@testset "Bridge readout: constant synaptic block" begin
    N, Nt = 4, 12
    u = zeros(3N, Nt);
    s = 0.37
    u[(2N + 1):3N, :] .= s
    res = _T.SimResult{Float64}(collect(1.0:Nt), u, "synthetic", nothing)
    mean_s, mod_dep = _T._mean_sbar(res, N)
    @info "constant-block readout" mean_s mod_dep
    @test mean_s ≈ s rtol=1e-12
    @test mod_dep ≈ 0.0 atol=1e-10
    a = 2.5
    R_eff = _T.compute_R_eff(res, N, a)
    Rv = R_eff isa Number ? [R_eff] : vec(R_eff)
    @test all(x -> isapprox(x, a * s; rtol=1e-10), Rv)
end

@testset "Bridge readout: graded synaptic rows (independent mean)" begin
    N, Nt = 4, 12
    u = zeros(3N, Nt);
    rows = (2N + 1):3N;
    vals = (0.1, 0.2, 0.3, 0.4)
    for (k, r) in enumerate(rows)
        ;
        u[r, :] .= vals[k];
    end
    res = _T.SimResult{Float64}(collect(1.0:Nt), u, "synthetic", nothing)
    mean_s, mod_dep = _T._mean_sbar(res, N)
    @info "graded-rows readout" mean_s mod_dep
    @test mean_s ≈ 0.25 rtol=1e-12
end
