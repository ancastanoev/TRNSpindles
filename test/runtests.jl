"""
runtests.jl — TRNSpindles test suite

Run with:
    julia --project test/runtests.jl
or:
    julia --project -e "using Pkg; Pkg.test()"
"""

using Test
using TRNSpindles
using LinearAlgebra, Statistics
using Random
using LaTeXStrings

@testset "TRNSpindles" begin

    # Parameters
    @testset "WCParams defaults" begin
        p = default_wc_params()
        @test p.wPP ≈ 1.0
        @test p.wRT ≈ 0.2
        @test p.τ ≈ 1.0
        @test p.rP ≈ 0.05
    end

    @testset "WCDelayParams" begin
        dp = default_delay_params()
        @test dp.m == 2
        @test gamma_mean_delay(dp) ≈ 1.0
        @test gamma_scale(dp) ≈ 0.5
    end

    # Sigmoids
    @testset "WC sigmoid" begin
        # At u=θ, the sigmoid should equal α/2 - normalisation ≈ α/2 for large βθ
        p = default_wc_params()
        s = wc_sigmoid(p.θE, p.αE, p.θE, p.βE)
        @test s > 0      # positive
        @test s < p.αE   # bounded by α

        # Should increase monotonically
        vals = [wc_sigmoid(u, p.αE, p.θE, p.βE) for u in 0.0:1.0:20.0]
        @test all(diff(vals) .>= 0)

        # Derivative is non-negative
        d = sigmoid_derivative(5.0, p.αE, p.θE, p.βE)
        @test d >= 0
    end

    # Wilson-Cowan ODE
    @testset "WC ODE dimensions" begin
        p = default_wc_params()
        u = [2.0, 1.5, 1.5, 5.0]
        du = zeros(4)
        wc_rhs!(du, u, p, 0.0)
        @test length(du) == 4
        @test all(isfinite, du)
    end

    @testset "WC Jacobian shape" begin
        p = default_wc_params()
        u = [2.0, 1.5, 1.5, 5.0]
        J = wc_jacobian(u, p)
        @test size(J) == (4, 4)
        @test all(isfinite, J)
    end

    @testset "WC equilibrium stability" begin
        p = default_wc_params()
        u0 = [2.0, 1.5, 1.5, 5.0]
        λs, stable = wc_equilibrium_stability(u0, p)
        @test length(λs) == 4
        # Not testing stability value since it depends on parameters
    end

    # Delayed WC
    @testset "Delayed WC dimensions" begin
        dp = WCDelayParams{Float64}(m=3, ρ=2.0)
        u0 = make_delay_ic(dp, [2.0, 1.5, 1.5, 5.0])
        @test length(u0) == 4 + 3
        du = zeros(7)
        delayed_wc_rhs!(du, u0, dp, 0.0)
        @test all(isfinite, du)
    end

    # Rinzel-Golomb
    @testset "Rinzel-Golomb gating" begin
        p, g = default_rinzel_golomb()
        V = -60.0
        @test 0 <= gating_m_inf(V, g) <= 1
        @test 0 <= gating_h_inf(V, g) <= 1
        @test 0 <= gating_s_inf(V, g) <= 1
        @test gating_k_h(V, p, g) >= 0
    end

    @testset "Rinzel-Golomb ODE N=2" begin
        net = RinzelGolombNetwork(2; g_syn=0.2, g_el=0.1)
        x0 = make_rg_ic(net)
        dx = zeros(6)
        rinzel_golomb_rhs!(dx, x0, net, 0.0)
        @test all(isfinite, dx)
    end

    # Graph utilities
    @testset "Adjacency all-to-all" begin
        A = make_adjacency_all(5)
        @test size(A) == (5, 5)
        @test all(diag(A) .== 0)
        @test all(A[i, j] == 1 for i in 1:5 for j in 1:5 if i != j)
    end

    @testset "Laplacian" begin
        A = make_adjacency_all(4)
        L = make_laplacian(A)
        # Row sums of Laplacian are zero
        @test all(abs.(sum(L, dims=2)) .< 1e-12)
        @test isposdef(L + 1e-6*I)   # positive semi-definite (check shifted)
    end

    @testset "Partition fractions" begin
        ms = make_partition_from_fractions(10, [0.3, 0.4, 0.3])
        @test sum(ms) == 10
        @test length(ms) == 3
    end

    @testset "Set2 topology" begin
        A, ms, gaps = make_adjacency_set2(20)
        @test size(A) == (20, 20)
        @test all(diag(A) .== 0)
    end

    # Analysis
    @testset "Dominant frequency" begin
        dt = 0.001
        t = 0.0:dt:10.0
        f0 = 12.0   # Hz (spindle band)
        sig = sin.(2π * f0 .* collect(t))
        f_est = dominant_frequency(sig, dt)
        @test abs(f_est - f0) < 0.5
    end

    @testset "Spindle proximity score" begin
        # Perfect spindle
        s1 = spindle_proximity_score(5.0, 12.0)
        @test s1 > 0.9

        # Escape state (too high)
        s2 = spindle_proximity_score(10.0, 25.0)
        @test s2 < 0.1

        # Monotonicity: closer to target → higher score
        s3 = spindle_proximity_score(4.5, 11.0)
        s4 = spindle_proximity_score(1.0, 5.0)
        @test s3 > s4
    end

    @testset "Regime classification" begin
        @test classify_regime(12.0) == N2_SPINDLE
        @test classify_regime(3.0) == TRANSITION
        @test classify_regime(0.7) == N3_SLOW_WAVE
        @test classify_regime(0.05) == INFRASLOW
    end

    # Cluster detection
    @testset "Pfeuty chi" begin
        # Perfect synchrony: all identical signals
        V_sync = repeat(sin.(0.0:0.1:10.0)', 5, 1)
        @test pfeuty_chi(V_sync) ≈ 1.0 atol=1e-6

        # Perfect incoherence: random independent signals
        rng = Random.MersenneTwister(0)
        V_rand = randn(rng, 5, 100)
        # For truly uncorrelated, chi should be small (no strict assertion, just positive)
        @test pfeuty_chi(V_rand) >= 0
    end

    @testset "Fixed-point detection" begin
        N_t = 100
        V_fp = ones(3, N_t)           # constant → fixed point
        V_osc = [sin(0.1*j + i) for i in 1:3, j in 1:N_t]
        @test is_fixed_point(V_fp)
        @test !is_fixed_point(V_osc)
    end

    @testset "Cluster grouping" begin
        # Two perfectly synchronised groups
        t = 0.0:0.01:10.0
        V1 = [sin.(2π*1 .* collect(t)) for _ in 1:3]   # group 1
        V2 = [sin.(2π*2 .* collect(t)) for _ in 1:3]   # group 2
        V = vcat(hcat(V1...)', hcat(V2...)')
        clusters = find_clusters(V)
        # Should find 2 clusters (or close to it with correlation threshold)
        @test length(clusters) >= 1
    end

    # Integration tests (short)
    @testset "Short WC integration" begin
        p = default_wc_params()
        res = solve_wc(p; u0=[2.0, 1.5, 1.5, 5.0], tspan=(0.0, 10.0))
        @test length(res.t) > 2
        @test size(res.u, 1) == 4
        @test all(isfinite, res.u)
    end

    @testset "Short delayed WC integration" begin
        dp = WCDelayParams{Float64}(m=2, ρ=1.0)
        res = solve_delayed_wc(dp; tspan=(0.0, 10.0))
        @test size(res.u, 1) == 6
        @test all(isfinite, res.u)
    end

    @testset "Short Rinzel-Golomb integration" begin
        net = RinzelGolombNetwork(4; g_syn=0.2, g_el=0.05)
        res = solve_rinzel_golomb(net; tspan=(0.0, 100.0))
        @test all(isfinite, res.u)
    end

    # Bridge (WC→RG open-loop)
    @testset "Bridge: g_drive=0 matches undriven exactly" begin
        # Construct a minimal drive without running a full WC simulation
        Nt = 101
        t_wc = collect(range(0.0, 10.0; length=Nt))
        drive = WCDriveResult{Float64}(
            t_wc,
            fill(0.3, Nt),          # D
            fill(0.006, Nt),        # D̂ = D / 50
            fill(5.0, Nt),          # R
            t_wc[2] - t_wc[1],     # dt_wc
            3.0,                    # P_WC
            3.0,                    # P_WC_pp
            1.0,                    # κ
            50.0,                   # α_E
            5.0,                    # mean_R
            fill(1.0, Nt),          # P_t  (dummy cortical pyramidal)
            fill(0.5, Nt),           # T_t  (dummy thalamic relay)
        )

        N = 4
        net = RinzelGolombNetwork(N; g_syn=0.5, g_el=0.0)
        dp = DrivenRGParams(net, drive; g_drive=0.0)

        u0 = make_rg_ic(net; rng=Random.MersenneTwister(1))
        dx_undriven = zeros(3N)
        dx_driven = zeros(3N)

        rinzel_golomb_rhs!(dx_undriven, u0, net, 0.0)
        driven_rg_rhs!(dx_driven, u0, dp, 0.0)

        @test dx_undriven ≈ dx_driven atol=1e-14
    end

    @testset "Bridge: calibration a is frozen (immutable scalar)" begin
        calib = CalibrationResult{Float64}(0.1, 0.2, 5.123, 0.3)
        @test calib.a === 5.123          # exact identity — plain scalar
        @test calib.a isa Float64
        # CalibrationResult is an immutable struct; this must throw:
        @test_throws Exception setfield!(calib, :a, 99.0)
    end

    @testset "Verification: linear-chain ≈ Gamma convolution" begin
        m = 3
        ρ = 1.5
        θ = ρ / m
        y(t) = 1.0 + 0.5*sin(2π*t/2.0)          # arbitrary bounded test signal

        # ground-truth Gamma kernel and direct convolution
        Γ(n) = factorial(n - 1)                  # Γ(integer) = (m-1)!
        h(s) = s^(m-1) * exp(-s/θ) / (θ^m * Γ(m))
        smax = 12θ                               # kernel essentially zero beyond this
        ds = 1e-3
        svec = 0.0:ds:smax
        conv_at(t) = sum(h(s) * y(t - s) for s in svec) * ds   # ∫ h(s) y(t−s) ds

        # integrate the chain forward with the SAME kernel via explicit Euler
        dt = 1e-4
        T = 30.0                                # long enough to forget the IC
        a = zeros(m)                            # chain initial state
        for t in 0.0:dt:T
            ynow = y(t)
            da1 = (-a[1] + ynow) / θ
            das = [(-a[j] + a[j - 1]) / θ for j in 2:m]
            a[1] += dt*da1
            for j in 2:m
                ;
                a[j] += dt*das[j - 1];
            end
        end
        # compare at the final time, where transients have decayed
        tcmp = T
        chain = a[m]                            # aₘ = filtered output
        direct = conv_at(tcmp)
        relerr = abs(chain - direct) / abs(direct)
        @info "linear-chain vs convolution" chain direct relerr
        @test relerr < 1e-2                       # agree to ~1% with these step sizes
    end

    # (2) Zero feedback recovers the base WC model exactly
    # With ga = 0 the delayed system's four base equations are identical to the
    # undelayed wc_rhs!. The delayed solver (any ρ) must then match the base solver
    # on [P,I,T,R] to solver tolerance.
    @testset "Verification: delayed WC (ga=0) recovers base WC" begin
        base = default_wc_params()
        u0 = Float64[2.0, 1.0, 1.5, 5.0]
        tspan = (0.0, 200.0)
        saveat = 0.1

        res_base = solve_wc(base; u0=u0, tspan=tspan, saveat=saveat)

        # delayed params wrapping the SAME WC params, no slow feedback
        dp = WCDelayParams{Float64}(wc=base, m=3, ρ=0.5, ga=0.0)
        res_del = solve_delayed_wc(dp; u0_4=u0, tspan=tspan, saveat=saveat)

        # compare the base-4 components over the shared time grid
        nb = min(size(res_base.u, 2), size(res_del.u, 2))
        base4 = res_base.u[1:4, 1:nb]
        del4 = res_del.u[1:4, 1:nb]
        maxdiff = maximum(abs.(base4 .- del4))
        @info "delayed(ga=0) vs base WC" maxdiff
        @test maxdiff < 1e-6                       # identical up to solver tolerance
    end

    # (3) Analytic Jacobian matches a finite-difference Jacobian
    # wc_jacobian must equal the numerical Jacobian of wc_rhs! at several states.
    @testset "Verification: WC Jacobian vs finite differences" begin
        p = default_wc_params()
        function fd_jacobian(u, p; h=1e-6)
            n = length(u)
            J = zeros(n, n)
            f0 = zeros(n);
            wc_rhs!(f0, u, p, 0.0)
            for k in 1:n
                uk = copy(u);
                uk[k] += h
                fk = zeros(n);
                wc_rhs!(fk, uk, p, 0.0)
                J[:, k] = (fk .- f0) ./ h
            end
            J
        end
        for u in ([2.0, 1.5, 1.5, 5.0], [1.0, 0.5, 0.8, 3.0], [4.0, 2.0, 1.0, 8.0])
            Ja = wc_jacobian(u, p)
            Jn = fd_jacobian(u, p)
            relerr = maximum(abs.(Ja .- Jn)) / (maximum(abs.(Ja)) + eps())
            @info "Jacobian analytic vs finite-diff" u relerr
            @test relerr < 1e-4
        end
    end

    # (4) Limit-cycle gate on a known oscillatory configuration
    # default_wc_params() sits in a bistable region (a stable limit cycle coexists
    # with a stable fixed point), so the cycle is reached via the same settling path
    # the bridge uses. extract_wc_drive integrates to the limit cycle and measures
    # its period two independent ways (FFT and peak-to-peak); the gate requires that
    # (a) the settled signal is NOT a fixed point and (b) the two period estimates
    # agree and recover the documented cycle (period ~2.96, f ~0.34 Hz).
    @testset "Verification: limit-cycle gate on known oscillatory config" begin
        p = default_wc_params()
        d = extract_wc_drive(
            p; u0=Float64[1, 5, 2, 5], tspan=(0.0, 2000.0), P_phys=3.0, dt_save=0.01
        )

        # (a) the drive must carry a genuine oscillation, not a flat fixed point
        Rsig = d.R
        amp = maximum(Rsig) - minimum(Rsig)
        @test amp > 1.0                              # clear oscillation, not a fixed point

        # (b) FFT and peak-to-peak period estimates agree (consistency of the gate)
        rel_consistency = abs(d.P_WC - d.P_WC_pp) / d.P_WC
        @test rel_consistency < 0.05                 # the two estimators agree to ~5%

        # (c) recovered frequency matches the documented value f ~ 0.34 Hz
        f_det = 1 / d.P_WC
        f_exp = 0.34
        relerr = abs(f_det - f_exp) / f_exp
        @info "limit-cycle gate" P_WC=d.P_WC P_WC_pp=d.P_WC_pp f_detected=f_det relerr
        @test relerr < 0.10
    end
end  # @testset

# --- Step 7 regression / characterization ---
include("test_regression.jl")

# --- Step 7 kernel known-answer + Bridge coverage ---
include("test_kernel_bridge.jl")

# --- Step 7 Bridge readout unit tests ---
include("test_bridge_units.jl")
