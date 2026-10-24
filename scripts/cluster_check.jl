# cluster_check.jl — fast pre-flight for the gap-junction cluster sweep.
#
# KEY POINT: a symmetric all-to-all network started from IDENTICAL initial
# conditions stays on the synchronized manifold forever (chi=1, 1 cluster) no
# matter what g_syn/g_el are. Clusters require HETEROGENEOUS initial conditions.
# The sweep (sweep_gsyn_gel) already randomizes ICs; this check does too, and
# shows a homogeneous control so the contrast is explicit.
#
#   julia --project scripts/cluster_check.jl

using TRNSpindles, Random, Printf, LinearAlgebra

function main()
    N = 10
    tspan = (0.0, 2000.0)
    A = let M = ones(Float64, N, N);
        M[diagind(M)] .= 0.0;
        M
    end  # all-to-all

    @printf(
        "%-22s %7s %7s %9s %9s %s\n", "case", "g_syn", "g_el", "clusters", "chi", "fixed_pt"
    )
    println("─"^70)

    # heterogeneous ICs (correct): expect clusters at g_el=0, sync as g_el rises
    for (gsyn, gel) in [(0.5, 0.0), (0.5, 0.1), (0.5, 0.4)]
        net = RinzelGolombNetwork(N; g_syn=gsyn, g_el=gel)
        # average over a few random ICs (cluster count can depend on basin)
        for seed in (1, 2, 3)
            rng = MersenneTwister(seed)
            u0 = make_rg_ic(net; rng=rng)          # <-- heterogeneous, per-neuron
            res = solve_rinzel_golomb(net; u0=u0, tspan=tspan, rng=rng)
            cr = detect_clusters(res, net)
            @printf(
                "%-22s %7.2f %7.2f %9d %9.4f %s\n",
                "hetero IC seed=$seed",
                gsyn,
                gel,
                cr.n_clusters,
                cr.pfeuty_chi,
                cr.is_fixed_point ? "yes" : "no"
            )
        end
    end

    # homogeneous IC control: must give chi=1, 1 cluster (sanity of mechanism)
    let net = RinzelGolombNetwork(N; g_syn=0.5, g_el=0.0)
        u0 = vcat(fill(-55.0, N), fill(0.5, N), fill(0.2, N))  # identical neurons
        res = solve_rinzel_golomb(net; u0=u0, tspan=tspan)
        cr = detect_clusters(res, net)
        @printf(
            "%-22s %7.2f %7.2f %9d %9.4f %s\n",
            "homogeneous CONTROL",
            0.5,
            0.0,
            cr.n_clusters,
            cr.pfeuty_chi,
            cr.is_fixed_point ? "yes" : "no"
        )
    end

    # prove the IC actually differs across neurons
    let net = RinzelGolombNetwork(N; g_syn=0.5, g_el=0.0)
        u0 = make_rg_ic(net; rng=MersenneTwister(1))
        V0 = u0[1:N]
        @printf(
            "\nIC spread check: V0 range = [%.2f, %.2f] mV, std = %.2f (should be >0)\n",
            minimum(V0),
            maximum(V0),
            std(V0)
        )
    end

    return println(
        """
Interpretation:
  hetero IC, g_el=0   → expect >1 cluster, chi < 1   (inhibition clusters)
  hetero IC, g_el up  → clusters merge toward 1, chi → 1 (gap junctions synchronize)
  homogeneous CONTROL → 1 cluster, chi = 1.0000  (symmetry never broken)
If hetero cases also give chi=1, the network synchronizes from any IC;
if so, inspect the synaptic kinetics (kr/kf) against trnlib.py.
"""
    )
end

using Statistics: std
main()
