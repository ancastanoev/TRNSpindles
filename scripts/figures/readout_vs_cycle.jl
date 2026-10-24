# readout_vs_cycle.jl — validation overlay: reconstructed R_eff(t) vs mean-field R_WC(t).
#
# Shows how the open-loop read-out R_eff = a·s̄(t) compares to the Wilson--Cowan
# reticular cycle R_WC(t) it is derived from.
#
# Output: vector PDF.
#
# Usage, from the TRNSpindles project root:
#     julia --project scripts/figures/readout_vs_cycle.jl

using TRNSpindles, CairoMakie, Random, SparseArrays, Statistics, Printf

const P_PHYS = 3.0
const N_CYC = 8
const N = 1000
const NCL = 50
const ADJSEED = 42
const SAVEAT = 0.01
const RELTOL = 1e-6
const ABSTOL = 1e-8

const G_SYN = 0.60
const G_EL = 0.06
const G_DRIVE = 1.8
const SEED = 0

const NAVY = RGBf(0.094, 0.212, 0.408)
const BRICK = RGBf(0.72, 0.18, 0.16)

movmean(x, w) = begin
    n = length(x)
    y = similar(x)
    h = w ÷ 2

    @inbounds for i in 1:n
        lo = max(1, i - h)
        hi = min(n, i + h)
        y[i] = sum(@view x[lo:hi]) / (hi - lo + 1)
    end

    y
end

nz(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x) + eps())

function main()

    # Wilson--Cowan drive
    drive = extract_wc_drive(
        default_wc_params();
        u0=Float64[1, 5, 2, 5],
        tspan=(0.0, 2000.0),
        P_phys=P_PHYS,
        dt_save=0.01,
    )

    Tdrive = drive.P_WC * drive.κ

    # network construction
    A = sparse(
        first(make_adjacency_random_clusters(MersenneTwister(ADJSEED), N, NCL; gap=0))
    )

    Pbase, Gb = default_rinzel_golomb(Float64)

    make_net(gsyn, gel) = RinzelGolombNetwork(
        RinzelGolombParams{Float64}(;
            g_syn=gsyn,
            g_el=gel,

            C=Pbase.C,
            g_Ca=Pbase.g_Ca,
            g_L=Pbase.g_L,
            V_Ca=Pbase.V_Ca,
            V_L=Pbase.V_L,
            V_syn=Pbase.V_syn,

            kr=Pbase.kr,
            kf=Pbase.kf,
            phi=Pbase.phi,
        ),
        Gb,
        A,
    )

    run_sbar(net, gd) = begin
        dp = DrivenRGParams(net, drive; g_drive=gd)
        u0 = make_rg_ic(net; rng=MersenneTwister(SEED))

        solve_driven_rg(
            dp;
            u0=u0,
            tspan=(0.0, N_CYC * Tdrive),
            saveat=SAVEAT,
            reltol=RELTOL,
            abstol=ABSTOL,
        )
    end

    # single synchronised, entrained operating point
    net = make_net(G_SYN, G_EL)
    res = run_sbar(net, G_DRIVE)

    nskip = round(Int, 0.5 * size(res.u, 2))
    tw = res.t[(nskip + 1):end]

    w = max(3, round(Int, (Tdrive / 12) / (tw[2] - tw[1])))

    sbar = movmean(vec(mean(res.u[(2N + 1):3N, (nskip + 1):end]; dims=1)), w)

    a = drive.mean_R / mean(sbar)
    Reff = a .* sbar

    interp_RWC(t) = begin
        tw_wc = mod(t / drive.κ, drive.P_WC)
        g = drive.t .- drive.t[1]

        k = clamp(searchsortedlast(g, tw_wc), 1, length(g) - 1)
        f = (tw_wc - g[k]) / (g[k + 1] - g[k])

        (1 - f) * drive.R[k] + f * drive.R[k + 1]
    end

    Rwc = interp_RWC.(tw)

    resid = sqrt(mean((Reff .- Rwc) .^ 2))

    sc = sbar .- mean(sbar)
    rc = Rwc .- mean(Rwc)

    corr = sum(sc .* rc) / (sqrt(sum(sc .^ 2)) * sqrt(sum(rc .^ 2)) + eps())

    fS = dominant_frequency(sbar, tw[2] - tw[1])
    Tenv = fS > 0 ? 1 / fS : Inf

    # lag alignment of R_eff relative to R_WC
    dt = tw[2] - tw[1]
    Lmax = round(Int, Tdrive / dt)
    n = length(Reff)

    af = Reff .- mean(Reff)
    bf = Rwc .- mean(Rwc)

    function corr_at(lag)
        if lag >= 0
            a_seg = af[(1 + lag):n]
            b_seg = bf[1:(n - lag)]
        else
            a_seg = af[1:(n + lag)]
            b_seg = bf[(1 - lag):n]
        end

        return sum(a_seg .* b_seg) / (sqrt(sum(a_seg .^ 2)) * sqrt(sum(b_seg .^ 2)) + eps())
    end

    lags = (-Lmax):Lmax
    cvals = [corr_at(l) for l in lags]
    bestlag = lags[argmax(cvals)]
    corr_alig = maximum(cvals)
    lag_time = bestlag * dt

    @printf(
        "g_syn=%.2f g_el=%.2f g_drive=%.2f  a=%.3f  corr=%.3f  corr_aligned=%.3f  lag=%.3f  period(s̄)=%.3f  rms=%.3f\n",
        G_SYN,
        G_EL,
        G_DRIVE,
        a,
        corr,
        corr_alig,
        lag_time,
        Tenv,
        resid,
    )

    # build aligned time series without boundary clamping
    if bestlag >= 0
        idx_wc = 1:(n - bestlag)
        idx_eff = (1 + bestlag):n
    else
        idx_wc = (1 - bestlag):n
        idx_eff = 1:(n + bestlag)
    end

    tw_a = tw[idx_wc]
    Rwc_a = Rwc[idx_wc]
    Reff_a = Reff[idx_eff]

    # Show several cycles so the agreement reads as systematic.
    avail = (tw_a[end] - tw_a[1]) / Tdrive
    NSHOW = min(6, floor(Int, avail))

    keep = tw_a .>= (tw_a[end] - NSHOW * Tdrive)

    t0 = (tw_a[keep] .- tw_a[keep][1]) ./ Tdrive
    Rwc_n = nz(Rwc_a[keep])
    Reff_n = nz(Reff_a[keep])

    # publication-style figure
    set_theme!(
        Theme(;
            fontsize=10,
            fonts=(; regular="CMU Serif", bold="CMU Serif Bold", italic="CMU Serif Italic"),
            figure_padding=(8, 10, 6, 6),
            Axis=(
                xgridvisible=false,
                ygridvisible=false,
                xminorgridvisible=false,
                yminorgridvisible=false,
                xminorticksvisible=false,
                yminorticksvisible=false,
                topspinevisible=false,
                rightspinevisible=false,
                spinewidth=0.8,
                xtickwidth=0.8,
                ytickwidth=0.8,
                xticksize=4,
                yticksize=4,
                xtickalign=1,
                ytickalign=1,
                xlabelpadding=7,
                ylabelpadding=9,
            ),
            Legend=(
                framevisible=false, labelsize=9, patchsize=(24, 8), rowgap=3, colgap=12
            ),
        ),
    )

    fig = Figure(; size=(680, 315))

    ax = Axis(
        fig[2, 1];
        xlabel=L"\mathrm{model\ time}\;t/T",
        ylabel=L"\mathrm{normalised\ reticular\ activity}",
        xticks=0:NSHOW,
        yticks=([0.0, 0.5, 1.0], [L"0", L"0.5", L"1"]),
    )

    xlims!(ax, 0, NSHOW)
    ylims!(ax, -0.04, 1.04)

    # Subtle cycle separators.
    if NSHOW > 1
        vlines!(ax, 1:(NSHOW - 1); color=(:black, 0.12), linewidth=0.8, linestyle=:dot)
    end

    lines!(ax, t0, Rwc_n; color=NAVY, linewidth=2.35, label=L"R_{\mathrm{WC}}(t)")

    lines!(
        ax,
        t0,
        Reff_n;
        color=BRICK,
        linewidth=2.25,
        linestyle=:dash,
        label=L"R_{\mathrm{eff}}(t)=a\,\bar{s}(t)",
    )

    # Legend outside the plotting area.
    Legend(
        fig[1, 1],
        ax;
        orientation=:horizontal,
        halign=:center,
        valign=:center,
        tellwidth=false,
        nbanks=1,
    )

    rowgap!(fig.layout, 4)

    # save
    mkpath(joinpath(@__DIR__, "..", "..", "results", "figures"))

    out = joinpath(@__DIR__, "..", "..", "results", "figures", "fig9_reff_vs_r.pdf")
    save(out, fig; pt_per_unit=1)

    @printf(
        "saved %s  (g_syn=%.2f, g_drive=%.2f, corr_aligned=%.2f, lag=%.2f, period=%.2f)\n",
        out,
        G_SYN,
        G_DRIVE,
        corr_alig,
        lag_time,
        Tenv,
    )

    set_theme!()

    return fig
end

main()
