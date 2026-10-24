# test/test_regression.jl — Step 7 regression / characterization
#  (1) WC-core regression: CI-portable; pins solve_wc summary stats vs committed golden.
#  (2) Publication characterization: local only (isfile-guarded); pins published
#      driven-sweep arrays vs sha256 manifest (excludes volatile timing keys).
using Test, TRNSpindles, Statistics, JLD2, SHA, Printf

const GOLDEN = joinpath(@__DIR__, "golden")

function _read_kv(path)
    d = Dict{String,String}()
    for ln in eachline(path)
        (isempty(ln) || startswith(ln, "#")) && continue
        p = split(ln; limit=2)
        length(p) == 2 && (d[p[1]] = p[2])
    end
    return d
end

@testset "WC-core regression (vs golden)" begin
    ref = _read_kv(joinpath(GOLDEN, "wc_core_reference.txt"))
    @test !isempty(ref)
    r = if haskey(ref, "saveat")
        solve_wc(default_wc_params(); tspan=(0.0, 200.0), saveat=parse(Float64, ref["saveat"]))
    else
        solve_wc(default_wc_params(); tspan=(0.0, 200.0))
    end
    u = r.u
    popdim = parse(Int, ref["popdim"])
    finals = popdim == 1 ? u[:, end] : u[end, :]
    means = popdim == 1 ? vec(mean(u; dims=2)) : vec(mean(u; dims=1))
    stds = popdim == 1 ? vec(std(u; dims=2)) : vec(std(u; dims=1))
    for i in eachindex(finals)
        @test finals[i] ≈ parse(Float64, ref["final_$i"]) rtol=1e-4
    end
    for i in eachindex(means)
        @test means[i] ≈ parse(Float64, ref["mean_$i"]) rtol=1e-3
        @test stds[i] ≈ parse(Float64, ref["std_$i"]) rtol=1e-3
    end
end

@testset "Publication sweep characterization (local, guarded)" begin
    pubfile = joinpath(
        @__DIR__, "..", "results", "data", "sweep_driven_publication_2026-06-10_03-46.jld2"
    )
    manifest = joinpath(GOLDEN, "sweep_publication_manifest.txt")
    if !isfile(pubfile) || !isfile(manifest)
        @info "publication file or manifest absent — skipping (expected on CI)"
        @test_skip true
    else
        d = load(pubfile)
        for ln in eachline(manifest)
            (isempty(ln) || startswith(ln, "#")) && continue
            f = split(ln)
            key, kind = f[1], f[2]
            @test haskey(d, key)
            v = d[key]
            if kind == "sha256"
                @test bytes2hex(sha256(reinterpret(UInt8, vec(float.(v))))) == f[3]
            elseif kind == "value"
                @test string(v) == f[3]
            end
        end
    end
end
