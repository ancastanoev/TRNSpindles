# setup.jl — one-time environment setup.
# Run with:  julia setup.jl
#
# Installs the exact dependency set recorded in Project.toml / Manifest.toml
# and precompiles it. Equivalent to `make install`.

using Pkg

Pkg.activate(@__DIR__)

println("Instantiating TRNSpindles environment from Project.toml / Manifest.toml...")
Pkg.instantiate()

println("\nPrecompiling...")
Pkg.precompile()

println("\nSetup complete. Run tests with:  julia --project test/runtests.jl")
