# Makefile for TRNSpindles
# Usage: make <target>

JULIA := julia --project
SCRIPTS := scripts

.PHONY: all install test clean baseline delay sweep docs

# ── meta targets ─────────────────────────────────────────────────────────────
all: install test baseline delay

install:
	$(JULIA) -e "using Pkg; Pkg.instantiate()"

test:
	$(JULIA) test/runtests.jl

# ── individual scripts ────────────────────────────────────────────────────────
baseline:
	$(JULIA) $(SCRIPTS)/wc_baseline.jl

delay:
	$(JULIA) $(SCRIPTS)/delay_sweep.jl

sweep:
	$(JULIA) $(SCRIPTS)/sweep_rg_driven.jl

# ── quick smoke tests (short time spans) ─────────────────────────────────────
smoke:
	WC_T=50 $(JULIA) $(SCRIPTS)/wc_baseline.jl
	$(JULIA) -e "\
	  using TRNSpindles; \
	  res = solve_wc(default_wc_params(); tspan=(0.0,20.0)); \
	  @assert all(isfinite, res.u) \"NaN/Inf in WC output!\" ; \
	  println(\"Smoke: WC OK\")"

# ── clean generated files ─────────────────────────────────────────────────────
clean:
	find results/data -name '*.jld2' ! -name 'sweep_driven_publication_2026-06-10_03-46.jld2' -delete
	rm -f results/figures/*.png results/figures/*.pdf results/*.png results/*.pdf
	rm -f *.so

docs:
	@echo "See docs/GUIDE.md for the full project build guide."
	@cat docs/GUIDE.md | head -60
