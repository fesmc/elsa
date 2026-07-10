# Regenerate the committed documentation figures.
#
# The figures on the docs site are the benchmark output figures, drawn by the
# plot functions in `analysis/`. This script points them at `docs/figures/` so
# the same code that backs `make validate` also produces the committed PNGs --
# there is no second copy of the plotting logic.
#
# It needs the benchmark output, so run the benchmarks first:
#
#     make check
#     julia --project=analysis docs/figures/make_figures.jl
#
# The rendered site (CI) never runs this: the PNGs are committed, so building the
# docs needs only Quarto.

const ANALYSIS = joinpath(@__DIR__, "..", "..", "analysis")
const HERE     = @__DIR__

include(joinpath(ANALYSIS, "plot_column.jl"))
include(joinpath(ANALYSIS, "plot_greenland.jl"))

plot_column(out    = joinpath(HERE, "column_nye.png"))
plot_greenland(out = joinpath(HERE, "greenland_isochrones.png"))

println("\n  documentation figures written to docs/figures/")
