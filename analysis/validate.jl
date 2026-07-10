# Quantitative validation of the benchmark output, and the figures.
#
# The division of labour: the Fortran benchmarks assert structural properties and
# exit nonzero (`make check`); this asserts the physics against closed-form
# answers and draws (`make validate`). Run `make check` first -- these read the
# NetCDF the benchmarks write.

using Printf

include(joinpath(@__DIR__, "plot_column.jl"))
include(joinpath(@__DIR__, "plot_greenland.jl"))

const H_COLUMN = 3000.0

n_fail = 0

function check(ok, name)
    global n_fail
    if ok
        println("   pass   $name")
    else
        println("   FAIL   $name")
        n_fail += 1
    end
end

println("\n== elsa validation ==\n")

println(" column (Nye)")
col = plot_column()

# elsa lands on the discrete recursion exactly -- but `dsum_iso` is written single
# precision, so the floor here is float32 epsilon at 3 km of ice, about 2e-4 m.
# The bit-level statement (4.9e-15 relative, on the in-memory doubles) is asserted
# in test_column.x; this only has to catch an isochrone moved by a real error,
# which would be metres.
check(col.err_discrete < 1e-6 * H_COLUMN, "matches the discrete solution")

# And approach Nye to first order in dt. At dt = 100 yr and an age of 20 kyr the
# leading term is (a*dt/2H) * (a*τ/H) * z ~ 5 m; much larger means the thinning is
# not first-order accurate.
check(col.err_nye < 0.01 * H_COLUMN, "approaches the Nye solution")

println("\n greenland (16 km)")
grl = plot_greenland()

check(grl.closure < 1e-9,     "columns close onto H_ice")
check(grl.crossings == 0,     "isochrones never cross")
check(grl.n_isochrones > 0,   "isochrones were laid down")

println()
if n_fail > 0
    println("  $n_fail check(s) FAILED\n")
    exit(1)
end
println("  all checks passed\n")
