# The Nye ice-divide benchmark: modelled isochrone depth against the analytic
# and the discrete solutions.
#
# elsa must reproduce the *discrete* recursion to roundoff -- that is a statement
# about the code -- and approach the *analytic* Nye profile as dt_coupling -> 0,
# which is a statement about the scheme.

using CairoMakie
using Printf

include(joinpath(@__DIR__, "elsa_analysis.jl"))

# These mirror src/test_column.f90; they are properties of the driver, not of
# elsa, so they are not carried in the output file.
const H_CONST = 3000.0    # [m]
const ACC     = 0.3       # [m/yr]

function plot_column(; path = joinpath(@__DIR__, "..", "output", "column", "elsa.nc"),
                       nml  = joinpath(@__DIR__, "..", "par", "test_column.nml"),
                       out  = joinpath(@__DIR__, "..", "plots", "column_nye.png"))

    dt = nml_get(nml, "column", "dt_coupling")
    r  = ElsaRun(path)

    it = length(r.time)
    ages, z = isochrone_height(r, it; i = 3, j = 3)

    z_exact = nye_discrete(ages, H_CONST, ACC, dt)
    z_nye   = nye(ages, H_CONST, ACC)

    err_discrete = maximum(abs.(z .- z_exact))
    err_nye      = maximum(abs.(z .- z_nye))

    fig = Figure(size = (900, 380))

    ax1 = Axis(fig[1, 1],
               title  = "Isochrone height above the bed",
               xlabel = "age (yr)", ylabel = "z (m)")

    τ = range(0, maximum(ages); length = 400)
    lines!(ax1, τ, nye(τ, H_CONST, ACC), color = :black,
           label = "Nye  H exp(-a τ/H)")
    lines!(ax1, τ, nye_discrete(τ, H_CONST, ACC, dt), color = :orangered,
           linestyle = :dash, label = @sprintf("discrete, dt = %g yr", dt))
    scatter!(ax1, ages, z, color = :dodgerblue, markersize = 9, label = "elsa")
    axislegend(ax1, position = :rt, framevisible = false)

    ax2 = Axis(fig[1, 2],
               title  = "elsa − reference",
               xlabel = "age (yr)", ylabel = "Δz (m)")
    scatterlines!(ax2, ages, z .- z_exact, color = :orangered,
                  label = @sprintf("vs discrete (max %.1e m)", err_discrete))
    scatterlines!(ax2, ages, z .- z_nye, color = :black,
                  label = @sprintf("vs Nye (max %.2f m)", err_nye))
    axislegend(ax2, position = :lb, framevisible = false)

    mkpath(dirname(out))
    save(out, fig)
    close(r)

    @printf("  column: max |elsa - discrete| = %.3e m\n", err_discrete)
    @printf("  column: max |elsa - Nye|      = %.3f m   (dt = %g yr)\n", err_nye, dt)
    println("  wrote $out")

    (; err_discrete, err_nye, dt)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    plot_column()
end
