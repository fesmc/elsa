# The 3D Greenland benchmark: isochrone depth below the surface, as a map and as
# a transect. The transect is the view the isochronal scheme exists to produce --
# the layer stack the model would compare against radiostratigraphy.

using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "elsa_analysis.jl"))

function plot_greenland(; path = joinpath(@__DIR__, "..", "output", "GRL-16KM", "elsa.nc"),
                          out  = joinpath(@__DIR__, "..", "plots", "greenland_isochrones.png"))

    r  = ElsaRun(path)
    it = length(r.time)
    t  = r.time[it]

    H    = field2d(r, "H_ice", it)
    dsum = field(r, "dsum_iso", it)

    ice = H .> 1.0

    # The oldest isochrone present: the one laid down at time_init.
    ks     = [k for k in isochrone_layers(r) if k >= 2 && r.layer_time[k] <= t]
    k_old  = first(ks)
    age    = t - r.layer_time[k_old]

    # Depth below the ice surface, which is what radar sees.
    depth = H .- dsum[:, :, k_old-1]
    depth[.!ice] .= NaN

    # Transect across the ice sheet at the row with the most ice.
    j_cut = argmax(vec(sum(ice, dims = 1)))

    fig = Figure(size = (1000, 430))

    ax1 = Axis(fig[1, 1], aspect = DataAspect(),
               title  = @sprintf("Depth of the %.0f yr isochrone", age),
               xlabel = "x (km)", ylabel = "y (km)")
    hm = heatmap!(ax1, r.x ./ 1e3, r.y ./ 1e3, depth,
                  colormap = :dense, nan_color = :transparent)
    Colorbar(fig[1, 2], hm, label = "depth below surface (m)")
    hlines!(ax1, [r.y[j_cut] / 1e3], color = :orangered, linestyle = :dash)

    ax2 = Axis(fig[1, 3],
               title  = @sprintf("Isochrones along y = %.0f km", r.y[j_cut] / 1e3),
               xlabel = "x (km)", ylabel = "z above bed (m)")

    xkm = r.x ./ 1e3
    surf = copy(H[:, j_cut]);  surf[.!ice[:, j_cut]] .= NaN
    for k in ks
        z = dsum[:, j_cut, k-1]
        z[.!ice[:, j_cut]] .= NaN
        lines!(ax2, xkm, z, color = t - r.layer_time[k],
               colorrange = (0, age), colormap = :viridis)
    end
    lines!(ax2, xkm, surf, color = :black, linewidth = 2)
    Colorbar(fig[1, 4], colorrange = (0, age), colormap = :viridis,
             label = "isochrone age (yr)")

    mkpath(dirname(out))
    save(out, fig)

    # Diagnostics worth asserting on: the layer stack must close onto H, and
    # isochrones must not cross.
    closure = maximum(abs.(dsum[:, :, r.n_top[it]][ice] .- H[ice]) ./ H[ice])
    crossings = 0
    for k in 2:r.n_top[it]
        crossings += count(dsum[:, :, k][ice] .< dsum[:, :, k-1][ice] .- 1e-9)
    end

    close(r)

    @printf("  greenland: %d isochrones, oldest %.0f yr\n", length(ks), age)
    @printf("  greenland: max column closure error = %.2e\n", closure)
    @printf("  greenland: isochrone crossings      = %d\n", crossings)
    println("  wrote $out")

    (; closure, crossings, n_isochrones = length(ks))
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    plot_greenland()
end
