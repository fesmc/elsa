# Shared helpers for reading elsa output.
#
# elsa's vertical axis is time. Layer `k` is bounded below by the isochrone laid
# down at `layer_time[k]`, so the height of that isochrone above the bed is
# `dsum_iso[:,:,k-1]`. The initialization layers carry a missing `layer_time`:
# they are not isochrones of anything.

using NCDatasets
using Printf

"""
    ElsaRun

One elsa output file, with the time axis loaded and the layer bookkeeping
resolved into a list of isochrones.
"""
struct ElsaRun
    x::Vector{Float64}          # [m]
    y::Vector{Float64}          # [m]
    time::Vector{Float64}       # [yr]
    layer_time::Vector{Float64} # [yr] NaN for the initialization layers
    n_top::Vector{Int}
    ds::NCDataset
end

function ElsaRun(path::AbstractString)
    ds = NCDataset(path, "r")
    lt = Array{Float64}(coalesce.(ds["layer_time"][:], NaN))
    # ncio writes the missing value as -9999, not as a _FillValue attribute.
    lt[lt .< -9998.0] .= NaN
    ElsaRun(Array{Float64}(ds["xc"][:]),
            Array{Float64}(ds["yc"][:]),
            Array{Float64}(ds["time"][:]),
            lt,
            Array{Int}(ds["n_top"][:]),
            ds)
end

Base.close(r::ElsaRun) = close(r.ds)

"Layer indices that have an isochrone at their base, i.e. not initialization layers."
isochrone_layers(r::ElsaRun) = findall(!isnan, r.layer_time)

"Field `name` at time index `it`, as (x, y, layer) or (x, y)."
field(r::ElsaRun, name, it) = Array{Float64}(r.ds[name][:, :, :, it])
field2d(r::ElsaRun, name, it) = Array{Float64}(r.ds[name][:, :, it])

"""
    isochrone_height(r, it) -> (ages, z)

Height above the bed of every isochrone present at time index `it`, and its age
at that time. The isochrone at the base of layer `k` sits at `dsum_iso[k-1]`.
"""
function isochrone_height(r::ElsaRun, it::Int; i=1, j=1)
    dsum = field(r, "dsum_iso", it)
    t = r.time[it]
    ks = [k for k in isochrone_layers(r) if k >= 2 && r.layer_time[k] <= t]
    ages = [t - r.layer_time[k] for k in ks]
    z = [dsum[i, j, k-1] for k in ks]
    ages, z
end

"Nye's steady ice-divide solution: an isochrone of age `τ` sits at `H exp(-a τ / H)`."
nye(τ, H, a) = H .* exp.(-a .* τ ./ H)

"""
elsa's *discrete* thinning: each coupling step multiplies every layer height by
`H/(H + a dt)`, so after `τ/dt` steps the isochrone is at `H (1 + a dt/H)^(-τ/dt)`.
This tends to `nye` as `dt -> 0`, and it is what elsa must reproduce exactly.
"""
nye_discrete(τ, H, a, dt) = H .* (1 + a * dt / H) .^ (-τ ./ dt)

"Read one scalar from a Fortran namelist group. Cheap, and enough for the benchmarks."
function nml_get(path::AbstractString, group::AbstractString, key::AbstractString)
    in_group = false
    for line in eachline(path)
        s = strip(line)
        startswith(s, "!") && continue
        if startswith(s, "&")
            in_group = (lowercase(s[2:end]) == lowercase(group))
            continue
        end
        s == "/" && (in_group = false)
        if in_group && occursin("=", s)
            k, v = split(s, "=", limit=2)
            if lowercase(strip(k)) == lowercase(key)
                v = strip(first(split(v, "!")))
                return parse(Float64, v)
            end
        end
    end
    error("$key not found in group &$group of $path")
end
