# elsa

**E**nglacial **L**ayer **S**imulation **A**rchitecture: an isochronal model for
ice-sheet layer tracing.

elsa advects a stack of isochronal layers through an ice-sheet domain. Its
vertical axis is time — each layer is bounded by two isochrones and never
exchanges mass with its neighbours — so vertical numerical diffusion is
eliminated by construction. It is driven by a host ice-sheet model, which
supplies horizontal velocity, ice thickness, and surface and basal mass
balance, and which elsa never modifies.

elsa builds as a static library, `libelsa.a`, and ships stand-alone benchmarks.

The method is described in [Born (2017)](https://doi.org/10.1017/jog.2016.111),
[Born and Robinson (2021)](https://doi.org/10.5194/tc-15-4539-2021) and
[Rieckh et al. (2024)](https://doi.org/10.5194/gmd-17-6987-2024). This is a
reimplementation of the Bergen ELSA v2.0
([git.app.uib.no/melt-team-bergen/elsa](https://git.app.uib.no/melt-team-bergen/elsa)),
not a fork. See [docs/DESIGN.md](docs/DESIGN.md) for where it departs from the
published scheme, and why.

Documentation, including a page per benchmark, is published at
[fesmc.github.io/elsa](https://fesmc.github.io/elsa/).

> **Status: under construction.** The library, its public API, NetCDF output,
> restart, both benchmarks and the Julia analysis are in place: `make check`
> passes serial, OpenMP and under bounds checking, and `make validate` passes.
> The Yelmox coupling is still to come.

## Install

elsa is a [configme](https://github.com/fesmc/configme) package. Its only
dependency is [fesm-utils](https://github.com/fesmc/fesm-utils), for the `ncio`
and `nml` modules. The advection is an explicit sub-stepped upwind scheme, so
there is no linear solver and no LIS — v2.0 needed one, and needed you to build
and locate it.

```bash
configme install elsa --only     # elsa + fesm-utils
```

To configure an existing checkout for your machine:

```bash
cd elsa
configme -m macbook -c gfortran    # writes the repo-root Makefile
```

## Build

```bash
make elsa-static     # libelsa/include/libelsa.a
make all             # the library and every benchmark
make usage           # all targets
```

Add `openmp=1` to thread the layer loop, or `debug=3` for bounds checking and
floating-point traps.

## Benchmarks

```bash
make check           # runs every benchmark; nonzero exit on any failure
```

**`test_physics.x`** exercises the advection and layer kernels: uniform,
rotational, convergent and divergent flow, plus the mass-balance bookkeeping.
The convergent and divergent cases run at a raw CFL of ~976 and stay
non-negative, bounded and mass-conserving with no clipping anywhere in the code.

**`test_interp.x`** asserts the exactness properties of the maps: conservative
remap conserves mass at a non-integer `grid_factor`, bilinear reproduces a linear
field at the faces for either staggering, and the layer-mean integral is exact
for a linear velocity profile.

**`test_column.x`** is the quantitative one. At an ice divide with no horizontal
flow, constant thickness and constant accumulation, elsa's isochrones must follow
Nye's `z = H exp(-a t/H)`. elsa never computes a vertical velocity — the thinning
emerges from adding accumulation to the top layer and renormalizing the column —
and it converges onto Nye at first order in the coupling period.

**`test_greenland.x`** runs the 3D ice sheet at 16 km, forced offline from a
Yelmo restart (`data/initmip-grl-16km/yelmo_restart.nc`, the same file `tracer`
uses). The restart carries no spun-up isochrone field, so this is a structural
check against a real velocity field, real staggering and a real margin, rather
than a comparison with a known answer. It writes `output/GRL-16KM/elsa.nc`.

## Using the library

```fortran
use elsa

type(elsa_class) :: els

call elsa_init(els,"elsa.nml","elsa",time,time_end,xc,yc,zeta_aa,H_ice,"acx_acy")

call elsa_update(els,time,H_ice,ux,uy,smb,bmb)   ! every host timestep

call elsa_end(els)
```

`zeta` must be a strictly ascending sigma axis with `zeta(1) = 0` at the bed and
`zeta(nz) = 1` at the surface, and `size(zeta) == size(ux,3)`. `stagger`
declares where the host's velocity samples sit: `"acx_acy"` for staggered
velocities (Yelmo), `"aa"` for cell-centred ones. Host fields may be single or
double precision — elsa converts at the boundary, so a host never casts.

`elsa_update` takes an absolute time, works out its own `dt`, decides internally
whether an update is due, and keeps its own previous-step ice thickness. The
host calls it unconditionally, once per timestep, and manages none of elsa's
bookkeeping. `time_end` is needed at init only to size the layer stack, which is
allocated once and never grown.

To restart, write a file and pass it back:

```fortran
call elsa_restart_write(els,"elsa_restart.nc")
...
call elsa_init(els,"elsa.nml","elsa",time,time_end,xc,yc,zeta_aa,H_ice,"acx_acy", &
               restart="elsa_restart.nc")
```

The layer stack and the isochrone schedule then come from the file rather than
from `time_end` and `layer_resolution`. A restarted run is bit-identical to the
run that never stopped, which `test_greenland.x` asserts.

Build with `openmp=1` to thread the layer loop. The layers never exchange mass,
so the result is bit-identical to the serial one at any thread count.

## Analysis

```bash
make validate        # runs check, then the Julia validation and figures
```

The division of labour: the Fortran benchmarks assert structural properties and
exit nonzero; Julia checks the physics against closed-form answers and draws.
`analysis/` is a Julia project using CairoMakie and NCDatasets, and the figures
land in `plots/`:

  - `column_nye.png` — modelled isochrone height against Nye's analytic solution
    and against elsa's own discrete recursion.
  - `greenland_isochrones.png` — the depth of the oldest isochrone as a map, and
    the full layer stack along a transect. This second panel is the view the
    isochronal scheme exists to produce: the layers a model compares against
    radiostratigraphy.

`analysis/elsa_analysis.jl` carries the shared readers. The key thing it knows is
that elsa's vertical axis is time: layer `k` is bounded below by the isochrone
laid down at `layer_time[k]`, so that isochrone sits at `dsum_iso[:,:,k-1]`, and
the initialization layers have no `layer_time` at all.

## License

GPL-3.0, following ELSA v2.0.
