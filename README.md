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

> **Status: under construction.** The library and its public API are in place
> and validated end-to-end against the Nye analytic solution: `make check`
> passes. NetCDF output, restart, the Greenland benchmark and the Julia analysis
> are still being written.

## Install

elsa is a [configme](https://github.com/fesmc/configme) package. It depends on
[fesm-utils](https://github.com/fesmc/fesm-utils) for its `ncio`, `nml`,
`staggering` and `coords` modules, and links the LIS solver that fesm-utils
vendors and builds — so there is no system-wide LIS or netCDF path to set by
hand.

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
make usage           # all targets
```

Add `debug=1` to any build for bounds checking and floating-point traps.

## Benchmarks

```bash
make check           # runs every benchmark; nonzero exit on any failure
```

`test_column.x` is the quantitative one. At an ice divide with no horizontal
flow, constant thickness and constant accumulation, elsa's isochrones must
follow Nye's `z = H exp(-a t/H)`. elsa never computes a vertical velocity — the
thinning emerges from adding accumulation to the top layer and renormalizing the
column — and it converges onto Nye at first order in the coupling period.

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

Build with `openmp=1` to thread the layer loop. The layers never exchange mass,
so the result is bit-identical to the serial one at any thread count.

Diagnostics and figures are Julia, under `analysis/`, plotted with CairoMakie:
the Fortran benchmarks assert and exit nonzero, Julia validates quantitatively
and draws.

## License

GPL-3.0, following ELSA v2.0.
