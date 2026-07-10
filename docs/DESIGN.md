# elsa design notes

This is a reimplementation of ELSA, not a fork. The method is that of Born
(2017), Born and Robinson (2021) and Rieckh et al. (2024, GMD 17, 6987–7000,
hereafter **R24**); the code is new.

This document records where the implementation departs from the published v2.0
scheme, and why. Every departure changes results relative to R24, so each one
is stated explicitly rather than absorbed silently.

## What elsa is

A stack of isochronal layers on a 2D grid. The vertical axis is *time*: each
layer is bounded by two isochrones and never exchanges mass with its
neighbours. Layers thin as ice flows toward the margin, and vertical motion
emerges from changes in the individual layer thicknesses rather than from a
vertical velocity. This eliminates vertical numerical diffusion by
construction, which is the whole point of the scheme.

elsa is not an ice-sheet model. It is driven by a host model, which supplies
horizontal velocity, ice thickness, and surface and basal mass balance, and
which elsa never modifies.

## The advection scheme

### What R24 does

Layer thickness `d` obeys the flux-form continuity equation (R24 Eq. B1)

    ∂d/∂t = -∂(ud)/∂x - ∂(vd)/∂y

discretized implicitly in time, upstream in space, with velocities at the cell
faces. For `u > 0`, `v > 0` this gives (R24 Eq. B5)

    d^t(i,j) = d^{t+1}(i,j) · [1 + (Δt/Δx)·u(i+½,j) + (Δt/Δy)·v(i,j+½)]
             - d^{t+1}(i-1,j) · (Δt/Δx)·u(i-½,j)
             - d^{t+1}(i,j-1) · (Δt/Δy)·v(i,j-½)

with three companion cases (Eqs. B6–B8) for the other sign combinations, and
the system `A·d^{t+1} = d^t` solved with LIS.

### The defect

Equations B5–B8 select the upwind direction for **both** x-faces from a single
sign test on `u(i+½,j)` (and both y-faces from `v(i,j+½)`). But the flux
through face `i-½` is upwinded according to `u(i-½,j)`, which may have the
opposite sign.

In a convergent cell, `u(i-½,j) > 0 > u(i+½,j)`, Eq. B6 is selected and the
main diagonal becomes `1 - (Δt/Δx)·u(i-½,j)`, which is **negative** for
`u(i-½,j) > Δx/Δt`. The matrix loses diagonal dominance and the solve produces
unbounded layer thicknesses.

This is the instability R24 documents in Sect. 3.5 — *"layer thickness can
become unrealistically large during one advection step... often at the
ice-sheet boundaries where velocities are large"* — and works around with a
maximum allowed thickness change per step, `d_max = 100 + 10·update_factor +
layer_resolution/10` (Eq. B9). The v2.0 code carries two further guards: a
fallback to the previous value whenever the solution is non-finite, and a skip
of all advection where `|u| < 0.1 m/yr`, commented *"using 0 is unstable"*.

Zero velocity is not unstable. Sign-changing velocity is.

### What elsa does instead

Upwind each face on its own velocity:

    F(i+½) = max(u(i+½),0)·d(i) + min(u(i+½),0)·d(i+1)
    F(i-½) = max(u(i-½),0)·d(i-1) + min(u(i-½),0)·d(i)

so that

    d^{t+1}(i) · [1 + (Δt/Δx)·(max(u(i+½),0) - min(u(i-½),0))]
      + d^{t+1}(i+1) · (Δt/Δx)·min(u(i+½),0)
      - d^{t+1}(i-1) · (Δt/Δx)·max(u(i-½),0)   =   d^t(i)

The diagonal is `≥ 1` and the off-diagonals are `≤ 0` for any velocity field.
The matrix is a strictly diagonally dominant M-matrix, so the scheme is
unconditionally stable, positivity-preserving, and exactly mass-conserving.
`u = 0` needs no special case.

**Consequence.** `d_max` (B9), the non-finite fallback, and the `0.1 m/yr`
velocity threshold are all deleted. None is ported. Being an M-matrix also
conditions the system far better, so the BiCG/Jacobi solve converges in fewer
iterations than in v2.0.

The solve stays implicit and stays on LIS: there is no CFL constraint, which
matters because the coupling period reaches 200 yr.

## Vertical interpolation of the host velocity

R24 states only that host velocities are *"linearly interpolated in the
vertical... onto the isochronal grid"*. The v2.0 code evaluates them at
`dsum(i,j,iz)`, the **upper interface** of layer `iz`.

In the flux form `∂d/∂t = -∇·(ū d)`, the velocity that transports layer `iz` is
its thickness-average

    ū(iz) = (1/d(iz)) · ∫ over [dsum(iz-1), dsum(iz)] of u(z) dz

Sampling at the upper interface instead biases every layer toward the faster
ice above it, systematically over-advecting the stack. Since `u(z)` is
piecewise-linear on the host's vertical grid, the integral is exact and cheap.

**Consequence.** elsa uses the layer-mean velocity. Deep layers advect more
slowly than in v2.0; the effect grows with layer thickness, so it is largest
for coarse `layer_resolution` and near the bed.

## Grids, staggering, and the host contract

### Vertical

elsa requires host velocities on `nz` levels co-located with a strictly
ascending sigma axis `zeta`, with `zeta(1) = 0` at the bed and `zeta(nz) = 1`
at the surface, and `size(zeta) == size(ux,3)`. This is checked at
`elsa_init`. It is exactly Yelmo's `zeta_aa` convention, and exactly the
contract `tracer` states, so the two packages can be compared directly on the
same host fields.

This resolves an inconsistency in v2.0, where the offline path read `nz+1`
values from `zeta_aa.txt` into an `nz`-element array while the coupled path
received Yelmo's `nz`.

### Horizontal

The scheme wants velocities **at cell faces** (R24 App. B). Rather than treat
that as two code paths, `stagger` declares where the host's velocity samples
live:

  - `"acx_acy"` — `ux` at `(x + dx/2, y)`, `uy` at `(x, y + dy/2)`. Yelmo.
  - `"aa"` — both at `(x, y)`, cell-centred.

This changes only the *coordinates* of the source samples. The same
interpolation then lands them on elsa's own acx/acy faces, whatever the two
grids' relative resolution or offset. There is no unstagger/restagger round
trip and no smoothing penalty for the common `grid_factor = 1` case.

Mass-like fields (`H_ice`, `smb`, `bmb`) are remapped conservatively
(area-weighted), so column mass balance is exact. v2.0's `regrid_xy` did
integer-factor box averaging only, required `grid_factor` to divide both `nx`
and `ny` (it warned but continued otherwise), and box-averaged the staggered
velocities as though they were cell-centred.

## Precision

`wp = dp` internally. Layer thicknesses are summed over a stack reaching O(10³)
layers and renormalized against host ice thickness every coupling period across
O(10⁵) yr. Run time is dominated by the linear solve, not array traffic, so the
wider type is nearly free.

Host fields are accepted as `real(sp)` or `real(dp)` through generic interfaces
and converted at the boundary, so a host never casts. Yelmo's `wp` is currently
single; this decouples elsa from that choice.

## Dependencies

  - **fesm-utils** — `ncio`, `nml`, `staggering`, `coords` (`interp1D`,
    `interp2D`, `conservative`). Nothing is vendored that fesm-utils provides.
  - **LIS** — vendored and built by fesm-utils, so elsa needs no system-wide
    install. Sources that `#include "lisf.h"` are named `.F90` so gfortran and
    ifx both preprocess them without a compiler-specific flag.

## Deliberate omissions

  - **The dye tracer** (`tracer_iso`) is not carried over. In v2.0 it was
    passed unallocated into the advection routine whenever `use_dye_tracer` was
    false. It will return as a tested feature with its own benchmark.
  - **`misc_1`**, a debug array written to output, is deleted.
  - **Restart** is new. v2.0 could not restart at all: `elsa_dealloc` freed
    seven of the nine arrays `elsa_init` allocated, so a second `elsa_init` on
    the same object aborted on an already-allocated array.
