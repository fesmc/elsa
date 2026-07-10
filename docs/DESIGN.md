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

### The fix: upwind each face on its own velocity

    F(i+½) = max(u(i+½),0)·d(i) + min(u(i+½),0)·d(i+1)
    F(i-½) = max(u(i-½),0)·d(i-1) + min(u(i-½),0)·d(i)

    ∂d(i)/∂t = -[F(i+½) - F(i-½)]/Δx   (and likewise in y)

Zero velocity needs no special case; sign reversal across a cell needs no
special case. `d_max` (B9), the non-finite fallback, and the `0.1 m/yr`
threshold are all deleted. None is ported.

Written implicitly this yields a matrix with diagonal `≥ 1` and off-diagonals
`≤ 0`. It is **column** diagonally dominant, not row diagonally dominant: each
off-diagonal entry in column `j` is the coefficient with which cell `j` feeds a
neighbour, so the off-diagonal column sum is exactly the total outflow from `j`,
and

    A(j,j) - Σ_{i≠j} |A(i,j)|  =  1     exactly, for every j and every u

A strongly convergent cell has inflow far exceeding outflow, so its *row* is not
dominant at all — the margin goes negative. It is easy to get this backwards.
Three consequences follow, all verified numerically at CFL ≫ 1:

  - `A` is a Z-matrix with strictly dominant columns and positive diagonal, so
    it is a nonsingular M-matrix and `A⁻¹ ≥ 0`: layer thicknesses stay
    non-negative for any velocity field and any timestep.
  - The interior column sums are exactly `1`, so `Σ d^{t+1} = Σ d^t` up to
    boundary fluxes: mass is conserved to roundoff.
  - `‖A⁻¹‖₁ = 1`, so `‖d^{t+1}‖₁ ≤ ‖d^t‖₁`: unconditionally stable in the mass
    norm, which is the norm that matters for a thickness field.

### The time discretization: explicit, sub-stepped

elsa evaluates the same fluxes **explicitly**, sub-stepping the coupling period
so that the outflow rate obeys

    Δt_sub · [ (max(u(i+½),0) - min(u(i-½),0))/Δx
             + (max(v(j+½),0) - min(v(j-½),0))/Δy ]  ≤  CFL   ≤ 1

for every cell. This is the exact positivity condition for the explicit update:
the coefficient multiplying `d(i,j)` is `1 - Δt·r(i,j)`, and everything else
entering the cell is non-negative. The substep count is computed per layer, per
coupling step, from the actual velocity field.

Explicit upwinding under this condition inherits the same three guarantees
listed above — non-negativity, exact mass conservation, L1 stability — and adds
none of the extra damping the implicit solve applies on top of the upwind
truncation error. For a model whose purpose is not smearing layers, that matters.

There is no linear solve, so **LIS is not a dependency of elsa.** R24 notes
*"Its only dependency is the Library of Iterative Solvers"*; this
implementation has none beyond fesm-utils.

The reasoning, at the R24 CTRL configuration (16 km, `CP = 10 yr`,
`u_max ≈ 12 km/yr`):

  - Implicit BiCG to R24's `-tol 1.0e-12` runs ~30–50 iterations, each a
    5-point matvec plus vector work: order 800 flops/cell. Plus `5·nx·ny`
    `lis_matrix_set_value` calls per layer per coupling step.
  - Explicit needs `⌈u_max·Δt/Δx⌉ = 8` substeps at ~15 flops each: order
    120 flops/cell. No assembly, no solver objects, no library calls.

The asymptotics are the same, not just the constant. Explicit cost grows as
`nl·N·(u_max Δt/Δx) ∝ nl/Δx³`. But Krylov iteration count on an
advection-dominated M-matrix grows like the characteristic path length,
`O(1/Δx)`, absent a strong preconditioner — so the implicit cost is *also*
`∝ nl/Δx³`, with a worse constant.

Assembling the `nl` layers into one block-diagonal system would have been worse
still: the blocks are decoupled, so a Krylov method on the whole thing iterates
until the *worst* layer converges while the rest pay full matvec cost, and its
Krylov vectors are `nl·N` long rather than `N`.

**Consequence.** A large coupling period no longer buys speed: cost is now
roughly proportional to simulated time regardless of `CP`, because the substep
count absorbs it. This is the intended trade. R24's own Figs. 7–8 show `CP = 200`
already carries 40–80 m RMSE, so the cheap-large-`CP` regime was not one worth
preserving. `grid_factor` and `layer_resolution` remain the effective levers for
ensemble cost, and Fig. 9 shows they always were the larger ones.

### Boundaries

Zero normal flux at the domain edge, so global mass is conserved exactly. v2.0
instead froze layer thickness at boundary cells. Ice reaching the elsa domain
boundary is a host-domain problem, not something elsa should paper over.

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
