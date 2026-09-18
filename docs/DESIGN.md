# elsa design notes

elsa is a reimplementation of ELSA rather than a fork. The method follows Born
(2017), Born and Robinson (2021) and Rieckh et al. (2024, GMD 17, 6987–7000,
hereafter **R24**); the code is new.

This document records the departures from the published v2.0 scheme and their
motivation. Every departure changes results relative to R24, so each is stated
explicitly rather than left implicit.

## What elsa is

elsa consists of a stack of isochronal layers on a 2D grid. The vertical axis
is *time*: each layer is bounded by two isochrones and does not exchange mass
with its neighbours. Layers thin as ice flows toward the margin, and vertical
motion emerges from the changes in individual layer thicknesses, rather than
from an imposed vertical velocity. This eliminates vertical numerical diffusion
by construction, which is the central motivation of the scheme.

elsa is not an ice-sheet model. It is driven by a host model, which supplies
the horizontal velocity, the ice thickness, and the surface and basal mass
balance. elsa does not modify these fields.

## The advection scheme

### What R24 does

The layer thickness $d$ obeys the flux-form continuity equation (R24 Eq. B1)

    ∂d/∂t = -∂(ud)/∂x - ∂(vd)/∂y

which is discretized implicitly in time, upstream in space, with velocities
defined at the cell faces. For $u > 0$, $v > 0$ this gives (R24 Eq. B5)

    d^t(i,j) = d^{t+1}(i,j) · [1 + (Δt/Δx)·u(i+½,j) + (Δt/Δy)·v(i,j+½)]
             - d^{t+1}(i-1,j) · (Δt/Δx)·u(i-½,j)
             - d^{t+1}(i,j-1) · (Δt/Δy)·v(i,j-½)

with three companion cases (Eqs. B6–B8) for the remaining sign combinations,
and the system $A \cdot d^{t+1} = d^t$ solved with LIS.

### The defect

Equations B5–B8 select the upwind direction for **both** x-faces from a single
sign test on $u(i+½,j)$, and for both y-faces from $v(i,j+½)$. However, the
flux through face $i-½$ is upwinded according to $u(i-½,j)$, which may have the
opposite sign.

In a convergent cell $u(i-½,j) > 0 > u(i+½,j)$, Eq. B6 is selected and the main
diagonal becomes $1 - (\Delta t/\Delta x) \cdot u(i-½,j)$, which is
**negative** for $u(i-½,j) > \Delta x/\Delta t$. The matrix then loses diagonal
dominance and the solve produces unbounded layer thicknesses.

This is the instability that R24 documents in Sect. 3.5 — *"layer thickness can
become unrealistically large during one advection step... often at the
ice-sheet boundaries where velocities are large"* — and works around with a
maximum allowed thickness change per step, $d_\mathrm{max} = 100 + 10 \cdot
\mathrm{update\_factor} + \mathrm{layer\_resolution}/10$ (Eq. B9). The v2.0
code carries two further guards: a fallback to the previous value whenever the
solution is non-finite, and a skip of all advection where $|u| < 0.1$ m yr$^{-1}$,
commented *"using 0 is unstable"*.

Zero velocity is not unstable. Sign-changing velocity is.

### The fix: upwind each face on its own velocity

    F(i+½) = max(u(i+½),0)·d(i) + min(u(i+½),0)·d(i+1)
    F(i-½) = max(u(i-½),0)·d(i-1) + min(u(i-½),0)·d(i)

    ∂d(i)/∂t = -[F(i+½) - F(i-½)]/Δx   (and likewise in y)

Zero velocity requires no special case, and neither does sign reversal across a
cell. $d_\mathrm{max}$ (B9), the non-finite fallback, and the 0.1 m yr$^{-1}$
threshold have all been removed. None is ported.

Written implicitly, this yields a matrix with diagonal $\geq 1$ and
off-diagonals $\leq 0$. It is **column** diagonally dominant, rather than row
diagonally dominant: each off-diagonal entry in column $j$ is the coefficient
with which cell $j$ feeds a neighbour, so the off-diagonal column sum is
exactly the total outflow from $j$, and

    A(j,j) - Σ_{i≠j} |A(i,j)|  =  1     exactly, for every j and every u

A strongly convergent cell has inflow that far exceeds outflow, so its *row*
is not dominant at all: the margin goes negative. It is easy to get this
backwards. Three consequences follow, all verified numerically at CFL $\gg 1$:

  - $A$ is a Z-matrix with strictly dominant columns and a positive diagonal,
    and is therefore a non-singular M-matrix with $A^{-1} \geq 0$. Layer
    thicknesses thus remain non-negative for any velocity field and any
    timestep.
  - The interior column sums are exactly 1, so $\sum d^{t+1} = \sum d^t$ up to
    boundary fluxes. Mass is conserved to roundoff.
  - $\|A^{-1}\|_1 = 1$, so $\|d^{t+1}\|_1 \leq \|d^t\|_1$. The scheme is
    unconditionally stable in the mass norm, which is the relevant norm for a
    thickness field.

### The time discretization: explicit, sub-stepped

elsa evaluates the same fluxes **explicitly**, sub-stepping the coupling
period so that the outflow rate obeys

    Δt_sub · [ (max(u(i+½),0) - min(u(i-½),0))/Δx
             + (max(v(j+½),0) - min(v(j-½),0))/Δy ]  ≤  CFL   ≤ 1

for every cell. This is the exact positivity condition for the explicit
update: the coefficient multiplying $d(i,j)$ is $1 - \Delta t \cdot r(i,j)$,
and everything else entering the cell is non-negative. The substep count is
computed per layer, per coupling step, from the actual velocity field.

Explicit upwinding under this condition inherits the same three guarantees as
the implicit form — non-negativity, exact mass conservation and L1 stability —
without the additional damping that the implicit solve applies on top of the
upwind truncation error. This is important for a model whose purpose is not to
smear layers.

There is no linear solve, and therefore **LIS is not a dependency of elsa.**
R24 notes *"Its only dependency is the Library of Iterative Solvers"*; this
implementation has no dependency beyond fesm-utils.

The reasoning, at the R24 CTRL configuration (16 km, CP = 10 yr,
$u_\mathrm{max} \approx 12$ km yr$^{-1}$):

  - Implicit BiCG to R24's `-tol 1.0e-12` runs roughly 30–50 iterations, each
    a 5-point matvec plus vector work: on the order of 800 flops per cell.
    In addition, there are $5 \cdot n_x \cdot n_y$ `lis_matrix_set_value`
    calls per layer per coupling step.
  - Explicit needs $\lceil u_\mathrm{max} \Delta t / \Delta x \rceil = 8$
    substeps at roughly 15 flops each: on the order of 120 flops per cell.
    No assembly, no solver objects and no library calls.

The asymptotics are the same, not just the constant. Explicit cost grows as
$n_l \cdot N \cdot (u_\mathrm{max} \Delta t/\Delta x) \propto n_l/\Delta x^3$.
However, the Krylov iteration count on an advection-dominated M-matrix grows
like the characteristic path length, $O(1/\Delta x)$, in the absence of a
strong preconditioner. The implicit cost is therefore *also*
$\propto n_l/\Delta x^3$, with a worse constant.

Assembling the $n_l$ layers into one block-diagonal system would have been
worse still: the blocks are decoupled, so a Krylov method applied to the whole
system iterates until the *worst* layer converges, while the rest pay the
full matvec cost. Its Krylov vectors are of length $n_l \cdot N$ rather than
$N$.

**Consequence.** A large coupling period no longer reduces cost: the cost is
now roughly proportional to simulated time regardless of CP, since the substep
count absorbs it. This is the intended trade-off. R24's Figs. 7–8 show that
CP = 200 already carries an RMSE of 40–80 m, so the cheap-CP regime was in
practice not one worth preserving. The parameters `grid_factor` and
`layer_resolution` remain the effective levers for ensemble cost, and R24's
Fig. 9 shows they always were the more important ones.

### Boundaries

Zero normal flux is imposed at the domain edge, so global mass is conserved
exactly. v2.0 instead froze the layer thickness at boundary cells. Ice
reaching the elsa domain boundary is a host-domain problem, and is not
something that elsa should correct for silently.

## The update sequence

`elsa_update` receives an absolute time and does nothing unless a coupling
period has elapsed. When it fires, the sequence is:

1. Map $H_\mathrm{ice}$, smb, bmb and the velocities onto elsa's grid.
2. Normalize the layer stack onto $H_\mathrm{ice, prev}$ — the host thickness
   on which the incoming velocities were computed.
3. Add smb$\cdot$dt to the top layer and bmb$\cdot$dt to the bottom, exhausting
   layers in turn where the balance is negative.
4. Accumulate `dsum`, then take the layer-mean velocity at every face.
5. Advect each layer independently.
6. Reseed any column that elsa has emptied but in which the host still has ice
   (see below).
7. Normalize onto $H_\mathrm{ice}$. Horizontal layer advection does not enforce
   the host's mass conservation, so this step applies the drift correction.
   Relative layer thicknesses are unchanged.
8. Lay down any isochrone whose time has been reached.

Steps 2 and 6 are where elsa's vertical motion originates. elsa never computes
a vertical velocity: adding accumulation on top and renormalizing the column
multiplies every layer height by $H/(H + a \cdot dt)$, which in the limit
$dt \to 0$ is Nye's uniform vertical strain. `test_column.x` asserts both the
discrete result (exact to $5 \times 10^{-15}$) and its first-order convergence
onto $z = H \exp(-a \cdot t / H)$.

elsa owns `H_ice_prev` and the isochrone schedule. v2.0 required the driver to
manage the previous ice thickness across two lines 180 apart in
`yelmox_elsa.f90`, and to wrap the call in its own `its_time` check.

Where `layer_resolution < dt_coupling`, more than one isochrone falls inside a
coupling period. v2.0 rejected such a configuration at init; elsa instead lays
down each isochrone in turn — the additional layers simply receive no
accumulation — and issues a warning at init.

### Emptied columns

Surface ablation removes ice from the top layer downward and stops when the
column runs out; it cannot remove ice that is not there. Where $-\mathrm{smb}
\cdot dt$ exceeds the whole column, the column is left at exactly zero and
the normalization has nothing to rescale: the column is then dead for the rest
of the run. v2.0 has the same hole — `normalise_d` zeroes a column whose sum
is not positive, silently and permanently.

On the 16 km Greenland benchmark, 85 of 7204 ice cells satisfy $-\mathrm{smb}
\cdot dt > H$ at `dt_coupling = 50 yr`. All are thin margin cells (3–244 m,
against a median ice thickness of 1718 m). Advection refills most of them, and
around nine per step do not recover.

elsa returns the host's ice to those columns, **in layer 1, at the bed**.
Ablation eats a column from the top, so the last ice standing before
exhaustion is the deepest. Placing the survivor there is the consistent
continuation of the process that destroyed the column, and is where basal
freeze-on already deposits ice. Reseeding at the top would instead claim that
the ice is young, which contradicts the process that removed it.

This is defined behaviour for an under-determined state — elsa cannot date ice
that it believes it has removed — rather than a correction to a bug. It is
**counted** rather than silent: `now%n_reseed` (last update) and
`now%n_reseed_total`, together with a note at first occurrence. A large or
growing per-step count means that `dt_coupling` is too long, or that the
forcing is inconsistent.

### The face column

The layer-mean velocity requires the column geometry *at each face*. elsa
takes this as the average of the two cell columns that the face separates, so
that `dsum_face(n_top) = H_face` holds by construction. This includes the
margin, where one neighbour is ice-free and the face thickness is half that
of the ice-covered cell. Mapping the host's $H$ directly onto a face grid
with its own conservative weights would not stay consistent with elsa's layer
stack.

## Vertical interpolation of the host velocity

R24 states only that host velocities are *"linearly interpolated in the
vertical... onto the isochronal grid"*. The v2.0 code evaluates them at
`dsum(i,j,iz)`, the **upper interface** of layer `iz`.

In the flux form $\partial d/\partial t = -\nabla \cdot (\bar u \, d)$, the
velocity that transports layer `iz` is its thickness-average

    ū(iz) = (1/d(iz)) · ∫ over [dsum(iz-1), dsum(iz)] of u(z) dz

Sampling at the upper interface instead biases every layer toward the faster
ice above it, and so systematically over-advects the stack. Since $u(z)$ is
piecewise-linear on the host's vertical grid, the integral is exact and
computationally cheap.

**Consequence.** elsa uses the layer-mean velocity. Deep layers advect more
slowly than in v2.0; the effect grows with layer thickness, so it is largest
for coarse `layer_resolution` and near the bed.

## Grids, staggering, and the host contract

### Vertical

elsa requires host velocities on `nz` levels co-located with a strictly
ascending sigma axis `zeta`, where `zeta(1) = 0` at the bed and `zeta(nz) = 1`
at the surface, and where `size(zeta) == size(ux,3)`. This is checked at
`elsa_init`. It matches Yelmo's `zeta_aa` convention and the contract that
`tracer` states, so that the two packages can be compared directly on the
same host fields.

This resolves an inconsistency in v2.0, where the offline path read `nz+1`
values from `zeta_aa.txt` into an `nz`-element array, while the coupled path
received Yelmo's `nz`.

### Horizontal

The scheme requires velocities **at cell faces** (R24 App. B). Rather than
treat this as two code paths, `stagger` declares where the host's velocity
samples live:

  - `"acx_acy"` — `ux` at $(x + dx/2, y)$, `uy` at $(x, y + dy/2)$. Yelmo.
  - `"aa"` — both at $(x, y)$, cell-centred.

This changes only the *coordinates* of the source samples. The same
interpolation then lands them on elsa's acx/acy faces, whatever the relative
resolution or offset of the two grids. There is no unstagger/restagger round
trip, and no smoothing penalty for the common case `grid_factor = 1`.

Mass-like fields (`H_ice`, `smb`, `bmb`) are remapped conservatively (area
weighted), so that column mass balance is exact. Velocities are bilinear.
v2.0's `regrid_xy` did integer-factor box averaging only, required
`grid_factor` to divide both `nx` and `ny` (it warned but continued
otherwise), and box-averaged the staggered velocities as though they were
cell-centred.

`grid_factor` is now a real value $\geq 1$ and is no longer required to be an
integer.

### Why not coords

Both grids are axis-aligned and share a plane by construction — elsa's grid
is a coarsening of the host's, over the same axes — so every horizontal map
is *separable*: the 2D weight is the outer product of a 1D weight along x and
a 1D weight along y, and it is exact. Both weight vectors are precomputed
once at `elsa_map_init` and applied as a fixed stencil, so no index search
takes place in the time loop.

fesm-utils' `coords` does provide this functionality. `conservative_weights`
even detects the same-Cartesian-system case and falls through to an analytic
separable overlap — the same algorithm. However, reaching it requires
constructing two `grid_class` objects with projection metadata and a
`weight_map_t`, none of which elsa's $(x, y)$ host contract carries. In
addition, `interp2D::interp_bilinear` re-searches its bracketing indices on
every call, which elsa would pay per level, per component, per coupling step.
The separable weights add roughly 200 lines of code and provide exactness
together with a search-free inner loop. fesm-utils remains elsa's only
dependency, for `nml` and `ncio`.

Since `stagger` is expressed as a coordinate offset rather than a code path,
fesm-utils' `staggering` module is not needed either: nothing downstream
branches on where the host's velocities were originally located.

## Precision

Internally, `wp = dp`. Layer thicknesses are summed over a stack of order
$10^3$ layers, and are renormalized against the host ice thickness every
coupling period across $O(10^5)$ yr. Run time is dominated by the linear
solve rather than by array traffic, so the wider type is nearly free.

Host fields are accepted as `real(sp)` or `real(dp)` through generic
interfaces and converted at the boundary, so the host never casts. Yelmo's
`wp` is currently single; this decouples elsa from that choice.

## Dependencies

  - **fesm-utils** — `ncio`, `nml`, `staggering`, `coords` (`interp1D`,
    `interp2D`, `conservative`). Nothing is vendored that fesm-utils already
    provides.
  - **LIS** — vendored and built by fesm-utils, so elsa needs no system-wide
    install. Sources that `#include "lisf.h"` are named `.F90` so that
    gfortran and ifx both preprocess them without a compiler-specific flag.

## Deliberate omissions

  - **The dye tracer** (`tracer_iso`) is not carried over. In v2.0 it was
    passed unallocated into the advection routine whenever `use_dye_tracer`
    was false. It will return as a tested feature with its own benchmark.
  - **`misc_1`**, a debug array written to output, is deleted.
  - **Restart** is new. v2.0 could not restart at all: `elsa_dealloc` freed
    seven of the nine arrays that `elsa_init` allocates, so a second
    `elsa_init` on the same object aborted on an already-allocated array.

## Restart

`elsa_restart_write` stores only what cannot be reconstructed: `d_iso`,
`H_ice_prev`, `time`, `n_top`, `i_add`, and the isochrone schedule
`time_add`. `dsum_iso` is derived from `d_iso`, and the velocities and mass
balance are remapped from the host on every update.

Two details matter, and the round-trip test would fail without either:

  - `d_iso` and `H_ice_prev` are written in **double** precision, unlike the
    diagnostic output of `elsa_write_step`, which is single. A restart that
    loses bits does not reproduce the run that it continues.
  - `time_add` travels in the restart file, rather than being regenerated from
    `layer_resolution`. Rebuilding it from the restart time would shift every
    subsequent isochrone.

The grid on which the restart was written must match the grid onto which it
is read — dimensions, axes and `zeta`. A mismatch is a hard error, rather
than something to interpolate away silently.

`test_greenland.x` asserts that a run that is stopped, written and restarted
is **bit-identical** to the run that never stopped. This is the only check
strong enough to catch a piece of state being reconstructed rather than
carried.
