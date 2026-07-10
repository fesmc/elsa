module elsa_interp
    ! Mapping the host ice-sheet model's fields onto elsa's grid.
    !
    ! Both grids are axis-aligned and share a plane by construction: elsa's grid
    ! is a coarsening of the host's, over the same x/y axes. Every horizontal map
    ! is therefore separable — the 2D weight is the outer product of a weight
    ! along x and a weight along y — and exact. Both maps are precomputed once
    ! at elsa_map_init and applied as a fixed stencil thereafter, so no index
    ! search happens inside the time loop.
    !
    ! Three maps are needed, and they differ only in the axes they connect:
    !
    !   scalars (H_ice, smb, bmb)   source aa nodes  -> elsa aa nodes
    !                               area-weighted (conservative) overlap
    !
    !   ux                          source velocity x-nodes -> elsa acx faces
    !   uy                          source velocity y-nodes -> elsa acy faces
    !                               bilinear
    !
    ! `stagger` enters only as an offset applied to the source velocity node
    ! coordinates: `"acx_acy"` puts ux at x_src + dx_src/2, `"aa"` leaves it at
    ! x_src. Nothing downstream branches on it. This is what lets any host grid
    ! -- staggered or centred, any resolution ratio, integer or not -- reach
    ! elsa's faces through one code path.
    !
    ! Vertically, the host supplies velocity on `nz` sigma levels. elsa needs the
    ! *thickness-average* of that velocity over each isochronal layer, because
    ! the layer thickness obeys d(d)/dt = -div(u_bar d). Since u is piecewise
    ! linear in z, the average integrates exactly.

    use elsa_precision, only : wp, H_ICE_MIN

    implicit none

    private

    ! A layer thinner than this is treated as a point: the layer-mean velocity
    ! degenerates to the value at the layer's own position rather than 0/0.
    real(wp), parameter :: DZ_MIN = 1.0e-8_wp

    type axis_map_class
        ! Sparse weights along one axis: target point t draws on source points
        ! idx(1:stencil,t) with weights wt(1:stencil,t), which sum to 1.
        integer :: n_tgt   = 0
        integer :: stencil = 0
        integer,  allocatable :: idx(:,:)
        real(wp), allocatable :: wt(:,:)
    end type axis_map_class

    type elsa_map_class
        ! elsa's own grid
        integer  :: nx, ny, nz
        real(wp) :: dx, dy
        real(wp), allocatable :: x(:), y(:)
        real(wp), allocatable :: zeta(:)

        ! source aa nodes -> elsa aa nodes, conservative
        type(axis_map_class) :: cons_x, cons_y

        ! source ux nodes -> elsa acx faces, bilinear
        type(axis_map_class) :: ux_x, ux_y

        ! source uy nodes -> elsa acy faces, bilinear
        type(axis_map_class) :: uy_x, uy_y
    end type elsa_map_class

    public :: elsa_map_class
    public :: elsa_map_init
    public :: elsa_map_end
    public :: map_scalar
    public :: map_velocity
    public :: interp_u_column
    public :: interp_u_to_layers

contains

    subroutine elsa_map_init(map,x_src,y_src,zeta,stagger,grid_factor)
        ! Define elsa's grid as a coarsening of the host's by `grid_factor` (a
        ! real >= 1, need not be an integer), and build every weight stencil.

        type(elsa_map_class), intent(inout) :: map
        real(wp),             intent(in)    :: x_src(:), y_src(:)
        real(wp),             intent(in)    :: zeta(:)
        character(len=*),     intent(in)    :: stagger
        real(wp),             intent(in)    :: grid_factor

        integer  :: i, j, nx_src, ny_src
        real(wp) :: dx_src, dy_src, x0, y0
        real(wp), allocatable :: xu_src(:), yu_src(:), xv_src(:), yv_src(:)
        real(wp), allocatable :: x_acx(:), y_acy(:)

        nx_src = size(x_src)
        ny_src = size(y_src)

        call check_axis(x_src,"x_src")
        call check_axis(y_src,"y_src")
        call check_zeta(zeta,size(zeta))

        if (grid_factor .lt. 1.0_wp) then
            write(*,*) "elsa_map_init:: Error: grid_factor must be >= 1, got ", grid_factor
            error stop 1
        end if

        dx_src = x_src(2) - x_src(1)
        dy_src = y_src(2) - y_src(1)

        map%nz = size(zeta)
        if (allocated(map%zeta)) deallocate(map%zeta)
        allocate(map%zeta(map%nz))
        map%zeta = zeta

        ! elsa's grid tiles the source domain from its outer edge.
        map%dx = grid_factor*dx_src
        map%dy = grid_factor*dy_src
        map%nx = max(1, floor(real(nx_src,wp)/grid_factor))
        map%ny = max(1, floor(real(ny_src,wp)/grid_factor))

        if (allocated(map%x)) deallocate(map%x)
        if (allocated(map%y)) deallocate(map%y)
        allocate(map%x(map%nx),map%y(map%ny))

        x0 = x_src(1) - 0.5_wp*dx_src        ! outer edge of the source domain
        y0 = y_src(1) - 0.5_wp*dy_src

        do i = 1, map%nx
            map%x(i) = x0 + (real(i,wp) - 0.5_wp)*map%dx
        end do
        do j = 1, map%ny
            map%y(j) = y0 + (real(j,wp) - 0.5_wp)*map%dy
        end do

        ! Where the host's velocity samples actually sit. This is the whole of
        ! the staggering treatment.
        allocate(xu_src(nx_src),yu_src(ny_src),xv_src(nx_src),yv_src(ny_src))

        select case (trim(stagger))
            case ("acx_acy")
                xu_src = x_src + 0.5_wp*dx_src
                yu_src = y_src
                xv_src = x_src
                yv_src = y_src + 0.5_wp*dy_src
            case ("aa")
                xu_src = x_src
                yu_src = y_src
                xv_src = x_src
                yv_src = y_src
            case default
                write(*,*) "elsa_map_init:: Error: unknown stagger '"//trim(stagger)//"'"
                write(*,*) "  expected 'acx_acy' (velocities on staggered faces) or 'aa' (cell-centred)"
                error stop 1
        end select

        ! elsa's own face coordinates.
        allocate(x_acx(map%nx),y_acy(map%ny))
        x_acx = map%x + 0.5_wp*map%dx
        y_acy = map%y + 0.5_wp*map%dy

        call axis_map_conservative(map%cons_x,x_src,dx_src,map%x,map%dx)
        call axis_map_conservative(map%cons_y,y_src,dy_src,map%y,map%dy)

        call axis_map_linear(map%ux_x,xu_src,x_acx)
        call axis_map_linear(map%ux_y,yu_src,map%y)

        call axis_map_linear(map%uy_x,xv_src,map%x)
        call axis_map_linear(map%uy_y,yv_src,y_acy)

    end subroutine elsa_map_init

    subroutine elsa_map_end(map)
        type(elsa_map_class), intent(inout) :: map

        if (allocated(map%x))    deallocate(map%x)
        if (allocated(map%y))    deallocate(map%y)
        if (allocated(map%zeta)) deallocate(map%zeta)

        call axis_map_end(map%cons_x)
        call axis_map_end(map%cons_y)
        call axis_map_end(map%ux_x)
        call axis_map_end(map%ux_y)
        call axis_map_end(map%uy_x)
        call axis_map_end(map%uy_y)

    end subroutine elsa_map_end

    ! ======================================================================
    ! Horizontal
    ! ======================================================================

    subroutine map_scalar(map,f_src,f)
        ! Area-weighted mean of a host field onto elsa's cells. Exact, and it
        ! degenerates to plain box averaging when grid_factor is an integer.
        type(elsa_map_class), intent(in)  :: map
        real(wp),             intent(in)  :: f_src(:,:)
        real(wp),             intent(out) :: f(:,:)

        call map_field(map%cons_x,map%cons_y,f_src,f)

    end subroutine map_scalar

    subroutine map_velocity(map,ux_src,uy_src,ux,uy)
        ! Host velocity onto elsa's acx/acy faces, level by level. The outermost
        ! face in each direction carries zero flux (see docs/DESIGN.md), so its
        ! velocity is never read; it is zeroed rather than extrapolated.
        type(elsa_map_class), intent(in)  :: map
        real(wp),             intent(in)  :: ux_src(:,:,:), uy_src(:,:,:)
        real(wp),             intent(out) :: ux(:,:,:), uy(:,:,:)

        integer :: k

        do k = 1, map%nz
            call map_field(map%ux_x,map%ux_y,ux_src(:,:,k),ux(:,:,k))
            call map_field(map%uy_x,map%uy_y,uy_src(:,:,k),uy(:,:,k))
        end do

        ux(map%nx,:,:) = 0.0_wp
        uy(:,map%ny,:) = 0.0_wp

    end subroutine map_velocity

    subroutine map_field(mx,my,f_src,f)
        ! Apply the separable stencil: f(i,j) = sum_ab wx(a,i) wy(b,j) f_src(..).
        type(axis_map_class), intent(in)  :: mx, my
        real(wp),             intent(in)  :: f_src(:,:)
        real(wp),             intent(out) :: f(:,:)

        integer  :: i, j, a, b
        real(wp) :: acc, wy

        do j = 1, my%n_tgt
        do i = 1, mx%n_tgt
            acc = 0.0_wp
            do b = 1, my%stencil
                wy = my%wt(b,j)
                if (wy .eq. 0.0_wp) cycle
                do a = 1, mx%stencil
                    acc = acc + mx%wt(a,i)*wy*f_src(mx%idx(a,i),my%idx(b,j))
                end do
            end do
            f(i,j) = acc
        end do
        end do

    end subroutine map_field

    ! ======================================================================
    ! Vertical
    ! ======================================================================

    subroutine interp_u_to_layers(u_layer,u_lev,zeta,H,dsum,n_top)
        ! Layer-mean velocity for every column of a face grid. H and dsum must
        ! already be evaluated on that face grid (the caller staggers them).
        real(wp), intent(out) :: u_layer(:,:,:)
        real(wp), intent(in)  :: u_lev(:,:,:)
        real(wp), intent(in)  :: zeta(:)
        real(wp), intent(in)  :: H(:,:)
        real(wp), intent(in)  :: dsum(:,:,:)
        integer,  intent(in)  :: n_top

        integer :: i, j

        do j = 1, size(H,2)
        do i = 1, size(H,1)
            call interp_u_column(u_layer(i,j,1:n_top),u_lev(i,j,:),zeta,H(i,j),dsum(i,j,1:n_top))
        end do
        end do

    end subroutine interp_u_to_layers

    subroutine interp_u_column(u_layer,u_lev,zeta,H,dsum)
        ! Thickness-average of a piecewise-linear velocity profile over each
        ! isochronal layer:
        !
        !     u_layer(k) = 1/d(k) * Integral over [dsum(k-1), dsum(k)] of u(z) dz
        !
        ! v2.0 instead sampled u at dsum(k), the layer's *upper* interface, which
        ! biases every layer toward the faster ice above it.
        !
        ! Exact for a linear profile, in which case the result is u at the layer
        ! midpoint. A single monotone walk over the source levels keeps this
        ! O(n_top + nz) rather than O(n_top * nz).

        real(wp), intent(out) :: u_layer(:)
        real(wp), intent(in)  :: u_lev(:), zeta(:)
        real(wp), intent(in)  :: H
        real(wp), intent(in)  :: dsum(:)

        integer  :: k, l, nz, n_top
        real(wp) :: z(size(zeta)), fcum(size(zeta))
        real(wp) :: za, zb, fa, fb, dz, t, u_at

        nz    = size(zeta)
        n_top = size(dsum)

        u_layer = 0.0_wp
        if (H .le. H_ICE_MIN) return

        ! Source level heights above the bed, and the running integral of u to each.
        z = zeta*H

        fcum(1) = 0.0_wp
        do l = 2, nz
            fcum(l) = fcum(l-1) + 0.5_wp*(u_lev(l-1)+u_lev(l))*(z(l)-z(l-1))
        end do

        l  = 1
        za = 0.0_wp
        fa = 0.0_wp

        do k = 1, n_top

            zb = min(dsum(k),z(nz))     ! dsum(n_top) == H up to roundoff

            do while (l .lt. nz-1 .and. z(l+1) .lt. zb)
                l = l + 1
            end do

            t     = (zb - z(l))/(z(l+1) - z(l))
            t     = min(max(t,0.0_wp),1.0_wp)
            u_at  = u_lev(l) + t*(u_lev(l+1) - u_lev(l))
            fb    = fcum(l) + 0.5_wp*(u_lev(l) + u_at)*(zb - z(l))

            dz = zb - za
            if (dz .gt. DZ_MIN) then
                u_layer(k) = (fb - fa)/dz
            else
                u_layer(k) = u_at       ! vanishing layer: use its own position
            end if

            za = zb
            fa = fb

        end do

    end subroutine interp_u_column

    ! ======================================================================
    ! Axis maps
    ! ======================================================================

    subroutine axis_map_conservative(m,xs,dxs,xt,dxt)
        ! Overlap length of each source cell with each target cell, normalized.
        ! The target cell's value is the area-weighted mean of the source cells
        ! it covers. Where a target cell overhangs the source domain, it is the
        ! mean over the overlapped part only.
        type(axis_map_class), intent(inout) :: m
        real(wp),             intent(in)    :: xs(:), dxs, xt(:), dxt

        integer  :: ns, nt, i, t, i0, i1, n
        real(wp) :: lo, hi, s_lo, s_hi, ov, wsum, x_edge

        ns = size(xs)
        nt = size(xt)

        call axis_map_alloc(m,nt,ceiling(dxt/dxs) + 1)

        x_edge = xs(1) - 0.5_wp*dxs

        do t = 1, nt

            lo = xt(t) - 0.5_wp*dxt
            hi = xt(t) + 0.5_wp*dxt

            i0 = max(1, floor((lo - x_edge)/dxs) + 1)
            i1 = min(ns, ceiling((hi - x_edge)/dxs))

            n    = 0
            wsum = 0.0_wp

            do i = i0, i1
                s_lo = xs(i) - 0.5_wp*dxs
                s_hi = xs(i) + 0.5_wp*dxs
                ov   = min(hi,s_hi) - max(lo,s_lo)
                if (ov .le. 0.0_wp) cycle

                n = n + 1
                if (n .gt. m%stencil) then
                    write(*,*) "axis_map_conservative:: Error: stencil overflow at target ", t
                    error stop 1
                end if
                m%idx(n,t) = i
                m%wt(n,t)  = ov
                wsum       = wsum + ov
            end do

            if (wsum .gt. 0.0_wp) then
                m%wt(1:n,t) = m%wt(1:n,t)/wsum
            else
                ! Target cell lies entirely outside the source domain.
                m%idx(1,t) = min(max(1, nint((xt(t)-x_edge)/dxs)), ns)
                m%wt(1,t)  = 1.0_wp
            end if

        end do

    end subroutine axis_map_conservative

    subroutine axis_map_linear(m,xs,xt)
        ! Linear interpolation weights, clamped to the end values outside the
        ! source axis. xs must be ascending; it need not be uniform.
        type(axis_map_class), intent(inout) :: m
        real(wp),             intent(in)    :: xs(:), xt(:)

        integer  :: ns, nt, t, i
        real(wp) :: a

        ns = size(xs)
        nt = size(xt)

        call axis_map_alloc(m,nt,2)

        do t = 1, nt

            if (xt(t) .le. xs(1)) then
                m%idx(:,t) = [1,min(2,ns)]
                m%wt(:,t)  = [1.0_wp,0.0_wp]
            else if (xt(t) .ge. xs(ns)) then
                m%idx(:,t) = [max(1,ns-1),ns]
                m%wt(:,t)  = [0.0_wp,1.0_wp]
            else
                i = bracket(xs,xt(t))
                a = (xt(t) - xs(i))/(xs(i+1) - xs(i))
                m%idx(:,t) = [i,i+1]
                m%wt(:,t)  = [1.0_wp-a,a]
            end if

        end do

    end subroutine axis_map_linear

    function bracket(xs,x) result(i)
        ! Largest i with xs(i) <= x, for xs(1) < x < xs(n). Binary search.
        real(wp), intent(in) :: xs(:), x
        integer :: i

        integer :: lo, hi, mid

        lo = 1
        hi = size(xs)

        do while (hi - lo .gt. 1)
            mid = (lo + hi)/2
            if (xs(mid) .le. x) then
                lo = mid
            else
                hi = mid
            end if
        end do

        i = lo

    end function bracket

    subroutine axis_map_alloc(m,n_tgt,stencil)
        type(axis_map_class), intent(inout) :: m
        integer,              intent(in)    :: n_tgt, stencil

        call axis_map_end(m)

        m%n_tgt   = n_tgt
        m%stencil = stencil
        allocate(m%idx(stencil,n_tgt))
        allocate(m%wt(stencil,n_tgt))

        m%idx = 1
        m%wt  = 0.0_wp

    end subroutine axis_map_alloc

    subroutine axis_map_end(m)
        type(axis_map_class), intent(inout) :: m

        if (allocated(m%idx)) deallocate(m%idx)
        if (allocated(m%wt))  deallocate(m%wt)
        m%n_tgt   = 0
        m%stencil = 0

    end subroutine axis_map_end

    ! ======================================================================
    ! Validation of the host contract
    ! ======================================================================

    subroutine check_axis(x,name)
        real(wp),         intent(in) :: x(:)
        character(len=*), intent(in) :: name

        integer  :: i
        real(wp) :: dx

        if (size(x) .lt. 2) then
            write(*,*) "elsa_map_init:: Error: "//name//" needs at least 2 points"
            error stop 1
        end if

        dx = x(2) - x(1)
        if (dx .le. 0.0_wp) then
            write(*,*) "elsa_map_init:: Error: "//name//" must be ascending"
            error stop 1
        end if

        do i = 2, size(x)
            if (abs((x(i)-x(i-1)) - dx) .gt. 1.0e-6_wp*dx) then
                write(*,*) "elsa_map_init:: Error: "//name//" must be uniformly spaced, breaks at i = ", i
                error stop 1
            end if
        end do

    end subroutine check_axis

    subroutine check_zeta(zeta,nz)
        ! elsa's host contract: velocities live on nz sigma levels, ascending,
        ! zeta(1) = 0 at the bed and zeta(nz) = 1 at the surface. This is Yelmo's
        ! zeta_aa convention and the one `tracer` states.
        real(wp), intent(in) :: zeta(:)
        integer,  intent(in) :: nz

        integer :: k

        if (nz .lt. 2) then
            write(*,*) "elsa_map_init:: Error: zeta needs at least 2 levels"
            error stop 1
        end if

        if (abs(zeta(1)) .gt. 1.0e-8_wp .or. abs(zeta(nz)-1.0_wp) .gt. 1.0e-8_wp) then
            write(*,*) "elsa_map_init:: Error: zeta must run 0 (bed) to 1 (surface)."
            write(*,*) "  got zeta(1) = ", zeta(1), " zeta(nz) = ", zeta(nz)
            error stop 1
        end if

        do k = 2, nz
            if (zeta(k) .le. zeta(k-1)) then
                write(*,*) "elsa_map_init:: Error: zeta must be strictly ascending, breaks at k = ", k
                error stop 1
            end if
        end do

    end subroutine check_zeta

end module elsa_interp
