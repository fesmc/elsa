module elsa_physics
    ! The isochronal layer kernels: array in, array out.
    !
    ! Nothing here holds state, opens a file, or knows about a derived type, so
    ! every routine is directly testable and safe to call from inside an OpenMP
    ! region. The layer loop lives in the caller; each routine here acts on a
    ! single layer or on an explicit layer range.
    !
    ! Index conventions, throughout:
    !
    !   d(i,j)      layer thickness [m], cell centre (aa node)
    !   ux(i,j)     velocity [m/yr] on the x face between (i,j) and (i+1,j),
    !               i.e. the acx node at x(i) + dx/2. Only i = 1..nx-1 is read.
    !   uy(i,j)     velocity [m/yr] on the y face between (i,j) and (i,j+1),
    !               i.e. the acy node at y(j) + dy/2. Only j = 1..ny-1 is read.
    !
    ! The domain edge carries zero normal flux, so global mass is conserved to
    ! roundoff. See docs/DESIGN.md.

    use elsa_precision, only : wp, H_ICE_MIN

    implicit none

    private

    public :: calc_n_substeps
    public :: advect_layer
    public :: normalize_layers
    public :: calc_dsum
    public :: apply_smb
    public :: apply_bmb

contains

    function calc_n_substeps(ux,uy,dx,dy,dt,cfl) result(n_sub)
        ! Substeps needed so that every cell satisfies the positivity condition
        ! of the explicit upwind update. The explicit step multiplies d(i,j) by
        ! (1 - dt*r(i,j)) and adds non-negative inflow, where r is the total
        ! outflow rate below. Requiring dt_sub*r <= cfl <= 1 keeps d >= 0.
        !
        ! This is the exact condition, not the looser max|u|/dx + max|v|/dy
        ! bound: a cell whose faces both carry inflow contributes r = 0.

        real(wp), intent(in) :: ux(:,:), uy(:,:)
        real(wp), intent(in) :: dx, dy, dt, cfl
        integer :: n_sub

        integer  :: i, j, nx, ny
        real(wp) :: uxp, uxm, uyp, uym, r, r_max

        nx = size(ux,1)
        ny = size(ux,2)

        r_max = 0.0_wp

        do j = 1, ny
        do i = 1, nx
            ! Face velocities bounding cell (i,j); the domain edge has no flux.
            uxp = 0.0_wp
            uxm = 0.0_wp
            uyp = 0.0_wp
            uym = 0.0_wp
            if (i .lt. nx) uxp = ux(i,j)
            if (i .gt. 1)  uxm = ux(i-1,j)
            if (j .lt. ny) uyp = uy(i,j)
            if (j .gt. 1)  uym = uy(i,j-1)

            r = (max(uxp,0.0_wp) - min(uxm,0.0_wp)) / dx &
              + (max(uyp,0.0_wp) - min(uym,0.0_wp)) / dy

            r_max = max(r_max,r)
        end do
        end do

        n_sub = 1
        if (r_max*dt .gt. 0.0_wp) n_sub = max(1, ceiling(dt*r_max/cfl))

    end function calc_n_substeps

    subroutine advect_layer(d,ux,uy,dx,dy,dt,cfl,fx,fy)
        ! Advect one isochronal layer over the interval dt, sub-stepping to hold
        ! the CFL condition. Thread-safe: the caller may run one layer per
        ! thread. fx/fy are scratch, passed in so that a threaded caller can
        ! allocate them once per thread rather than per call.

        real(wp), intent(inout) :: d(:,:)
        real(wp), intent(in)    :: ux(:,:), uy(:,:)
        real(wp), intent(in)    :: dx, dy, dt, cfl
        real(wp), intent(inout) :: fx(:,:), fy(:,:)

        integer  :: k, n_sub
        real(wp) :: dt_sub

        n_sub  = calc_n_substeps(ux,uy,dx,dy,dt,cfl)
        dt_sub = dt / real(n_sub,wp)

        do k = 1, n_sub
            call advect_substep(d,ux,uy,dx,dy,dt_sub,fx,fy)
        end do

    end subroutine advect_layer

    subroutine advect_substep(d,ux,uy,dx,dy,dt,fx,fy)
        ! One explicit upwind step of the flux-form continuity equation
        !
        !     dd/dt = -d(u d)/dx - d(v d)/dy
        !
        ! Each face is upwinded on its own velocity. This is where v2.0 chose
        ! the upwind direction for both x faces from the sign of a single one,
        ! which is what made convergent cells blow up.

        real(wp), intent(inout) :: d(:,:)
        real(wp), intent(in)    :: ux(:,:), uy(:,:)
        real(wp), intent(in)    :: dx, dy, dt
        real(wp), intent(inout) :: fx(:,:), fy(:,:)

        integer  :: i, j, nx, ny
        real(wp) :: fxm, fym

        nx = size(d,1)
        ny = size(d,2)

        ! Zero the edge fluxes once; the loops below never write them.
        fx(nx,:) = 0.0_wp
        fy(:,ny) = 0.0_wp

        do j = 1, ny
        do i = 1, nx-1
            fx(i,j) = max(ux(i,j),0.0_wp)*d(i,j) + min(ux(i,j),0.0_wp)*d(i+1,j)
        end do
        end do

        do j = 1, ny-1
        do i = 1, nx
            fy(i,j) = max(uy(i,j),0.0_wp)*d(i,j) + min(uy(i,j),0.0_wp)*d(i,j+1)
        end do
        end do

        do j = 1, ny
        do i = 1, nx
            fxm = 0.0_wp
            fym = 0.0_wp
            if (i .gt. 1) fxm = fx(i-1,j)
            if (j .gt. 1) fym = fy(i,j-1)

            d(i,j) = d(i,j) - (dt/dx)*(fx(i,j)-fxm) - (dt/dy)*(fy(i,j)-fym)
        end do
        end do

    end subroutine advect_substep

    subroutine normalize_layers(d,H_ice,n_top)
        ! Rescale each column of layers so that it sums to the host's ice
        ! thickness. This is elsa's drift control: advection and the mass
        ! balance terms move layers around, and this ties the stack back to the
        ! host state without changing the layers' relative thicknesses.

        real(wp), intent(inout) :: d(:,:,:)
        real(wp), intent(in)    :: H_ice(:,:)
        integer,  intent(in)    :: n_top

        integer  :: i, j, nx, ny
        real(wp) :: d_sum

        nx = size(d,1)
        ny = size(d,2)

        do j = 1, ny
        do i = 1, nx
            d_sum = sum(d(i,j,1:n_top))

            if (H_ice(i,j) .gt. H_ICE_MIN .and. d_sum .gt. 0.0_wp) then
                d(i,j,1:n_top) = d(i,j,1:n_top) * (H_ice(i,j)/d_sum)
            else
                d(i,j,1:n_top) = 0.0_wp
            end if
        end do
        end do

    end subroutine normalize_layers

    subroutine calc_dsum(dsum,d,n_top)
        ! Height of each layer's upper interface above the bed [m]. dsum(:,:,n_top)
        ! is the ice thickness of elsa's own layered column.

        real(wp), intent(out) :: dsum(:,:,:)
        real(wp), intent(in)  :: d(:,:,:)
        integer,  intent(in)  :: n_top

        integer :: k

        dsum(:,:,1) = d(:,:,1)
        do k = 2, n_top
            dsum(:,:,k) = dsum(:,:,k-1) + d(:,:,k)
        end do

    end subroutine calc_dsum

    subroutine apply_smb(d,dm,n_top)
        ! Surface mass balance over one coupling period, as a thickness change
        ! dm [m]. Accumulation lands on the topmost layer. Ablation is taken
        ! from the top downward, exhausting each layer before moving deeper.

        real(wp), intent(inout) :: d(:,:,:)
        real(wp), intent(in)    :: dm(:,:)
        integer,  intent(in)    :: n_top

        integer  :: i, j, k, nx, ny
        real(wp) :: to_melt, take

        nx = size(d,1)
        ny = size(d,2)

        do j = 1, ny
        do i = 1, nx

            if (dm(i,j) .ge. 0.0_wp) then
                d(i,j,n_top) = d(i,j,n_top) + dm(i,j)
            else
                to_melt = -dm(i,j)
                do k = n_top, 1, -1
                    take       = min(to_melt,d(i,j,k))
                    d(i,j,k)   = d(i,j,k) - take
                    to_melt    = to_melt - take
                    if (to_melt .le. 0.0_wp) exit
                end do
            end if

        end do
        end do

    end subroutine apply_smb

    subroutine apply_bmb(d,dm,n_top,allow_pos_bmb)
        ! Basal mass balance over one coupling period, as a thickness change
        ! dm [m]. Freeze-on thickens the bottom (initialization) layer, and only
        ! if allow_pos_bmb; no new layer is created. Melt is taken from the
        ! bottom upward.

        real(wp), intent(inout) :: d(:,:,:)
        real(wp), intent(in)    :: dm(:,:)
        integer,  intent(in)    :: n_top
        logical,  intent(in)    :: allow_pos_bmb

        integer  :: i, j, k, nx, ny
        real(wp) :: to_melt, take

        nx = size(d,1)
        ny = size(d,2)

        do j = 1, ny
        do i = 1, nx

            if (dm(i,j) .ge. 0.0_wp) then
                if (allow_pos_bmb) d(i,j,1) = d(i,j,1) + dm(i,j)
            else
                to_melt = -dm(i,j)
                do k = 1, n_top
                    take       = min(to_melt,d(i,j,k))
                    d(i,j,k)   = d(i,j,k) - take
                    to_melt    = to_melt - take
                    if (to_melt .le. 0.0_wp) exit
                end do
            end if

        end do
        end do

    end subroutine apply_bmb

end module elsa_physics
