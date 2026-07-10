program test_physics
    ! Benchmarks for the elsa_physics kernels. Each check asserts; any failure
    ! exits nonzero so `make check` fails.
    !
    ! The advection cases test the three properties the scheme is supposed to
    ! guarantee unconditionally — non-negativity, exact mass conservation, and
    ! (for uniform flow) exact transport of the first moment — and specifically
    ! exercise the convergent and divergent cells that made ELSA v2.0 unstable.

    use elsa_precision, only : wp
    use elsa_physics

    implicit none

    integer :: n_fail

    n_fail = 0

    write(*,*) ""
    write(*,*) "== elsa_physics =="

    call test_uniform_translation(n_fail)
    call test_solid_body_rotation(n_fail)
    call test_sign_reversal(n_fail)
    call test_layer_bookkeeping(n_fail)

    write(*,*) ""
    if (n_fail .gt. 0) then
        write(*,'(a,i0,a)') "  ", n_fail, " check(s) FAILED"
        write(*,*) ""
        error stop 1
    end if

    write(*,*) "  all checks passed"
    write(*,*) ""

contains

    subroutine check(ok,name,n_fail)
        logical,          intent(in)    :: ok
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: n_fail

        if (ok) then
            write(*,'(a,a)') "   pass   ", name
        else
            write(*,'(a,a)') "   FAIL   ", name
            n_fail = n_fail + 1
        end if

    end subroutine check

    subroutine check_close(val,ref,tol,name,n_fail)
        ! Relative comparison, falling back to absolute when ref is ~0.
        real(wp),         intent(in)    :: val, ref, tol
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: n_fail

        real(wp) :: err

        err = abs(val-ref)
        if (abs(ref) .gt. 1.0e-30_wp) err = err/abs(ref)

        if (err .le. tol) then
            write(*,'(a,a,a,es9.2,a,es8.1,a)') "   pass   ", name, "  (err ", err, " <= ", tol, ")"
        else
            write(*,'(a,a,a,es9.2,a,es8.1,a)') "   FAIL   ", name, "  (err ", err, " >  ", tol, ")"
            write(*,'(a,es22.14,a,es22.14)')   "            got ", val, "  want ", ref
            n_fail = n_fail + 1
        end if

    end subroutine check_close

    subroutine gaussian_blob(d,x,y,x0,y0,sigma)
        real(wp), intent(out) :: d(:,:)
        real(wp), intent(in)  :: x(:), y(:), x0, y0, sigma

        integer :: i, j

        do j = 1, size(d,2)
        do i = 1, size(d,1)
            d(i,j) = exp( -((x(i)-x0)**2 + (y(j)-y0)**2) / (2.0_wp*sigma**2) )
        end do
        end do

    end subroutine gaussian_blob

    function centroid_x(d,x) result(xbar)
        real(wp), intent(in) :: d(:,:), x(:)
        real(wp) :: xbar

        integer :: i
        real(wp) :: num

        num = 0.0_wp
        do i = 1, size(d,1)
            num = num + x(i)*sum(d(i,:))
        end do
        xbar = num/sum(d)

    end function centroid_x

    subroutine test_uniform_translation(n_fail)
        ! A blob in a uniform flow. Upwind advection transports the first moment
        ! of d exactly, whatever the numerical diffusion does to its shape, so
        ! the centroid displacement must equal u*dt to roundoff.
        integer, intent(inout) :: n_fail

        integer,  parameter :: nx = 161, ny = 41
        real(wp), parameter :: dx = 1000.0_wp, dy = 1000.0_wp
        real(wp), parameter :: uu = 100.0_wp        ! [m/yr]
        real(wp), parameter :: dt = 300.0_wp        ! [yr]
        real(wp), parameter :: cfl = 0.9_wp

        real(wp) :: d(nx,ny), ux(nx,ny), uy(nx,ny), fx(nx,ny), fy(nx,ny)
        real(wp) :: x(nx), y(ny)
        real(wp) :: m0, m1, xbar0, xbar1
        integer  :: i, j, n_sub

        write(*,*) ""
        write(*,*) " uniform translation"

        do i = 1, nx
            x(i) = real(i-1,wp)*dx
        end do
        do j = 1, ny
            y(j) = real(j-1,wp)*dy
        end do

        ! sigma = 4 cells, centred 40 cells in: the tails are ~1e-21 at the
        ! edge, so the zero-flux boundary never truncates the blob.
        call gaussian_blob(d,x,y,x(40),y(21),4.0_wp*dx)

        ux = uu
        uy = 0.0_wp

        m0    = sum(d)
        xbar0 = centroid_x(d,x)

        n_sub = calc_n_substeps(ux,uy,dx,dy,dt,cfl)
        write(*,'(a,i0)') "   n_sub = ", n_sub

        call advect_layer(d,ux,uy,dx,dy,dt,cfl,fx,fy)

        m1    = sum(d)
        xbar1 = centroid_x(d,x)

        call check_close(m1,m0,1.0e-13_wp,"mass conserved            ",n_fail)
        call check(minval(d) .ge. 0.0_wp,"non-negative              ",n_fail)
        call check_close(xbar1-xbar0,uu*dt,1.0e-10_wp,"centroid transported exact",n_fail)

    end subroutine test_uniform_translation

    subroutine test_solid_body_rotation(n_fail)
        ! Rigid rotation. The flow is non-uniform and strongly sheared, so the
        ! first moment is not preserved, but mass and positivity must be.
        integer, intent(inout) :: n_fail

        integer,  parameter :: nx = 101, ny = 101
        real(wp), parameter :: dx = 1000.0_wp, dy = 1000.0_wp
        real(wp), parameter :: t_rev = 2000.0_wp    ! [yr] one full revolution
        real(wp), parameter :: cfl = 0.9_wp
        real(wp), parameter :: pi = 3.14159265358979323846_wp

        real(wp) :: d(nx,ny), ux(nx,ny), uy(nx,ny), fx(nx,ny), fy(nx,ny)
        real(wp) :: x(nx), y(ny)
        real(wp) :: m0, m1, peak0, peak1, omega, xc, yc
        integer  :: i, j, n_sub

        write(*,*) ""
        write(*,*) " solid-body rotation"

        do i = 1, nx
            x(i) = real(i-1,wp)*dx
        end do
        do j = 1, ny
            y(j) = real(j-1,wp)*dy
        end do

        xc = x(51)
        yc = y(51)

        ! Blob offset from the axis, well inside the domain.
        call gaussian_blob(d,x,y,x(31),y(51),5.0_wp*dx)

        omega = 2.0_wp*pi/t_rev

        ! ux lives on the x face at x(i)+dx/2; uy on the y face at y(j)+dy/2.
        do j = 1, ny
        do i = 1, nx
            ux(i,j) = -omega*(y(j) - yc)
            uy(i,j) =  omega*(x(i) - xc)
        end do
        end do

        m0    = sum(d)
        peak0 = maxval(d)

        n_sub = calc_n_substeps(ux,uy,dx,dy,t_rev,cfl)
        write(*,'(a,i0)') "   n_sub = ", n_sub

        call advect_layer(d,ux,uy,dx,dy,t_rev,cfl,fx,fy)

        m1    = sum(d)
        peak1 = maxval(d)

        call check_close(m1,m0,1.0e-12_wp,"mass conserved            ",n_fail)
        call check(minval(d) .ge. 0.0_wp,"non-negative              ",n_fail)

        ! First-order upwind is heavily diffusive over a full revolution; this
        ! is reported, not asserted, and is the number to beat if the scheme is
        ! ever raised to higher order.
        write(*,'(a,f6.3)') "   peak amplitude retained after one revolution: ", peak1/peak0

    end subroutine test_solid_body_rotation

    subroutine test_sign_reversal(n_fail)
        ! The regression test for the v2.0 defect.
        !
        ! Convergent: u(i-1/2) > 0 > u(i+1/2). v2.0 selected its upwind branch
        ! from the sign of u(i+1/2) alone, which put a *negative* entry on the
        ! main diagonal and made the solve diverge. It was masked with a cap on
        ! the layer thickness change per step (R24 Eq. B9), a non-finite
        ! fallback, and a |u| < 0.1 m/yr cutoff.
        !
        ! Divergent is the mirror case. Both are run at a CFL that would be far
        ! outside the explicit stability limit if it were not sub-stepped.
        integer, intent(inout) :: n_fail

        integer,  parameter :: nx = 41, ny = 41
        real(wp), parameter :: dx = 1000.0_wp, dy = 1000.0_wp
        real(wp), parameter :: dt = 100.0_wp
        real(wp), parameter :: cfl = 0.9_wp
        real(wp), parameter :: u0 = 5000.0_wp       ! [m/yr] at the domain edge

        real(wp) :: d(nx,ny), ux(nx,ny), uy(nx,ny), fx(nx,ny), fy(nx,ny)
        real(wp) :: x(nx), y(ny), xc, yc, ell
        real(wp) :: m0, m1, cfl_raw
        integer  :: i, j, s, n_sub

        do i = 1, nx
            x(i) = real(i-1,wp)*dx
        end do
        do j = 1, ny
            y(j) = real(j-1,wp)*dy
        end do

        xc  = x(21)
        yc  = y(21)
        ell = 20.0_wp*dx

        do s = 1, 2      ! s=1 convergent, s=2 divergent

            if (s .eq. 1) then
                write(*,*) ""
                write(*,*) " convergent flow (v2.0 regression)"
            else
                write(*,*) ""
                write(*,*) " divergent flow (v2.0 regression)"
            end if

            d = 100.0_wp

            do j = 1, ny
            do i = 1, nx
                ! Face velocities pointing toward (s=1) / away from (s=2) centre.
                ux(i,j) = u0*(xc - (x(i)+0.5_wp*dx))/ell
                uy(i,j) = u0*(yc - (y(j)+0.5_wp*dy))/ell
                if (s .eq. 2) then
                    ux(i,j) = -ux(i,j)
                    uy(i,j) = -uy(i,j)
                end if
            end do
            end do

            m0 = sum(d)

            n_sub   = calc_n_substeps(ux,uy,dx,dy,dt,cfl)
            cfl_raw = real(n_sub,wp)*cfl
            write(*,'(a,f7.1,a,i0)') "   unsubstepped CFL ~ ", cfl_raw, " ; n_sub = ", n_sub

            call advect_layer(d,ux,uy,dx,dy,dt,cfl,fx,fy)

            m1 = sum(d)

            call check(all(d .eq. d),                  "all finite (no NaN)       ",n_fail)
            call check(minval(d) .ge. 0.0_wp,          "non-negative              ",n_fail)
            call check_close(m1,m0,1.0e-12_wp,         "mass conserved            ",n_fail)
            ! No cap on layer thickness change is applied anywhere; if the
            ! scheme were still the v2.0 one, d would be unbounded here.
            call check(maxval(d) .lt. 1.0e6_wp,        "bounded (no cap applied)  ",n_fail)

        end do

    end subroutine test_sign_reversal

    subroutine test_layer_bookkeeping(n_fail)
        ! normalize_layers, apply_smb, apply_bmb.
        integer, intent(inout) :: n_fail

        integer,  parameter :: nx = 3, ny = 2, nl = 4

        real(wp) :: d(nx,ny,nl), H_ice(nx,ny), dm(nx,ny)
        integer  :: n_reseed
        real(wp) :: m_before

        write(*,*) ""
        write(*,*) " layer bookkeeping"

        ! -- normalize_layers ------------------------------------------------
        d     = 1.0_wp
        H_ice = 800.0_wp
        H_ice(3,2) = 0.0_wp            ! ice-free column

        call normalize_layers(d,H_ice,nl)

        call check_close(sum(d(1,1,:)),800.0_wp,1.0e-14_wp,"column normalized to H    ",n_fail)
        call check(all(d(3,2,:) .eq. 0.0_wp),              "ice-free column zeroed    ",n_fail)

        ! -- apply_smb: accumulation lands on the top layer -------------------
        d  = 10.0_wp
        dm = 5.0_wp
        call apply_smb(d,dm,nl)

        call check_close(d(1,1,nl),15.0_wp,1.0e-14_wp,"smb+ added to top layer   ",n_fail)
        call check_close(d(1,1,1), 10.0_wp,1.0e-14_wp,"smb+ leaves base alone    ",n_fail)

        ! -- apply_smb: ablation eats the top layers downward -----------------
        d  = 10.0_wp
        dm = -25.0_wp                  ! two full layers plus half of the third
        m_before = sum(d(1,1,:))
        call apply_smb(d,dm,nl)

        call check_close(d(1,1,nl),  0.0_wp,1.0e-14_wp,"smb- empties top layer    ",n_fail)
        call check_close(d(1,1,nl-1),0.0_wp,1.0e-14_wp,"smb- empties next layer   ",n_fail)
        call check_close(d(1,1,nl-2),5.0_wp,1.0e-14_wp,"smb- partial third layer  ",n_fail)
        call check_close(sum(d(1,1,:)),m_before-25.0_wp,1.0e-14_wp,"smb- removes exactly |dm| ",n_fail)

        ! -- apply_smb: ablation cannot remove more than is there -------------
        d  = 1.0_wp
        dm = -999.0_wp
        call apply_smb(d,dm,nl)
        call check(minval(d) .ge. 0.0_wp,             "smb- cannot go negative   ",n_fail)
        call check_close(sum(d(1,1,:)),0.0_wp,1.0e-14_wp,"smb- exhausts the column  ",n_fail)

        ! -- reseed_empty_columns ---------------------------------------------
        ! A column the surface mass balance ate through, where the host still has
        ! ice, must come back at the bed rather than stay dead forever.
        d          = 10.0_wp
        H_ice      = 800.0_wp
        H_ice(3,2) = 0.0_wp                          ! genuinely ice-free
        dm         = -999.0_wp
        call apply_smb(d,dm,nl)                      ! annihilates every column

        call reseed_empty_columns(d,H_ice,nl,n_reseed)

        call check(n_reseed .eq. nx*ny - 1,           "reseeds every iced column ",n_fail)
        call check_close(d(1,1,1),800.0_wp,1.0e-14_wp,"reseed lands at the bed   ",n_fail)
        call check(all(d(1,1,2:nl) .eq. 0.0_wp),      "reseed leaves layers above",n_fail)
        call check(all(d(3,2,:) .eq. 0.0_wp),         "ice-free column not seeded",n_fail)

        ! A healthy column must not be touched.
        d = 10.0_wp
        call reseed_empty_columns(d,H_ice,nl,n_reseed)
        call check(n_reseed .eq. 0,                   "healthy column untouched  ",n_fail)

        ! -- apply_bmb: freeze-on gated by allow_pos_bmb ----------------------
        d  = 10.0_wp
        dm = 5.0_wp
        call apply_bmb(d,dm,nl,allow_pos_bmb=.false.)
        call check_close(d(1,1,1),10.0_wp,1.0e-14_wp,"bmb+ ignored when disabled",n_fail)

        call apply_bmb(d,dm,nl,allow_pos_bmb=.true.)
        call check_close(d(1,1,1),15.0_wp,1.0e-14_wp,"bmb+ thickens base layer  ",n_fail)

        ! -- apply_bmb: melt eats the bottom layers upward --------------------
        d  = 10.0_wp
        dm = -25.0_wp
        call apply_bmb(d,dm,nl,allow_pos_bmb=.true.)

        call check_close(d(1,1,1),0.0_wp,1.0e-14_wp,"bmb- empties base layer   ",n_fail)
        call check_close(d(1,1,2),0.0_wp,1.0e-14_wp,"bmb- empties next layer   ",n_fail)
        call check_close(d(1,1,3),5.0_wp,1.0e-14_wp,"bmb- partial third layer  ",n_fail)
        call check(minval(d) .ge. 0.0_wp,            "bmb- cannot go negative   ",n_fail)

    end subroutine test_layer_bookkeeping

end program test_physics
