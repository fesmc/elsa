program test_interp
    ! Benchmarks for elsa_interp.
    !
    ! Each map has an exactness property it must satisfy, and those are what is
    ! asserted here rather than a tolerance pulled from the air:
    !
    !   conservative   partition of unity; total mass unchanged when elsa's grid
    !                  tiles the host's exactly, integer ratio or not
    !   bilinear       reproduces a linear field exactly
    !   layer-mean     reproduces a linear velocity profile exactly, giving the
    !                  value at each layer's midpoint; and the thickness-weighted
    !                  sum equals the exact integral of the profile

    use elsa_precision, only : wp
    use elsa_interp

    implicit none

    integer :: n_fail

    n_fail = 0

    write(*,*) ""
    write(*,*) "== elsa_interp =="

    call test_conservative(n_fail)
    call test_bilinear_velocity(n_fail)
    call test_layer_mean(n_fail)

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

    subroutine make_axes(x,y,n,d)
        real(wp), allocatable, intent(out) :: x(:), y(:)
        integer,  intent(in) :: n
        real(wp), intent(in) :: d

        integer :: i

        allocate(x(n),y(n))
        do i = 1, n
            x(i) = (real(i,wp) - 0.5_wp)*d
            y(i) = (real(i,wp) - 0.5_wp)*d
        end do

    end subroutine make_axes

    subroutine test_conservative(n_fail)
        ! Host grid: 10 x 10 cells of size 1, spanning [0,10] x [0,10].
        ! grid_factor 2.0 and 2.5 both tile that domain exactly (5 and 4 cells),
        ! so total mass must be unchanged in both.
        integer, intent(inout) :: n_fail

        integer,  parameter :: ns = 10
        real(wp), parameter :: ds = 1.0_wp

        type(elsa_map_class)  :: map
        real(wp), allocatable :: x(:), y(:), zeta(:)
        real(wp)              :: f_src(ns,ns)
        real(wp), allocatable :: f(:,:)
        real(wp)              :: m_src, m_tgt
        integer               :: i, j, gf_case
        real(wp)              :: gf

        write(*,*) ""
        write(*,*) " conservative remap"

        call make_axes(x,y,ns,ds)
        allocate(zeta(3))
        zeta = [0.0_wp,0.5_wp,1.0_wp]

        do j = 1, ns
        do i = 1, ns
            f_src(i,j) = real(i,wp) + 10.0_wp*real(j,wp)
        end do
        end do

        m_src = sum(f_src)*ds*ds

        do gf_case = 1, 3

            select case (gf_case)
                case (1) ; gf = 1.0_wp
                case (2) ; gf = 2.0_wp
                case (3) ; gf = 2.5_wp
            end select

            call elsa_map_init(map,x,y,zeta,"aa",gf)

            if (allocated(f)) deallocate(f)
            allocate(f(map%nx,map%ny))

            ! Partition of unity: a constant field must map to that constant.
            call map_scalar(map,spread(spread(7.0_wp,1,ns),2,ns),f)
            call check(all(abs(f-7.0_wp) .lt. 1.0e-13_wp),"constant preserved        ",n_fail)

            call map_scalar(map,f_src,f)
            m_tgt = sum(f)*map%dx*map%dy

            write(*,'(a,f4.1,a,i0,a,i0)') "   grid_factor ", gf, " -> ", map%nx, " x ", map%ny
            call check_close(m_tgt,m_src,1.0e-13_wp,"total mass conserved      ",n_fail)

            if (gf_case .eq. 1) then
                call check(all(abs(f-f_src) .lt. 1.0e-13_wp),"grid_factor=1 is identity ",n_fail)
            end if

            if (gf_case .eq. 2) then
                ! Integer ratio must degenerate to a plain 2x2 box average.
                call check_close(f(1,1),0.25_wp*(f_src(1,1)+f_src(2,1)+f_src(1,2)+f_src(2,2)), &
                                 1.0e-14_wp,"integer ratio = box mean  ",n_fail)
            end if

            call elsa_map_end(map)

        end do

    end subroutine test_conservative

    subroutine test_bilinear_velocity(n_fail)
        ! A field linear in x and y must be reproduced exactly at elsa's faces,
        ! for either staggering of the source velocity. The outermost face is
        ! zeroed by convention, so it is excluded from the comparison.
        integer, intent(inout) :: n_fail

        integer,  parameter :: ns = 10, nz = 4
        real(wp), parameter :: ds = 1.0_wp
        real(wp), parameter :: a = 3.0_wp, b = -2.0_wp, c = 0.5_wp

        type(elsa_map_class)  :: map
        real(wp), allocatable :: x(:), y(:), zeta(:)
        real(wp)              :: ux_src(ns,ns,nz), uy_src(ns,ns,nz)
        real(wp), allocatable :: ux(:,:,:), uy(:,:,:)
        real(wp)              :: xu, yv, err_ux, err_uy, want
        integer               :: i, j, k, s
        character(len=8)      :: stag

        write(*,*) ""
        write(*,*) " bilinear velocity onto faces"

        call make_axes(x,y,ns,ds)
        allocate(zeta(nz))
        do k = 1, nz
            zeta(k) = real(k-1,wp)/real(nz-1,wp)
        end do

        do s = 1, 2

            if (s .eq. 1) then
                stag = "acx_acy"
            else
                stag = "aa"
            end if

            ! Sample the same linear function at wherever the host's velocity
            ! nodes actually are for this staggering.
            do k = 1, nz
            do j = 1, ns
            do i = 1, ns
                if (s .eq. 1) then
                    xu = x(i) + 0.5_wp*ds
                    yv = y(j) + 0.5_wp*ds
                else
                    xu = x(i)
                    yv = y(j)
                end if
                ux_src(i,j,k) = a + b*xu   + c*y(j)
                uy_src(i,j,k) = a + b*x(i) + c*yv
            end do
            end do
            end do

            call elsa_map_init(map,x,y,zeta,trim(stag),1.0_wp)

            if (allocated(ux)) deallocate(ux,uy)
            allocate(ux(map%nx,map%ny,nz),uy(map%nx,map%ny,nz))

            call map_velocity(map,ux_src,uy_src,ux,uy)

            err_ux = 0.0_wp
            do k = 1, nz
            do j = 1, map%ny
            do i = 1, map%nx-1                     ! outermost face is zeroed
                want   = a + b*(map%x(i)+0.5_wp*map%dx) + c*map%y(j)
                err_ux = max(err_ux,abs(ux(i,j,k)-want))
            end do
            end do
            end do

            err_uy = 0.0_wp
            do k = 1, nz
            do j = 1, map%ny-1
            do i = 1, map%nx
                want   = a + b*map%x(i) + c*(map%y(j)+0.5_wp*map%dy)
                err_uy = max(err_uy,abs(uy(i,j,k)-want))
            end do
            end do
            end do

            write(*,'(a,a)') "   stagger = ", trim(stag)
            call check(err_ux .lt. 1.0e-12_wp,"ux linear field exact     ",n_fail)
            call check(err_uy .lt. 1.0e-12_wp,"uy linear field exact     ",n_fail)
            call check(all(ux(map%nx,:,:) .eq. 0.0_wp),"outer x face zeroed       ",n_fail)
            call check(all(uy(:,map%ny,:) .eq. 0.0_wp),"outer y face zeroed       ",n_fail)

            call elsa_map_end(map)

        end do

    end subroutine test_bilinear_velocity

    subroutine test_layer_mean(n_fail)
        integer, intent(inout) :: n_fail

        integer,  parameter :: nz = 12, nl = 40
        real(wp), parameter :: H = 3000.0_wp

        real(wp) :: zeta(nz), u_lev(nz), z(nz)
        real(wp) :: dsum(nl), u_layer(nl), d_layer
        real(wp) :: u0, slope, want, err, integ_exact, integ_layers, z_mid
        integer  :: k

        write(*,*) ""
        write(*,*) " layer-mean velocity"

        do k = 1, nz
            zeta(k) = real(k-1,wp)/real(nz-1,wp)
        end do
        z = zeta*H

        d_layer = H/real(nl,wp)
        do k = 1, nl
            dsum(k) = real(k,wp)*d_layer
        end do

        ! -- constant profile -------------------------------------------------
        u_lev = 42.0_wp
        call interp_u_column(u_layer,u_lev,zeta,H,dsum)
        call check(all(abs(u_layer-42.0_wp) .lt. 1.0e-12_wp),"constant profile exact    ",n_fail)

        ! -- linear profile: layer mean = value at the layer midpoint ---------
        u0    = 5.0_wp
        slope = 0.03_wp
        u_lev = u0 + slope*z

        call interp_u_column(u_layer,u_lev,zeta,H,dsum)

        err = 0.0_wp
        do k = 1, nl
            z_mid = (real(k,wp) - 0.5_wp)*d_layer
            want  = u0 + slope*z_mid
            err   = max(err,abs(u_layer(k)-want))
        end do
        call check(err .lt. 1.0e-10_wp,"linear profile -> midpoint",n_fail)

        ! The v2.0 choice sampled the layer's upper interface. Report the bias so
        ! the sign and size of the change from v2.0 are on the record.
        write(*,'(a,f8.4,a)') "   v2.0 (layer-top) bias vs layer-mean: +", slope*0.5_wp*d_layer, " m/yr"

        ! -- arbitrary profile: thickness-weighted sum == exact integral -------
        do k = 1, nz
            u_lev(k) = 20.0_wp*sin(3.0_wp*zeta(k)) + 4.0_wp*zeta(k)**2
        end do

        call interp_u_column(u_layer,u_lev,zeta,H,dsum)

        ! elsa integrates the piecewise-linear interpolant of u exactly, so the
        ! reference is the trapezoid rule over the source levels.
        integ_exact = 0.0_wp
        do k = 2, nz
            integ_exact = integ_exact + 0.5_wp*(u_lev(k-1)+u_lev(k))*(z(k)-z(k-1))
        end do

        integ_layers = sum(u_layer)*d_layer

        call check_close(integ_layers,integ_exact,1.0e-12_wp,"sum(u_k d_k) = integral   ",n_fail)

        ! -- ice-free column ---------------------------------------------------
        call interp_u_column(u_layer,u_lev,zeta,0.0_wp,dsum)
        call check(all(u_layer .eq. 0.0_wp),"ice-free column -> zero   ",n_fail)

    end subroutine test_layer_mean

end program test_interp
