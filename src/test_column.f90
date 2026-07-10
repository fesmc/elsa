program test_column
    ! The Nye analytic ice-divide benchmark: elsa's end-to-end quantitative test.
    !
    ! At an ice divide with no horizontal flow, constant ice thickness H and
    ! constant accumulation a, an isochrone laid down at the surface sinks under
    ! uniform vertical strain. Nye's steady solution puts it, after an elapsed
    ! time tau, at a height above the bed of
    !
    !     z(tau) = H * exp(-a*tau/H)
    !
    ! elsa never computes a vertical velocity. Each coupling step it adds a*dt to
    ! the top layer and renormalizes the column onto H, so every layer height is
    ! multiplied by r = H/(H + a*dt). After n steps,
    !
    !     z = H * r**n = H * (1 + a*dt/H)**(-tau/dt)     -> H*exp(-a*tau/H)
    !
    ! as dt -> 0. The vertical thinning is thus an emergent property of the layer
    ! bookkeeping, not something imposed. Both statements are checked: the
    ! discrete result to roundoff, and its first-order convergence onto Nye.

    use elsa

    implicit none

    integer,  parameter :: NX = 5, NY = 5, NZ = 6
    integer,  parameter :: N_INIT = 10
    real(wp), parameter :: H_CONST = 3000.0_wp      ! [m]
    real(wp), parameter :: ACC     = 0.3_wp         ! [m/yr]
    real(wp), parameter :: DX      = 1000.0_wp      ! [m]
    real(wp), parameter :: TIME_0  = 0.0_wp         ! [yr]
    real(wp), parameter :: TIME_1  = 20000.0_wp     ! [yr]

    integer :: n_fail

    n_fail = 0

    write(*,*) ""
    write(*,*) "== elsa column (Nye) =="

    call test_discrete_exact(n_fail)
    call test_nye_convergence(n_fail)

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

    subroutine run_divide(els,group,dt)
        ! Drive elsa with a steady, horizontally uniform ice divide.
        type(elsa_class), intent(inout) :: els
        character(len=*), intent(in)    :: group
        real(wp),         intent(in)    :: dt

        real(wp) :: x(NX), y(NY), zeta(NZ)
        real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)
        real(wp) :: time
        integer  :: i, n, n_steps

        do i = 1, NX
            x(i) = (real(i,wp) - 0.5_wp)*DX
        end do
        do i = 1, NY
            y(i) = (real(i,wp) - 0.5_wp)*DX
        end do
        do i = 1, NZ
            zeta(i) = real(i-1,wp)/real(NZ-1,wp)
        end do

        H_ice = H_CONST
        smb   = ACC
        bmb   = 0.0_wp
        ux    = 0.0_wp
        uy    = 0.0_wp

        call elsa_init(els,"par/test_column.nml",group,TIME_0,TIME_1,x,y,zeta,H_ice,"aa")

        n_steps = nint((TIME_1-TIME_0)/dt)
        do n = 1, n_steps
            time = TIME_0 + real(n,wp)*dt
            call elsa_update(els,time,H_ice,ux,uy,smb,bmb)
        end do

    end subroutine run_divide

    subroutine test_discrete_exact(n_fail)
        ! Every isochrone must sit exactly where the discrete recursion puts it.
        ! This exercises init, the smb accounting, the normalization and the
        ! layer insertion together; an error in any of them moves an isochrone.
        integer, intent(inout) :: n_fail

        real(wp), parameter :: DT = 100.0_wp

        type(elsa_class) :: els
        real(wp) :: r, t_create, z_want, z_got, err, err_top
        integer  :: jj, n_add, n_steps

        write(*,*) ""
        write(*,*) " discrete layer thinning"

        call run_divide(els,"column",DT)

        r     = H_CONST/(H_CONST + ACC*DT)
        n_add = size(els%par%time_add)

        ! Isochrone jj = 0 is laid down at time_init and sits at dsum(N_INIT).
        ! Isochrone jj > 0 is laid down at time_add(jj) and sits at dsum(N_INIT+jj).
        err = 0.0_wp
        do jj = 0, n_add
            if (jj .eq. 0) then
                t_create = TIME_0
            else
                t_create = els%par%time_add(jj)
            end if

            n_steps = nint((TIME_1 - t_create)/DT)
            z_want  = H_CONST * r**n_steps
            z_got   = els%now%dsum_iso(3,3,N_INIT+jj)

            err = max(err,abs(z_got-z_want)/z_want)
        end do

        write(*,'(a,i0,a,i0)')     "   isochrones: ", n_add, "   layers: ", els%now%n_top
        write(*,'(a,es9.2)')       "   max relative isochrone error: ", err

        call check(err .lt. 1.0e-11_wp,           "isochrone heights exact   ",n_fail)

        err_top = abs(els%now%dsum_iso(3,3,els%now%n_top) - H_CONST)/H_CONST
        call check(err_top .lt. 1.0e-13_wp,       "column sums to H          ",n_fail)

        call check(minval(els%now%d_iso) .ge. 0.0_wp,"layers non-negative       ",n_fail)

        ! Horizontally uniform forcing must stay horizontally uniform.
        call check(maxval(abs(els%now%dsum_iso(1,1,1:els%now%n_top) &
                            - els%now%dsum_iso(4,2,1:els%now%n_top))) .lt. 1.0e-12_wp, &
                                                  "column-to-column identical",n_fail)

        call elsa_end(els)

        ! elsa_end must free everything, so a second init on the same object works.
        ! v2.0 leaked two arrays here and aborted on the second allocate.
        call run_divide(els,"column",DT)
        call check(els%now%n_top .eq. els%par%n_layers,"re-init after end succeeds",n_fail)
        call elsa_end(els)

    end subroutine test_discrete_exact

    subroutine test_nye_convergence(n_fail)
        ! The discrete solution converges onto Nye at first order in dt, so
        ! halving dt must halve the error.
        integer, intent(inout) :: n_fail

        type(elsa_class) :: els
        real(wp) :: dt(3), err(3), z_nye, z_got, ratio
        integer  :: k
        character(len=16) :: group(3)

        write(*,*) ""
        write(*,*) " convergence onto the Nye solution"

        dt    = [200.0_wp,100.0_wp,50.0_wp]
        group = ["column_dt200    ","column_dt100    ","column_dt50     "]

        z_nye = H_CONST*exp(-ACC*(TIME_1-TIME_0)/H_CONST)

        do k = 1, 3
            call run_divide(els,trim(group(k)),dt(k))
            z_got  = els%now%dsum_iso(3,3,N_INIT)
            err(k) = abs(z_got - z_nye)
            call elsa_end(els)
            write(*,'(a,f6.1,a,f9.4,a,f8.4,a)') "   dt = ", dt(k), " yr   z = ", z_got, &
                                                " m   error = ", err(k), " m"
        end do

        write(*,'(a,f9.4,a)') "   Nye analytic z = ", z_nye, " m"

        do k = 1, 2
            ratio = err(k)/err(k+1)
            write(*,'(a,f6.1,a,f6.1,a,f6.3)') "   error ratio dt=", dt(k), " / dt=", dt(k+1), " : ", ratio
            call check(abs(ratio-2.0_wp) .lt. 0.1_wp,"first-order convergence   ",n_fail)
        end do

    end subroutine test_nye_convergence

end program test_column
