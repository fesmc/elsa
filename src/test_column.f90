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

    character(len=*), parameter :: FILE_OUT = "output/column/elsa.nc"

    integer :: n_fail

    n_fail = 0

    write(*,*) ""
    write(*,*) "== elsa column (Nye) =="

    call test_discrete_exact(n_fail)
    call test_nye_convergence(n_fail)
    call test_time_mean(n_fail)
    call test_single_precision(n_fail)
    call test_layer_file(n_fail)
    call test_stagger(n_fail)

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

    subroutine set_divide(x,y,zeta,H_ice,smb,bmb,ux,uy)
        ! The host fields of a steady, horizontally uniform ice divide.
        real(wp), intent(out) :: x(NX), y(NY), zeta(NZ)
        real(wp), intent(out) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp), intent(out) :: ux(NX,NY,NZ), uy(NX,NY,NZ)

        integer :: i

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

    end subroutine set_divide

    subroutine run_divide(els,group,dt,write_output)
        ! Drive elsa with a steady, horizontally uniform ice divide.
        type(elsa_class),  intent(inout) :: els
        character(len=*),  intent(in)    :: group
        real(wp),          intent(in)    :: dt
        logical, optional, intent(in)    :: write_output

        real(wp) :: x(NX), y(NY), zeta(NZ)
        real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)
        real(wp) :: time
        integer  :: n, n_steps, n_out
        logical  :: do_write

        call set_divide(x,y,zeta,H_ice,smb,bmb,ux,uy)

        do_write = .false.
        if (present(write_output)) do_write = write_output

        call elsa_init(els,"par/test_column.nml",group,TIME_0,TIME_1,x,y,zeta,H_ice,"aa")

        if (do_write) then
            call elsa_write_init(els,FILE_OUT,TIME_0)
            call elsa_write_step(els,FILE_OUT,TIME_0,1)
            n_out = 1
        end if

        n_steps = nint((TIME_1-TIME_0)/dt)
        do n = 1, n_steps
            time = TIME_0 + real(n,wp)*dt
            call elsa_update(els,time,H_ice,ux,uy,smb,bmb)

            if (do_write .and. mod(n,20) .eq. 0) then
                n_out = n_out + 1
                call elsa_write_step(els,FILE_OUT,time,n_out)
            end if
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

        call run_divide(els,"column",DT,write_output=.true.)

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

    subroutine test_time_mean(n_fail)
        ! elsa applies the host's mean forcing over the coupling period, not its
        ! value at the final instant. An accumulation that alternates about its
        ! mean from one host step to the next must therefore give the layers of
        ! the constant one -- the instantaneous value would be off by a third
        ! here. The run is also stopped between two updates: the restart must
        ! carry the part of the period already integrated.
        integer, intent(inout) :: n_fail

        real(wp), parameter :: DT       = 100.0_wp
        integer,  parameter :: N_HOST   = 4                 ! host steps per coupling period
        real(wp), parameter :: DT_HOST  = DT/real(N_HOST,wp)
        real(wp), parameter :: TIME_MID = 10050.0_wp        ! between two updates

        character(len=*), parameter :: FILE_RST = "output/column/elsa_restart.nc"

        type(elsa_class) :: els
        real(wp) :: x(NX), y(NY), zeta(NZ)
        real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)
        real(wp), allocatable :: d_const(:,:,:), d_ref(:,:,:)
        real(wp) :: err
        integer  :: n, n_steps, n_mid

        write(*,*) ""
        write(*,*) " time-mean forcing"

        call set_divide(x,y,zeta,H_ice,smb,bmb,ux,uy)

        n_steps = N_HOST*nint((TIME_1-TIME_0)/DT)
        n_mid   = nint((TIME_MID-TIME_0)/DT_HOST)

        ! -- constant accumulation ------------------------------------------------
        call elsa_init(els,"par/test_column.nml","column",TIME_0,TIME_1,x,y,zeta,H_ice,"aa")
        do n = 1, n_steps
            call elsa_update(els,TIME_0+real(n,wp)*DT_HOST,H_ice,ux,uy,smb,bmb)
        end do
        allocate(d_const(size(els%now%d_iso,1),size(els%now%d_iso,2),size(els%now%d_iso,3)))
        allocate(d_ref,mold=d_const)
        d_const = els%now%d_iso
        call elsa_end(els)

        ! -- alternating about the same mean, in one run ------------------------
        call elsa_init(els,"par/test_column.nml","column",TIME_0,TIME_1,x,y,zeta,H_ice,"aa")
        do n = 1, n_steps
            call elsa_update(els,TIME_0+real(n,wp)*DT_HOST,H_ice,ux,uy,smb_alt(n),bmb)
        end do
        d_ref = els%now%d_iso

        err = maxval(abs(d_ref - d_const))/H_CONST
        write(*,'(a,es9.2)') "   max relative difference to constant forcing: ", err
        call check(err .lt. 1.0e-12_wp,               "mean forcing is applied   ",n_fail)
        call check(maxval(abs(els%now%smb - ACC)) .lt. 1.0e-12_wp, &
                                                      "smb diagnostic is the mean",n_fail)
        call elsa_end(els)

        ! -- the same, stopped between two updates --------------------------------
        call elsa_init(els,"par/test_column.nml","column",TIME_0,TIME_1,x,y,zeta,H_ice,"aa")
        do n = 1, n_mid
            call elsa_update(els,TIME_0+real(n,wp)*DT_HOST,H_ice,ux,uy,smb_alt(n),bmb)
        end do
        call check(els%now%time .lt. TIME_MID,        "stopped between updates   ",n_fail)
        call elsa_restart_write(els,FILE_RST)
        call elsa_end(els)

        call elsa_init(els,"par/test_column.nml","column",TIME_MID,TIME_1,x,y,zeta,H_ice,"aa", &
                       restart=FILE_RST)
        do n = n_mid+1, n_steps
            call elsa_update(els,TIME_0+real(n,wp)*DT_HOST,H_ice,ux,uy,smb_alt(n),bmb)
        end do

        call check(maxval(abs(els%now%d_iso - d_ref)) .eq. 0.0_wp, &
                                                      "mid-period restart exact  ",n_fail)

        call elsa_end(els)
        deallocate(d_const,d_ref)

    end subroutine test_time_mean

    function smb_alt(n) result(smb)
        ! Accumulation at host step n: ACC +/- a third, mean ACC over any even
        ! number of steps.
        integer, intent(in) :: n
        real(wp) :: smb(NX,NY)

        if (mod(n,2) .eq. 0) then
            smb = ACC*(1.0_wp + 1.0_wp/3.0_wp)
        else
            smb = ACC*(1.0_wp - 1.0_wp/3.0_wp)
        end if

    end function smb_alt

    subroutine test_single_precision(n_fail)
        ! A single-precision host must get the double-precision answer for the
        ! same (single-precision) numbers. The host is called four times per
        ! coupling period, so three calls in four are not due and must not
        ! advance elsa.
        integer, intent(inout) :: n_fail

        real(wp), parameter :: DT = 100.0_wp

        type(elsa_class) :: els_dp, els_sp
        real(wp) :: x(NX), y(NY), zeta(NZ)
        real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)
        real(wp) :: time, time_prev
        integer  :: n, n_steps
        logical  :: gated

        write(*,*) ""
        write(*,*) " single-precision interface"

        call set_divide(x,y,zeta,H_ice,smb,bmb,ux,uy)

        ! The reference sees the accumulation the single-precision host has.
        smb = real(real(smb,sp),wp)

        n_steps = nint((TIME_1-TIME_0)/DT)

        call elsa_init(els_dp,"par/test_column.nml","column",TIME_0,TIME_1,x,y,zeta,H_ice,"aa")
        do n = 1, 4*n_steps
            call elsa_update(els_dp,TIME_0+real(n,wp)*0.25_wp*DT,H_ice,ux,uy,smb,bmb)
        end do

        call elsa_init(els_sp,"par/test_column.nml","column",real(TIME_0,sp),real(TIME_1,sp), &
                       real(x,sp),real(y,sp),real(zeta,sp),real(H_ice,sp),"aa")

        gated = .true.
        do n = 1, 4*n_steps
            time      = TIME_0 + real(n,wp)*0.25_wp*DT
            time_prev = els_sp%now%time
            call elsa_update(els_sp,real(time,sp),real(H_ice,sp), &
                             real(ux,sp),real(uy,sp),real(smb,sp),real(bmb,sp))
            if (mod(n,4) .ne. 0 .and. els_sp%now%time .ne. time_prev) gated = .false.
        end do

        call check(gated,                             "calls not due do not step ",n_fail)
        call check(els_sp%now%n_top .eq. els_dp%now%n_top,"same n_top as double      ",n_fail)
        call check(maxval(abs(els_sp%now%d_iso - els_dp%now%d_iso)) .eq. 0.0_wp, &
                                                      "bit-identical to double   ",n_fail)

        call elsa_end(els_dp)
        call elsa_end(els_sp)

    end subroutine test_single_precision

    subroutine test_layer_file(n_fail)
        ! An explicit, irregular isochrone list: the schedule is the file's, the
        ! isochrones sit where the discrete recursion puts them, and a restart
        ! whose first segment ended mid-list picks up the rest of the list.
        !
        ! The second entry falls between updates. Its layer is laid down at the
        ! next update, and must be stamped with that time, not the scheduled one.
        integer, intent(inout) :: n_fail

        real(wp), parameter :: DT       = 100.0_wp
        real(wp), parameter :: TIME_MID = 5000.0_wp
        real(wp), parameter :: T_ISO(4)  = [1000.0_wp,2550.0_wp,7000.0_wp,15000.0_wp]
        real(wp), parameter :: T_LAID(4) = [1000.0_wp,2600.0_wp,7000.0_wp,15000.0_wp]

        character(len=*), parameter :: FILE_LAYERS = "output/column/layers.txt"
        character(len=*), parameter :: FILE_RST    = "output/column/elsa_restart.nc"

        type(elsa_class) :: els
        real(wp) :: x(NX), y(NY), zeta(NZ)
        real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)
        real(wp), allocatable :: d_ref(:,:,:)
        real(wp) :: r, err
        integer  :: n, jj, unit

        write(*,*) ""
        write(*,*) " isochrones from a layer file"

        open(newunit=unit,file=FILE_LAYERS,status="replace",action="write")
        do jj = 1, size(T_ISO)
            write(unit,*) T_ISO(jj)
        end do
        close(unit)

        call set_divide(x,y,zeta,H_ice,smb,bmb,ux,uy)

        call elsa_init(els,"par/test_column.nml","column_file",TIME_0,TIME_1,x,y,zeta,H_ice,"aa")

        call check(size(els%par%time_add) .eq. size(T_ISO),"schedule has file's length",n_fail)
        if (size(els%par%time_add) .eq. size(T_ISO)) then
            call check(all(els%par%time_add .eq. T_ISO),   "schedule has file's times ",n_fail)
        end if

        do n = 1, nint((TIME_1-TIME_0)/DT)
            call elsa_update(els,TIME_0+real(n,wp)*DT,H_ice,ux,uy,smb,bmb)
        end do

        r   = H_CONST/(H_CONST + ACC*DT)
        err = 0.0_wp
        do jj = 1, size(T_ISO)
            err = max(err,abs(els%now%dsum_iso(3,3,N_INIT+jj) &
                              /(H_CONST*r**nint((TIME_1-T_LAID(jj))/DT)) - 1.0_wp))
        end do
        call check(err .lt. 1.0e-11_wp,               "isochrone heights exact   ",n_fail)
        call check(all(els%now%t_dep(N_INIT+2:N_INIT+1+size(T_ISO)) .eq. T_LAID), &
                                                      "t_dep is the update time  ",n_fail)

        allocate(d_ref(size(els%now%d_iso,1),size(els%now%d_iso,2),size(els%now%d_iso,3)))
        d_ref = els%now%d_iso
        call elsa_end(els)

        ! -- first segment knows only its own end, mid-list -----------------------
        call elsa_init(els,"par/test_column.nml","column_file",TIME_0,TIME_MID,x,y,zeta,H_ice,"aa")
        call check(size(els%par%time_add) .eq. 2,     "segment skips later times ",n_fail)
        do n = 1, nint((TIME_MID-TIME_0)/DT)
            call elsa_update(els,TIME_0+real(n,wp)*DT,H_ice,ux,uy,smb,bmb)
        end do
        call elsa_restart_write(els,FILE_RST)
        call elsa_end(els)

        call elsa_init(els,"par/test_column.nml","column_file",TIME_MID,TIME_1,x,y,zeta,H_ice,"aa", &
                       restart=FILE_RST)
        call check(size(els%par%time_add) .eq. size(T_ISO),"restart extends the list  ",n_fail)
        do n = 1, nint((TIME_1-TIME_MID)/DT)
            call elsa_update(els,TIME_MID+real(n,wp)*DT,H_ice,ux,uy,smb,bmb)
        end do

        if (size(els%now%d_iso,3) .eq. size(d_ref,3)) then
            call check(maxval(abs(els%now%d_iso - d_ref)) .eq. 0.0_wp, &
                                                      "extended is bit-identical ",n_fail)
        else
            call check(.false.,                       "extended is bit-identical ",n_fail)
        end if
        call check(all(els%now%t_dep(N_INIT+2:N_INIT+1+size(T_ISO)) .eq. T_LAID), &
                                                      "extended keeps t_dep      ",n_fail)

        call elsa_end(els)
        deallocate(d_ref)

    end subroutine test_layer_file

    subroutine test_stagger(n_fail)
        ! End to end with flow: a uniform velocity is the same field whether the
        ! host declares it on staggered faces or at cell centres, so the two
        ! declarations must give the same layers.
        integer, intent(inout) :: n_fail

        real(wp), parameter :: DT = 100.0_wp
        real(wp), parameter :: U0 = 10.0_wp         ! [m/yr]

        type(elsa_class) :: els_ac, els_aa
        real(wp) :: x(NX), y(NY), zeta(NZ)
        real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
        real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)
        real(wp) :: time, err
        integer  :: n

        write(*,*) ""
        write(*,*) " staggering, with flow"

        call set_divide(x,y,zeta,H_ice,smb,bmb,ux,uy)
        ux =  U0
        uy = -0.5_wp*U0

        call elsa_init(els_ac,"par/test_column.nml","column",TIME_0,TIME_1,x,y,zeta,H_ice,"acx_acy")
        call elsa_init(els_aa,"par/test_column.nml","column",TIME_0,TIME_1,x,y,zeta,H_ice,"aa")

        do n = 1, nint((TIME_1-TIME_0)/DT)
            time = TIME_0 + real(n,wp)*DT
            call elsa_update(els_ac,time,H_ice,ux,uy,smb,bmb)
            call elsa_update(els_aa,time,H_ice,ux,uy,smb,bmb)
        end do

        err = maxval(abs(els_aa%now%d_iso - els_ac%now%d_iso))/H_CONST
        write(*,'(a,es9.2)') "   max relative difference aa vs acx_acy: ", err

        call check(maxval(abs(els_ac%now%ux_iso(1:NX-1,:,1:els_ac%now%n_top) - U0)) .lt. 1.0e-9_wp, &
                                                      "flow reaches the layers   ",n_fail)
        call check(err .lt. 1.0e-10_wp,               "aa and acx_acy agree      ",n_fail)
        call check(abs(els_aa%now%dsum_iso(3,3,els_aa%now%n_top) - H_CONST) .lt. 1.0e-9_wp, &
                                                      "column sums to H          ",n_fail)

        call elsa_end(els_ac)
        call elsa_end(els_aa)

    end subroutine test_stagger

end program test_column
