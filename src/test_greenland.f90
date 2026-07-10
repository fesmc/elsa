program test_greenland
    ! The 3D Greenland benchmark: elsa forced offline from a 16 km Yelmo restart.
    !
    ! The restart carries no spun-up isochrone field, so there is no known answer
    ! to compare against. This is a structural check -- layers stay non-negative
    ! and ordered, columns sum to the host's ice thickness, no cell needs
    ! clipping -- run against a real velocity field with real staggering, real
    ! non-uniform sigma levels, and a real ice margin.
    !
    ! Yelmo writes ux/uy on staggered acx/acy nodes. elsa's advection wants face
    ! velocities, so they are handed over as-is with stagger = "acx_acy". The
    ! `tracer` package reads the same file and destaggers to aa nodes, because it
    ! wants cell-centred velocities; the contrast is the point of the convention.
    !
    ! The same driver serves an offline transient run: the forcing is read at an
    ! explicit index along the restart's time dimension.

    use elsa
    use elsa_physics, only : calc_n_substeps
    use ncio

    implicit none

    integer,  parameter :: NX = 106, NY = 181, NZ = 10
    integer,  parameter :: TIME_INDEX = 1
    real(wp), parameter :: TIME_0 = 0.0_wp
    real(wp), parameter :: TIME_1 = 2000.0_wp
    real(wp), parameter :: H_MIN  = 1.0e-6_wp

    character(len=*), parameter :: FILE_RESTART = "data/initmip-grl-16km/yelmo_restart.nc"
    character(len=*), parameter :: FILE_OUT     = "output/GRL-16KM/elsa.nc"

    real(wp) :: xc(NX), yc(NY), zeta(NZ)
    real(wp) :: H_ice(NX,NY), smb(NX,NY), bmb(NX,NY)
    real(wp) :: ux(NX,NY,NZ), uy(NX,NY,NZ)

    integer :: n_fail

    n_fail = 0

    write(*,*) ""
    write(*,*) "== elsa Greenland 16 km =="

    call load_forcing()

    call run_case("greenland",       .true., n_fail)
    call run_case("greenland_coarse",.false.,n_fail)

    call check_output(n_fail)

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

    subroutine load_forcing()

        call nc_read(FILE_RESTART,"xc",  xc)
        call nc_read(FILE_RESTART,"yc",  yc)
        call nc_read(FILE_RESTART,"zeta",zeta)

        ! The restart's axes are in km; elsa works in metres.
        xc = xc*1.0e3_wp
        yc = yc*1.0e3_wp

        call nc_read(FILE_RESTART,"H_ice",H_ice,start=[1,1,TIME_INDEX],  count=[NX,NY,1])
        call nc_read(FILE_RESTART,"smb",  smb,  start=[1,1,TIME_INDEX],  count=[NX,NY,1])
        call nc_read(FILE_RESTART,"ux",   ux,   start=[1,1,1,TIME_INDEX],count=[NX,NY,NZ,1])
        call nc_read(FILE_RESTART,"uy",   uy,   start=[1,1,1,TIME_INDEX],count=[NX,NY,NZ,1])

        ! The restart carries no basal mass balance.
        bmb = 0.0_wp

        write(*,*) ""
        write(*,'(a,i0,a,i0,a,i0)') "   host grid : ", NX, " x ", NY, " x ", NZ
        write(*,'(a,f8.1,a)')       "   host dx   : ", xc(2)-xc(1), " m"
        write(*,'(a,f10.1,a)')      "   max |ux|  : ", maxval(abs(ux)), " m/yr"
        write(*,'(a,f10.1,a)')      "   max |uy|  : ", maxval(abs(uy)), " m/yr"
        write(*,'(a,f10.1,a)')      "   max H_ice : ", maxval(H_ice), " m"
        write(*,'(a,f10.3,a)')      "   smb range : ", minval(smb), " ..."
        write(*,'(a,f10.3,a)')      "               ", maxval(smb), " m/yr"

    end subroutine load_forcing

    subroutine run_case(group,write_output,n_fail)
        character(len=*), intent(in)    :: group
        logical,          intent(in)    :: write_output
        integer,          intent(inout) :: n_fail

        type(elsa_class) :: els
        real(wp) :: time, vol_src, vol_els, err, e_ij
        integer  :: n, n_steps, n_out, i, j, k, nx_e, ny_e, n_top, n_bad, i_bad, j_bad
        integer  :: t0, t1, rate, t_sum, n_sub, n_sub_tot, n_sub_max, n_sub_bed, n_sub_srf
        logical  :: ok

        write(*,*) ""
        write(*,*) " case: "//group

        call elsa_init(els,"par/test_greenland.nml",group,TIME_0,TIME_1,xc,yc,zeta,H_ice,"acx_acy")

        nx_e = els%map%nx
        ny_e = els%map%ny

        if (write_output) then
            call elsa_write_init(els,FILE_OUT,TIME_0)
            call elsa_write_step(els,FILE_OUT,TIME_0,1)
            n_out = 1
        end if

        n_steps = nint((TIME_1-TIME_0)/els%par%dt_coupling)

        ! Time elsa_update only. The netCDF writes are serial and would otherwise
        ! be charged to the compute, hiding how the solver actually scales.
        call system_clock(count_rate=rate)
        t_sum = 0

        do n = 1, n_steps
            time = TIME_0 + real(n,wp)*els%par%dt_coupling

            call system_clock(t0)
            call elsa_update(els,time,H_ice,ux,uy,smb,bmb)
            call system_clock(t1)
            t_sum = t_sum + (t1-t0)

            if (write_output .and. mod(n,10) .eq. 0) then
                n_out = n_out + 1
                call elsa_write_step(els,FILE_OUT,time,n_out)
            end if
        end do

        write(*,'(a,f8.2,a)') "   wall time (elsa_update) : ", real(t_sum,wp)/real(rate,wp), " s"

        n_top = els%now%n_top

        ! Per-layer substep counts. A layer's cost is proportional to its own
        ! substep count, and no layer can be split across threads, so the
        ! parallel floor is max(n_sub)/sum(n_sub) of the serial time.
        n_sub_tot = 0
        n_sub_max = 0
        do k = 1, n_top
            n_sub = calc_n_substeps(els%now%ux_iso(:,:,k),els%now%uy_iso(:,:,k), &
                                    els%map%dx,els%map%dy,els%par%dt_coupling,els%par%cfl)
            n_sub_tot = n_sub_tot + n_sub
            n_sub_max = max(n_sub_max,n_sub)
            if (k .eq. 1)     n_sub_bed = n_sub
            if (k .eq. n_top) n_sub_srf = n_sub
        end do
        write(*,'(a,i0,a,i0,a,i0)') "   substeps/step : bed layer ", n_sub_bed, &
                    ", surface layer ", n_sub_srf, ", total ", n_sub_tot
        write(*,'(a,f6.2,a)')       "   layer-parallel speedup ceiling : ", &
                    real(n_sub_tot,wp)/real(n_sub_max,wp), " x"

        ! -- finiteness and positivity ---------------------------------------
        call check(all(els%now%d_iso .eq. els%now%d_iso),   "all finite (no NaN)       ",n_fail)
        call check(minval(els%now%d_iso) .ge. 0.0_wp,       "layers non-negative       ",n_fail)
        call check(n_top .eq. els%par%n_layers,             "all isochrones laid down  ",n_fail)

        ! -- every column sums to the host's ice thickness ---------------------
        err = 0.0_wp
        n_bad = 0
        do j = 1, ny_e
        do i = 1, nx_e
            if (els%now%H_ice(i,j) .gt. H_MIN) then
                e_ij = abs(els%now%dsum_iso(i,j,n_top) - els%now%H_ice(i,j))/els%now%H_ice(i,j)
                if (e_ij .gt. 1.0e-9_wp) then
                    n_bad = n_bad + 1
                    if (e_ij .gt. err) then
                        i_bad = i
                        j_bad = j
                    end if
                end if
                err = max(err,e_ij)
            end if
        end do
        end do
        write(*,'(a,es9.2,a,i0,a)') "   max column closure error : ", err, "  (", n_bad, " cells)"
        if (n_bad .gt. 0) then
            write(*,'(a,i0,a,i0,a)')  "     worst cell (", i_bad, ",", j_bad, "):"
            write(*,'(a,f10.3,a)')    "       H_ice (elsa grid) = ", els%now%H_ice(i_bad,j_bad), " m"
            write(*,'(a,f10.3,a)')    "       smb              = ", els%now%smb(i_bad,j_bad), " m/yr"
            write(*,'(a,f10.3,a)')    "       smb * dt         = ", els%now%smb(i_bad,j_bad)*els%par%dt_coupling, " m"
            write(*,'(a,f10.3,a)')    "       elsa column      = ", els%now%dsum_iso(i_bad,j_bad,n_top), " m"
        end if
        call check(err .lt. 1.0e-12_wp,                     "columns sum to H_ice      ",n_fail)

        ! -- reseeding: exercised here, and rare -------------------------------
        !    85 of the 7204 ice cells have -smb*dt > H at dt_coupling = 50 yr, so
        !    the mass balance can annihilate their column in one step. Advection
        !    refills most of them; the rest are reseeded at the bed.
        !
        !    The forcing here is steady, so the handful that need reseeding need
        !    it again every step: the column is normalized back to H = 8 m, then
        !    ablated by 141 m. The total is therefore n_reseed * n_steps, and it
        !    is the *per-step* count that says whether this is rare.
        write(*,'(a,i0,a,i0,a)') "   columns reseeded at the bed : ", els%now%n_reseed, &
                                 " per step (", els%now%n_reseed_total, " events total)"
        call check(els%now%n_reseed_total .gt. 0,            "reseed path exercised     ",n_fail)
        call check(els%now%n_reseed .lt. nx_e*ny_e/100,      "reseeding stays rare      ",n_fail)

        ! -- ice-free columns carry no layers ---------------------------------
        ok = .true.
        do j = 1, ny_e
        do i = 1, nx_e
            if (els%now%H_ice(i,j) .le. H_MIN) then
                if (any(els%now%d_iso(i,j,1:n_top) .ne. 0.0_wp)) ok = .false.
            end if
        end do
        end do
        call check(ok,                                      "ice-free columns empty    ",n_fail)

        ! -- isochrones ordered: younger ice sits higher ----------------------
        ok = .true.
        do k = 2, n_top
            if (any(els%now%dsum_iso(:,:,k) .lt. els%now%dsum_iso(:,:,k-1) - 1.0e-9_wp)) ok = .false.
        end do
        call check(ok,                                      "isochrones monotone in z  ",n_fail)

        ! -- no isochrone escapes the ice column ------------------------------
        call check(maxval(els%now%dsum_iso(:,:,1:n_top)) .le. maxval(els%now%H_ice)*(1.0_wp+1.0e-12_wp), &
                                                            "isochrones inside the ice ",n_fail)

        ! -- the conservative map preserves ice volume, exactly, at gf = 1 -----
        !    Cell areas are identical here, so the sums are compared directly.
        !    (Deriving the source area from xc(2)-xc(1) would not do: the axis is
        !    single precision, and that spacing differs from the endpoint-derived
        !    dx by 1e-7 relative -- see axis_spacing in elsa_interp.)
        if (abs(els%par%grid_factor - 1.0_wp) .lt. 1.0e-12_wp) then
            vol_src = sum(H_ice)
            vol_els = sum(els%now%H_ice)
            call check(abs(vol_els-vol_src) .le. 1.0e-13_wp*vol_src,"grid_factor=1 volume exact",n_fail)
            call check(maxval(abs(els%now%H_ice-H_ice)) .eq. 0.0_wp,"grid_factor=1 H identity  ",n_fail)
        else
            write(*,'(a,i0,a,i0,a,f6.1,a)') "   coarsened to ", nx_e, " x ", ny_e, &
                                            " (dx = ", els%map%dx*1.0e-3_wp, " km)"
        end if

        write(*,'(a,f9.2,a)') "   deepest isochrone below surface : ", &
            maxval(els%now%H_ice) - minval(els%now%dsum_iso(:,:,els%par%n_layers_init), &
                                           mask=els%now%H_ice .gt. H_MIN), " m"

        call elsa_end(els)

    end subroutine run_case

    subroutine check_output(n_fail)
        ! Read the file back: elsa_io must round-trip.
        integer, intent(inout) :: n_fail

        real(sp), allocatable :: d_iso(:,:,:,:)
        real(wp), allocatable :: layer_time(:)
        integer :: nl

        write(*,*) ""
        write(*,*) " output round-trip"

        nl = nc_size(FILE_OUT,"layer")
        allocate(layer_time(nl))
        call nc_read(FILE_OUT,"layer_time",layer_time)

        call check(nl .eq. 20,                              "layer dimension written   ",n_fail)
        call check(nc_size(FILE_OUT,"time") .eq. 5,         "5 time slices written     ",n_fail)

        ! layer_time is missing for the initialization layers and increasing after.
        call check(all(layer_time(11:nl) - layer_time(10:nl-1) .gt. 0.0_wp), &
                                                            "layer_time increasing     ",n_fail)

        allocate(d_iso(NX,NY,nl,1))
        call nc_read(FILE_OUT,"d_iso",d_iso,start=[1,1,1,5],count=[NX,NY,nl,1])
        call check(minval(d_iso) .ge. 0.0_sp,               "d_iso round-trips         ",n_fail)

        write(*,'(a,a)') "   wrote ", FILE_OUT

    end subroutine check_output

end program test_greenland
