module elsa
    ! Public facade of the Englacial Layer Simulation Architecture.
    !
    ! A host ice-sheet model should `use elsa` and nothing else: every type and
    ! procedure it needs is re-exported here. The internal modules (elsa_defs,
    ! elsa_interp, elsa_physics) are implementation detail.
    !
    !   call elsa_init(els,"elsa.nml","elsa",time,time_end,x,y,zeta,H_ice,stagger)
    !   call elsa_update(els,time,H_ice,ux,uy,smb,bmb)     ! every host timestep
    !   call elsa_end(els)
    !
    ! elsa_update takes an absolute time, works out its own dt, decides for
    ! itself whether an update is due, and keeps its own previous-step ice
    ! thickness. The host calls it unconditionally, once per step, and manages
    ! none of elsa's bookkeeping. Host fields may be single or double precision.
    !
    ! Method: Born (2017); Born and Robinson (2021); Rieckh et al. (2024),
    ! Geosci. Model Dev. 17, 6987-7000. See docs/DESIGN.md for where this
    ! implementation deliberately departs from the published v2.0 scheme.

    use elsa_precision, only : sp, dp, wp
    use elsa_defs
    use elsa_physics
    use elsa_interp
    use elsa_io
    use nml, only : nml_read

    implicit none

    private

    public :: sp, dp, wp
    public :: elsa_class, elsa_param_class, elsa_state_class, elsa_map_class
    public :: elsa_init, elsa_update, elsa_end
    public :: elsa_write_init, elsa_write_step
    public :: elsa_restart_write
    public :: elsa_version

    character(len=*), parameter :: elsa_version = "3.0.0-dev"

    ! Fractional slack when comparing model times, which arrive as accumulated
    ! sums and need not land exactly on a multiple of dt_coupling.
    real(wp), parameter :: TIME_TOL = 1.0e-6_wp

    interface elsa_init
        module procedure elsa_init_dp
        module procedure elsa_init_sp
    end interface

    interface elsa_update
        module procedure elsa_update_dp
        module procedure elsa_update_sp
    end interface

contains

    ! ======================================================================
    ! Public interface
    ! ======================================================================

    subroutine elsa_init_dp(els,filename,group,time,time_end,x,y,zeta,H_ice,stagger,restart)
        ! `x`, `y` are the host's cell-centre axes; `zeta` its sigma levels,
        ! ascending, 0 at the bed and 1 at the surface, with size(zeta) equal to
        ! size(ux,3). `stagger` says where the host's velocity samples sit:
        ! "acx_acy" (staggered faces, as Yelmo) or "aa" (cell-centred).
        !
        ! `time_end` sizes the layer stack, which is allocated once and never
        ! grown: a 100 kyr run at 200 yr resolution holds ~500 layers, and
        ! reallocating that is not something to do inside a time loop.
        !
        ! With `restart`, the layer stack and the isochrone schedule come from the
        ! file rather than from `time_end` and `layer_resolution`, and `time` is
        ! taken from the file too. The remaining namelist parameters (dt_coupling,
        ! cfl, grid_factor, allow_pos_bmb) are still read, so they may be changed
        ! across a restart -- but a changed grid_factor will fail the grid check.

        type(elsa_class),           intent(inout) :: els
        character(len=*),           intent(in)    :: filename, group
        real(dp),                   intent(in)    :: time, time_end
        real(dp),                   intent(in)    :: x(:), y(:), zeta(:)
        real(dp),                   intent(in)    :: H_ice(:,:)
        character(len=*),           intent(in)    :: stagger
        character(len=*), optional, intent(in)    :: restart

        integer :: k
        logical :: is_restart

        is_restart = present(restart)
        if (is_restart) then
            if (trim(restart) .eq. "None" .or. len_trim(restart) .eq. 0) is_restart = .false.
        end if

        call elsa_par_load(els%par,filename,group)

        if (is_restart) then
            call elsa_restart_read_par(els%par,restart)
        else
            call elsa_build_time_add(els%par,time,time_end)
            els%par%n_layers = els%par%n_layers_init + size(els%par%time_add) + 1
        end if

        call elsa_map_init(els%map,x,y,zeta,stagger,els%par%grid_factor)

        call elsa_alloc(els%now,els%map%nx,els%map%ny,els%map%nz,els%par%n_layers)

        call map_scalar(els%map,H_ice,els%now%H_ice)

        els%now%n_reseed = 0

        if (is_restart) then

            call elsa_restart_read_state(els,restart)

            if (abs(els%now%time - time) .gt. TIME_TOL*max(1.0_wp,abs(time))) then
                write(*,'(a)')        " elsa:: Warning: restart time differs from the host's."
                write(*,'(a,f14.2)')  "   restart : ", els%now%time
                write(*,'(a,f14.2)')  "   host    : ", time
                write(*,'(a)')        "   Using the restart time; the first update will span the gap."
            end if

        else

            els%now%H_ice_prev = els%now%H_ice

            ! Initialization layers: equal thickness, filling the column.
            els%now%n_top = els%par%n_layers_init
            do k = 1, els%now%n_top
                els%now%d_iso(:,:,k) = els%now%H_ice/real(els%now%n_top,wp)
            end do
            call normalize_layers(els%now%d_iso,els%now%H_ice,els%now%n_top)

            ! One empty layer on top, to receive accumulation until the first
            ! isochrone is laid down. It is laid down now, at time_init.
            els%now%n_top = els%now%n_top + 1
            els%now%t_dep(els%now%n_top) = time

            els%now%time           = time
            els%now%i_add          = 1
            els%now%n_reseed_total = 0

        end if

        call calc_dsum(els%now%dsum_iso,els%now%d_iso,els%now%n_top)

        call elsa_print_summary(els,is_restart)

    end subroutine elsa_init_dp

    subroutine elsa_update_dp(els,time,H_ice,ux,uy,smb,bmb)
        ! Advance elsa to `time`, if an update is due. All fields are on the
        ! host's grid, in the units elsa's namelist documents: H_ice [m],
        ! ux/uy [m/yr], smb/bmb [m/yr].

        type(elsa_class), intent(inout) :: els
        real(dp),         intent(in)    :: time
        real(dp),         intent(in)    :: H_ice(:,:)
        real(dp),         intent(in)    :: ux(:,:,:), uy(:,:,:)
        real(dp),         intent(in)    :: smb(:,:), bmb(:,:)

        real(wp) :: dt

        dt = time - els%now%time

        if (dt .le. 0.0_wp) return
        if (dt .lt. els%par%dt_coupling*(1.0_wp - TIME_TOL)) return

        ! -- host state onto elsa's grid
        call map_scalar(els%map,H_ice,els%now%H_ice)
        call map_scalar(els%map,smb,els%now%smb)
        call map_scalar(els%map,bmb,els%now%bmb)
        call map_velocity(els%map,ux,uy,els%now%ux_lev,els%now%uy_lev)

        ! -- tie the stack to the geometry the host's velocities were computed on
        call normalize_layers(els%now%d_iso,els%now%H_ice_prev,els%now%n_top)

        ! -- surface and basal mass exchange over the coupling period
        call apply_smb(els%now%d_iso,els%now%smb*dt,els%now%n_top)
        call apply_bmb(els%now%d_iso,els%now%bmb*dt,els%now%n_top,els%par%allow_pos_bmb)

        call calc_dsum(els%now%dsum_iso,els%now%d_iso,els%now%n_top)

        ! -- transport
        call elsa_interp_layer_velocities(els)
        call elsa_advect_layers(els,dt)

        ! -- restore any column the mass balance emptied that the host still has
        !    ice in, before the normalization tries to rescale nothing
        call reseed_empty_columns(els%now%d_iso,els%now%H_ice,els%now%n_top,els%now%n_reseed)

        if (els%now%n_reseed .gt. 0) then
            els%now%n_reseed_total = els%now%n_reseed_total + els%now%n_reseed
            if (els%now%n_reseed_total .eq. els%now%n_reseed) then
                write(*,'(a,i0,a,f10.1)') " elsa:: Note: reseeded ", els%now%n_reseed, &
                        " emptied column(s) at the bed, first at time ", time
                write(*,'(a)')            "   Expected in small numbers at thin margins where -smb*dt exceeds"
                write(*,'(a)')            "   the ice thickness. A large or growing count means dt_coupling is"
                write(*,'(a)')            "   too long, or the forcing is inconsistent. See now%n_reseed_total."
            end if
        end if

        ! -- drift correction: horizontal layer advection is not the host's mass
        !    conservation, so renormalize onto the host's ice thickness. The
        !    layers' relative thicknesses are unchanged.
        call normalize_layers(els%now%d_iso,els%now%H_ice,els%now%n_top)

        call elsa_add_due_layers(els,time)

        call calc_dsum(els%now%dsum_iso,els%now%d_iso,els%now%n_top)

        els%now%H_ice_prev = els%now%H_ice
        els%now%time       = time

    end subroutine elsa_update_dp

    subroutine elsa_end(els)
        type(elsa_class), intent(inout) :: els

        call elsa_dealloc(els%now)
        call elsa_map_end(els%map)
        if (allocated(els%par%time_add)) deallocate(els%par%time_add)

    end subroutine elsa_end

    ! ======================================================================
    ! Single-precision entry points. A host never casts.
    ! ======================================================================

    subroutine elsa_init_sp(els,filename,group,time,time_end,x,y,zeta,H_ice,stagger,restart)
        type(elsa_class),           intent(inout) :: els
        character(len=*),           intent(in)    :: filename, group
        real(sp),                   intent(in)    :: time, time_end
        real(sp),                   intent(in)    :: x(:), y(:), zeta(:)
        real(sp),                   intent(in)    :: H_ice(:,:)
        character(len=*),           intent(in)    :: stagger
        character(len=*), optional, intent(in)    :: restart

        call elsa_init_dp(els,filename,group,real(time,dp),real(time_end,dp), &
                          real(x,dp),real(y,dp),real(zeta,dp),real(H_ice,dp),stagger,restart)

    end subroutine elsa_init_sp

    subroutine elsa_update_sp(els,time,H_ice,ux,uy,smb,bmb)
        type(elsa_class), intent(inout) :: els
        real(sp),         intent(in)    :: time
        real(sp),         intent(in)    :: H_ice(:,:)
        real(sp),         intent(in)    :: ux(:,:,:), uy(:,:,:)
        real(sp),         intent(in)    :: smb(:,:), bmb(:,:)

        call elsa_update_dp(els,real(time,dp),real(H_ice,dp), &
                            real(ux,dp),real(uy,dp),real(smb,dp),real(bmb,dp))

    end subroutine elsa_update_sp

    ! ======================================================================
    ! Internals
    ! ======================================================================

    subroutine elsa_interp_layer_velocities(els)
        ! Layer-mean host velocity at each of elsa's faces.
        !
        ! The face column's geometry is the average of the two cell columns it
        ! separates, so dsum_face(n_top) is the face ice thickness by
        ! construction -- including at the margin, where one neighbour is
        ! ice-free and the face thickness is half the ice-covered cell's.
        !
        ! The outermost face in each direction carries zero flux and is never
        ! read; it is left at zero.

        type(elsa_class), intent(inout) :: els

        integer :: i, j, nx, ny, n_top
        real(wp), allocatable :: dsum_f(:)

        nx    = els%map%nx
        ny    = els%map%ny
        n_top = els%now%n_top

        !$omp parallel private(i,j,dsum_f)
        allocate(dsum_f(n_top))

        !$omp do
        do j = 1, ny
        do i = 1, nx-1
            dsum_f = 0.5_wp*(els%now%dsum_iso(i,j,1:n_top) + els%now%dsum_iso(i+1,j,1:n_top))
            call interp_u_column(els%now%ux_iso(i,j,1:n_top),els%now%ux_lev(i,j,:), &
                                 els%map%zeta,dsum_f(n_top),dsum_f)
        end do
        end do
        !$omp end do

        !$omp do
        do j = 1, ny-1
        do i = 1, nx
            dsum_f = 0.5_wp*(els%now%dsum_iso(i,j,1:n_top) + els%now%dsum_iso(i,j+1,1:n_top))
            call interp_u_column(els%now%uy_iso(i,j,1:n_top),els%now%uy_lev(i,j,:), &
                                 els%map%zeta,dsum_f(n_top),dsum_f)
        end do
        end do
        !$omp end do

        deallocate(dsum_f)
        !$omp end parallel

        els%now%ux_iso(nx,:,:) = 0.0_wp
        els%now%uy_iso(:,ny,:) = 0.0_wp

    end subroutine elsa_interp_layer_velocities

    subroutine elsa_advect_layers(els,dt)
        ! One thread per layer. The layers are independent -- they never exchange
        ! mass, which is the whole point of the isochronal scheme -- so each
        ! writes a disjoint slice and no synchronization is needed. v2.0 wrapped
        ! this loop in an OMP CRITICAL section, and called the LIS solver's
        ! global initialize/finalize inside every thread.
        !
        ! Dynamic scheduling: a layer's cost is its substep count, which is set by
        ! its own maximum outflow rate anywhere on the grid. On Greenland the
        ! layers turn out near-uniform (16 substeps at the bed against 17 at the
        ! surface, because the fast outlets slide, so every layer sees a similar
        ! grid maximum), and dynamic buys nothing. It is cheap insurance for a
        ! domain with a largely frozen bed, where the deep layers are genuinely
        ! slow and a static split would hand one thread all the fast ones.

        type(elsa_class), intent(inout) :: els
        real(wp),         intent(in)    :: dt

        integer :: iz, nx, ny
        real(wp), allocatable :: fx(:,:), fy(:,:)

        nx = els%map%nx
        ny = els%map%ny

        !$omp parallel private(iz,fx,fy)
        allocate(fx(nx,ny),fy(nx,ny))

        !$omp do schedule(dynamic)
        do iz = 1, els%now%n_top
            call advect_layer(els%now%d_iso(:,:,iz),els%now%ux_iso(:,:,iz), &
                              els%now%uy_iso(:,:,iz),els%map%dx,els%map%dy, &
                              dt,els%par%cfl,fx,fy)
        end do
        !$omp end do

        deallocate(fx,fy)
        !$omp end parallel

    end subroutine elsa_advect_layers

    subroutine elsa_add_due_layers(els,time)
        ! Lay down every isochrone whose time has been reached. More than one may
        ! fall inside a coupling period if layer_resolution < dt_coupling; v2.0
        ! rejected that configuration at init, this one just handles it (and
        ! warns at init, since the extra layers carry no accumulation).

        type(elsa_class), intent(inout) :: els
        real(wp),         intent(in)    :: time

        do while (els%now%i_add .le. size(els%par%time_add))

            if (els%par%time_add(els%now%i_add) .gt. time + TIME_TOL*abs(time) + TIME_TOL) exit

            if (els%now%n_top .ge. els%par%n_layers) then
                write(*,*) "elsa_add_due_layers:: Error: layer stack exhausted at time ", time
                write(*,*) "  n_layers = ", els%par%n_layers
                error stop 1
            end if

            els%now%n_top                    = els%now%n_top + 1
            els%now%d_iso(:,:,els%now%n_top) = 0.0_wp
            els%now%t_dep(els%now%n_top)     = els%par%time_add(els%now%i_add)
            els%now%i_add                    = els%now%i_add + 1

        end do

    end subroutine elsa_add_due_layers

    subroutine elsa_par_load(par,filename,group)
        type(elsa_param_class), intent(inout) :: par
        character(len=*),       intent(in)    :: filename, group

        call nml_read(filename,group,"n_layers_init",   par%n_layers_init)
        call nml_read(filename,group,"layer_resolution",par%layer_resolution)
        call nml_read(filename,group,"layer_file",      par%layer_file)
        call nml_read(filename,group,"grid_factor",     par%grid_factor)
        call nml_read(filename,group,"dt_coupling",     par%dt_coupling)
        call nml_read(filename,group,"cfl",             par%cfl)
        call nml_read(filename,group,"allow_pos_bmb",   par%allow_pos_bmb)

        if (par%n_layers_init .lt. 1) then
            write(*,*) "elsa_par_load:: Error: n_layers_init must be >= 1, got ", par%n_layers_init
            error stop 1
        end if

        if (par%dt_coupling .le. 0.0_wp) then
            write(*,*) "elsa_par_load:: Error: dt_coupling must be > 0, got ", par%dt_coupling
            error stop 1
        end if

        if (par%cfl .le. 0.0_wp .or. par%cfl .gt. 1.0_wp) then
            write(*,*) "elsa_par_load:: Error: cfl must be in (0,1], got ", par%cfl
            error stop 1
        end if

    end subroutine elsa_par_load

    subroutine elsa_build_time_add(par,time_init,time_end)
        ! The isochrone times, from either a regular resolution or an explicit list.
        type(elsa_param_class), intent(inout) :: par
        real(wp),               intent(in)    :: time_init, time_end

        integer  :: k, n
        real(wp) :: t, dt_min

        logical :: use_file, use_res

        use_file = trim(par%layer_file) .ne. "None"
        use_res  = par%layer_resolution .gt. 0.0_wp

        if (use_file .eqv. use_res) then
            write(*,*) "elsa_build_time_add:: Error: set exactly one of layer_file and layer_resolution."
            write(*,*) "  layer_file       = '"//trim(par%layer_file)//"'  (use 'None' to disable)"
            write(*,*) "  layer_resolution = ", par%layer_resolution, " (use 0 to disable)"
            error stop 1
        end if

        if (allocated(par%time_add)) deallocate(par%time_add)

        if (use_file) then

            call elsa_read_layer_file(par%time_add,par%layer_file)

        else

            n = max(0, ceiling((time_end-time_init)/par%layer_resolution) - 1)
            allocate(par%time_add(n))
            do k = 1, n
                par%time_add(k) = time_init + real(k,wp)*par%layer_resolution
            end do

        end if

        n = size(par%time_add)

        do k = 1, n
            t = par%time_add(k)
            if (t .le. time_init .or. t .ge. time_end) then
                write(*,*) "elsa_build_time_add:: Error: isochrone time outside the simulation window."
                write(*,*) "  time_add(",k,") = ", t
                write(*,*) "  time_init, time_end = ", time_init, time_end
                error stop 1
            end if
            if (k .gt. 1) then
                if (t .le. par%time_add(k-1)) then
                    write(*,*) "elsa_build_time_add:: Error: isochrone times must be strictly increasing."
                    write(*,*) "  breaks at k = ", k
                    error stop 1
                end if
            end if
        end do

        if (n .gt. 1) then
            dt_min = minval(par%time_add(2:n) - par%time_add(1:n-1))
            if (dt_min .lt. par%dt_coupling) then
                write(*,*) "elsa_build_time_add:: Warning: isochrones closer together (", dt_min, &
                           " yr) than the coupling period (", par%dt_coupling, " yr)."
                write(*,*) "  Layers laid down within one coupling period receive no accumulation."
            end if
        end if

    end subroutine elsa_build_time_add

    subroutine elsa_read_layer_file(time_add,filename)
        ! One isochrone time per line, ascending.
        real(wp), allocatable, intent(out) :: time_add(:)
        character(len=*),      intent(in)  :: filename

        integer :: unit, io, n, k
        logical :: exists
        character(len=256) :: line

        inquire(file=trim(filename),exist=exists)
        if (.not. exists) then
            write(*,*) "elsa_read_layer_file:: Error: layer_file not found: '"//trim(filename)//"'"
            error stop 1
        end if

        n = 0
        open(newunit=unit,file=trim(filename),status="old",action="read")
        do
            read(unit,'(a)',iostat=io) line
            if (io .ne. 0) exit
            if (len_trim(line) .eq. 0) cycle
            n = n + 1
        end do
        close(unit)

        allocate(time_add(n))

        k = 0
        open(newunit=unit,file=trim(filename),status="old",action="read")
        do
            read(unit,'(a)',iostat=io) line
            if (io .ne. 0) exit
            if (len_trim(line) .eq. 0) cycle
            k = k + 1
            read(line,*) time_add(k)
        end do
        close(unit)

    end subroutine elsa_read_layer_file

    subroutine elsa_print_summary(els,is_restart)
        type(elsa_class), intent(in) :: els
        logical,          intent(in) :: is_restart

        write(*,*) ""
        if (is_restart) then
            write(*,'(a,f0.2,a)') " elsa "//elsa_version//" (restarted at time ", els%now%time, ")"
        else
            write(*,'(a)')        " elsa "//elsa_version
        end if
        write(*,'(a,i0,a,i0)')  "   grid          : ", els%map%nx, " x ", els%map%ny
        write(*,'(a,f10.1,a)')  "   dx            : ", els%map%dx, " m"
        write(*,'(a,f10.1,a)')  "   grid_factor   : ", els%par%grid_factor, ""
        write(*,'(a,i0)')       "   host levels   : ", els%map%nz
        write(*,'(a,i0)')       "   layers total  : ", els%par%n_layers
        write(*,'(a,i0)')       "   layers init   : ", els%par%n_layers_init
        write(*,'(a,i0)')       "   isochrones    : ", size(els%par%time_add)
        write(*,'(a,f10.1,a)')  "   dt_coupling   : ", els%par%dt_coupling, " yr"
        write(*,'(a,f10.3)')    "   cfl           : ", els%par%cfl
        write(*,*) ""

    end subroutine elsa_print_summary

end module elsa
