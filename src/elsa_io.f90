module elsa_io
    ! NetCDF output for elsa, via fesm-utils' ncio.
    !
    ! The file carries `layer_time`: the model time at which each layer was laid
    ! down. That is the only thing that makes the output interpretable, because
    ! elsa's vertical axis is time -- layer k is bounded below by the isochrone
    ! of age `layer_time(k)`. The initialization layers carry a missing value,
    ! since they are not isochrones of anything.

    use elsa_precision, only : sp, wp, MV
    use elsa_defs
    use ncio

    implicit none

    private

    public :: elsa_write_init
    public :: elsa_write_step
    public :: elsa_restart_write
    public :: elsa_restart_read_par
    public :: elsa_restart_read_state

contains

    ! ======================================================================
    ! Restart
    ! ======================================================================
    !
    ! The file carries the full state object, so it doubles as a single-slice
    ! diagnostic snapshot. But only a subset is read back to continue a run:
    ! d_iso, H_ice_prev, and the scalar bookkeeping. dsum_iso is derived from
    ! d_iso, and the velocities and mass balance are remapped from the host every
    ! update, so restoring them would have no effect -- see elsa_restart_read_state.
    !
    ! d_iso and H_ice_prev are written in double precision: they are read back, and
    ! a restart that loses bits does not reproduce the run it continues. The other
    ! fields are diagnostic only and never read back, so they are written single
    ! precision like the output in elsa_write_step.
    !
    ! The isochrone schedule (time_add) travels in the restart rather than being
    ! regenerated from layer_resolution. Rebuilding it from the restart time would
    ! shift every subsequent isochrone.

    subroutine elsa_restart_write(els,filename)
        type(elsa_class), intent(in) :: els
        character(len=*), intent(in) :: filename

        integer :: n_add

        n_add = size(els%par%time_add)

        call nc_create(filename)

        call nc_write_dim(filename,"xc",   x=els%map%x,units="m")
        call nc_write_dim(filename,"yc",   x=els%map%y,units="m")
        call nc_write_dim(filename,"zeta", x=els%map%zeta,units="1")
        call nc_write_dim(filename,"layer",x=1.0_wp,dx=1.0_wp,nx=els%par%n_layers,units="1")
        call nc_write_dim(filename,"one",  x=1.0_wp,dx=1.0_wp,nx=1,units="1")

        call nc_write(filename,"time",         els%now%time,          dim1="one")
        call nc_write(filename,"n_top",        els%now%n_top,         dim1="one")
        call nc_write(filename,"i_add",        els%now%i_add,         dim1="one")
        call nc_write(filename,"n_reseed",     els%now%n_reseed,      dim1="one")
        call nc_write(filename,"n_reseed_total",els%now%n_reseed_total,dim1="one")
        call nc_write(filename,"n_layers_init",els%par%n_layers_init, dim1="one")

        if (n_add .gt. 0) then
            call nc_write_dim(filename,"isochrone",x=1.0_wp,dx=1.0_wp,nx=n_add,units="1")
            call nc_write(filename,"time_add",els%par%time_add,dim1="isochrone",units="years")
        end if

        ! -- read back to continue the run: double precision
        call nc_write(filename,"d_iso",els%now%d_iso, &
                      dim1="xc",dim2="yc",dim3="layer",units="m")
        call nc_write(filename,"H_ice_prev",els%now%H_ice_prev, &
                      dim1="xc",dim2="yc",units="m")

        ! -- diagnostic only, never read back: single precision
        call nc_write(filename,"dsum_iso",real(els%now%dsum_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",units="m", &
                      long_name="Height of the layer top above the bed")
        call nc_write(filename,"ux_iso",real(els%now%ux_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",units="m/yr", &
                      long_name="Layer-mean velocity (x), acx nodes")
        call nc_write(filename,"uy_iso",real(els%now%uy_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",units="m/yr", &
                      long_name="Layer-mean velocity (y), acy nodes")
        call nc_write(filename,"H_ice",real(els%now%H_ice,sp), &
                      dim1="xc",dim2="yc",units="m", &
                      long_name="Ice thickness on the elsa grid")
        call nc_write(filename,"smb",real(els%now%smb,sp), &
                      dim1="xc",dim2="yc",units="m/yr")
        call nc_write(filename,"bmb",real(els%now%bmb,sp), &
                      dim1="xc",dim2="yc",units="m/yr")
        call nc_write(filename,"ux_lev",real(els%now%ux_lev,sp), &
                      dim1="xc",dim2="yc",dim3="zeta",units="m/yr", &
                      long_name="Host velocity (x), acx nodes, host sigma levels")
        call nc_write(filename,"uy_lev",real(els%now%uy_lev,sp), &
                      dim1="xc",dim2="yc",dim3="zeta",units="m/yr", &
                      long_name="Host velocity (y), acy nodes, host sigma levels")

    end subroutine elsa_restart_write

    subroutine elsa_restart_read_par(par,filename)
        ! The layer structure, which must be known before the state is allocated.
        type(elsa_param_class), intent(inout) :: par
        character(len=*),       intent(in)    :: filename

        integer :: n_add
        logical :: exists

        inquire(file=trim(filename),exist=exists)
        if (.not. exists) then
            write(*,*) "elsa_restart_read:: Error: restart file not found: '"//trim(filename)//"'"
            error stop 1
        end if

        call nc_read(filename,"n_layers_init",par%n_layers_init,start=[1],count=[1])

        par%n_layers = nc_size(filename,"layer")

        if (allocated(par%time_add)) deallocate(par%time_add)
        n_add = par%n_layers - par%n_layers_init - 1
        allocate(par%time_add(n_add))
        if (n_add .gt. 0) call nc_read(filename,"time_add",par%time_add)

    end subroutine elsa_restart_read_par

    subroutine elsa_restart_read_state(els,filename)
        ! The state, after the grid has been built and the arrays allocated. The
        ! grid the restart was written on must be the grid we are about to use:
        ! a changed grid_factor or a changed host domain is a hard error, not
        ! something to interpolate away silently.
        type(elsa_class), intent(inout) :: els
        character(len=*), intent(in)    :: filename

        real(wp), allocatable :: x_chk(:), y_chk(:), zeta_chk(:)

        if (nc_size(filename,"xc") .ne. els%map%nx .or. &
            nc_size(filename,"yc") .ne. els%map%ny) then
            write(*,*) "elsa_restart_read:: Error: restart grid does not match the current grid."
            write(*,*) "  restart : ", nc_size(filename,"xc"), " x ", nc_size(filename,"yc")
            write(*,*) "  current : ", els%map%nx, " x ", els%map%ny
            error stop 1
        end if

        if (nc_size(filename,"zeta") .ne. els%map%nz) then
            write(*,*) "elsa_restart_read:: Error: restart has ", nc_size(filename,"zeta"), &
                       " host levels, current grid has ", els%map%nz
            error stop 1
        end if

        allocate(x_chk(els%map%nx),y_chk(els%map%ny),zeta_chk(els%map%nz))
        call nc_read(filename,"xc",  x_chk)
        call nc_read(filename,"yc",  y_chk)
        call nc_read(filename,"zeta",zeta_chk)

        if (maxval(abs(x_chk-els%map%x)) .gt. 1.0e-6_wp*els%map%dx .or. &
            maxval(abs(y_chk-els%map%y)) .gt. 1.0e-6_wp*els%map%dy) then
            write(*,*) "elsa_restart_read:: Error: restart axes do not match the current grid."
            error stop 1
        end if

        if (maxval(abs(zeta_chk-els%map%zeta)) .gt. 1.0e-8_wp) then
            write(*,*) "elsa_restart_read:: Error: restart zeta does not match the host's."
            error stop 1
        end if

        call nc_read(filename,"time",          els%now%time,          start=[1],count=[1])
        call nc_read(filename,"n_top",         els%now%n_top,         start=[1],count=[1])
        call nc_read(filename,"i_add",         els%now%i_add,         start=[1],count=[1])
        call nc_read(filename,"n_reseed_total",els%now%n_reseed_total,start=[1],count=[1])

        call nc_read(filename,"d_iso",     els%now%d_iso)
        call nc_read(filename,"H_ice_prev",els%now%H_ice_prev)

    end subroutine elsa_restart_read_state

    subroutine elsa_write_init(els,filename,time)
        ! Create the output file and write everything that does not change.
        type(elsa_class), intent(in) :: els
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time

        integer :: k, n_init
        real(wp), allocatable :: layer_time(:), layer_index(:)

        call nc_create(filename)

        call nc_write_dim(filename,"xc",   x=els%map%x,units="m")
        call nc_write_dim(filename,"yc",   x=els%map%y,units="m")
        call nc_write_dim(filename,"layer",x=1.0_wp,dx=1.0_wp,nx=els%par%n_layers,units="1")
        call nc_write_dim(filename,"time", x=time,dx=1.0_wp,nx=1,units="years",unlimited=.TRUE.)

        ! The time at which each layer was created. Layer n_layers_init+1 is laid
        ! down at time_init; each subsequent layer at its entry in time_add.
        n_init = els%par%n_layers_init
        allocate(layer_time(els%par%n_layers))

        layer_time(1:n_init) = MV
        layer_time(n_init+1) = time
        do k = 1, size(els%par%time_add)
            layer_time(n_init+1+k) = els%par%time_add(k)
        end do

        call nc_write(filename,"layer_time",layer_time,dim1="layer", &
                      units="years",long_name="Time at which the layer was laid down", &
                      missing_value=MV)

        deallocate(layer_time)

    end subroutine elsa_write_init

    subroutine elsa_write_step(els,filename,time,n)
        ! Append one time slice. `n` is the 1-based index along the time axis.
        !
        ! The 3D fields are written single precision: they dominate the file size
        ! and nothing downstream needs more than 7 digits of a layer thickness.
        type(elsa_class), intent(in) :: els
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        integer,          intent(in) :: n

        integer :: nx, ny, nl

        nx = els%map%nx
        ny = els%map%ny
        nl = els%par%n_layers

        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1])

        call nc_write(filename,"n_top",els%now%n_top,dim1="time",start=[n],count=[1], &
                      long_name="Index of the topmost active layer")

        call nc_write(filename,"H_ice",real(els%now%H_ice,sp), &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[nx,ny,1], &
                      units="m",long_name="Ice thickness on the elsa grid")

        call nc_write(filename,"d_iso",real(els%now%d_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",dim4="time", &
                      start=[1,1,1,n],count=[nx,ny,nl,1], &
                      units="m",long_name="Layer thickness")

        call nc_write(filename,"dsum_iso",real(els%now%dsum_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",dim4="time", &
                      start=[1,1,1,n],count=[nx,ny,nl,1], &
                      units="m",long_name="Height of the layer top above the bed")

        call nc_write(filename,"ux_iso",real(els%now%ux_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",dim4="time", &
                      start=[1,1,1,n],count=[nx,ny,nl,1], &
                      units="m/yr",long_name="Layer-mean velocity (x), acx nodes")

        call nc_write(filename,"uy_iso",real(els%now%uy_iso,sp), &
                      dim1="xc",dim2="yc",dim3="layer",dim4="time", &
                      start=[1,1,1,n],count=[nx,ny,nl,1], &
                      units="m/yr",long_name="Layer-mean velocity (y), acy nodes")

    end subroutine elsa_write_step

end module elsa_io
