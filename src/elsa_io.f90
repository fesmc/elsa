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

contains

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
