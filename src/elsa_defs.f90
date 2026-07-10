module elsa_defs
    ! The elsa state object and its parameters.
    !
    ! elsa owns everything it needs to advance itself: its own grid, its own
    ! mapping weights, its own previous-step ice thickness, and its own notion
    ! of when an update is due. A host passes fields and an absolute time; it
    ! never manages elsa's bookkeeping.

    use elsa_precision, only : wp
    use elsa_interp,    only : elsa_map_class

    implicit none

    type elsa_param_class

        ! -- read from the namelist
        integer            :: n_layers_init    ! [1]  layers present at t = time_init
        real(wp)           :: layer_resolution ! [yr] add an isochrone this often; 0 => use layer_file
        character(len=512) :: layer_file       ! [-]  one isochrone time per line; "None" => use layer_resolution
        real(wp)           :: grid_factor      ! [1]  coarsen the host grid by this factor (real, >= 1)
        real(wp)           :: dt_coupling      ! [yr] elsa update interval
        real(wp)           :: cfl              ! [1]  advection sub-step target, <= 1
        logical            :: allow_pos_bmb    ! [-]  let freeze-on thicken the bottom layer

        ! -- derived at init
        integer               :: n_layers      ! [1]  total layers allocated = n_layers_init + size(time_add) + 1
        real(wp), allocatable :: time_add(:)   ! [yr] times at which a new isochrone is created

    end type elsa_param_class

    type elsa_state_class

        real(wp) :: time     ! [yr] time of the last completed update
        integer  :: n_top    ! [1]  index of the topmost active layer
        integer  :: i_add    ! [1]  next entry of par%time_add still to be applied

        ! Layers, on elsa's grid. Layer 1 is at the bed.
        real(wp), allocatable :: d_iso(:,:,:)     ! [m] layer thickness, aa nodes
        real(wp), allocatable :: dsum_iso(:,:,:)  ! [m] height of each layer's top above the bed, aa nodes

        ! Layer-mean host velocity, on elsa's faces.
        real(wp), allocatable :: ux_iso(:,:,:)    ! [m/yr] acx nodes
        real(wp), allocatable :: uy_iso(:,:,:)    ! [m/yr] acy nodes

        ! Host fields mapped onto elsa's grid, aa nodes.
        real(wp), allocatable :: H_ice(:,:)       ! [m]    this update
        real(wp), allocatable :: H_ice_prev(:,:)  ! [m]    previous update
        real(wp), allocatable :: smb(:,:)         ! [m/yr]
        real(wp), allocatable :: bmb(:,:)         ! [m/yr]

        ! Host velocity mapped horizontally but not yet vertically. Small: nz levels.
        real(wp), allocatable :: ux_lev(:,:,:)    ! [m/yr] acx nodes, host sigma levels
        real(wp), allocatable :: uy_lev(:,:,:)    ! [m/yr] acy nodes, host sigma levels

    end type elsa_state_class

    type elsa_class
        type(elsa_param_class) :: par
        type(elsa_state_class) :: now
        type(elsa_map_class)   :: map
    end type elsa_class

    public

contains

    subroutine elsa_alloc(now,nx,ny,nz,n_layers)
        type(elsa_state_class), intent(inout) :: now
        integer,                intent(in)    :: nx, ny, nz, n_layers

        call elsa_dealloc(now)

        allocate(now%d_iso(nx,ny,n_layers))
        allocate(now%dsum_iso(nx,ny,n_layers))
        allocate(now%ux_iso(nx,ny,n_layers))
        allocate(now%uy_iso(nx,ny,n_layers))

        allocate(now%H_ice(nx,ny))
        allocate(now%H_ice_prev(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%bmb(nx,ny))

        allocate(now%ux_lev(nx,ny,nz))
        allocate(now%uy_lev(nx,ny,nz))

        ! Everything is set explicitly below, but a layer above n_top must never
        ! carry stale memory: v2.0 added its first layer over whatever was in the
        ! freshly allocated array and then accumulated smb into it.
        now%d_iso      = 0.0_wp
        now%dsum_iso   = 0.0_wp
        now%ux_iso     = 0.0_wp
        now%uy_iso     = 0.0_wp
        now%H_ice      = 0.0_wp
        now%H_ice_prev = 0.0_wp
        now%smb        = 0.0_wp
        now%bmb        = 0.0_wp
        now%ux_lev     = 0.0_wp
        now%uy_lev     = 0.0_wp

    end subroutine elsa_alloc

    subroutine elsa_dealloc(now)
        ! Frees every array elsa_alloc allocates. v2.0's equivalent freed seven
        ! of its nine, so a second init aborted on an already-allocated array,
        ! which is why it could not restart.
        type(elsa_state_class), intent(inout) :: now

        if (allocated(now%d_iso))      deallocate(now%d_iso)
        if (allocated(now%dsum_iso))   deallocate(now%dsum_iso)
        if (allocated(now%ux_iso))     deallocate(now%ux_iso)
        if (allocated(now%uy_iso))     deallocate(now%uy_iso)
        if (allocated(now%H_ice))      deallocate(now%H_ice)
        if (allocated(now%H_ice_prev)) deallocate(now%H_ice_prev)
        if (allocated(now%smb))        deallocate(now%smb)
        if (allocated(now%bmb))        deallocate(now%bmb)
        if (allocated(now%ux_lev))     deallocate(now%ux_lev)
        if (allocated(now%uy_lev))     deallocate(now%uy_lev)

    end subroutine elsa_dealloc

end module elsa_defs
