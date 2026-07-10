module elsa
    ! Public facade of the Englacial Layer Simulation Architecture.
    !
    ! A host ice-sheet model should `use elsa` and nothing else: every type and
    ! procedure it needs is re-exported here. The internal modules (elsa_defs,
    ! elsa_interp, elsa_physics, elsa_io) are implementation detail.
    !
    ! Method: Born (2017); Born and Robinson (2021); Rieckh et al. (2024),
    ! Geosci. Model Dev. 17, 6987-7000. See docs/DESIGN.md for where this
    ! implementation deliberately departs from the published v2.0 scheme.

    use elsa_precision, only : sp, dp, wp

    implicit none

    private

    public :: sp, dp, wp
    public :: elsa_version

    character(len=*), parameter :: elsa_version = "3.0.0-dev"

end module elsa
