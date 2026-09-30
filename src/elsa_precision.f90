module elsa_precision

    implicit none

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.d0)

    ! elsa works internally in double precision. Layer thicknesses are summed
    ! over a stack that reaches O(10^3) layers and renormalized against the
    ! host ice thickness every coupling period across O(10^5) yr, and a restart
    ! must reproduce the continuous run bit for bit. The wider type costs only
    ! memory. Host fields are accepted in either precision through the
    ! generic interfaces in `elsa` and converted at the boundary.
    integer, parameter :: wp = dp

    ! Missing value aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_wp
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT

    ! Below this ice thickness a column is treated as ice-free: no layers, no
    ! advection, zero velocity. Guards the division in the layer normalization.
    real(wp), parameter :: H_ICE_MIN = 1e-6_wp

end module elsa_precision
