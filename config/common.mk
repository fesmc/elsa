# Shared build configuration for elsa (dependency wiring).
#
# Loaded after the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# fesm-utils provides the ncio, nml, staggering and coords modules that elsa
# uses. elsa carries no netCDF calls of its own, so INC_NC is not needed at
# compile time; only LIB_NC at link time, via ncio inside libfesmutils.
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

# LIS (Library of Iterative Solvers) solves elsa's implicit advection system.
# It is vendored and built by fesm-utils, so elsa needs no system-wide install.
# The lisf.h header is pulled in by an #include, so INC_LIS is a compile flag;
# sources that need it are named .F90 so both gfortran and ifx preprocess them
# without a compiler-specific -cpp / -fpp flag.
LISROOT = ${FESMUTILSROOT}/lis/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib -llis

# elsa needs the fesm-utils and LIS includes in its compile flags.
FFLAGS += $(INC_FESMUTILS) $(INC_LIS)

ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    LISROOT = ${FESMUTILSROOT}/lis/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib -llis

    FFLAGS += $(FFLAGS_OPENMP)
endif

# elsa build outputs, for downstream linking by a consumer (e.g. yelmox, which
# clones this checkout at yelmox/elsa alongside yelmo and FastIsostasy).
ELSAROOT = ${CURDIR}
INC_ELSA = -I${ELSAROOT}/libelsa/include
LIB_ELSA = -L${ELSAROOT}/libelsa/include -lelsa

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_LIS) $(LFLAGS_EXTRA)
