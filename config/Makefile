.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libelsa/include
bindir = libelsa/bin

# Command-line options at make call
debug  ?= 0
openmp ?= 0

## COMPILER CONFIGURATION ##
# Assembled by `configme`. The fragment defines the compiler (FC, FFLAGS,
# FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and any
# host-specific overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC and LFLAGS_EXTRA set there. Defines the FESMUTILS and
# LIS paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

# debug=3: as debug=1, but without trapping underflow.
#
# elsa advects layer thickness with an explicit upwind scheme. A cell with
# outflow through one face and no inflow drains geometrically, so its thickness
# decays through the denormal range on its way to zero. That is IEEE-defined and
# harmless, but gfortran's DFLAGS_DEBUG traps it and aborts. The substitution
# below drops that one trap and leaves -fcheck=all and the invalid / zero /
# overflow traps in place. It is a no-op for the Intel fragments, whose -fpe0
# does not trap underflow to begin with.
comma := ,
DFLAGS_DEBUG_ELSA = $(subst $(comma)underflow,,$(DFLAGS_DEBUG))
ifeq ($(debug), 3)
	DFLAGS = $(DFLAGS_DEBUG_ELSA)
endif

###############################################
##
## List of elsa rules and source files
##
###############################################

$(objdir)/elsa_precision.o: $(srcdir)/elsa_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/elsa_physics.o: $(srcdir)/elsa_physics.f90 $(objdir)/elsa_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/elsa_interp.o: $(srcdir)/elsa_interp.f90 $(objdir)/elsa_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/elsa.o: $(srcdir)/elsa.f90 $(objdir)/elsa_precision.o $(objdir)/elsa_physics.o \
				  $(objdir)/elsa_interp.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

obj_elsa =	$(objdir)/elsa_precision.o \
			$(objdir)/elsa_physics.o \
			$(objdir)/elsa_interp.o \
			$(objdir)/elsa.o

###############################################
##
## Compilation of complete programs
##
###############################################

# Static library compilation
elsa-static: $(obj_elsa)
	ar rc $(objdir)/libelsa.a $(obj_elsa)
	ranlib $(objdir)/libelsa.a
	@echo " "
	@echo "    $(objdir)/libelsa.a is ready."
	@echo " "

# Benchmark: the elsa_physics kernels. Advection of a blob under uniform,
# rotational, convergent and divergent flow, plus the layer bookkeeping.
physics: elsa-static
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_physics.x $(srcdir)/test_physics.f90 \
		$(LIB_ELSA) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_physics.x is ready."
	@echo " "

# Benchmark: the elsa_interp maps. Conservative horizontal remap, bilinear
# velocity onto faces, and the layer-mean vertical integral.
interp: elsa-static
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_interp.x $(srcdir)/test_interp.f90 \
		$(LIB_ELSA) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_interp.x is ready."
	@echo " "

.PHONY: all usage clean check run-physics run-interp
all: physics interp

###############################################
##
## Running the benchmarks
##
###############################################

run-physics: physics
	./$(bindir)/test_physics.x

run-interp: interp
	./$(bindir)/test_interp.x

# Run every benchmark. Any of them failing exits nonzero and fails here.
check: run-physics run-interp
	@echo " "
	@echo "    All benchmarks passed."
	@echo " "

usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make elsa-static  : compiles the static library libelsa.a"
	@echo " make physics      : compiles test_physics.x, the kernel benchmarks"
	@echo " make interp       : compiles test_interp.x, the interpolation benchmarks"
	@echo " make all          : compiles everything"
	@echo ""
	@echo " make run-physics  : runs the kernel benchmarks"
	@echo " make run-interp   : runs the interpolation benchmarks"
	@echo " make check        : runs all benchmarks; nonzero exit on any failure"
	@echo ""
	@echo " make clean        : cleans object files and programs"
	@echo ""
	@echo " debug=1 : bounds checking and floating-point traps (traps underflow)"
	@echo " debug=2 : optimized, with profiling output"
	@echo " debug=3 : as debug=1, without the underflow trap -- use this one"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM $(bindir)/*.x.dSYM
