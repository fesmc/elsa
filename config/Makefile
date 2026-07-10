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

###############################################
##
## List of elsa rules and source files
##
###############################################

$(objdir)/elsa_precision.o: $(srcdir)/elsa_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/elsa.o: $(srcdir)/elsa.f90 $(objdir)/elsa_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

obj_elsa =	$(objdir)/elsa_precision.o \
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

.PHONY: all usage clean check
all: elsa-static

# Run every benchmark. Any of them failing exits nonzero and fails here.
check: all
	@echo " "
	@echo "    No benchmarks yet."
	@echo " "

usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make elsa-static  : compiles the static library libelsa.a"
	@echo " make all          : compiles everything"
	@echo " make check        : runs all benchmarks; nonzero exit on any failure"
	@echo ""
	@echo " make clean        : cleans object files and programs"
	@echo " Add debug=1 to any build for bounds checking and floating-point traps."
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM $(bindir)/*.x.dSYM
