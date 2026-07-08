# Makefile for ATLAS12 + SYNTHE
#
# Usage:
#   make              Build both atlas12c.exe and synthe.exe
#   make atlas        Build atlas12c.exe only
#   make synthe       Build synthe.exe only
#
# Compile order:
#   mod_mklinelist -> atlas12_modules -> {atlas12c, synthe_module -> synthe}
# (ATLAS12 now reads diatomic line lists via mod_mklinelist, the same
#  reader SYNTHE uses, so mklinelist must be built first for both.)
#
# Executables are installed to ../bin/
# All builds start clean (recompile everything) and leave the source
# directory free of .o/.mod files.
#
# Override optimization:
#   make OPT=-O3
# for testing:
#  make OPT="-O0 -g -finit-real=snan -ffpe-trap=invalid -fbacktrace"
#
# NOTE: -march=native tunes for the build machine's CPU.  Binaries in
# ../bin/ may emit illegal instructions if copied to a different machine
# (e.g. an older Mac or a cluster node).  Drop it for portable builds:
#   make OPT="-O3 -mtune=native"
#
#---------------------------------------------------------------------

FC       = gfortran
OPT     ?= -O3
FFLAGS   = $(OPT) -ftree-vectorize -march=native -ffpe-summary=none
BINDIR   = ../bin

# Object lists.  COMMON_OBJ holds the shared modules linked into both
# executables; they are compiled in dependency order by the recipes
# below (mod_mklinelist before atlas12_modules, which USEs it).
COMMON_OBJ = mod_mklinelist.o atlas12_modules.o
ATLAS_OBJ  = $(COMMON_OBJ) atlas12c.o
SYNTHE_OBJ = $(COMMON_OBJ) synthe_module.o synthe.o

all: clean
	$(FC) $(FFLAGS) -c mod_mklinelist.f90
	$(FC) $(FFLAGS) -c atlas12_modules.f90
	$(FC) $(FFLAGS) -c atlas12c.f90
	$(FC) $(FFLAGS) -o atlas12c.exe $(ATLAS_OBJ)
	$(FC) $(FFLAGS) -c synthe_module.f90
	$(FC) $(FFLAGS) -c synthe.f90
	$(FC) $(FFLAGS) -o synthe.exe $(SYNTHE_OBJ)
	mv atlas12c.exe synthe.exe $(BINDIR)/
	rm -f *.o *.mod

atlas: clean
	$(FC) $(FFLAGS) -c mod_mklinelist.f90
	$(FC) $(FFLAGS) -c atlas12_modules.f90
	$(FC) $(FFLAGS) -c atlas12c.f90
	$(FC) $(FFLAGS) -o atlas12c.exe $(ATLAS_OBJ)
	mv atlas12c.exe $(BINDIR)/
	rm -f *.o *.mod

synthe: clean
	$(FC) $(FFLAGS) -c mod_mklinelist.f90
	$(FC) $(FFLAGS) -c atlas12_modules.f90
	$(FC) $(FFLAGS) -c synthe_module.f90
	$(FC) $(FFLAGS) -c synthe.f90
	$(FC) $(FFLAGS) -o synthe.exe $(SYNTHE_OBJ)
	mv synthe.exe $(BINDIR)/
	rm -f *.o *.mod

# clean removes build artifacts only.  It must NOT delete the installed
# executables: 'make atlas' and 'make synthe' would otherwise wipe the
# other program's binary out of ../bin/.
clean:
	rm -f *.o *.mod

.PHONY: all atlas synthe clean
