# Makefile for ATLAS12 + SYNTHE
#
# Usage:
#   make              Build both atlas12c.exe and synthe.exe
#   make atlas        Build atlas12c.exe only
#   make synthe       Build synthe.exe only
#
# Executables are installed to ../bin/
# All builds start clean (recompile everything).
#
# Override optimization:
#   make OPT=-O3
# for testing:
#  make OPT="-O0 -g -finit-real=snan -ffpe-trap=invalid -fbacktrace"
#
#---------------------------------------------------------------------

FC       = gfortran
OPT     ?= -O3
FFLAGS   = $(OPT) -ftree-vectorize -march=native -ffpe-summary=none
BINDIR   = ../bin

all: clean
	$(FC) $(FFLAGS) -c atlas12_modules.f90
	$(FC) $(FFLAGS) -c atlas12c.f90
	$(FC) $(FFLAGS) -o atlas12c.exe atlas12_modules.o atlas12c.o
	$(FC) $(FFLAGS) -c synthe_module.f90
	$(FC) $(FFLAGS) -c synthe.f90
	$(FC) $(FFLAGS) -o synthe.exe atlas12_modules.o synthe_module.o synthe.o
	mv atlas12c.exe synthe.exe $(BINDIR)/
	rm -f *.o *.mod

atlas: clean
	$(FC) $(FFLAGS) -c atlas12_modules.f90
	$(FC) $(FFLAGS) -c atlas12c.f90
	$(FC) $(FFLAGS) -o atlas12c.exe atlas12_modules.o atlas12c.o
	mv atlas12c.exe $(BINDIR)/
	rm -f *.o *.mod

synthe: clean
	$(FC) $(FFLAGS) -c atlas12_modules.f90
	$(FC) $(FFLAGS) -c synthe_module.f90
	$(FC) $(FFLAGS) -c synthe.f90
	$(FC) $(FFLAGS) -o synthe.exe atlas12_modules.o synthe_module.o synthe.o
	mv synthe.exe $(BINDIR)/
	rm -f *.o *.mod

clean:
	rm -f *.o *.mod

.PHONY: all atlas synthe clean
