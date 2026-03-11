# Makefile for ATLAS12 + SYNTHE
#
# Usage:
#   make              Build both atlas12c.exe and synthe_spectrv.exe
#   make atlas        Build atlas12c.exe only
#   make synthe       Build synthe_spectrv.exe only
#   make clean        Remove .o and .mod files (keeps executables)
#
# Override optimization:
#   make OPT=-O3

FC       = gfortran
OPT     ?= -O3
FFLAGS   = $(OPT) -ftree-vectorize -march=native -ffpe-summary=none

# Shared
MODOBJ   = atlas12_modules.o

# ATLAS12
ATLAS    = atlas12c.exe
ATLSRC   = atlas12c.f90
ATLOBJ   = atlas12c.o

# SYNTHE
SYNTHE   = synthe_spectrv.exe
SYNSRC   = synthe_module.f90 synthe_spectrv.f90
SYNMOBJ  = synthe_module.o
SYNPOBJ  = synthe_spectrv.o

all: $(ATLAS) $(SYNTHE)

atlas: $(ATLAS)

synthe: $(SYNTHE)

$(ATLAS): $(MODOBJ) $(ATLOBJ)
	$(FC) $(FFLAGS) -o $(ATLAS) $(MODOBJ) $(ATLOBJ)

$(SYNTHE): $(MODOBJ) $(SYNMOBJ) $(SYNPOBJ)
	$(FC) $(FFLAGS) -o $(SYNTHE) $(MODOBJ) $(SYNMOBJ) $(SYNPOBJ)

$(MODOBJ): atlas12_modules.f90
	$(FC) $(FFLAGS) -c atlas12_modules.f90

$(ATLOBJ): $(ATLSRC) $(MODOBJ)
	$(FC) $(FFLAGS) -c $(ATLSRC)

$(SYNMOBJ): synthe_module.f90 $(MODOBJ)
	$(FC) $(FFLAGS) -c synthe_module.f90

$(SYNPOBJ): synthe_spectrv.f90 $(SYNMOBJ) $(MODOBJ)
	$(FC) $(FFLAGS) -c synthe_spectrv.f90

clean:
	rm -f *.o *.mod

.PHONY: all atlas synthe clean
