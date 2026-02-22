# Makefile for ATLAS12
#
# Usage:
#   make                  Build with gfortran (default)
#   make FC=ifort         Build with ifort
#   make FC=ifx           Build with ifx
#   make clean            Remove build artifacts
#
# Override optimization:
#   make OPT=-O3

# Default to gfortran (override make's built-in FC)
ifndef FC
  FC = gfortran
endif
ifeq ($(FC),f77)
  FC = gfortran
endif

OPT     ?= -O3
TARGET   = atlas12c.exe
SRC      = atlas12c.f90

# Compiler-specific flags
ifeq ($(FC),gfortran)
  FFLAGS = $(OPT) -ftree-vectorize -march=native
else
  # ifort, ifx, and others
  FFLAGS = $(OPT) -ipo
endif

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) *.o *.mod

.PHONY: clean
