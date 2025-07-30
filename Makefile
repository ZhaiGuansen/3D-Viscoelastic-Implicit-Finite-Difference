# --- Makefile for Fortran Project ---

FC ?= gfortran
FFLAGS = -O2 -Wall -fopenmp -m64 -ffree-line-length-none
# FFLAGS = -O2 -ffree-line-length-none
LDFLAGS = -llapack -lblas

SRC = src/DIA_SM.f90 \
      src/INI.f90 \
      src/matrix_construction.f90 \
      src/PML.f90 \
      src/PDC.f90 \
      src/cube_of_five_meshes.f90 \
      src/make_cube.f90 \
      src/partial_main.f90 \
      src/main.f90

OBJ = $(SRC:.f90=.o)

EXE = main

ifeq ($(OS),Windows_NT)
    EXE := $(EXE).exe
    ENABLE_LARGE_ADDRESS := yes
endif

.PHONY: all clean run set-large

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
ifeq ($(ENABLE_LARGE_ADDRESS),yes)
	@echo "Setting LARGEADDRESSAWARE flag for Windows..."
	@where editbin >nul 2>&1 && editbin /LARGEADDRESSAWARE $(EXE) || echo "editbin not found, please install Visual Studio Build Tools."
endif

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f src/*.o *.mod $(EXE)

run: $(EXE)
	@echo "Running..."
	@stdbuf -o0 -e0 ./$(EXE)
# 	@GFORTRAN_STACKSIZE=4G OMP_STACKSIZE=4G  ulimit -v unlimited && ./$(EXE)


set-env:
ifeq ($(OS),Windows_NT)
	@set OMP_WAIT_POLICY=active
	@set GFORTRAN_STACKSIZE=4G
	@set OMP_STACKSIZE=4G
else
	@export OMP_WAIT_POLICY=active
	@export GFORTRAN_STACKSIZE=4G
	@export OMP_STACKSIZE=4G
endif
