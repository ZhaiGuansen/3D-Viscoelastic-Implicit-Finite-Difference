# --- Cross-platform Makefile for Fortran Project ---

# 编译器与编译选项
FC ?= gfortran
FFLAGS = -O2 -Wall -fopenmp -m64 -ffree-line-length-none
LDFLAGS = -llapack -lblas

# 源文件与目标文件
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

# 系统平台识别
UNAME_S := $(shell uname -s)

ifeq ($(OS),Windows_NT)
    # Windows 原生（非 WSL）
    EXE := $(EXE).exe
    ENABLE_LARGE_ADDRESS := yes
    IS_WINDOWS := true
else ifeq ($(UNAME_S),Linux)
    # Linux or WSL
    ifneq (,$(findstring Microsoft,$(shell uname -r)))
        # WSL
        IS_WSL := true
    else
        IS_LINUX := true
    endif
else ifeq ($(UNAME_S),Darwin)
    # macOS
    IS_MAC := true
endif

# stdbuf 适配
ifeq ($(IS_LINUX),true)
    STDBUF := stdbuf -o0 -e0
else ifeq ($(IS_WSL),true)
    STDBUF := stdbuf -o0 -e0
else
    STDBUF :=
endif

.PHONY: all clean run set-env

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
	$(STDBUF) ./$(EXE)

set-env:
ifeq ($(IS_WINDOWS),true)
	@set OMP_WAIT_POLICY=active
	@set GFORTRAN_STACKSIZE=4G
	@set OMP_STACKSIZE=4G
else
	@export OMP_WAIT_POLICY=active
	@export GFORTRAN_STACKSIZE=4G
	@export OMP_STACKSIZE=4G
endif
