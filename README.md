#     3D Viscoelastic Implicit Finite-Difference
This repository contains a Fortran-based geophysical modeling software originally developed for Windows/Visual Studio, now ported to Linux/macOS with Makefile support. The software performs forward modeling and outputs results that can be processed with MATLAB and Python into SEGY format.


## Table of Contents
- [Installation](#installation)
- [Compilation & Execution](#compilation--execution)
- [Code Description](#code-description)
  - [Parameter Settings](#parameter-settings)
  - [Model Data Settings](#model-data-settings)

## Installation

### Linux (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install build-essential gfortran make
sudo apt install liblapack-dev libblas-dev
```

### macOS 
```bash
brew update
brew install gcc make  # gcc includes gfortran
brew install lapack
```
## Compilation & Execution

### Clone the repository

### Compile the code
```bash
make clean      # Clean previous builds
make FC=gfortran  # Compile with gfortran
```
### Run the program
```bash
make run
```
## Code Description

### Parameter Settings

./src/main.f90 line18-28

### Model Data Settings

./src/partial_main.f90 line43-63

## Initial version:Compile and run this Fortran program using Visual Studio 2022
Installation of the library mkl_lapack95 and parallel computing with openmp

## run.ps1
After compiling with VS, run _run.ps1_ to automatically run the program(mainly used to automatically continue running the program after an interruption.)

## Modified on July 8, 2025.
Uploaded the latest code modifications and improved the Fortran code comments.

