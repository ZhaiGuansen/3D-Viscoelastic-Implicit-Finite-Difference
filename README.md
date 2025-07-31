# 3D Viscoelastic Implicit Finite-Difference

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

### 1.Clone the repository

### 2.Compile the code

```bash
make clean      # Clean previous builds
make FC=gfortran  # Compile with gfortran
```

### 3.Run the program

```bash
make run
```

## Code Description

### Parameter Settings

./src/main.f90 line18-28

### Model Data Settings

./src/partial_main.f90 line43-63

### Working Example

<table>
  <tr>
    <th>Parameter</th>
    <th>Value</th>
    <th>Units</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>Grid dimensions</td>
    <td>43×43×43</td>
    <td>-</td>
    <td>Excluding 15-layer PML</td>
  </tr>
  <tr>
    <td>Source position</td>
    <td>(22,22,22)</td>
    <td>grid</td>
    <td>-</td>
  </tr>
    <tr>
    <td>Spatial step</td>
    <td>1.0</td>
    <td>m</td>
    <td>-</td>
  </tr>
    <tr>
    <td>Time step</td>
    <td>2×10⁻⁴</td>
    <td>s</td>
    <td>-</td>
  </tr>
    <tr>
    <td>P-wave velocity</td>
    <td>3500</td>
    <td>m/s</td>
    <td>Surrounding rock</td>
  </tr>
    <tr>
    <td>S-wave velocity</td>
    <td>2100</td>
    <td>m/s</td>
    <td>Surrounding rock</td>
  </tr>
    <tr>
    <td>Q<sub>p</sub></td>
    <td>200</td>
    <td>-</td>
    <td>Quality factor (P-wave)</td>
  </tr>
    <tr>
    <td>Q<sub>s</sub></td>
    <td>100</td>
    <td>-</td>
    <td>Quality factor (S-wave)</td>
  </tr>
</table>

## Initial version:Compile and run this Fortran program using Visual Studio 2022

Installation of the library mkl_lapack95 and parallel computing with openmp

## run.ps1

After compiling with VS, run _run.ps1_ to automatically run the program(mainly used to automatically continue running the program after an interruption.)

## Modified on July 8, 2025.

Uploaded the latest code modifications and improved the Fortran code comments.

