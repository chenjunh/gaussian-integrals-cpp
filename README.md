# Gaussian Integrals in C++

This project implements Hermite Gaussian coefficient recursion, an essential component in evaluating quantum chemistry integrals.

## Features

- Hermite coefficient recurrence using Obara–Saika relation.
- Ready for extension to overlap, kinetic, and nuclear attraction integrals.

## Build Instructions

### Prerequisites
- C++17 compiler (GCC/Clang/MSVC)
- CMake >= 3.10
- Armadillo (optional, not actively used)

### Build

```bash
mkdir build && cd build
cmake ..
make
./hermite_test
