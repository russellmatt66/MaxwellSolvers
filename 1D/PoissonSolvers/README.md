# Overview
Numerically solve Poisson's equation in 1D with periodic boundary conditions.

# Directory Structure
fft/ 
- Solve Poisson's equation using the Fast Fourier Transform

finitedifferences/
- Solve Poisson's equation using matrices based on a finite difference stencil

iterative/
- Solve Poisson's equation using an iterative procedure

test/
- Validate the codes

src/
- Driver and CMakeLists.txt

build/
- Build the project in here with cmake