# StarLord

StarLord is a stripped-down version of Castro meant to serve as a
mini-app for exploring GPU offloading.

Castro is available at:

https://github.com/AMReX-Astro/Castro

# Running on GPUs

StarLord can run some algorithms on GPUs. Some algorithms are ported using CUDA
Fortran (supported only on the PGI Fortran compiler), and others are ported
using OpenACC. Below are instructions for compiling for each:

## Compiling on Titan (OLCF)

Nothing works.

## Compiling CUDA Fortran on summitdev (OLCF)

First, swap the `xl` module for `pgi`. Only the PGI Fortran compiler
supports CUDA Fortran. Then load the CUDA 9 module. The latest
versions of these modules tested on summitdev are:

- pgi/17.10
- cuda/9.0.69

To compile with CUDA 8, use the following modules:

- pgi/17.9
- cuda/8.0.61-1

Then after setting the modules, in the GNUmakefile, set
`USE_CUDA=TRUE` and type `make`. The code should now compile and run.

*NOTE*: If using NVIDIA's profiling tool `nvprof` on StarLord with
CUDA 8, it will likely encounter an error on summitdev. This is an MPI
bug in CUDA 8, which is fixed in CUDA 9. To work around this problem,
compile the code with CUDA 8, but before running, swap the CUDA 8
module with CUDA 9. Now the code should run with `nvprof`.
