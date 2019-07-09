[![DOI](https://zenodo.org/badge/92557777.svg)](https://zenodo.org/badge/latestdoi/92557777)

# mini-Castro

mini-Castro is a stripped-down version of Castro meant to serve as a
mini-app for exploring GPU offloading.

Castro is available at:

https://github.com/AMReX-Astro/Castro

## Compiling

mini-Castro depends on the AMReX library and uses its build system. Compiling
is done in the Exec/ directory using `make`. By default a serial build will
be generated. MPI can be enabled with `make USE_MPI=TRUE`, which will require
valid MPI C/C++/Fortran compilers in your PATH. CUDA can be enabled by adding
`USE_CUDA=TRUE`, which will require both `nvcc` and a CUDA Fortran compiler
(either `pgfortran` or `xlf`). Parallel builds with `make -j` are acceptable.
The compiler used is specified with COMP (PGI or IBM respectively).

Below are instructions for compiling on various systems:

### Compiling on bender

```
module load gcc/7.3
make CUDA_ARCH=60 COMPILE_CUDA_PATH=/usr/local/cuda-10.0 USE_MPI=FALSE USE_CUDA=TRUE COMP=PGI -j 4
```

### Compiling on groot

```
module load gcc/7.3
make CUDA_ARCH=70 COMPILE_CUDA_PATH=/usr/local/cuda-10.0 USE_MPI=FALSE USE_CUDA=TRUE COMP=PGI -j 4
```

### Compiling on Summitdev (OLCF)

```
# IBM's xlf compiler is supported, but we have found
# that it results in much lower performance than pgfortran
module load pgi/18.7
# Note that CUDA 9.2 is not currently supported due to a compiler issue.
module load cuda/9.0.69
make USE_CUDA=TRUE USE_MPI=TRUE COMP=PGI -j 4
```

## Running

mini-Castro has a few options to control the amount of work done
in the problem (determined by the spatial resolution of the
simulation) and the load balancing across MPI ranks. Run the program
after compiling to get a summary of the available options.

### Running mini-Castro on Summit/Summitdev

The following command requests an interactive job for 30 minutes on one node.

`bsub -P [project ID] -XF -nnodes 1 -W 30 -Is $SHELL`

Then launch mini-Castro using a `jsrun` command similar to the following,
which uses 4 GPUs on a single node:

`jsrun -n 4 -a 1 -g 1 ./mini-Castro3d.pgi.MPI.CUDA.ex`

To run on multiple nodes, request them using the `-nnodes` option to `bsub`.
For 1 MPI task per GPU, and 4 nodes with 4 GPUs per node, launch
mini-Castro via a jsrun command like:

`jsrun -n 16 -a 1 -g 1 -r 4 ./mini-Castro3d.pgi.MPI.CUDA.ex`

Where the '-r' option specifies the number of 1 MPI task/1 GPU
pairings (i.e. resource sets) per node.

## History

mini-Castro was originally called StarLord.  The name change reflects
its heritage.
