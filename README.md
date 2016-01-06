CloverLeaf3D_ref
==============

The reference version of CloverLeaf3D v1.2 - OpenACC port

This repo is forked from the 2D CloverLeaf reference version at https://github.com/UK-MAC/CloverLeaf_ref from version 1.1.
The code has since been updated to mimic the style of CloverLeaf2D v1.3.
This adds the tiling scheme and a number of optimisations.

This branch contains the port to OpenACC, based off of the Kernels directives in https://github.com/UK-MAC/CloverLeaf_OpenACC.
With this version we have support for 1 GPU per MPI task.

# Build Instructions

These remain the same as CloverLeaf:
In many case just typing make in the required software directory will work. This is the case if the mpif90 and mpicc wrappers are available on the system. This is true even for the Serial and OpenMP versions.

If the MPI compilers have different names then the build process needs to notified of this by defining two environment variables, `MPI_COMPILER` and `C_MPI_COMPILER`.

For example on some Intel systems:
```
make MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc
```
Or on Cray systems:
```
make MPI_COMPILER=ftn C_MPI_COMPILER=cc
```

