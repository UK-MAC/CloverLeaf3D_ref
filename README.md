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

# Tiles

The use of tiles for CloverLeaf3D has been shown to have potential performance gains.
By increasing the tile count you decrease the size of the tile, increasing the cache efficiency when processing the tile.
This can have serious gains for certain cache expensive kernels.
However this comes at the cost of additional halo communications costs.
A tradeoff must be made between these two factors to achieve an overall performance gain.

Running the problem with a a single tile (per MPI rank), the default behaviour, will mimic the original untiled version of the code.

Currently the OpenACC implementation only supports a single Tile per MPI rank.
This is because the GPU implementation makes use of the loop level threading.

This capability is currently in development.
