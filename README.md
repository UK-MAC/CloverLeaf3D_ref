CloverLeaf3D_ref
==============

The reference version of CloverLeaf3D v1.2

This repo is forked from the 2D CloverLeaf reference version at https://github.com/UK-MAC/CloverLeaf_ref from version 1.1.
The code has since been updated to mimic the style of CloverLeaf2D v1.3.
This adds the tiling scheme and a number of optimisations.


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

Tiles can be controlled in the following way:

### Tile Count

Using the input deck you can specify the number of tiles you want to include in the problem.
`tiles_per_chunk` Specifies how many tiles each MPI rank should have.
`tiles_per_problem` Specifies how many tiles are to be used across all of the MPI ranks combined. This will automatically calculate the `tiles_per_chunk` at runtime based on the number of MPI ranks. If the total number of tiles is not divisable by the number of ranks, the floor of that value will be taken and usued.


### Decomposition 

In CloverLeaf3D the default tile decomposition is 3D with the closest match to the problem dimensions.
This can be controlled at runtime to force the use of either a 3D, 2D or 1D decomposition:
* `1d_tile_decomposition`
* `2d_tile_decomposition`
* `3d_tile_decomposition`
