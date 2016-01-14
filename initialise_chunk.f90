!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Driver for chunk initialisation.
!>  @author Wayne Gaudin
!>  @details Invokes the user specified chunk initialisation kernel.

SUBROUTINE initialise_chunk(tile)

  USE clover_module
  USE initialise_chunk_kernel_module

  IMPLICIT NONE

  INTEGER :: tile

  REAL(KIND=8) :: xmin,ymin,zmin,dx,dy,dz

  dx=(grid%xmax-grid%xmin)/float(grid%x_cells)
  dy=(grid%ymax-grid%ymin)/float(grid%y_cells)
  dz=(grid%zmax-grid%zmin)/float(grid%z_cells)

  !write(*,*) "Left:", chunk%tiles(tile)%t_left
  !write(*,*) "Bottom:", chunk%tiles(tile)%t_left
  !write(*,*) "Back:", chunk%tiles(tile)%t_left

  xmin=grid%xmin+dx*float(chunk%tiles(tile)%t_left-1)

  ymin=grid%ymin+dy*float(chunk%tiles(tile)%t_bottom-1)

  zmin=grid%zmin+dz*float(chunk%tiles(tile)%t_back-1)

  !write(*,*) xmin, ymin, zmin

  CALL initialise_chunk_kernel(chunk%tiles(tile)%t_xmin,    &
                               chunk%tiles(tile)%t_xmax,    &
                               chunk%tiles(tile)%t_ymin,    &
                               chunk%tiles(tile)%t_ymax,    &
                               chunk%tiles(tile)%t_zmin,    &
                               chunk%tiles(tile)%t_zmax,    &
                               xmin,ymin,zmin,dx,dy,dz,      &
                               chunk%tiles(tile)%field%vertexx,  &
                               chunk%tiles(tile)%field%vertexdx, &
                               chunk%tiles(tile)%field%vertexy,  &
                               chunk%tiles(tile)%field%vertexdy, &
                               chunk%tiles(tile)%field%vertexz,  &
                               chunk%tiles(tile)%field%vertexdz, &
                               chunk%tiles(tile)%field%cellx,    &
                               chunk%tiles(tile)%field%celldx,   &
                               chunk%tiles(tile)%field%celly,    &
                               chunk%tiles(tile)%field%celldy,   &
                               chunk%tiles(tile)%field%cellz,    &
                               chunk%tiles(tile)%field%celldz,   &
                               chunk%tiles(tile)%field%volume,   &
                               chunk%tiles(tile)%field%xarea,    &
                               chunk%tiles(tile)%field%yarea,    &
                               chunk%tiles(tile)%field%zarea     )



END SUBROUTINE initialise_chunk
