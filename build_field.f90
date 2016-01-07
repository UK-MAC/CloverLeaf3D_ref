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

!>  @brief  Allocates the data for each mesh chunk
!>  @author Wayne Gaudin
!>  @details The data fields for the mesh chunk are allocated based on the mesh
!>  size.

MODULE build_field_module

CONTAINS

SUBROUTINE build_field(tile)

   USE clover_module

   IMPLICIT NONE

   INTEGER :: tile
   INTEGER :: z

   !write(*,*) "Tile ", tile, chunk%tiles(tile)%t_xmin, chunk%tiles(tile)%t_xmax, chunk%tiles(tile)%t_ymin, chunk%tiles(tile)%t_ymax, &
   !         chunk%tiles(tile)%t_zmin, chunk%tiles(tile)%t_zmax



   ALLOCATE(chunk%tiles(tile)%field%density0  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%density1  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%energy0   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%energy1   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%pressure  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%viscosity (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%soundspeed(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))

   ALLOCATE(chunk%tiles(tile)%field%xvel0(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%xvel1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%yvel0(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%yvel1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%zvel0(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%zvel1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))


   ALLOCATE(chunk%tiles(tile)%field%vol_flux_x (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%mass_flux_x(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vol_flux_y (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%mass_flux_y(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vol_flux_z (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%mass_flux_z(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))

   ALLOCATE(chunk%tiles(tile)%field%work_array1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array2(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array3(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array4(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array5(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array6(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array7(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                            chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                            chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))

   ALLOCATE(chunk%tiles(tile)%field%cellx   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2))
   ALLOCATE(chunk%tiles(tile)%field%celly   (chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2))
   ALLOCATE(chunk%tiles(tile)%field%cellz   (chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vertexx (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexy (chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexz (chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%celldx  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2))
   ALLOCATE(chunk%tiles(tile)%field%celldy  (chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2))
   ALLOCATE(chunk%tiles(tile)%field%celldz  (chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vertexdx(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexdy(chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexdz(chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%volume  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%xarea   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%yarea   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%zarea   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))

   ! Zeroing isn't strictly neccessary but it ensures physical pages
   ! are allocated. This prevents first touch overheads in the main code
   ! cycle which can skew timings in the first step

   !$OMP PARALLEL

!$OMP DO
   DO z=chunk%tiles(tile)%t_zmin-2, chunk%tiles(tile)%t_zmax+3
     chunk%tiles(tile)%field%work_array1(:,:,z)=0.0
     chunk%tiles(tile)%field%work_array2(:,:,z)=0.0
     chunk%tiles(tile)%field%work_array3(:,:,z)=0.0
     chunk%tiles(tile)%field%work_array4(:,:,z)=0.0
     chunk%tiles(tile)%field%work_array5(:,:,z)=0.0
     chunk%tiles(tile)%field%work_array6(:,:,z)=0.0
     chunk%tiles(tile)%field%work_array7(:,:,z)=0.0
   END DO
!$OMP END DO

!$OMP DO
   DO z=chunk%tiles(tile)%t_zmin-2, chunk%tiles(tile)%t_zmax+2
   chunk%tiles(tile)%field%density0(:,:,z)=0.0
   chunk%tiles(tile)%field%density1(:,:,z)=0.0
   chunk%tiles(tile)%field%energy0(:,:,z)=0.0
   chunk%tiles(tile)%field%energy1(:,:,z)=0.0
   chunk%tiles(tile)%field%pressure(:,:,z)=0.0
   chunk%tiles(tile)%field%viscosity(:,:,z)=0.0
   chunk%tiles(tile)%field%soundspeed(:,:,z)=0.0

   chunk%tiles(tile)%field%volume(:,:,z)=0.0
   chunk%tiles(tile)%field%xarea(:,:,z)=0.0
   chunk%tiles(tile)%field%yarea(:,:,z)=0.0
   END DO
!$OMP END DO

!$OMP DO
   DO z=chunk%tiles(tile)%t_zmin-2, chunk%tiles(tile)%t_zmax+3
   chunk%tiles(tile)%field%xvel0(:,:,z)=0.0
   chunk%tiles(tile)%field%xvel1(:,:,z)=0.0
   chunk%tiles(tile)%field%yvel0(:,:,z)=0.0
   chunk%tiles(tile)%field%yvel1(:,:,z)=0.0
   chunk%tiles(tile)%field%zvel0(:,:,z)=0.0
   chunk%tiles(tile)%field%zvel1(:,:,z)=0.0
   END DO
!$OMP END DO


!$OMP DO
   DO z=chunk%tiles(tile)%t_zmin-2, chunk%tiles(tile)%t_zmax+2
   chunk%tiles(tile)%field%vol_flux_x(:,:,z)=0.0
   chunk%tiles(tile)%field%mass_flux_x(:,:,z)=0.0
   chunk%tiles(tile)%field%vol_flux_y(:,:,z)=0.0
   chunk%tiles(tile)%field%mass_flux_y(:,:,z)=0.0
   END DO
!$OMP END DO
   
!$OMP DO
   DO z=chunk%tiles(tile)%t_zmin-2, chunk%tiles(tile)%t_zmax+3
   chunk%tiles(tile)%field%vol_flux_z(:,:,z)=0.0
   chunk%tiles(tile)%field%mass_flux_z(:,:,z)=0.0
   
   chunk%tiles(tile)%field%zarea(:,:,z)=0.0
   END DO
!$OMP END DO

!$OMP END PARALLEL

   chunk%tiles(tile)%field%cellx=0.0
   chunk%tiles(tile)%field%celly=0.0
   chunk%tiles(tile)%field%cellz=0.0
   chunk%tiles(tile)%field%vertexx=0.0
   chunk%tiles(tile)%field%vertexy=0.0
   chunk%tiles(tile)%field%vertexz=0.0
   chunk%tiles(tile)%field%celldx=0.0
   chunk%tiles(tile)%field%celldy=0.0
   chunk%tiles(tile)%field%celldz=0.0
   chunk%tiles(tile)%field%vertexdx=0.0
   chunk%tiles(tile)%field%vertexdy=0.0
   chunk%tiles(tile)%field%vertexdz=0.0



  
END SUBROUTINE build_field

END MODULE build_field_module
