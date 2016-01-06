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

!>  @brief Fortran kernel to update the external halo cells in a chunk.
!>  @author Wayne Gaudin
!>  @details Updates halo cells for the required fields at the required depth
!>  for any halo cells that lie on an external boundary. The location and type
!>  of data governs how this is carried out. External boundaries are always
!>  reflective.

MODULE update_tile_halo_kernel_module

    USE data_module 

CONTAINS

    SUBROUTINE update_tile_halo_left_kernel(x_min,x_max,y_min,y_max,z_min,z_max,  &
        density0,                                                   &
        energy0,                                                    &
        pressure,                                                   &
        viscosity,                                                  &
        soundspeed,                                                 &
        density1,                                                   &
        energy1,                                                    &
        xvel0,                                                      &
        yvel0,                                                      &
        zvel0,                                                      &
        xvel1,                                                      &
        yvel1,                                                      &
        zvel1,                                                      &
        vol_flux_x,                                                 &
        vol_flux_y,                                                 &
        vol_flux_z,                                                 &
        mass_flux_x,                                                &
        mass_flux_y,                                                &
        mass_flux_z,                                                &
        left_xmin, left_xmax, left_ymin, left_ymax,  left_zmin, left_zmax,               &
        left_density0,                                                   &
        left_energy0,                                                    &
        left_pressure,                                                   &
        left_viscosity,                                                  &
        left_soundspeed,                                                 &
        left_density1,                                                   &
        left_energy1,                                                    &
        left_xvel0,                                                      &
        left_yvel0,                                                      &
        left_zvel0,                                                      &
        left_xvel1,                                                      &
        left_yvel1,                                                      &
        left_zvel1,                                                      &
        left_vol_flux_x,                                                 &
        left_vol_flux_y,                                                 &
        left_vol_flux_z,                                                 &
        left_mass_flux_x,                                                &
        left_mass_flux_y,                                                &
        left_mass_flux_z,                                                &
        fields,                                                     &
        depth)

        IMPLICIT NONE

        INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density0,energy0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure,viscosity,soundspeed
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1,energy1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel0,yvel0,zvel0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel1,yvel1,zvel1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x,mass_flux_x
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y,mass_flux_y
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z,mass_flux_z
        INTEGER :: left_xmin, left_xmax, left_ymin, left_ymax, left_zmin, left_zmax
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+2,left_ymin-2:left_ymax+2,left_zmin-2:left_zmax+2) :: left_density0,left_energy0
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+2,left_ymin-2:left_ymax+2,left_zmin-2:left_zmax+2) :: left_pressure,left_viscosity,left_soundspeed
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+2,left_ymin-2:left_ymax+2,left_zmin-2:left_zmax+2) :: left_density1,left_energy1
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+3,left_ymin-2:left_ymax+3,left_zmin-2:left_zmax+3) :: left_xvel0,left_yvel0,left_zvel0
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+3,left_ymin-2:left_ymax+3,left_zmin-2:left_zmax+3) :: left_xvel1,left_yvel1,left_zvel1
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+3,left_ymin-2:left_ymax+2,left_zmin-2:left_zmax+2) :: left_vol_flux_x,left_mass_flux_x
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+2,left_ymin-2:left_ymax+3,left_zmin-2:left_zmax+2) :: left_vol_flux_y,left_mass_flux_y
        REAL(KIND=8), DIMENSION(left_xmin-2:left_xmax+2,left_ymin-2:left_ymax+2,left_zmin-2:left_zmax+3) :: left_vol_flux_z,left_mass_flux_z

        INTEGER :: fields(:),depth

        INTEGER :: j,k,l


!$OMP PARALLEL PRIVATE(k,j)

        ! Density 0
   

        IF(fields(FIELD_DENSITY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        density0(x_min-j,k,l)=left_density0(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            END DO
            !$OMP END DO
        ENDIF

        ! Density 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        density1(x_min-j,k,l)=left_density1(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

   
        ! Energy 0
        IF(fields(FIELD_ENERGY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        energy0(x_min-j,k,l)=left_energy0(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        energy1(x_min-j,k,l)=left_energy1(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

  
        ! Pressure
        IF(fields(FIELD_PRESSURE).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        pressure(x_min-j,k,l)=left_pressure(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Viscocity
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        viscosity(x_min-j,k,l)=left_viscosity(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Soundspeed
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        soundspeed(x_min-j,k,l)=left_soundspeed(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF


        ! XVEL 0
        IF(fields(FIELD_XVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        xvel0(x_min-j,k,l)=left_xvel0(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 1
        IF(fields(FIELD_XVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        xvel1(x_min-j,k,l)=left_xvel1(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 0
        IF(fields(FIELD_YVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        yvel0(x_min-j,k,l)=left_yvel0(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 1
        IF(fields(FIELD_YVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        yvel1(x_min-j,k,l)=left_yvel1(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 0
        IF(fields(FIELD_ZVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        zvel0(x_min-j,k,l)=left_zvel0(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 1
        IF(fields(FIELD_ZVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        zvel1(x_min-j,k,l)=left_zvel1(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_X
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        vol_flux_x(x_min-j,k,l)=left_vol_flux_x(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_X
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        mass_flux_x(x_min-j,k,l)=left_mass_flux_x(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Y
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        vol_flux_y(x_min-j,k,l)=left_vol_flux_y(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Y
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        mass_flux_y(x_min-j,k,l)=left_mass_flux_y(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Z
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        vol_flux_z(x_min-j,k,l)=left_vol_flux_z(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Z
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        mass_flux_z(x_min-j,k,l)=left_mass_flux_z(left_xmax+1-j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

!$OMP END PARALLEL

    END SUBROUTINE update_tile_halo_left_kernel




    SUBROUTINE update_tile_halo_right_kernel(x_min,x_max,y_min,y_max,z_min,z_max,  &
        density0,                                                   &
        energy0,                                                    &
        pressure,                                                   &
        viscosity,                                                  &
        soundspeed,                                                 &
        density1,                                                   &
        energy1,                                                    &
        xvel0,                                                      &
        yvel0,                                                      &
        zvel0,                                                      &
        xvel1,                                                      &
        yvel1,                                                      &
        zvel1,                                                      &
        vol_flux_x,                                                 &
        vol_flux_y,                                                 &
        vol_flux_z,                                                 &
        mass_flux_x,                                                &
        mass_flux_y,                                                &
        mass_flux_z,                                                &
        right_xmin, right_xmax, right_ymin, right_ymax,  right_zmin, right_zmax,               &
        right_density0,                                                   &
        right_energy0,                                                    &
        right_pressure,                                                   &
        right_viscosity,                                                  &
        right_soundspeed,                                                 &
        right_density1,                                                   &
        right_energy1,                                                    &
        right_xvel0,                                                      &
        right_yvel0,                                                      &
        right_zvel0,                                                      &
        right_xvel1,                                                      &
        right_yvel1,                                                      &
        right_zvel1,                                                      &
        right_vol_flux_x,                                                 &
        right_vol_flux_y,                                                 &
        right_vol_flux_z,                                                 &
        right_mass_flux_x,                                                &
        right_mass_flux_y,                                                &
        right_mass_flux_z,                                                &
        fields,                                                     &
        depth)

        IMPLICIT NONE

        INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density0,energy0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure,viscosity,soundspeed
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1,energy1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel0,yvel0,zvel0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel1,yvel1,zvel1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x,mass_flux_x
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y,mass_flux_y
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z,mass_flux_z
        INTEGER :: right_xmin, right_xmax, right_ymin, right_ymax, right_zmin, right_zmax
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+2,right_ymin-2:right_ymax+2,right_zmin-2:right_zmax+2) :: right_density0,right_energy0
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+2,right_ymin-2:right_ymax+2,right_zmin-2:right_zmax+2) :: right_pressure,right_viscosity,right_soundspeed
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+2,right_ymin-2:right_ymax+2,right_zmin-2:right_zmax+2) :: right_density1,right_energy1
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+3,right_ymin-2:right_ymax+3,right_zmin-2:right_zmax+3) :: right_xvel0,right_yvel0,right_zvel0
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+3,right_ymin-2:right_ymax+3,right_zmin-2:right_zmax+3) :: right_xvel1,right_yvel1,right_zvel1
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+3,right_ymin-2:right_ymax+2,right_zmin-2:right_zmax+2) :: right_vol_flux_x,right_mass_flux_x
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+2,right_ymin-2:right_ymax+3,right_zmin-2:right_zmax+2) :: right_vol_flux_y,right_mass_flux_y
        REAL(KIND=8), DIMENSION(right_xmin-2:right_xmax+2,right_ymin-2:right_ymax+2,right_zmin-2:right_zmax+3) :: right_vol_flux_z,right_mass_flux_z

        INTEGER :: fields(:),depth

        INTEGER :: j,k,l


!$OMP PARALLEL PRIVATE(k,j)

        ! Density 0
        IF(fields(FIELD_DENSITY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        density0(x_max+j,k,l)=right_density0(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Density 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        density1(x_max+j,k,l)=right_density1(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

   
        ! Energy 0
        IF(fields(FIELD_ENERGY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        energy0(x_max+j,k,l)=right_energy0(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        energy1(x_max+j,k,l)=right_energy1(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

  
        ! Pressure
        IF(fields(FIELD_PRESSURE).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        pressure(x_max+j,k,l)=right_pressure(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Viscocity
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        viscosity(x_max+j,k,l)=right_viscosity(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Soundspeed
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        soundspeed(x_max+j,k,l)=right_soundspeed(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF


        ! XVEL 0
        IF(fields(FIELD_XVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        xvel0(x_max+1+j,k,l)=right_xvel0(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 1
        IF(fields(FIELD_XVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        xvel1(x_max+1+j,k,l)=right_xvel1(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 0
        IF(fields(FIELD_YVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        yvel0(x_max+1+j,k,l)=right_yvel0(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 1
        IF(fields(FIELD_YVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        yvel1(x_max+1+j,k,l)=right_yvel1(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 0
        IF(fields(FIELD_ZVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        yvel0(x_max+1+j,k,l)=right_yvel0(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 1
        IF(fields(FIELD_ZVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        zvel1(x_max+1+j,k,l)=right_zvel1(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_X
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        vol_flux_x(x_max+1+j,k,l)=right_vol_flux_x(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_X
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        mass_flux_x(x_max+1+j,k,l)=right_mass_flux_x(right_xmin+1-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Y
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        vol_flux_y(x_max+j,k,l)=right_vol_flux_y(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Y
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=y_min-depth,y_max+1+depth
                    DO j=1,depth
                        mass_flux_y(x_max+j,k,l)=right_mass_flux_y(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Z
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        vol_flux_z(x_max+j,k,l)=right_vol_flux_z(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Z
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=y_min-depth,y_max+depth
                    DO j=1,depth
                        mass_flux_z(x_max+j,k,l)=right_mass_flux_z(right_xmin-1+j,k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

!$OMP END PARALLEL

    END SUBROUTINE update_tile_halo_right_kernel


    ! Top and bottom only do xmin -> xmax
    ! This is because the corner ghosts will get communicated in the left right communication

    SUBROUTINE update_tile_halo_top_kernel(x_min,x_max,y_min,y_max,z_min,z_max,  &
        density0,                                                   &
        energy0,                                                    &
        pressure,                                                   &
        viscosity,                                                  &
        soundspeed,                                                 &
        density1,                                                   &
        energy1,                                                    &
        xvel0,                                                      &
        yvel0,                                                      &
        zvel0,                                                      &
        xvel1,                                                      &
        yvel1,                                                      &
        zvel1,                                                      &
        vol_flux_x,                                                 &
        vol_flux_y,                                                 &
        vol_flux_z,                                                 &
        mass_flux_x,                                                &
        mass_flux_y,                                                &
        mass_flux_z,                                                &
        top_xmin, top_xmax, top_ymin, top_ymax,  top_zmin, top_zmax,               &
        top_density0,                                                   &
        top_energy0,                                                    &
        top_pressure,                                                   &
        top_viscosity,                                                  &
        top_soundspeed,                                                 &
        top_density1,                                                   &
        top_energy1,                                                    &
        top_xvel0,                                                      &
        top_yvel0,                                                      &
        top_zvel0,                                                      &
        top_xvel1,                                                      &
        top_yvel1,                                                      &
        top_zvel1,                                                      &
        top_vol_flux_x,                                                 &
        top_vol_flux_y,                                                 &
        top_vol_flux_z,                                                 &
        top_mass_flux_x,                                                &
        top_mass_flux_y,                                                &
        top_mass_flux_z,                                                &
        fields,                                                     &
        depth)

        IMPLICIT NONE

        INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density0,energy0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure,viscosity,soundspeed
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1,energy1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel0,yvel0,zvel0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel1,yvel1,zvel1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x,mass_flux_x
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y,mass_flux_y
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z,mass_flux_z
        INTEGER :: top_xmin, top_xmax, top_ymin, top_ymax, top_zmin, top_zmax
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+2,top_ymin-2:top_ymax+2,top_zmin-2:top_zmax+2) :: top_density0,top_energy0
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+2,top_ymin-2:top_ymax+2,top_zmin-2:top_zmax+2) :: top_pressure,top_viscosity,top_soundspeed
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+2,top_ymin-2:top_ymax+2,top_zmin-2:top_zmax+2) :: top_density1,top_energy1
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+3,top_ymin-2:top_ymax+3,top_zmin-2:top_zmax+3) :: top_xvel0,top_yvel0,top_zvel0
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+3,top_ymin-2:top_ymax+3,top_zmin-2:top_zmax+3) :: top_xvel1,top_yvel1,top_zvel1
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+3,top_ymin-2:top_ymax+2,top_zmin-2:top_zmax+2) :: top_vol_flux_x,top_mass_flux_x
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+2,top_ymin-2:top_ymax+3,top_zmin-2:top_zmax+2) :: top_vol_flux_y,top_mass_flux_y
        REAL(KIND=8), DIMENSION(top_xmin-2:top_xmax+2,top_ymin-2:top_ymax+2,top_zmin-2:top_zmax+3) :: top_vol_flux_z,top_mass_flux_z

        INTEGER :: fields(:),depth

        INTEGER :: j,k,l


!$OMP PARALLEL PRIVATE(k,j)

        ! Density 0
        IF(fields(FIELD_DENSITY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        density0(j,y_max+k,l)=top_density0(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Density 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        density1(j,y_max+k,l)=top_density1(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

   
        ! Energy 0
        IF(fields(FIELD_ENERGY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        energy0(j,y_max+k,l)=top_energy0(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        energy1(j,y_max+k,l)=top_energy1(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

  
        ! Pressure
        IF(fields(FIELD_PRESSURE).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        pressure(j,y_max+k,l)=top_pressure(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Viscocity
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        viscosity(j,y_max+k,l)=top_viscosity(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Soundspeed
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        soundspeed(j,y_max+k,l)=top_soundspeed(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF


        ! XVEL 0
        IF(fields(FIELD_XVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel0(j,y_max+1+k,l)=top_xvel0(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 1
        IF(fields(FIELD_XVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel1(j,y_max+1+k,l)=top_xvel1(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 0
        IF(fields(FIELD_YVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel0(j,y_max+1+k,l)=top_yvel0(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 1
        IF(fields(FIELD_YVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel1(j,y_max+1+k,l)=top_yvel1(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 0
        IF(fields(FIELD_ZVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel0(j,y_max+1+k,l)=top_zvel0(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 1
        IF(fields(FIELD_ZVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel1(j,y_max+1+k,l)=top_zvel1(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_X
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        vol_flux_x(j,y_max+k,l)=top_vol_flux_x(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_X
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        mass_flux_x(j,y_max+k,l)=top_mass_flux_x(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Y
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_y(j,y_max+1+k,l)=top_vol_flux_y(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Y
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_y(j,y_max+1+k,l)=top_mass_flux_y(j,top_ymin+1-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Z
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_z(j,y_max+k,l)=top_vol_flux_z(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Z
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_z(j,y_max+k,l)=top_mass_flux_z(j,top_ymin-1+k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

!$OMP END PARALLEL

    END SUBROUTINE update_tile_halo_top_kernel








    SUBROUTINE update_tile_halo_bottom_kernel(x_min,x_max,y_min,y_max,z_min,z_max,  &
        density0,                                                   &
        energy0,                                                    &
        pressure,                                                   &
        viscosity,                                                  &
        soundspeed,                                                 &
        density1,                                                   &
        energy1,                                                    &
        xvel0,                                                      &
        yvel0,                                                      &
        zvel0,                                                      &
        xvel1,                                                      &
        yvel1,                                                      &
        zvel1,                                                      &
        vol_flux_x,                                                 &
        vol_flux_y,                                                 &
        vol_flux_z,                                                 &
        mass_flux_x,                                                &
        mass_flux_y,                                                &
        mass_flux_z,                                                &
        bottom_xmin, bottom_xmax, bottom_ymin, bottom_ymax,  bottom_zmin, bottom_zmax,               &
        bottom_density0,                                                   &
        bottom_energy0,                                                    &
        bottom_pressure,                                                   &
        bottom_viscosity,                                                  &
        bottom_soundspeed,                                                 &
        bottom_density1,                                                   &
        bottom_energy1,                                                    &
        bottom_xvel0,                                                      &
        bottom_yvel0,                                                      &
        bottom_zvel0,                                                      &
        bottom_xvel1,                                                      &
        bottom_yvel1,                                                      &
        bottom_zvel1,                                                      &
        bottom_vol_flux_x,                                                 &
        bottom_vol_flux_y,                                                 &
        bottom_vol_flux_z,                                                 &
        bottom_mass_flux_x,                                                &
        bottom_mass_flux_y,                                                &
        bottom_mass_flux_z,                                                &
        fields,                                                     &
        depth)

        IMPLICIT NONE

        INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density0,energy0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure,viscosity,soundspeed
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1,energy1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel0,yvel0,zvel0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel1,yvel1,zvel1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x,mass_flux_x
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y,mass_flux_y
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z,mass_flux_z
        INTEGER :: bottom_xmin, bottom_xmax, bottom_ymin, bottom_ymax, bottom_zmin, bottom_zmax
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+2,bottom_ymin-2:bottom_ymax+2,bottom_zmin-2:bottom_zmax+2) :: bottom_density0,bottom_energy0
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+2,bottom_ymin-2:bottom_ymax+2,bottom_zmin-2:bottom_zmax+2) :: bottom_pressure,bottom_viscosity,bottom_soundspeed
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+2,bottom_ymin-2:bottom_ymax+2,bottom_zmin-2:bottom_zmax+2) :: bottom_density1,bottom_energy1
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+3,bottom_ymin-2:bottom_ymax+3,bottom_zmin-2:bottom_zmax+3) :: bottom_xvel0,bottom_yvel0,bottom_zvel0
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+3,bottom_ymin-2:bottom_ymax+3,bottom_zmin-2:bottom_zmax+3) :: bottom_xvel1,bottom_yvel1,bottom_zvel1
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+3,bottom_ymin-2:bottom_ymax+2,bottom_zmin-2:bottom_zmax+2) :: bottom_vol_flux_x,bottom_mass_flux_x
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+2,bottom_ymin-2:bottom_ymax+3,bottom_zmin-2:bottom_zmax+2) :: bottom_vol_flux_y,bottom_mass_flux_y
        REAL(KIND=8), DIMENSION(bottom_xmin-2:bottom_xmax+2,bottom_ymin-2:bottom_ymax+2,bottom_zmin-2:bottom_zmax+3) :: bottom_vol_flux_z,bottom_mass_flux_z

        INTEGER :: fields(:),depth

        INTEGER :: j,k,l


!$OMP PARALLEL PRIVATE(k,j)

        ! Density 0
        IF(fields(FIELD_DENSITY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        density0(j,y_min-k,l)=bottom_density0(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Density 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        density1(j,y_min-k,l)=bottom_density1(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

   
        ! Energy 0
        IF(fields(FIELD_ENERGY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        energy0(j,y_min-k,l)=bottom_energy0(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        energy1(j,y_min-k,l)=bottom_energy1(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

  
        ! Pressure
        IF(fields(FIELD_PRESSURE).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        pressure(j,y_min-k,l)=bottom_pressure(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Viscocity
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        viscosity(j,y_min-k,l)=bottom_viscosity(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Soundspeed
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        soundspeed(j,y_min-k,l)=bottom_soundspeed(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF


        ! XVEL 0
        IF(fields(FIELD_XVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel0(j,y_min-k,l)=bottom_xvel0(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 1
        IF(fields(FIELD_XVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel1(j,y_min-k,l)=bottom_xvel1(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 0
        IF(fields(FIELD_YVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel0(j,y_min-k,l)=bottom_yvel0(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 1
        IF(fields(FIELD_YVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel1(j,y_min-k,l)=bottom_yvel1(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 0
        IF(fields(FIELD_ZVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel0(j,y_min-k,l)=bottom_zvel0(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 1
        IF(fields(FIELD_ZVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel1(j,y_min-k,l)=bottom_zvel1(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_X
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        vol_flux_x(j,y_min-k,l)=bottom_vol_flux_x(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_X
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+1+depth
                        mass_flux_x(j,y_min-k,l)=bottom_mass_flux_x(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Y
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_y(j,y_min-k,l)=bottom_vol_flux_y(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Y
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_y(j,y_min-k,l)=bottom_mass_flux_y(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Z
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_z(j,y_min-k,l)=bottom_vol_flux_z(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Z
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=z_min-depth,z_max+1+depth
                DO k=1,depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_z(j,y_min-k,l)=bottom_mass_flux_z(j,bottom_ymax+1-k,l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

!$OMP END PARALLEL

    END SUBROUTINE update_tile_halo_bottom_kernel





    SUBROUTINE update_tile_halo_front_kernel(x_min,x_max,y_min,y_max,z_min,z_max,  &
        density0,                                                   &
        energy0,                                                    &
        pressure,                                                   &
        viscosity,                                                  &
        soundspeed,                                                 &
        density1,                                                   &
        energy1,                                                    &
        xvel0,                                                      &
        yvel0,                                                      &
        zvel0,                                                      &
        xvel1,                                                      &
        yvel1,                                                      &
        zvel1,                                                      &
        vol_flux_x,                                                 &
        vol_flux_y,                                                 &
        vol_flux_z,                                                 &
        mass_flux_x,                                                &
        mass_flux_y,                                                &
        mass_flux_z,                                                &
        front_xmin, front_xmax, front_ymin, front_ymax,  front_zmin, front_zmax,               &
        front_density0,                                                   &
        front_energy0,                                                    &
        front_pressure,                                                   &
        front_viscosity,                                                  &
        front_soundspeed,                                                 &
        front_density1,                                                   &
        front_energy1,                                                    &
        front_xvel0,                                                      &
        front_yvel0,                                                      &
        front_zvel0,                                                      &
        front_xvel1,                                                      &
        front_yvel1,                                                      &
        front_zvel1,                                                      &
        front_vol_flux_x,                                                 &
        front_vol_flux_y,                                                 &
        front_vol_flux_z,                                                 &
        front_mass_flux_x,                                                &
        front_mass_flux_y,                                                &
        front_mass_flux_z,                                                &
        fields,                                                     &
        depth)

        IMPLICIT NONE

        INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density0,energy0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure,viscosity,soundspeed
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1,energy1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel0,yvel0,zvel0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel1,yvel1,zvel1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x,mass_flux_x
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y,mass_flux_y
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z,mass_flux_z
        INTEGER :: front_xmin, front_xmax, front_ymin, front_ymax, front_zmin, front_zmax
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+2,front_ymin-2:front_ymax+2,front_zmin-2:front_zmax+2) :: front_density0,front_energy0
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+2,front_ymin-2:front_ymax+2,front_zmin-2:front_zmax+2) :: front_pressure,front_viscosity,front_soundspeed
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+2,front_ymin-2:front_ymax+2,front_zmin-2:front_zmax+2) :: front_density1,front_energy1
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+3,front_ymin-2:front_ymax+3,front_zmin-2:front_zmax+3) :: front_xvel0,front_yvel0,front_zvel0
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+3,front_ymin-2:front_ymax+3,front_zmin-2:front_zmax+3) :: front_xvel1,front_yvel1,front_zvel1
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+3,front_ymin-2:front_ymax+2,front_zmin-2:front_zmax+2) :: front_vol_flux_x,front_mass_flux_x
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+2,front_ymin-2:front_ymax+3,front_zmin-2:front_zmax+2) :: front_vol_flux_y,front_mass_flux_y
        REAL(KIND=8), DIMENSION(front_xmin-2:front_xmax+2,front_ymin-2:front_ymax+2,front_zmin-2:front_zmax+3) :: front_vol_flux_z,front_mass_flux_z

        INTEGER :: fields(:),depth

        INTEGER :: j,k,l

!$OMP PARALLEL PRIVATE(k,j)

        ! Density 0
        IF(fields(FIELD_DENSITY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        density0(j,k,z_max+l)=front_density0(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Density 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        density1(j,k,z_max+l)=front_density1(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 0
        IF(fields(FIELD_ENERGY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        energy0(j,k,z_max+l)=front_energy0(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 1
        IF(fields(FIELD_ENERGY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        energy1(j,k,z_max+l)=front_energy1(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Pressure
        IF(fields(FIELD_PRESSURE).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        pressure(j,k,z_max+l)=front_pressure(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Viscosity
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        viscosity(j,k,z_max+l)=front_viscosity(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Soundspeed
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        soundspeed(j,k,z_max+l)=front_soundspeed(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 0
        IF(fields(FIELD_XVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel0(j,k,z_max+1+l)=front_xvel0(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 1
        IF(fields(FIELD_XVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel1(j,k,z_max+1+l)=front_xvel1(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 0
        IF(fields(FIELD_YVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel0(j,k,z_max+1+l)=front_yvel0(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 1
        IF(fields(FIELD_YVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel1(j,k,z_max+1+l)=front_yvel1(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 0
        IF(fields(FIELD_ZVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel0(j,k,z_max+1+l)=front_zvel0(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 1
        IF(fields(FIELD_ZVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel1(j,k,z_max+1+l)=front_zvel1(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_X
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+1+depth
                        vol_flux_x(j,k,z_max+l)=front_vol_flux_x(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_X
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+1+depth
                        mass_flux_x(j,k,z_max+l)=front_mass_flux_x(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Y
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_y(j,k,z_max+l)=front_vol_flux_y(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Y
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_y(j,k,z_max+l)=front_mass_flux_y(j,k,front_zmin-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Z
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_z(j,k,z_max+1+l)=front_vol_flux_z(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Z
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_z(j,k,z_max+1+l)=front_mass_flux_z(j,k,front_zmin+1-1+l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

!$OMP END PARALLEL

    END SUBROUTINE update_tile_halo_front_kernel



    SUBROUTINE update_tile_halo_back_kernel(x_min,x_max,y_min,y_max,z_min,z_max,  &
        density0,                                                   &
        energy0,                                                    &
        pressure,                                                   &
        viscosity,                                                  &
        soundspeed,                                                 &
        density1,                                                   &
        energy1,                                                    &
        xvel0,                                                      &
        yvel0,                                                      &
        zvel0,                                                      &
        xvel1,                                                      &
        yvel1,                                                      &
        zvel1,                                                      &
        vol_flux_x,                                                 &
        vol_flux_y,                                                 &
        vol_flux_z,                                                 &
        mass_flux_x,                                                &
        mass_flux_y,                                                &
        mass_flux_z,                                                &
        back_xmin, back_xmax, back_ymin, back_ymax,  back_zmin, back_zmax,               &
        back_density0,                                                   &
        back_energy0,                                                    &
        back_pressure,                                                   &
        back_viscosity,                                                  &
        back_soundspeed,                                                 &
        back_density1,                                                   &
        back_energy1,                                                    &
        back_xvel0,                                                      &
        back_yvel0,                                                      &
        back_zvel0,                                                      &
        back_xvel1,                                                      &
        back_yvel1,                                                      &
        back_zvel1,                                                      &
        back_vol_flux_x,                                                 &
        back_vol_flux_y,                                                 &
        back_vol_flux_z,                                                 &
        back_mass_flux_x,                                                &
        back_mass_flux_y,                                                &
        back_mass_flux_z,                                                &
        fields,                                                     &
        depth)

        IMPLICIT NONE

        INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density0,energy0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure,viscosity,soundspeed
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1,energy1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel0,yvel0,zvel0
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: xvel1,yvel1,zvel1
        REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x,mass_flux_x
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y,mass_flux_y
        REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z,mass_flux_z
        INTEGER :: back_xmin, back_xmax, back_ymin, back_ymax, back_zmin, back_zmax
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+2,back_ymin-2:back_ymax+2,back_zmin-2:back_zmax+2) :: back_density0,back_energy0
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+2,back_ymin-2:back_ymax+2,back_zmin-2:back_zmax+2) :: back_pressure,back_viscosity,back_soundspeed
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+2,back_ymin-2:back_ymax+2,back_zmin-2:back_zmax+2) :: back_density1,back_energy1
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+3,back_ymin-2:back_ymax+3,back_zmin-2:back_zmax+3) :: back_xvel0,back_yvel0,back_zvel0
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+3,back_ymin-2:back_ymax+3,back_zmin-2:back_zmax+3) :: back_xvel1,back_yvel1,back_zvel1
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+3,back_ymin-2:back_ymax+2,back_zmin-2:back_zmax+2) :: back_vol_flux_x,back_mass_flux_x
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+2,back_ymin-2:back_ymax+3,back_zmin-2:back_zmax+2) :: back_vol_flux_y,back_mass_flux_y
        REAL(KIND=8), DIMENSION(back_xmin-2:back_xmax+2,back_ymin-2:back_ymax+2,back_zmin-2:back_zmax+3) :: back_vol_flux_z,back_mass_flux_z

        INTEGER :: fields(:),depth

        INTEGER :: j,k,l

!$OMP PARALLEL PRIVATE(k,j)

        ! Density 0
        IF(fields(FIELD_DENSITY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        density0(j,k,z_min-l)=back_density0(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Density 1
        IF(fields(FIELD_DENSITY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        density1(j,k,z_min-l)=back_density1(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 0
        IF(fields(FIELD_ENERGY0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        energy0(j,k,z_min-l)=back_energy0(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Energy 1
        IF(fields(FIELD_ENERGY1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        energy1(j,k,z_min-l)=back_energy1(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Pressure
        IF(fields(FIELD_PRESSURE).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        pressure(j,k,z_min-l)=back_pressure(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Viscosity
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        viscosity(j,k,z_min-l)=back_viscosity(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! Soundspeed
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        soundspeed(j,k,z_min-l)=back_soundspeed(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 0
        IF(fields(FIELD_XVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel0(j,k,z_min-l)=back_xvel0(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! XVEL 1
        IF(fields(FIELD_XVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        xvel1(j,k,z_min-l)=back_xvel1(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 0
        IF(fields(FIELD_YVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel0(j,k,z_min-l)=back_yvel0(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! YVEL 1
        IF(fields(FIELD_YVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        yvel1(j,k,z_min-l)=back_yvel1(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 0
        IF(fields(FIELD_ZVEL0).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel0(j,k,z_min-l)=back_zvel0(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! ZVEL 1
        IF(fields(FIELD_ZVEL1).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+1+depth
                        zvel1(j,k,z_min-l)=back_zvel1(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_X
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+1+depth
                        vol_flux_x(j,k,z_min-l)=back_vol_flux_x(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_X
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+1+depth
                        mass_flux_x(j,k,z_min-l)=back_mass_flux_x(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Y
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_y(j,k,z_min-l)=back_vol_flux_y(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Y
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+1+depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_y(j,k,z_min-l)=back_mass_flux_y(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! VOL_FLUX_Z
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        vol_flux_z(j,k,z_min-l)=back_vol_flux_z(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

        ! MASS_FLUX_Z
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN
            !$OMP DO COLLAPSE(2)
            DO l=1,depth
                DO k=y_min-depth, y_max+depth
                    DO j=x_min-depth, x_max+depth
                        mass_flux_z(j,k,z_min-l)=back_mass_flux_z(j,k,back_zmax+1-l)
                    ENDDO
                ENDDO
            ENDDO
            !$OMP END DO
        ENDIF

!$OMP END PARALLEL

    END SUBROUTINE update_tile_halo_back_kernel




END  MODULE update_tile_halo_kernel_module
