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

!>  @brief Fortran ideal gas kernel.
!>  @author Wayne Gaudin
!>  @details Calculates the pressure and sound speed for the mesh chunk using
!>  the ideal gas equation of state, with a fixed gamma of 1.4.

MODULE ideal_gas_kernel_module

CONTAINS

SUBROUTINE ideal_gas_kernel(x_min,x_max,y_min,y_max,z_min,z_max,    &
                            density,                                &
                            energy,                                 &
                            pressure,                               &
                            soundspeed                              )                              

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: pressure
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: soundspeed

  INTEGER :: j,k,l

  REAL(KIND=8) :: sound_speed_squared,v,pressurebyenergy,pressurebyvolume

!$ACC DATA &
!$ACC PRESENT(density, energy, pressure, soundspeed)

!$ACC KERNELS

!$ACC LOOP INDEPENDENT
  DO l=z_min,z_max
!$ACC LOOP INDEPENDENT
    DO k=y_min,y_max
!$ACC LOOP INDEPENDENT PRIVATE(v,pressurebyenergy,pressurebyvolume,sound_speed_squared)
      DO j=x_min,x_max
        v=1.0_8/density(j,k,l)
        pressure(j,k,l)=(1.4_8-1.0_8)*density(j,k,l)*energy(j,k,l)
        pressurebyenergy=(1.4_8-1.0_8)*density(j,k,l)
        pressurebyvolume=-density(j,k,l)*pressure(j,k,l)
        sound_speed_squared=v*v*(pressure(j,k,l)*pressurebyenergy-pressurebyvolume)
        soundspeed(j,k,l)=SQRT(sound_speed_squared)
      ENDDO
    ENDDO
  ENDDO

!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE ideal_gas_kernel

END MODULE ideal_gas_kernel_module
