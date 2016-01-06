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

!>  @brief standalone driver for the acceleration kernels
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode


PROGRAM accelerate_driver

  USE set_data_module
  USE accelerate_kernel_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  INTEGER :: numargs,iargc,i
  CHARACTER (LEN=20)  :: command_line,temp

  INTEGER :: x_size,y_size,z_size

  REAL(KIND=8) :: kernel_time,timer,acceleration_time

  LOGICAL :: use_fortran_kernels
  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max,its,iteration
  REAL(KIND=8) :: dt
  REAL(KIND=8),ALLOCATABLE :: xarea(:,:,:),yarea(:,:,:),zarea(:,:,:),volume(:,:,:)
  REAL(KIND=8),ALLOCATABLE :: celldx(:),celldy(:),celldz(:)
  REAL(KIND=8),ALLOCATABLE :: density0(:,:,:),energy0(:,:,:),pressure(:,:,:),soundspeed(:,:,:),viscosity(:,:,:)
  REAL(KIND=8),ALLOCATABLE :: xvel0(:,:,:),yvel0(:,:,:),zvel0(:,:,:),xvel1(:,:,:),yvel1(:,:,:),zvel1(:,:,:),work_array1(:,:,:)

!$OMP PARALLEL
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
!$    WRITE(*,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
!$  ENDIF
!$OMP END PARALLEL

  x_size=100
  y_size=100
  z_size=100
  its=1
  use_fortran_kernels=.TRUE.

  numargs = iargc()

  DO i=1,numargs,2
    CALL GETARG(i,command_line)
    SELECT CASE (command_line)
      CASE("-help")
        WRITE(*,*) "Usage -nx 100 -ny 100 -nz 100 -its 10 -kernel fortran"
        stop
      CASE("-nx")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") x_size
      CASE("-ny")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") y_size
      CASE("-nz")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") z_size
      CASE("-its")
        CALL GETARG(i+1,temp)
        READ(UNIT=temp,FMT="(I20)") its
      CASE("-kernel")
        CALL GETARG(i+1,temp)
        IF(temp.EQ."fortran") THEN
          use_fortran_kernels=.TRUE.
        ENDIF
    END SELECT
  ENDDO

  x_min=1
  y_min=1
  z_min=1
  x_max=x_size
  y_max=y_size
  z_max=z_size

  WRITE(*,*) "Accelerate Kernel"
  WRITE(*,*) "Mesh size ",x_size,y_size,z_size
  WRITE(*,*) "Iterations ",its

  kernel_time=timer()

  CALL set_data(x_min,x_max,y_min,y_max,z_min,z_max, &
                xarea=xarea,                         &
                yarea=yarea,                         &
                zarea=zarea,                         &
                celldx=celldx,                       &
                celldy=celldy,                       &
                celldz=celldz,                       &
                volume=volume,                       &
                density0=density0,                   &
                energy0=energy0,                     &
                pressure=pressure,                   &
                soundspeed=soundspeed,               &
                viscosity=viscosity,                 &
                xvel0=xvel0,                         &
                xvel1=xvel1,                         &
                yvel0=yvel0,                         &
                yvel1=yvel1,                         &
                zvel0=zvel0,                         &
                zvel1=zvel1,                         &
                work_array1=work_array1,             &
                dt=dt                                )

  WRITE(*,*) "Setup time ",timer()-kernel_time

  WRITE(*,*) "Data initialised"

  IF(use_fortran_kernels) THEN
    WRITE(*,*) "Running Fortran kernel"
  ENDIF

  kernel_time=timer()

  IF(use_fortran_kernels) THEN
    DO iteration=1,its
      CALL accelerate_kernel(x_min,                  &
                             x_max,                  &
                             y_min,                  &
                             y_max,                  &
                             z_min,                  &
                             z_max,                  &
                             dt,                     &
                             xarea,                  &
                             yarea,                  &
                             zarea,                  &
                             volume,                 &
                             density0,               &
                             pressure,               &
                             viscosity,              &
                             xvel0,                  &
                             yvel0,                  &
                             zvel0,                  &
                             xvel1,                  &
                             yvel1,                  &
                             zvel1                   )
    ENDDO
  ENDIF

  acceleration_time=(timer()-kernel_time)

  WRITE(*,*) "Accelerate time ",acceleration_time 
  WRITE(*,*) "X vel ",SUM(xvel1)
  WRITE(*,*) "Y vel ",SUM(yvel1)
  WRITE(*,*) "Z vel ",SUM(zvel1)

  ! Answers need checking
  DEALLOCATE(xarea)
  DEALLOCATE(yarea)
  DEALLOCATE(zarea)
  DEALLOCATE(celldx)
  DEALLOCATE(celldy)
  DEALLOCATE(celldz)
  DEALLOCATE(volume)
  DEALLOCATE(density0)
  DEALLOCATE(energy0)
  DEALLOCATE(pressure)
  DEALLOCATE(soundspeed)
  DEALLOCATE(viscosity)
  DEALLOCATE(xvel0)
  DEALLOCATE(yvel0)
  DEALLOCATE(zvel0)
  DEALLOCATE(xvel1)
  DEALLOCATE(yvel1)
  DEALLOCATE(zvel1)
  DEALLOCATE(work_array1)

END PROGRAM accelerate_driver

