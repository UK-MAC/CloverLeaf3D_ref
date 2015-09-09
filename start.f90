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

!>  @brief Main set up routine
!>  @author Wayne Gaudin
!>  @details Invokes the mesh decomposer and sets up chunk connectivity. It then
!>  allocates the communication buffers and call the chunk initialisation and
!>  generation routines. It calls the equation of state to calculate initial
!>  pressure before priming the halo cells and writing an initial field summary.

SUBROUTINE start

  USE clover_module
  USE parse_module
  USE update_halo_module
  USE ideal_gas_module
  USE build_field_module

  IMPLICIT NONE

  INTEGER :: t

  INTEGER :: x_cells,y_cells,z_cells
  INTEGER :: right,left,top,bottom,back,front

  INTEGER :: fields(NUM_FIELDS) !, chunk_task_responsible_for 

  LOGICAL :: profiler_off

  IF(parallel%boss)THEN
     WRITE(g_out,*) 'Setting up initial geometry'
     WRITE(g_out,*)
  ENDIF

  time  = 0.0
  step  = 0
  dtold = dtinit
  dt    = dtinit

  CALL clover_barrier

  CALL clover_get_num_chunks(number_of_chunks)

  ALLOCATE(chunk%tiles(1:tiles_per_chunk))



  CALL clover_decompose(grid%x_cells,grid%y_cells,grid%z_cells,left,right,bottom,top,back,front)

      
    chunk%task = parallel%task

    !chunk_task_responsible_for = parallel%task+1

    x_cells = right -left  +1
    y_cells = top   -bottom+1
    z_cells = front -back  +1

    write(*,*) parallel%task, left, right, bottom, top, back, front

    CALL clover_decompose_tile(x_cells,y_cells,z_cells)
      
    DO t=1, tiles_per_chunk
      CALL build_field(t)
    END DO




  CALL clover_barrier

  CALL clover_allocate_buffers()

  IF(parallel%boss)THEN
     WRITE(g_out,*) 'Generating chunks'
  ENDIF
!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_chunk
      CALL initialise_chunk(t)
      CALL generate_chunk(t)
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  advect_x=.TRUE.

  CALL clover_barrier

  ! Do no profile the start up costs otherwise the total times will not add up
  ! at the end
  profiler_off=profiler_on
  profiler_on=.FALSE.

!$OMP PARALLEL
!$OMP DO
  DO t = 1, tiles_per_chunk
    CALL ideal_gas(t,.FALSE.)
  END DO
!$OMP END DO
!$OMP END PARALLEL

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY0)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_PRESSURE)=1
  fields(FIELD_VISCOSITY)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  fields(FIELD_ZVEL0)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1

  CALL update_halo(fields,2)

  IF(parallel%boss)THEN
     WRITE(g_out,*)
     WRITE(g_out,*) 'Problem initialised and generated'
  ENDIF

  CALL field_summary()

  IF(visit_frequency.NE.0) CALL visit()

  CALL clover_barrier

  profiler_on=profiler_off

END SUBROUTINE start
