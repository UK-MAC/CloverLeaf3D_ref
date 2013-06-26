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

!>  @brief Communication Utilities
!>  @author Wayne Gaudin
!>  @details Contains all utilities required to run CloverLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no 
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications CloverLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE clover_module

  USE data_module
  USE definitions_module
  USE MPI

  IMPLICIT NONE

CONTAINS

SUBROUTINE clover_barrier

  INTEGER :: err

  CALL MPI_BARRIER(MPI_COMM_WORLD,err)

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

  INTEGER :: ierr,err

  CALL MPI_ABORT(MPI_COMM_WORLD,ierr,err)

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)
  CALL MPI_FINALIZE(err)

END SUBROUTINE clover_finalize

SUBROUTINE clover_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size

  rank=0
  size=1

  CALL MPI_INIT(err) 

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err) 
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err) 

  parallel%parallel=.TRUE.
  parallel%task=rank

  IF(rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=size

END SUBROUTINE clover_init_comms

SUBROUTINE clover_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE clover_get_num_chunks

SUBROUTINE clover_decompose(x_cells,y_cells,z_cells,left,right,bottom,top,back,front)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Doesn't always return the best split if there are few factors
  ! All factors need to be stored and the best picked. But its ok for now

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,z_cells,left(:),right(:),top(:),bottom(:),back(:),front(:)
  INTEGER :: c,delta_x,delta_y,delta_z

  REAL(KIND=8) :: mesh_ratio,factor_x,factor_y,factor_z
  INTEGER  :: chunk_x,chunk_y,chunk_z,mod_x,mod_y,split_found

  INTEGER  :: cx,cy,cz,chunk,add_x,add_y,add_x_prev,add_y_prev

  ! 2D Decomposition of the mesh

  mesh_ratio=real(x_cells)/real(y_cells)

  chunk_x=number_of_chunks
  chunk_y=1

  split_found=0 ! Used to detect 1D decomposition
  DO c=1,number_of_chunks
    IF (MOD(number_of_chunks,c).EQ.0) THEN
      factor_x=number_of_chunks/real(c)
      factor_y=c
      !Compare the factor ratio with the mesh ratio
      IF(factor_x/factor_y.LE.mesh_ratio) THEN
        chunk_y=c
        chunk_x=number_of_chunks/c
        split_found=1
        EXIT
      ENDIF
    ENDIF
  ENDDO

  IF(split_found.EQ.0.OR.chunk_y.EQ.number_of_chunks) THEN ! Prime number or 1D decomp detected
    IF(mesh_ratio.GE.1.0) THEN
      chunk_x=number_of_chunks
      chunk_y=1
    ELSE
      chunk_x=1
      chunk_y=number_of_chunks
    ENDIF
  ENDIF

  chunk_z=1

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  ! Set up chunk mesh ranges and chunk connectivity

  back(1)=1
  front(1)=z_cells
  add_x_prev=0
  add_y_prev=0
  chunk=1
  DO cz=1,chunk_z
    DO cy=1,chunk_y
      DO cx=1,chunk_x
        add_x=0
        add_y=0
        IF(cx.LE.mod_x)add_x=1
        IF(cy.LE.mod_y)add_y=1
        left(chunk)=(cx-1)*delta_x+1+add_x_prev
        right(chunk)=left(chunk)+delta_x-1+add_x
        bottom(chunk)=(cy-1)*delta_y+1+add_y_prev
        top(chunk)=bottom(chunk)+delta_y-1+add_y
        chunks(chunk)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
        chunks(chunk)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
        chunks(chunk)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
        chunks(chunk)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx
        IF(cx.EQ.1)chunks(chunk)%chunk_neighbours(chunk_left)=external_face
        IF(cx.EQ.chunk_x)chunks(chunk)%chunk_neighbours(chunk_right)=external_face
        IF(cy.EQ.1)chunks(chunk)%chunk_neighbours(chunk_bottom)=external_face
        IF(cy.EQ.chunk_y)chunks(chunk)%chunk_neighbours(chunk_top)=external_face
        IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
        chunks(chunk)%chunk_neighbours(chunk_back)=external_face
        chunks(chunk)%chunk_neighbours(chunk_front)=external_face
        chunk=chunk+1
      ENDDO
      add_x_prev=0
      IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
    ENDDO
  ENDDO

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",mesh_ratio
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," by ",chunk_z," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE clover_decompose

SUBROUTINE clover_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk
  
  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN
    !IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%left_snd_buffer(10*2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%left_rcv_buffer(10*2*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%right_snd_buffer(10*2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%right_rcv_buffer(10*2*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%bottom_snd_buffer(10*2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%bottom_rcv_buffer(10*2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%top_snd_buffer(10*2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%top_rcv_buffer(10*2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_back).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%back_snd_buffer(10*2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%back_rcv_buffer(10*2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_front).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%front_snd_buffer(10*2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%front_rcv_buffer(10*2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
  ENDIF

END SUBROUTINE clover_allocate_buffers

SUBROUTINE clover_exchange(fields,depth)

    IMPLICIT NONE

    INTEGER      :: fields(:),depth, chunk
    INTEGER      :: left_right_offset(15),bottom_top_offset(15)
    INTEGER      :: request(4)
    INTEGER      :: message_count,err
    INTEGER      :: status(MPI_STATUS_SIZE,4)
    INTEGER      :: end_pack_index_left_right, end_pack_index_bottom_top,field

    ! Assuming 1 patch per task, this will be changed

    request=0
    message_count=0

    chunk = parallel%task+1 

    end_pack_index_left_right=0
    end_pack_index_bottom_top=0
    DO field=1,15
      IF(fields(field).EQ.1) THEN
        left_right_offset(field)=end_pack_index_left_right
        bottom_top_offset(field)=end_pack_index_bottom_top
        end_pack_index_left_right=end_pack_index_left_right+depth*(chunks(chunk)%field%y_max+5)
        end_pack_index_bottom_top=end_pack_index_bottom_top+depth*(chunks(chunk)%field%x_max+5)
      ENDIF
    ENDDO

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ! do left exchanges
      CALL clover_pack_left(chunk, fields, depth, left_right_offset)

      !send and recv messagse to the left
      CALL clover_send_recv_message_left(chunks(chunk)%left_snd_buffer,                      &
                                         chunks(chunk)%left_rcv_buffer,                      &
                                         chunk,end_pack_index_left_right,                    &
                                         1, 2,                                               &
                                         request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ! do right exchanges
      CALL clover_pack_right(chunk, fields, depth, left_right_offset)

      !send message to the right
      CALL clover_send_recv_message_right(chunks(chunk)%right_snd_buffer,                     &
                                          chunks(chunk)%right_rcv_buffer,                     &
                                          chunk,end_pack_index_left_right,                    &
                                          2, 1,                                               &
                                          request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    !make a call to wait / sync
    CALL MPI_WAITALL(message_count,request,status,err)

    !unpack in left direction
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      CALL clover_unpack_left(fields, chunk, depth,                      &
                              chunks(chunk)%left_rcv_buffer,             &
                              left_right_offset)                  
    ENDIF


    !unpack in right direction
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      CALL clover_unpack_right(fields, chunk, depth,                     &
                               chunks(chunk)%right_rcv_buffer,           &
                               left_right_offset)
    ENDIF

    message_count = 0
    request = 0

    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ! do bottom exchanges
      CALL clover_pack_bottom(chunk, fields, depth, bottom_top_offset)

      !send message downwards
      CALL clover_send_recv_message_bottom(chunks(chunk)%bottom_snd_buffer,                     &
                                           chunks(chunk)%bottom_rcv_buffer,                     &
                                           chunk,end_pack_index_bottom_top,                     &
                                           3, 4,                                                &
                                           request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ! do top exchanges
      CALL clover_pack_top(chunk, fields, depth, bottom_top_offset)

      !send message upwards
      CALL clover_send_recv_message_top(chunks(chunk)%top_snd_buffer,                           &
                                        chunks(chunk)%top_rcv_buffer,                           &
                                        chunk,end_pack_index_bottom_top,                        &
                                        4, 3,                                                   &
                                        request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    !need to make a call to wait / sync
    CALL MPI_WAITALL(message_count,request,status,err)

    !unpack in top direction
    IF( chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face ) THEN
      CALL clover_unpack_top(fields, chunk, depth,                       &
                             chunks(chunk)%top_rcv_buffer,               &
                             bottom_top_offset)
    ENDIF

    !unpack in bottom direction
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      CALL clover_unpack_bottom(fields, chunk, depth,                   &
                               chunks(chunk)%bottom_rcv_buffer,         &
                               bottom_top_offset)
    ENDIF

END SUBROUTINE clover_exchange

SUBROUTINE clover_pack_left(chunk, fields, depth, left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER      :: fields(:),depth, chunk
  INTEGER      :: left_right_offset(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%z_min,                    &
                                    chunks(chunk)%field%z_max,                    &
                                    chunks(chunk)%field%density0,                 &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density0,                 &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%z_min,                    &
                                    chunks(chunk)%field%z_max,                    &
                                    chunks(chunk)%field%density1,                 &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density1,                 &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%energy0,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy0,                  &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%energy1,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy1,                  &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%pressure,                 &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%pressure,                 &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%viscosity,                &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%viscosity,                &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%soundspeed,               &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%soundspeed,               &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%xvel0,                    &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, VERTEX_DATA,                           &
                                    left_right_offset(FIELD_XVEL0))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%xvel0,                    &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%xvel1,                    &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, VERTEX_DATA,                           &
                                    left_right_offset(FIELD_XVEL1))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%xvel1,                    &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%yvel0,                    &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, VERTEX_DATA,                           &
                                    left_right_offset(FIELD_YVEL0))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%yvel0,                    &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%yvel1,                    &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, VERTEX_DATA,                           &
                                    left_right_offset(FIELD_YVEL1))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%yvel1,                    &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%vol_flux_x,               &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, X_FACE_DATA,                           &
                                    left_right_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vol_flux_x,               &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, X_FACE_DATA,                           &
                                      left_right_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%vol_flux_y,               &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, Y_FACE_DATA,                           &
                                    left_right_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vol_flux_y,               &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, Y_FACE_DATA,                           &
                                      left_right_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%mass_flux_x,              &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, X_FACE_DATA,                           &
                                    left_right_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%mass_flux_x,              &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, X_FACE_DATA,                           &
                                      left_right_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%mass_flux_y,              &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, Y_FACE_DATA,                           &
                                    left_right_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_pack_message_left_c(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%mass_flux_y,              &
                                      chunks(chunk)%left_snd_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, Y_FACE_DATA,                           &
                                      left_right_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_pack_left

SUBROUTINE clover_send_recv_message_left(left_snd_buffer, left_rcv_buffer,      &
                                         chunk, total_size,                     &
                                         tag_send, tag_recv,                    &
                                         req_send, req_recv)

  REAL(KIND=8)    :: left_snd_buffer(:), left_rcv_buffer(:)
  INTEGER         :: left_task
  INTEGER         :: chunk
  INTEGER         :: total_size, tag_send, tag_recv, err
  INTEGER         :: req_send, req_recv

  left_task =chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task

  CALL MPI_ISEND(left_snd_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_send &
                ,MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(left_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_recv &
                ,MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE clover_send_recv_message_left

SUBROUTINE clover_unpack_left(fields, chunk, depth,                         &
                              left_rcv_buffer,                              &
                              left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, depth
  INTEGER         :: left_right_offset(:)
  REAL(KIND=8)    :: left_rcv_buffer(:)


  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density0,                 &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density0,                 &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density1,                 &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density1,                 &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy0,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy0,                  &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy1,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy1,                  &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%pressure,                 &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%pressure,                 &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%viscosity,                &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%viscosity,                &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%soundspeed,               &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%soundspeed,               &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        left_right_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%xvel0,                    &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_XVEL0))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%xvel0,                    &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        left_right_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%xvel1,                    &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_XVEL1))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%xvel1,                    &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        left_right_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%yvel0,                    &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_YVEL0))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%yvel0,                    &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        left_right_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%yvel1,                    &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      left_right_offset(FIELD_YVEL1))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%yvel1,                    &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        left_right_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vol_flux_x,               &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, X_FACE_DATA,                           &
                                      left_right_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vol_flux_x,               &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, X_FACE_DATA,                           &
                                        left_right_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vol_flux_y,               &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, Y_FACE_DATA,                           &
                                      left_right_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vol_flux_y,               &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, Y_FACE_DATA,                           &
                                        left_right_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%mass_flux_x,              &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, X_FACE_DATA,                           &
                                      left_right_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%mass_flux_x,              &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, X_FACE_DATA,                           &
                                        left_right_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%mass_flux_y,              &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, Y_FACE_DATA,                           &
                                      left_right_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_unpack_message_left_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%mass_flux_y,              &
                                        chunks(chunk)%left_rcv_buffer,                &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, Y_FACE_DATA,                           &
                                        left_right_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_unpack_left

SUBROUTINE clover_pack_right(chunk, fields, depth, left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER        :: chunk, fields(:), depth, tot_packr, left_right_offset(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density0,                 &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density0,                 &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density1,                 &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density1,                 &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy0,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy0,                  &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy1,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy1,                  &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%pressure,                 &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%pressure,                 &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%viscosity,                &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%viscosity,                &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%soundspeed,               &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%soundspeed,               &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%xvel0,                    &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     left_right_offset(FIELD_XVEL0))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%xvel0,                    &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%xvel1,                    &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     left_right_offset(FIELD_XVEL1))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%xvel1,                    &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%yvel0,                    &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     left_right_offset(FIELD_YVEL0))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%yvel0,                    &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%yvel1,                    &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     left_right_offset(FIELD_YVEL1))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%yvel1,                    &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vol_flux_x,               &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, X_FACE_DATA,                           &
                                     left_right_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vol_flux_x,               &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, X_FACE_DATA,                           &
                                       left_right_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vol_flux_y,               &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, Y_FACE_DATA,                           &
                                     left_right_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vol_flux_y,               &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, Y_FACE_DATA,                           &
                                       left_right_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%mass_flux_x,              &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, X_FACE_DATA,                           &
                                     left_right_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%mass_flux_x,              &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, X_FACE_DATA,                           &
                                       left_right_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%mass_flux_y,              &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, Y_FACE_DATA,                           &
                                     left_right_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_pack_message_right_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%mass_flux_y,              &
                                       chunks(chunk)%right_snd_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, Y_FACE_DATA,                           &
                                       left_right_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_pack_right

SUBROUTINE clover_send_recv_message_right(right_snd_buffer, right_rcv_buffer,   &
                                          chunk, total_size,                    &
                                          tag_send, tag_recv,                   &
                                          req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: right_snd_buffer(:), right_rcv_buffer(:)
  INTEGER      :: right_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  right_task=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task

  CALL MPI_ISEND(right_snd_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_send, &
                 MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(right_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_recv, &
                 MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE clover_send_recv_message_right

SUBROUTINE clover_unpack_right(fields, chunk, depth,                          &
                               right_rcv_buffer,                              &
                               left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, total_in_right_buff, depth, left_right_offset(:)
  REAL(KIND=8)    :: right_rcv_buffer(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density0,                 &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%density0,                 &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density1,                 &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%density1,                 &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy0,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%energy0,                  &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy1,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%energy1,                  &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%pressure,                 &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%pressure,                 &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%viscosity,                &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%viscosity,                &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%soundspeed,               &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%soundspeed,               &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, CELL_DATA,                             &
                                         left_right_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%xvel0,                    &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_XVEL0))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%xvel0,                    &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, VERTEX_DATA,                           &
                                         left_right_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%xvel1,                    &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_XVEL1))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%xvel1,                    &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, VERTEX_DATA,                           &
                                         left_right_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%yvel0,                    &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_YVEL0))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%yvel0,                    &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, VERTEX_DATA,                           &
                                         left_right_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%yvel1,                    &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       left_right_offset(FIELD_YVEL1))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%yvel1,                    &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, VERTEX_DATA,                           &
                                         left_right_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vol_flux_x,               &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, X_FACE_DATA,                           &
                                       left_right_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%vol_flux_x,               &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, X_FACE_DATA,                           &
                                         left_right_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vol_flux_y,               &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, Y_FACE_DATA,                           &
                                       left_right_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%vol_flux_y,               &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, Y_FACE_DATA,                           &
                                         left_right_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%mass_flux_x,              &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, X_FACE_DATA,                           &
                                       left_right_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%mass_flux_x,              &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, X_FACE_DATA,                           &
                                         left_right_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%mass_flux_y,              &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, Y_FACE_DATA,                           &
                                       left_right_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_unpack_message_right_c(chunks(chunk)%field%x_min,                    &
                                         chunks(chunk)%field%x_max,                    &
                                         chunks(chunk)%field%y_min,                    &
                                         chunks(chunk)%field%y_max,                    &
                                         chunks(chunk)%field%mass_flux_y,              &
                                         chunks(chunk)%right_rcv_buffer,               &
                                         CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                         depth, Y_FACE_DATA,                           &
                                         left_right_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_unpack_right

SUBROUTINE clover_pack_top(chunk, fields, depth, bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER        :: chunk, fields(:), depth, bottom_top_offset(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%density0,                 &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density0,                 &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%density1,                 &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density1,                 &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%energy0,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy0,                  &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%energy1,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy1,                  &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%pressure,                 &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%pressure,                 &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%viscosity,                &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%viscosity,                &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%soundspeed,               &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%soundspeed,               &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%xvel0,                    &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, VERTEX_DATA,                           &
                                   bottom_top_offset(FIELD_XVEL0))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%xvel0,                    &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%xvel1,                    &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, VERTEX_DATA,                           &
                                   bottom_top_offset(FIELD_XVEL1))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%xvel1,                    &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%yvel0,                    &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, VERTEX_DATA,                           &
                                   bottom_top_offset(FIELD_YVEL0))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%yvel0,                    &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%yvel1,                    &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, VERTEX_DATA,                           &
                                   bottom_top_offset(FIELD_YVEL1))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%yvel1,                    &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%vol_flux_x,               &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, X_FACE_DATA,                           &
                                   bottom_top_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vol_flux_x,               &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, X_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%vol_flux_y,               &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, Y_FACE_DATA,                           &
                                   bottom_top_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vol_flux_y,               &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, Y_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%mass_flux_x,              &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, X_FACE_DATA,                           &
                                   bottom_top_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%mass_flux_x,              &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, X_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%mass_flux_y,              &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, Y_FACE_DATA,                           &
                                   bottom_top_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_pack_message_top_c(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%mass_flux_y,              &
                                     chunks(chunk)%top_snd_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, Y_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_pack_top

SUBROUTINE clover_send_recv_message_top(top_snd_buffer, top_rcv_buffer,     &
                                        chunk, total_size,                  &
                                        tag_send, tag_recv,                 &
                                        req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: top_snd_buffer(:), top_rcv_buffer(:)
  INTEGER      :: top_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

    top_task=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task

  CALL MPI_ISEND(top_snd_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_send, &
                 MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(top_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_recv, &
                 MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE clover_send_recv_message_top

SUBROUTINE clover_unpack_top(fields, chunk, depth,                        &
                             top_rcv_buffer,                              &
                             bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, total_in_top_buff, depth, bottom_top_offset(:)
  REAL(KIND=8)    :: top_rcv_buffer(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density0,                 &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density0,                 &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density1,                 &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density1,                 &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy0,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy0,                  &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy1,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy1,                  &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%pressure,                 &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%pressure,                 &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%viscosity,                &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%viscosity,                &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%soundspeed,               &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%soundspeed,               &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       bottom_top_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%xvel0,                    &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_XVEL0))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%xvel0,                    &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       bottom_top_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%xvel1,                    &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_XVEL1))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%xvel1,                    &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       bottom_top_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%yvel0,                    &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_YVEL0))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%yvel0,                    &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       bottom_top_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%yvel1,                    &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, VERTEX_DATA,                           &
                                     bottom_top_offset(FIELD_YVEL1))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%yvel1,                    &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, VERTEX_DATA,                           &
                                       bottom_top_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vol_flux_x,               &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, X_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vol_flux_x,               &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, X_FACE_DATA,                           &
                                       bottom_top_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vol_flux_y,               &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, Y_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vol_flux_y,               &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, Y_FACE_DATA,                           &
                                       bottom_top_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%mass_flux_x,              &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, X_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%mass_flux_x,              &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, X_FACE_DATA,                           &
                                       bottom_top_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%mass_flux_y,              &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, Y_FACE_DATA,                           &
                                     bottom_top_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_unpack_message_top_c(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%mass_flux_y,              &
                                       chunks(chunk)%top_rcv_buffer,                 &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, Y_FACE_DATA,                           &
                                       bottom_top_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_unpack_top

SUBROUTINE clover_pack_bottom(chunk, fields, depth, bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER        :: chunk, fields(:), depth, tot_packb, bottom_top_offset(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density0,                 &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density0,                 &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density1,                 &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density1,                 &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy0,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy0,                  &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy1,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy1,                  &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%pressure,                 &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%pressure,                 &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%viscosity,                &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%viscosity,                &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%soundspeed,               &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%soundspeed,               &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%xvel0,                    &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      bottom_top_offset(FIELD_XVEL0))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%xvel0     ,               &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%xvel1,                    &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      bottom_top_offset(FIELD_XVEL1))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%xvel1,                    &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%yvel0,                    &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      bottom_top_offset(FIELD_YVEL0))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%yvel0,                    &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%yvel1,                    &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, VERTEX_DATA,                           &
                                      bottom_top_offset(FIELD_YVEL1))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%yvel1,                    &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vol_flux_x,               &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, X_FACE_DATA,                           &
                                      bottom_top_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vol_flux_x,               &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, X_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vol_flux_y,               &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, Y_FACE_DATA,                           &
                                      bottom_top_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vol_flux_y,               &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, Y_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%mass_flux_x,              &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, X_FACE_DATA,                           &
                                      bottom_top_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%mass_flux_x,              &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, X_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%mass_flux_y,              &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, Y_FACE_DATA,                           &
                                      bottom_top_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_pack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%mass_flux_y,              &
                                        chunks(chunk)%bottom_snd_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, Y_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_pack_bottom

SUBROUTINE clover_send_recv_message_bottom(bottom_snd_buffer, bottom_rcv_buffer,        &
                                           chunk, total_size,                           &   
                                           tag_send, tag_recv,                          &
                                           req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: bottom_snd_buffer(:), bottom_rcv_buffer(:)
  INTEGER      :: bottom_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  bottom_task=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task

  CALL MPI_ISEND(bottom_snd_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_send &
                ,MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(bottom_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_recv &
                ,MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE clover_send_recv_message_bottom

SUBROUTINE clover_unpack_bottom(fields, chunk, depth,                        &
                             bottom_rcv_buffer,                              &
                             bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, depth, bottom_top_offset(:)
  REAL(KIND=8)    :: bottom_rcv_buffer(:)

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density0,                 &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_DENSITY0))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%density0,                 &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_DENSITY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density1,                 &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_DENSITY1))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%density1,                 &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_DENSITY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy0,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_ENERGY0))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%energy0,                  &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy1,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_ENERGY1))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%energy1,                  &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%pressure,                 &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_PRESSURE))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%pressure,                 &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_PRESSURE))
    ENDIF
  ENDIF
  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%viscosity,                &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_VISCOSITY))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%viscosity,                &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_VISCOSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%soundspeed,               &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_SOUNDSPEED))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%soundspeed,               &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, CELL_DATA,                             &
                                          bottom_top_offset(FIELD_SOUNDSPEED))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%xvel0,                    &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_XVEL0))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%xvel0,                    &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, VERTEX_DATA,                           &
                                          bottom_top_offset(FIELD_XVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%xvel1,                    &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_XVEL1))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%xvel1,                    &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, VERTEX_DATA,                           &
                                          bottom_top_offset(FIELD_XVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%yvel0,                    &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_YVEL0))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%yvel0,                    &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, VERTEX_DATA,                           &
                                          bottom_top_offset(FIELD_YVEL0))
    ENDIF
  ENDIF
  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%yvel1,                    &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, VERTEX_DATA,                           &
                                        bottom_top_offset(FIELD_YVEL1))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%yvel1,                    &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, VERTEX_DATA,                           &
                                          bottom_top_offset(FIELD_YVEL1))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vol_flux_x,               &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, X_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_VOL_FLUX_X))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%vol_flux_x,               &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, X_FACE_DATA,                           &
                                          bottom_top_offset(FIELD_VOL_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vol_flux_y,               &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, Y_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_VOL_FLUX_Y))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%vol_flux_y,               &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, Y_FACE_DATA,                           &
                                          bottom_top_offset(FIELD_VOL_FLUX_Y))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%mass_flux_x,              &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, X_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_MASS_FLUX_X))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%mass_flux_x,              &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, X_FACE_DATA,                           &
                                          bottom_top_offset(FIELD_MASS_FLUX_X))
    ENDIF
  ENDIF
  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL clover_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%mass_flux_y,              &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, Y_FACE_DATA,                           &
                                        bottom_top_offset(FIELD_MASS_FLUX_Y))
    ELSE
      CALL clover_unpack_message_bottom_c(chunks(chunk)%field%x_min,                    &
                                          chunks(chunk)%field%x_max,                    &
                                          chunks(chunk)%field%y_min,                    &
                                          chunks(chunk)%field%y_max,                    &
                                          chunks(chunk)%field%mass_flux_y,              &
                                          chunks(chunk)%bottom_rcv_buffer,              &
                                          CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                          depth, Y_FACE_DATA,                           &
                                          bottom_top_offset(FIELD_MASS_FLUX_Y))
    ENDIF
  ENDIF

END SUBROUTINE clover_unpack_bottom

SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_REDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)

  value=total

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  INTEGER :: err

  minimum=value

  CALL MPI_ALLREDUCE(value,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,err)

  value=minimum

END SUBROUTINE clover_min

SUBROUTINE clover_max(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: maximum

  INTEGER :: err

  maximum=value

  CALL MPI_ALLREDUCE(value,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)

  value=maximum

END SUBROUTINE clover_max

SUBROUTINE clover_allgather(value,values)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: values(parallel%max_task)

  INTEGER :: err

  values(1)=value ! Just to ensure it will work in serial

  CALL MPI_ALLGATHER(value,1,MPI_DOUBLE_PRECISION,values,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,err)

END SUBROUTINE clover_allgather

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  INTEGER :: err

  maximum=error

  CALL MPI_ALLREDUCE(error,maximum,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)

  error=maximum

END SUBROUTINE clover_check_error


END MODULE clover_module
