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
!>  @author Wayne Gaudin, Ollie Perks
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
        ! Picks split with minimal surface area to volume ratio

        IMPLICIT NONE

        INTEGER :: x_cells,y_cells,z_cells,left,right,top,bottom,back,front
        INTEGER :: c,delta_x,delta_y,delta_z

        REAL(KIND=8) :: surface,volume,best_metric,current_metric
        INTEGER  :: chunk_x,chunk_y,chunk_z,mod_x,mod_y,mod_z

        INTEGER  :: cx,cy,cz,chunk_id,add_x,add_y,add_z,add_x_prev,add_y_prev,add_z_prev
        INTEGER  :: div1,j,current_x,current_y,current_z

        ! 3D Decomposition of the mesh

        current_x = 1
        current_y = 1
        current_z = number_of_chunks

        ! Initialise metric
        surface = (((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*2) &
            + (((1.0*x_cells)/current_x)*((1.0*z_cells)/current_z)*2) &
            + (((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)*2)
        volume  = ((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)
        best_metric = surface/volume
        chunk_x=current_x
        chunk_y=current_y
        chunk_z=current_z

        DO c=1,number_of_chunks

            ! If doesn't evenly divide loop
            IF(MOD(number_of_chunks,c).NE.0) CYCLE

            current_x=c

            div1 = number_of_chunks/c

            DO j=1,div1
                IF(MOD(div1,j).NE.0) CYCLE
                current_y = j

                IF(MOD(number_of_chunks,(c*j)).NE.0) CYCLE

                current_z = number_of_chunks/(c*j)

                surface = (((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*2) &
                    + (((1.0*x_cells)/current_x)*((1.0*z_cells)/current_z)*2) &
                    + (((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)*2)
                volume  = ((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)

                current_metric = surface/volume

                IF(current_metric < best_metric) THEN
                    chunk_x=current_x
                    chunk_y=current_y
                    chunk_z=current_z
                    best_metric=current_metric
                ENDIF

            ENDDO

        ENDDO

        ! Set up chunk mesh ranges and chunk connectivity

        delta_x=x_cells/chunk_x
        delta_y=y_cells/chunk_y
        delta_z=z_cells/chunk_z
        mod_x=MOD(x_cells,chunk_x)
        mod_y=MOD(y_cells,chunk_y)
        mod_z=MOD(z_cells,chunk_z)
        add_x_prev=0
        add_y_prev=0
        add_z_prev=0
        chunk_id=1
        DO cz=1,chunk_z
            DO cy=1,chunk_y
                DO cx=1,chunk_x
                    add_x=0
                    add_y=0
                    add_z=0
                    IF(cx.LE.mod_x)add_x=1
                    IF(cy.LE.mod_y)add_y=1
                    IF(cz.LE.mod_z)add_z=1
                    IF (chunk_id .EQ. parallel%task+1) THEN
                        ! Mesh chunks
                        left=(cx-1)*delta_x+1+add_x_prev
                        right=left+delta_x-1+add_x
                        bottom=(cy-1)*delta_y+1+add_y_prev
                        top=bottom+delta_y-1+add_y
                        back=(cz-1)*delta_z+1+add_z_prev
                        front=back+delta_z-1+add_z
                        ! Chunk connectivity
                        chunk%chunk_neighbours(chunk_left)=  chunk_id-1
                        chunk%chunk_neighbours(chunk_right)= chunk_id+1
                        chunk%chunk_neighbours(chunk_bottom)=chunk_id-chunk_x
                        chunk%chunk_neighbours(chunk_top)=   chunk_id+chunk_x
                        chunk%chunk_neighbours(chunk_back)=  chunk_id-chunk_x*chunk_y
                        chunk%chunk_neighbours(chunk_front)= chunk_id+chunk_x*chunk_y
                        IF(cx.EQ.1)chunk%chunk_neighbours(chunk_left)=external_face
                        IF(cx.EQ.chunk_x)chunk%chunk_neighbours(chunk_right)=external_face
                        IF(cy.EQ.1)chunk%chunk_neighbours(chunk_bottom)=external_face
                        IF(cy.EQ.chunk_y)chunk%chunk_neighbours(chunk_top)=external_face
                        IF(cz.EQ.1)chunk%chunk_neighbours(chunk_back)=external_face
                        IF(cz.EQ.chunk_z)chunk%chunk_neighbours(chunk_front)=external_face

                        chunk%left    = left
                        chunk%bottom  = bottom
                        chunk%right   = right
                        chunk%top     = top
                        chunk%back    = back
                        chunk%front   = front
                        chunk%left_boundary   = 1
                        chunk%bottom_boundary = 1
                        chunk%back_boundary   = 1
                        chunk%right_boundary  = grid%x_cells
                        chunk%top_boundary    = grid%y_cells
                        chunk%front_boundary  = grid%z_cells
                        chunk%x_min = 1
                        chunk%y_min = 1
                        chunk%z_min = 1
                        chunk%x_max = right-left+1
                        chunk%y_max = top-bottom+1
                        chunk%z_max = front-back+1

                    ENDIF
                    IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
                    chunk_id=chunk_id+1
                ENDDO
                add_x_prev=0
                IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
            ENDDO
            add_x_prev=0
            add_y_prev=0
            IF(cz.LE.mod_z)add_z_prev=add_z_prev+1
        ENDDO

        IF(parallel%boss)THEN
            WRITE(g_out,*)
            WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," by ",chunk_z," chunks"
            WRITE(g_out,*)
        ENDIF

    END SUBROUTINE clover_decompose

    SUBROUTINE clover_decompose_tile(x_cells,y_cells,z_cells)

        ! This decomposes the mesh into a number of chunks.
        ! The number of chunks may be a multiple of the number of mpi tasks
        ! Picks split with minimal surface area to volume ratio

        IMPLICIT NONE

        INTEGER :: x_cells,y_cells,z_cells,left,right,top,bottom,back,front
        INTEGER :: t,delta_x,delta_y,delta_z

        REAL(KIND=8) :: surface,volume,best_metric,current_metric
        INTEGER  :: chunk_x,chunk_y,chunk_z,mod_x,mod_y,mod_z

        INTEGER  :: cx,cy,cz,tile_id,add_x,add_y,add_z,add_x_prev,add_y_prev,add_z_prev
        INTEGER  :: div1,j,current_x,current_y,current_z

        INTEGER :: x_bound, y_bound, z_bound ! Bounds to force 1D and 2D depcompositions

        ! 3D Decomposition of the mesh

        current_x = 1
        current_y = 1
        current_z = tiles_per_chunk

        ! Initialise metric
        surface = (((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*2) &
            + (((1.0*x_cells)/current_x)*((1.0*z_cells)/current_z)*2) &
            + (((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)*2)
        volume  = ((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)
        best_metric = surface/volume
        chunk_x=current_x
        chunk_y=current_y
        chunk_z=current_z

        ! Set up bounds for decompositions

        IF(tile_3d)THEN
          x_bound = tiles_per_chunk
          y_bound = tiles_per_chunk
        ENDIF

        IF(tile_2d)THEN
          x_bound = 1
          y_bound = tiles_per_chunk
        ENDIF

        IF(tile_1d)THEN
          x_bound = 1
          y_bound = 1
        ENDIF


        DO t=1,x_bound

            ! If doesn't evenly divide loop
            IF(MOD(tiles_per_chunk,t).NE.0) CYCLE

            current_x=t

            div1 = tiles_per_chunk/t

            div1 = MIN(y_bound, div1)

            DO j=1,div1
                IF(MOD(div1,j).NE.0) CYCLE
                current_y = j

                IF(MOD(tiles_per_chunk,(t*j)).NE.0) CYCLE

                current_z = tiles_per_chunk/(t*j)

                surface = (((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*2) &
                    + (((1.0*x_cells)/current_x)*((1.0*z_cells)/current_z)*2) &
                    + (((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)*2)
                volume  = ((1.0*x_cells)/current_x)*((1.0*y_cells)/current_y)*((1.0*z_cells)/current_z)

                current_metric = surface/volume

                IF(current_metric < best_metric) THEN
                    chunk_x=current_x
                    chunk_y=current_y
                    chunk_z=current_z
                    best_metric=current_metric
                ENDIF

            ENDDO

        ENDDO

        ! Set up chunk mesh ranges and chunk connectivity

        delta_x=x_cells/chunk_x
        delta_y=y_cells/chunk_y
        delta_z=z_cells/chunk_z
        mod_x=MOD(x_cells,chunk_x)
        mod_y=MOD(y_cells,chunk_y)
        mod_z=MOD(z_cells,chunk_z)
        add_x_prev=0
        add_y_prev=0
        add_z_prev=0
        tile_id=1
        DO cz=1,chunk_z
            DO cy=1,chunk_y
                DO cx=1,chunk_x
                    add_x=0
                    add_y=0
                    add_z=0
                    IF(cx.LE.mod_x)add_x=1
                    IF(cy.LE.mod_y)add_y=1
                    IF(cz.LE.mod_z)add_z=1

                    ! Mesh chunks
                    left=(chunk%left-1) + (cx-1)*delta_x+1+add_x_prev
                    right=left+delta_x-1+add_x
                    bottom=(chunk%bottom-1) + (cy-1)*delta_y+1+add_y_prev
                    top=bottom+delta_y-1+add_y
                    back=(chunk%back-1) + (cz-1)*delta_z+1+add_z_prev
                    front=back+delta_z-1+add_z
                    ! Chunk connectivity
                    chunk%tiles(tile_id)%tile_neighbours(tile_left)=  tile_id-1
                    chunk%tiles(tile_id)%tile_neighbours(tile_right)= tile_id+1
                    chunk%tiles(tile_id)%tile_neighbours(tile_bottom)=tile_id-chunk_x
                    chunk%tiles(tile_id)%tile_neighbours(tile_top)=   tile_id+chunk_x
                    chunk%tiles(tile_id)%tile_neighbours(tile_back)=  tile_id-chunk_x*chunk_y
                    chunk%tiles(tile_id)%tile_neighbours(tile_front)= tile_id+chunk_x*chunk_y

                    chunk%tiles(tile_id)%external_tile_mask=0
                    IF(cx.EQ.1)THEN
                        chunk%tiles(tile_id)%tile_neighbours(tile_left)=external_tile
                        chunk%tiles(tile_id)%external_tile_mask(tile_left)=1
                    ENDIF
                    IF(cx.EQ.chunk_x)THEN
                        chunk%tiles(tile_id)%tile_neighbours(tile_right)=external_tile
                        chunk%tiles(tile_id)%external_tile_mask(tile_right)=1
                    ENDIF
                    IF(cy.EQ.1)THEN
                        chunk%tiles(tile_id)%tile_neighbours(tile_bottom)=external_tile
                        chunk%tiles(tile_id)%external_tile_mask(tile_bottom)=1
                    ENDIF
                    IF(cy.EQ.chunk_y)THEN
                        chunk%tiles(tile_id)%tile_neighbours(tile_top)=external_tile
                        chunk%tiles(tile_id)%external_tile_mask(tile_top)=1
                    ENDIF
                    IF(cz.EQ.1)THEN
                        chunk%tiles(tile_id)%tile_neighbours(tile_back)=external_tile
                        chunk%tiles(tile_id)%external_tile_mask(tile_back)=1
                    ENDIF
                    IF(cz.EQ.chunk_z)THEN
                        chunk%tiles(tile_id)%tile_neighbours(tile_front)=external_tile
                        chunk%tiles(tile_id)%external_tile_mask(tile_front)=1
                    ENDIF


                    ! Save local tile boundaries
                    chunk%tiles(tile_id)%t_xmin=1
                    chunk%tiles(tile_id)%t_xmax=right - left + 1
                    chunk%tiles(tile_id)%t_ymin=1
                    chunk%tiles(tile_id)%t_ymax=top - bottom + 1
                    chunk%tiles(tile_id)%t_zmin=1
                    chunk%tiles(tile_id)%t_zmax=front - back + 1

                    ! Save global tile boundaries
                    chunk%tiles(tile_id)%t_left=left
                    chunk%tiles(tile_id)%t_right=right
                    chunk%tiles(tile_id)%t_top=top
                    chunk%tiles(tile_id)%t_bottom=bottom
                    chunk%tiles(tile_id)%t_front=front
                    chunk%tiles(tile_id)%t_back=back



                    IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
                    tile_id=tile_id+1
                ENDDO
                add_x_prev=0
                IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
            ENDDO
            add_x_prev=0
            add_y_prev=0
            IF(cz.LE.mod_z)add_z_prev=add_z_prev+1
        ENDDO

        IF(parallel%boss)THEN
            WRITE(g_out,*)
            WRITE(g_out,*)"Decomposing the Chunk into ",chunk_x," by ",chunk_y," by ",chunk_z," tiles"
            WRITE(g_out,*)
        ENDIF

    END SUBROUTINE clover_decompose_tile



    SUBROUTINE clover_allocate_buffers()

        IMPLICIT NONE

  
        ! Unallocated buffers for external boundaries caused issues on some systems so they are now
        !  all allocated
        ALLOCATE(chunk%left_snd_buffer(19*2*(chunk%y_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%left_rcv_buffer(19*2*(chunk%y_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%right_snd_buffer(19*2*(chunk%y_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%right_rcv_buffer(19*2*(chunk%y_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%bottom_snd_buffer(19*2*(chunk%x_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%bottom_rcv_buffer(19*2*(chunk%x_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%top_snd_buffer(19*2*(chunk%x_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%top_rcv_buffer(19*2*(chunk%x_max+5)*(chunk%z_max+5)))
        ALLOCATE(chunk%back_snd_buffer(19*2*(chunk%x_max+5)*(chunk%y_max+5)))
        ALLOCATE(chunk%back_rcv_buffer(19*2*(chunk%x_max+5)*(chunk%y_max+5)))
        ALLOCATE(chunk%front_snd_buffer(19*2*(chunk%x_max+5)*(chunk%y_max+5)))
        ALLOCATE(chunk%front_rcv_buffer(19*2*(chunk%x_max+5)*(chunk%y_max+5)))


    END SUBROUTINE clover_allocate_buffers

    SUBROUTINE clover_exchange(fields,depth)

        IMPLICIT NONE

        INTEGER      :: fields(:),depth, tile
        INTEGER      :: left_right_offset(19),bottom_top_offset(19),back_front_offset(19)
        INTEGER      :: request(4)
        INTEGER      :: message_count,err
        INTEGER      :: status(MPI_STATUS_SIZE,4)
        INTEGER      :: end_pack_index_left_right, end_pack_index_bottom_top,end_pack_index_back_front,field

        ! Assuming 1 patch per task, this will be changed

        request=0
        message_count=0


        end_pack_index_left_right=0
        end_pack_index_bottom_top=0
        end_pack_index_back_front=0
        DO field=1,19
            IF(fields(field).EQ.1) THEN
                left_right_offset(field)=end_pack_index_left_right
                bottom_top_offset(field)=end_pack_index_bottom_top
                back_front_offset(field)=end_pack_index_back_front
                end_pack_index_left_right=end_pack_index_left_right+depth*(chunk%y_max+5)*(chunk%z_max+5)
                end_pack_index_bottom_top=end_pack_index_bottom_top+depth*(chunk%x_max+5)*(chunk%z_max+5)
                end_pack_index_back_front=end_pack_index_back_front+depth*(chunk%x_max+5)*(chunk%y_max+5)
            ENDIF
        ENDDO

        IF(chunk%chunk_neighbours(chunk_left).NE.external_face) THEN
                  ! do left exchanges
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_LEFT).EQ.1) THEN
                    CALL clover_pack_left(tile, fields, depth, left_right_offset)
                ENDIF
            ENDDO
            

            !send and recv messagse to the left
            CALL clover_send_recv_message_left(chunk%left_snd_buffer,                      &
                chunk%left_rcv_buffer,                      &
                end_pack_index_left_right,                    &
                1, 2,                                               &
                request(message_count+1), request(message_count+2))
            message_count = message_count + 2
        ENDIF

        IF(chunk%chunk_neighbours(chunk_right).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_RIGHT).EQ.1) THEN
                    ! do right exchanges
                    CALL clover_pack_right(tile, fields, depth, left_right_offset)
                ENDIF
            ENDDO
            

            !send message to the right
            CALL clover_send_recv_message_right(chunk%right_snd_buffer,                     &
                chunk%right_rcv_buffer,                     &
                end_pack_index_left_right,                    &
                2, 1,                                               &
                request(message_count+1), request(message_count+2))
            message_count = message_count + 2
        ENDIF

        !make a call to wait / sync
        CALL MPI_WAITALL(message_count,request,status,err)

        !unpack in left direction
        IF(chunk%chunk_neighbours(chunk_left).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_LEFT).EQ.1) THEN
                    CALL clover_unpack_left(fields, tile, depth,                      &
                        chunk%left_rcv_buffer,             &
                        left_right_offset)
                ENDIF
            ENDDO
        
        ENDIF


        !unpack in right direction
        IF(chunk%chunk_neighbours(chunk_right).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_RIGHT).EQ.1) THEN
                    CALL clover_unpack_right(fields, tile, depth,                     &
                        chunk%right_rcv_buffer,           &
                        left_right_offset)
                ENDIF
            ENDDO
        
        ENDIF

        message_count = 0
        request = 0

        IF(chunk%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_BOTTOM).EQ.1) THEN
                    ! do bottom exchanges
                    CALL clover_pack_bottom(tile, fields, depth, bottom_top_offset)
                ENDIF
            ENDDO
            

            !send message downwards
            CALL clover_send_recv_message_bottom(chunk%bottom_snd_buffer,                     &
                chunk%bottom_rcv_buffer,                     &
                end_pack_index_bottom_top,                     &
                3, 4,                                                &
                request(message_count+1), request(message_count+2))
            message_count = message_count + 2
        ENDIF

        IF(chunk%chunk_neighbours(chunk_top).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_TOP).EQ.1) THEN
                    ! do top exchanges
                    CALL clover_pack_top(tile, fields, depth, bottom_top_offset)
                ENDIF
            ENDDO
            

            !send message upwards
            CALL clover_send_recv_message_top(chunk%top_snd_buffer,                           &
                chunk%top_rcv_buffer,                           &
                end_pack_index_bottom_top,                        &
                4, 3,                                                   &
                request(message_count+1), request(message_count+2))
            message_count = message_count + 2
        ENDIF

        !need to make a call to wait / sync
        CALL MPI_WAITALL(message_count,request,status,err)

        !unpack in top direction
        IF( chunk%chunk_neighbours(chunk_top).NE.external_face ) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_TOP).EQ.1) THEN
                    CALL clover_unpack_top(fields, tile, depth,                       &
                        chunk%top_rcv_buffer,               &
                        bottom_top_offset)
                ENDIF
            ENDDO
        
        ENDIF

        !unpack in bottom direction
        IF(chunk%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_BOTTOM).EQ.1) THEN
                    CALL clover_unpack_bottom(fields, tile, depth,                   &
                        chunk%bottom_rcv_buffer,         &
                        bottom_top_offset)
                ENDIF
            ENDDO
        
        ENDIF

        message_count = 0
        request = 0

        IF(chunk%chunk_neighbours(chunk_back).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_BACK).EQ.1) THEN
                    ! do back exchanges
                    CALL clover_pack_back(tile, fields, depth, back_front_offset)
                ENDIF
            ENDDO
            

            !send message downwards
            CALL clover_send_recv_message_back(chunk%back_snd_buffer,                        &
                chunk%back_rcv_buffer,                      &
                end_pack_index_back_front,                    &
                5, 6,                                               &
                request(message_count+1), request(message_count+2))
            message_count = message_count + 2
        ENDIF

        IF(chunk%chunk_neighbours(chunk_front).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_FRONT).EQ.1) THEN
                    ! do top exchanges
                    CALL clover_pack_front(tile, fields, depth, back_front_offset)
                ENDIF
            ENDDO
            

            !send message upwards
            CALL clover_send_recv_message_front(chunk%front_snd_buffer,                       &
                chunk%front_rcv_buffer,                         &
                end_pack_index_back_front,                        &
                6, 5,                                                   &
                request(message_count+1), request(message_count+2))
            message_count = message_count + 2
        ENDIF

        !need to make a call to wait / sync
        CALL MPI_WAITALL(message_count,request,status,err)

        !unpack in front direction
        IF( chunk%chunk_neighbours(chunk_front).NE.external_face ) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_FRONT).EQ.1) THEN
                    CALL clover_unpack_front(fields, tile, depth,                       &
                        chunk%front_rcv_buffer,               &
                        back_front_offset)
                ENDIF
            ENDDO
        
        ENDIF

        !unpack in back direction
        IF(chunk%chunk_neighbours(chunk_back).NE.external_face) THEN
            
            DO tile=1,tiles_per_chunk
                IF(chunk%tiles(tile)%external_tile_mask(TILE_BACK).EQ.1) THEN
                    CALL clover_unpack_back(fields, tile, depth,                   &
                        chunk%back_rcv_buffer,         &
                        back_front_offset)
                ENDIF
            ENDDO
        
        ENDIF

    END SUBROUTINE clover_exchange

    SUBROUTINE clover_pack_left(tile, fields, depth, left_right_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER      :: fields(:),depth, tile
        INTEGER      :: left_right_offset(:)
        INTEGER         :: y_offset, z_offset

        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_PRESSURE),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_VISCOSITY),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_SOUNDSPEED),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%left_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF

    END SUBROUTINE clover_pack_left

    SUBROUTINE clover_send_recv_message_left(left_snd_buffer, left_rcv_buffer,      &
        total_size,                     &
        tag_send, tag_recv,                    &
        req_send, req_recv)

        REAL(KIND=8)    :: left_snd_buffer(:), left_rcv_buffer(:)
        INTEGER         :: left_task
        INTEGER         :: total_size, tag_send, tag_recv, err
        INTEGER         :: req_send, req_recv

        left_task =chunk%chunk_neighbours(chunk_left) - 1

        CALL MPI_ISEND(left_snd_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_send &
            ,MPI_COMM_WORLD,req_send,err)

        CALL MPI_IRECV(left_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_recv &
            ,MPI_COMM_WORLD,req_recv,err)

    END SUBROUTINE clover_send_recv_message_left

    SUBROUTINE clover_unpack_left(fields, tile, depth,                         &
        left_rcv_buffer,                              &
        left_right_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER         :: fields(:), tile, depth
        INTEGER         :: left_right_offset(:)
        REAL(KIND=8)    :: left_rcv_buffer(:)
        INTEGER         :: y_offset, z_offset

        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)


        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_PRESSURE),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_VISCOSITY),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_SOUNDSPEED),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_left(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%left_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF

    END SUBROUTINE clover_unpack_left

    SUBROUTINE clover_pack_right(tile, fields, depth, left_right_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER        :: tile, fields(:), depth, tot_packr, left_right_offset(:)
        INTEGER         :: y_offset, z_offset

        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_PRESSURE),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_VISCOSITY),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_SOUNDSPEED),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%right_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF

    END SUBROUTINE clover_pack_right

    SUBROUTINE clover_send_recv_message_right(right_snd_buffer, right_rcv_buffer,   &
        total_size,                    &
        tag_send, tag_recv,                   &
        req_send, req_recv)

        IMPLICIT NONE

        REAL(KIND=8) :: right_snd_buffer(:), right_rcv_buffer(:)
        INTEGER      :: right_task
        INTEGER      :: total_size, tag_send, tag_recv, err
        INTEGER      :: req_send, req_recv

        right_task=chunk%chunk_neighbours(chunk_right) - 1

        CALL MPI_ISEND(right_snd_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_send, &
            MPI_COMM_WORLD,req_send,err)

        CALL MPI_IRECV(right_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_recv, &
            MPI_COMM_WORLD,req_recv,err)

    END SUBROUTINE clover_send_recv_message_right

    SUBROUTINE clover_unpack_right(fields, tile, depth,                          &
        right_rcv_buffer,                              &
        left_right_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER         :: fields(:), tile, total_in_right_buff, depth, left_right_offset(:)
        REAL(KIND=8)    :: right_rcv_buffer(:)
        INTEGER         :: y_offset, z_offset

        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_DENSITY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_ENERGY1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_PRESSURE),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_VISCOSITY),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    left_right_offset(FIELD_SOUNDSPEED),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_XVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_YVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL0),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    left_right_offset(FIELD_ZVEL1),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_VOL_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_X),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Y),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_right(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%right_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    left_right_offset(FIELD_MASS_FLUX_Z),            &
                    y_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_unpack_right

    SUBROUTINE clover_pack_top(tile, fields, depth, bottom_top_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER        :: tile, fields(:), depth, bottom_top_offset(:)
        INTEGER         :: x_offset, z_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_PRESSURE),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_VISCOSITY),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_SOUNDSPEED),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%top_snd_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_pack_top

    SUBROUTINE clover_send_recv_message_top(top_snd_buffer, top_rcv_buffer,     &
        total_size,                  &
        tag_send, tag_recv,                 &
        req_send, req_recv)

        IMPLICIT NONE

        REAL(KIND=8) :: top_snd_buffer(:), top_rcv_buffer(:)
        INTEGER      :: top_task
        INTEGER      :: total_size, tag_send, tag_recv, err
        INTEGER      :: req_send, req_recv

        top_task=chunk%chunk_neighbours(chunk_top) - 1

        CALL MPI_ISEND(top_snd_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_send, &
            MPI_COMM_WORLD,req_send,err)

        CALL MPI_IRECV(top_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_recv, &
            MPI_COMM_WORLD,req_recv,err)

    END SUBROUTINE clover_send_recv_message_top

    SUBROUTINE clover_unpack_top(fields, tile, depth,                        &
        top_rcv_buffer,                              &
        bottom_top_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER         :: fields(:), tile, total_in_top_buff, depth, bottom_top_offset(:)
        REAL(KIND=8)    :: top_rcv_buffer(:)
        INTEGER         :: x_offset, z_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)


        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_PRESSURE),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_VISCOSITY),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_SOUNDSPEED),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_top(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%top_rcv_buffer,                 &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_unpack_top

    SUBROUTINE clover_pack_bottom(tile, fields, depth, bottom_top_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER        :: tile, fields(:), depth, tot_packb, bottom_top_offset(:)
        INTEGER         :: x_offset, z_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_PRESSURE),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_VISCOSITY),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_SOUNDSPEED),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%bottom_snd_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_pack_bottom

    SUBROUTINE clover_send_recv_message_bottom(bottom_snd_buffer, bottom_rcv_buffer,        &
        total_size,                           &
        tag_send, tag_recv,                          &
        req_send, req_recv)

        IMPLICIT NONE

        REAL(KIND=8) :: bottom_snd_buffer(:), bottom_rcv_buffer(:)
        INTEGER      :: bottom_task
        INTEGER      :: total_size, tag_send, tag_recv, err
        INTEGER      :: req_send, req_recv

        bottom_task=chunk%chunk_neighbours(chunk_bottom) - 1

        CALL MPI_ISEND(bottom_snd_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_send &
            ,MPI_COMM_WORLD,req_send,err)

        CALL MPI_IRECV(bottom_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_recv &
            ,MPI_COMM_WORLD,req_recv,err)

    END SUBROUTINE clover_send_recv_message_bottom

    SUBROUTINE clover_unpack_bottom(fields, tile, depth,                        &
        bottom_rcv_buffer,                              &
        bottom_top_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER         :: fields(:), tile, depth, bottom_top_offset(:)
        REAL(KIND=8)    :: bottom_rcv_buffer(:)
        INTEGER         :: x_offset, z_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        z_offset = (chunk%tiles(tile)%t_back - chunk%back)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_DENSITY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_ENERGY1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_PRESSURE),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_VISCOSITY),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    bottom_top_offset(FIELD_SOUNDSPEED),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_XVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_YVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL0),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    bottom_top_offset(FIELD_ZVEL1),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_bottom(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%bottom_rcv_buffer,              &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    bottom_top_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, z_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_unpack_bottom

    SUBROUTINE clover_pack_back(tile, fields, depth, back_front_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER        :: tile, fields(:), depth, back_front_offset(:)
        INTEGER         :: x_offset, y_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_PRESSURE),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_VISCOSITY),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_SOUNDSPEED),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%back_snd_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_pack_back

    SUBROUTINE clover_send_recv_message_back(back_snd_buffer, back_rcv_buffer,     &
        total_size,                  &
        tag_send, tag_recv,                 &
        req_send, req_recv)

        IMPLICIT NONE

        REAL(KIND=8) :: back_snd_buffer(:), back_rcv_buffer(:)
        INTEGER      :: back_task
        INTEGER      :: total_size, tag_send, tag_recv, err
        INTEGER      :: req_send, req_recv

        back_task=chunk%chunk_neighbours(chunk_back)-1

        CALL MPI_ISEND(back_snd_buffer,total_size,MPI_DOUBLE_PRECISION,back_task,tag_send, &
            MPI_COMM_WORLD,req_send,err)

        CALL MPI_IRECV(back_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,back_task,tag_recv, &
            MPI_COMM_WORLD,req_recv,err)

    END SUBROUTINE clover_send_recv_message_back

    SUBROUTINE clover_unpack_back(fields, tile, depth,                        &
        back_rcv_buffer,                             &
        back_front_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER         :: fields(:), tile, depth, back_front_offset(:)
        REAL(KIND=8)    :: back_rcv_buffer(:)
        INTEGER         :: x_offset, y_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_PRESSURE),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_VISCOSITY),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_SOUNDSPEED),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_back(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%back_rcv_buffer,                &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_unpack_back

    SUBROUTINE clover_pack_front(tile, fields, depth, back_front_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER        :: tile, fields(:), depth, tot_packb, back_front_offset(:)
        INTEGER         :: x_offset, y_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_PRESSURE),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_VISCOSITY),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_SOUNDSPEED),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_pack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%front_snd_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
    END SUBROUTINE clover_pack_front

    SUBROUTINE clover_send_recv_message_front(front_snd_buffer, front_rcv_buffer,        &
        total_size,                         &
        tag_send, tag_recv,                        &
        req_send, req_recv)

        IMPLICIT NONE

        REAL(KIND=8) :: front_snd_buffer(:), front_rcv_buffer(:)
        INTEGER      :: front_task
        INTEGER      :: total_size, tag_send, tag_recv, err
        INTEGER      :: req_send, req_recv

        front_task=chunk%chunk_neighbours(chunk_front)-1

        CALL MPI_ISEND(front_snd_buffer,total_size,MPI_DOUBLE_PRECISION,front_task,tag_send &
            ,MPI_COMM_WORLD,req_send,err)

        CALL MPI_IRECV(front_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,front_task,tag_recv &
            ,MPI_COMM_WORLD,req_recv,err)

    END SUBROUTINE clover_send_recv_message_front

    SUBROUTINE clover_unpack_front(fields, tile, depth,                        &
        front_rcv_buffer,                            &
        back_front_offset)

        USE pack_kernel_module

        IMPLICIT NONE

        INTEGER         :: fields(:), tile, depth, back_front_offset(:)
        REAL(KIND=8)    :: front_rcv_buffer(:)
        INTEGER         :: x_offset, y_offset

        x_offset = (chunk%tiles(tile)%t_left - chunk%left)
        y_offset = (chunk%tiles(tile)%t_bottom - chunk%bottom)

        IF(fields(FIELD_DENSITY0).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density0,                 &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_DENSITY1).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%density1,                 &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_DENSITY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY0).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy0,                  &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ENERGY1).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%energy1,                  &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_ENERGY1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_PRESSURE).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%pressure,                 &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_PRESSURE),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VISCOSITY).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%viscosity,                &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_VISCOSITY),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%soundspeed,               &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, CELL_DATA,                             &
                    back_front_offset(FIELD_SOUNDSPEED),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL0).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel0,                    &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_XVEL1).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%xvel1,                    &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_XVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL0).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel0,                    &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_YVEL1).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%yvel1,                    &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_YVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL0).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel0,                    &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL0),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_ZVEL1).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%zvel1,                    &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA, &
                    depth, VERTEX_DATA,                           &
                    back_front_offset(FIELD_ZVEL1),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_x,               &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_y,               &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_VOL_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%vol_flux_z,               &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_VOL_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_x,              &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, X_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_X),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_y,              &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Y_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Y),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF
        IF(fields(FIELD_MASS_FLUX_Z).EQ.1) THEN

                CALL clover_unpack_message_front(chunk%tiles(tile)%t_xmin,                    &
                    chunk%tiles(tile)%t_xmax,                    &
                    chunk%tiles(tile)%t_ymin,                    &
                    chunk%tiles(tile)%t_ymax,                    &
                    chunk%tiles(tile)%t_zmin,                    &
                    chunk%tiles(tile)%t_zmax,                    &
                    chunk%tiles(tile)%field%mass_flux_z,              &
                    chunk%front_rcv_buffer,               &
                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,Z_FACE_DATA,&
                    depth, Z_FACE_DATA,                           &
                    back_front_offset(FIELD_MASS_FLUX_Z),            &
                    x_offset, y_offset, chunk%x_max, chunk%y_max, chunk%z_max)
        ENDIF

    END SUBROUTINE clover_unpack_front

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
