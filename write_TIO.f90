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

!>  @brief Generates graphics output files.
!>  @author Wayne Gaudin
!>  @details The field data over all mesh chunks is written to a .vtk files and
!>  the .visit file is written that defines the time for each set of vtk files.
!>  The ideal gas and viscosity routines are invoked to make sure this data is
!>  up to data with the current energy, density and velocity.

MODULE write_TIO_module

  USE clover_module
  USE update_halo_module
  USE viscosity_module
  USE ideal_gas_module
  USE TyphonIO
  USE MPI

  IMPLICIT NONE


  INTEGER :: err
  INTEGER(kind=tio_sizek)                           :: counter=0
  INTEGER(kind=tio_errk)                            :: terr
  INTEGER(kind=tio_filek)                           :: file_id
  INTEGER(kind=tio_objk)                            :: state_id
  INTEGER(kind=tio_objk)                            :: mesh_id
  INTEGER(kind=tio_objk)                            :: quant_id , quant_id1, quant_id2, quant_id3 
  INTEGER(kind=tio_objk)                            :: quant_id4, quant_id5, quant_id6


  INTEGER(kind=tio_sizek)                           :: order
  INTEGER(kind=tio_sizek)                           :: numchunks
  INTEGER(kind=tio_sizek)                           :: index_timer=0
  INTEGER(kind=tio_sizek)                           :: tio_dump_freq
  CHARACTER(len=tio_strlen_f)                       :: path
 
  !File information
  CHARACTER(len=tio_strlen_f)                       :: code_name
  CHARACTER(len=tio_strlen_f)                       :: version_name
  CHARACTER(len=tio_strlen_f)                       :: title_name
  CHARACTER(len=tio_strlen_f)                       :: datec,tempc
  CHARACTER(len=tio_strlen_f)                       :: tiofname
  INTEGER(kind=tio_sizek)                           :: mpiinfo=MPI_INFO_NULL

  !State information
  REAL(kind=tio_timek)                              :: state_time
  INTEGER(kind=tio_stepk)                           :: state_num
  CHARACTER(len=tio_strlen_f)                       :: state_name
  INTEGER(kind=tio_stepk)                           :: step_num=0_tio_stepk

  !!$  Data structure to store variables required in this module
  TYPE :: clover_struct
     REAL(kind=8)                                   :: sim_time
     REAL(kind=8)                                   :: start_time
     REAL(kind=8)                                   :: dx,dy,dz
     INTEGER(kind=tio_sizek), dimension(3)          :: coords
     INTEGER(kind=tio_sizek), dimension(3)          :: dimsize
     INTEGER(kind=tio_sizek)                        :: nxmax, nymax, nzmax
     INTEGER(kind=tio_sizek)                        :: ndum_fx, ndum_fy, ndum_fz
     INTEGER(kind=tio_sizek)                        :: ngdix, ngdiy, ngdiz=1
     INTEGER(kind=tio_sizek)                        :: nxt, nyt, nzt
     INTEGER(kind=tio_sizek)                        :: ndims=1
     INTEGER(kind=tio_sizek)                        :: mpi_world

  END TYPE clover_struct

  save

contains


  SUBROUTINE start_TIO(file_id)

    IMPLICIT NONE

    integer(kind=TIO_FILEK), intent(out)       :: file_id
    INTEGER(kind=tio_sizek)                    :: counter=0

    INTEGER                                    :: fields(NUM_FIELDS)
    LOGICAL                                    :: fileExists =.false.


    !Set file information
    code_name         = "CloverLeaf3D"
    version_name      = "CloverLeaf3D-TyphonIO"
    CALL DATE_AND_TIME(date=tempc)
    datec             = TRIM(tempc(7:8)//'/'//tempc(5:6)//'/'//tempc(1:4))
    title_name        = 'TyphonIO porting'

    tiofname          = 'output.h5'

    
    IF(parallel%boss) THEN
      INQUIRE(FILE=tiofname,EXIST=fileExists)
    ENDIF
 
    IF(fileExists .EQV. .FALSE.) THEN
      !On the first instance, need to create a TIO file
      terr = TIO_Create_f( filename = trim(adjustl(tiofname)),   &
           fileID    = file_id, ACCESS = TIO_ACC_REPLACE_F,      &
           codename  = code_name, version = version_name,        &
           date      = datec, title = title_name,                &
           comm      = MPI_COMM_WORLD, info = MPIINFO,           &
           rank      = parallel%task)
      !write(*,*) "File name", trim(adjustl(tiofname))
      !write(*,*) "Task is", parallel%task
      IF(terr /=  TIO_SUCCESS_F) WRITE(*,*) "Error from TIO ", terr
    ELSE
        WRITE(g_out,*) 'Error File alread exists'
        WRITE(0    ,*) 'Error File alread exists'
        RETURN
    ENDIF  

  END SUBROUTINE start_TIO 



  SUBROUTINE write_TIO_Results(file_id)

    IMPLICIT NONE
    integer(kind=TIO_FILEK), intent(in) :: file_id
    CHARACTER(len=tio_strlen_f)   :: cname

    INTEGER(kind=TIO_SIZEK)                    :: n1, n2, n3
    INTEGER(kind=TIO_SIZEK)                    :: c=0
    INTEGER                                    :: xx, yy, zz

    INTEGER(kind=TIO_SIZEK), DIMENSION(6)      :: nn
    
    INTEGER(kind=TIO_SIZEK)                    :: il, ih, jl, jh, kl, kh

    !Number of points along x/y/z
    INTEGER(kind=tio_sizek)                    :: i,j,k
    !size of cells
    REAL(kind=4)                               :: delta_x, delta_y, delta_z 
    REAL(kind=4),DIMENSION(:), ALLOCATABLE     :: xcoords,ycoords,zcoords
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: ps_data
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: er_data
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: ds_data
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: vs_data
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: vsdatax
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: vsdatay
    REAL(kind=4),DIMENSION(:,:,:), ALLOCATABLE :: vsdataz

    !State info
    WRITE(cname,FMT="(I0)") step
    state_name="State_" // trim(adjustl(cname))
    state_time = time
    step_num = step + 1_TIO_STEPK
    state_num = step
 
    !Create State
    terr = TIO_Create_State_f( fileID = file_id,                 &
         name = state_name, stateID = state_id,                  &
         step = state_num, time = state_time,                    &
         units = "-"    )
    !If file is alaredy there, trying to create a state with same name will cause an error
    !following handles that error
    IF (terr /= TIO_SUCCESS_F) THEN
        WRITE(*,*) "Error from TIO - State already exists. Returning... ", terr
        terr = TIO_Close_F(fileID = file_id)
        IF (terr /= TIO_SUCCESS_F)  WRITE(*,*) "Error from TIO - File close ", terr
        RETURN
    ENDIF

    !Create the mesh
    xx = grid%x_cells+1
    yy = grid%y_cells+1
    zz = grid%z_cells+1


    terr = TIO_Create_Mesh_f( fileID       = file_id,                    &
                            stateID        = state_id,                   &
                            name           = "Mesh",                     &
                            meshID         = mesh_id,                    &
                            meshtype       = TIO_MESH_QUAD_COLINEAR_F,   &
                            coordtype      = TIO_COORD_CARTESIAN_F,      &
                            isAMR          = .FALSE.,                    &
                            group          = 'Entire region',            &
                            order          = 1,                          &
                            graph_datatype = TIO_DATATYPE_NULL_F,        &
                            coord_datatype = TIO_REAL4_F,                &
                            ndims          = TIO_3D_F,                   &
                            n1             = xx,                         &
                            n2             = yy,                         &
                            n3             = zz,                         &
                            n4             = TIO_NULL_F,                 &
                            nchunks        = number_of_chunks,           &
                            iunits         = "cm",                       &
                            junits         = "cm",                       &
                            kunits         = "cm",                       &
                            ilabel         = "x",                        &
                            jlabel         = "y",                        &
                            klabel         = "z"                          )



    IF (terr /= TIO_SUCCESS_F) write(*,*) "Error from TIO -Create mesh", terr
    ! Can use overload and function for cleaner structure (for cell+node)

    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "pressure",                      &
                              quantID  = quant_id,                        &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_CELL_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")


    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "energy",                        &
                              quantID  = quant_id1,                       &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_CELL_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")


    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "density",                       &
                              quantID  = quant_id2,                       &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_CELL_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")

    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "viscosity",                     &
                              quantID  = quant_id3,                       &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_CELL_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")


    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "x_vel",                         &
                              quantID  = quant_id4,                       &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_NODE_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")

    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "y_vel",                         &
                              quantID  = quant_id5,                       &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_NODE_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")

    terr = TIO_Create_Quant_f(fileID   = file_id,                         &
                              meshID   = mesh_id,                         &
                              name     = "z_vel",                         &
                              quantID  = quant_id6,                       &
                              datatype = TIO_REAL4_F,                     &
                              centring = TIO_CENTRE_NODE_F,               &
                              nghosts  = 0,                               &
                              ismixed  = .false.,                         &
                              units    = "unknown")

    DO c=1,number_of_chunks


    nn(1) = chunks(1)%field%left
    nn(2) = chunks(1)%field%right+1
    nn(3) = chunks(1)%field%bottom
    nn(4) = chunks(1)%field%top+1
    nn(5) = chunks(1)%field%back
    nn(6) = chunks(1)%field%front+1

    CALL MPI_Bcast( nn ,6, MPI_INTEGER, c-1, MPI_COMM_WORLD, err) 

    terr = TIO_Set_Quad_Chunk_f( fileID   = file_id,     & 
                                 meshID   = mesh_id,     &
                                 idx      = c,           &
                                 ndims    = TIO_3D_F,    &
                                 il       = nn(1),       &
                                 ih       = nn(2),       &
                                 jl       = nn(3),       &
                                 jh       = nn(4),       &
                                 kl       = nn(5),       &
                                 kh       = nn(6),       &
                                 nmixcell = 0,           &
                                 nmixcomp = 0            )

    END DO

    
   IF (parallel%boss) THEN    
    !Work out the mesh coordinates - 

      delta_x = (grid%xmax-grid%xmin)/grid%x_cells
      delta_y = (grid%ymax-grid%ymin)/grid%y_cells
      delta_z = (grid%zmax-grid%zmin)/grid%z_cells
      
      IF(.NOT. allocated(xcoords)) ALLOCATE (xcoords(xx))
      IF(.NOT. allocated(ycoords)) ALLOCATE (ycoords(yy))
      IF(.NOT. allocated(zcoords)) ALLOCATE (zcoords(zz))
     
      DO i=1,xx
        xcoords(i) = ((i-1)*delta_x)
      END DO
      DO i=1,yy
        ycoords(i) = ((i-1)*delta_y)
      END DO
      DO i=1,zz
        zcoords(i) = ((i-1)*delta_z)
      END DO


     terr = TIO_Write_QuadMesh_All_f( fileID = file_id,             &
                                      meshID     = mesh_id,         &
                                      datatype   = TIO_REAL4_F,     &
                                      icoords    = xcoords,         &
                                      jcoords    = ycoords,         &
                                      kcoords    = zcoords)

 
      IF (terr /= TIO_SUCCESS_F) write(*,*) "Error from TIO - write chunk", terr

      deallocate (xcoords,ycoords,zcoords)

    END IF
 
      c= parallel%task+1

      il = chunks(1)%field%x_min
      ih = chunks(1)%field%x_max
      jl = chunks(1)%field%y_min
      jh = chunks(1)%field%y_max
      kl = chunks(1)%field%z_min
      kh = chunks(1)%field%z_max

      IF(.NOT. allocated(ps_data)) ALLOCATE (ps_data(il:ih,jl:jh,kl:kh))

      ps_data(il:ih,jl:jh,kl:kh) = chunks(1)%field%pressure(il:ih,jl:jh,kl:kh)
      
      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id,                 &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = ps_data(il:ih,jl:jh,kl:kh) )

      deallocate (ps_data)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                             quantID  = quant_id)


      IF(.NOT. allocated(er_data)) ALLOCATE (er_data(il:ih,jl:jh,kl:kh))
      er_data(il:ih,jl:jh,kl:kh) = chunks(1)%field%energy0(il:ih,jl:jh,kl:kh)
           

      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id1,                &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = er_data(il:ih,jl:jh,kl:kh) )

      deallocate (er_data)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                             quantID  = quant_id1)

      IF(.NOT. allocated(ds_data)) ALLOCATE (ds_data(il:ih,jl:jh,kl:kh))
      ds_data(il:ih,jl:jh,kl:kh) = chunks(1)%field%density0(il:ih,jl:jh,kl:kh)
      
      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id2,                &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = ds_data(il:ih,jl:jh,kl:kh) )

      deallocate (ds_data)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                             quantID  = quant_id2)

      IF(.NOT. allocated(vs_data)) ALLOCATE (vs_data(il:ih,jl:jh,kl:kh))
      vs_data(il:ih,jl:jh,kl:kh) = chunks(1)%field%viscosity(il:ih,jl:jh,kl:kh)
      
      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id3,                &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = vs_data(il:ih,jl:jh,kl:kh) )

      deallocate (vs_data)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                             quantID  = quant_id3)


      il = chunks(1)%field%x_min
      ih = chunks(1)%field%x_max+1
      jl = chunks(1)%field%y_min
      jh = chunks(1)%field%y_max+1
      kl = chunks(1)%field%z_min
      kh = chunks(1)%field%z_max+1


      IF(.NOT. allocated(vsdatax)) ALLOCATE (vsdatax(il:ih,jl:jh,kl:kh))
      vsdatax(il:ih,jl:jh,kl:kh) = chunks(1)%field%xvel0(il:ih,jl:jh,kl:kh)
      
      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id4,                &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = vsdatax(il:ih,jl:jh,kl:kh) )

      deallocate (vsdatax)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                             quantID  = quant_id4)


      IF(.NOT. allocated(vsdatay)) ALLOCATE (vsdatay(il:ih,jl:jh,kl:kh))
      vsdatay(il:ih,jl:jh,kl:kh) = chunks(1)%field%yvel0(il:ih,jl:jh,kl:kh)
      
      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id5,                &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = vsdatay(il:ih,jl:jh,kl:kh) )

      deallocate (vsdatay)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                               quantID  = quant_id5)


      IF(.NOT. allocated(vsdataz)) ALLOCATE (vsdataz(il:ih,jl:jh,kl:kh))
      vsdataz(il:ih,jl:jh,kl:kh) = chunks(1)%field%zvel0(il:ih,jl:jh,kl:kh)
      
      terr = TIO_Write_QuadQuant_Chunk_f(fileID   = file_id,                  &
                                         quantID  = quant_id6,                &
                                         idx      = c,                        &
                                         xfer     = TIO_XFER_COLLECTIVE_F,    &
                                         datatype = TIO_REAL4_F,              &
                                         qdat     = vsdataz(il:ih,jl:jh,kl:kh) )

      deallocate (vsdataz)

      terr = TIO_Close_Quant_f(fileID   = file_id,                  &
                               quantID  = quant_id6)

      terr = TIO_Close_Mesh_f(fileID   = file_id,                  &
                              meshID  = mesh_id)

      terr = TIO_Close_State_f(fileID   = file_id,                  &
                               stateID  = state_id)


    IF (parallel%boss) THEN
      WRITE(*,*) "State is ", state_name
    ENDIF
             
  END SUBROUTINE write_TIO_Results 

  SUBROUTINE Close_TIO(file_id)

    implicit none

    integer(kind=TIO_FILEK), intent(in) :: file_id

    terr = TIO_Close_f(file_id)
    IF (terr /= TIO_SUCCESS_F) WRITE(*,*) "Error from TIO - File close", terr 

  END SUBROUTINE Close_TIO
                          
END MODULE write_TIO_module
