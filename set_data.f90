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

!>  @brief allocate and sets data for standalone mode
!>  @author Wayne Gaudin
!>  @details Calls user requested kernel in standalone mode

MODULE set_data_module

CONTAINS

SUBROUTINE set_data(x_min,x_max,y_min,y_max,z_min,z_max,     &
                    vertexdx,                                &
                    vertexdy,                                &
                    vertexdz,                                &
                    celldx,                                  &
                    celldy,                                  &
                    celldz,                                  &
                    xarea,                                   &
                    yarea,                                   &
                    zarea,                                   &
                    volume,                                  &
                    density0,                                &
                    density1,                                &
                    energy0,                                 &
                    energy1,                                 &
                    viscosity,                               &
                    pressure,                                &
                    soundspeed,                              &
                    xvel0,                                   &
                    xvel1,                                   &
                    yvel0,                                   &
                    yvel1,                                   &
                    zvel0,                                   &
                    zvel1,                                   &
                    vol_flux_x,                              &
                    vol_flux_y,                              &
                    vol_flux_z,                              &
                    mass_flux_x,                             &
                    mass_flux_y,                             &
                    mass_flux_z,                             &
                    work_array1,                             &
                    work_array2,                             &
                    work_array3,                             &
                    work_array4,                             &
                    work_array5,                             &
                    work_array6,                             &
                    work_array7,                             &
                    work_array8,                             &
                    dt                                       )

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
  REAL(KIND=8),OPTIONAL :: dt
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vertexdx(:),vertexdy(:),vertexdz(:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: celldx(:),celldy(:),celldz(:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: xarea(:,:,:),yarea(:,:,:),zarea(:,:,:),volume(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: density0(:,:,:),density1(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: energy0(:,:,:),energy1(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: pressure(:,:,:),viscosity(:,:,:),soundspeed(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: xvel0(:,:,:),yvel0(:,:,:),zvel0(:,:,:),xvel1(:,:,:),yvel1(:,:,:),zvel1(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_x(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_y(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: vol_flux_z(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_x(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_y(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: mass_flux_z(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: work_array1(:,:,:),work_array2(:,:,:),work_array3(:,:,:),work_array4(:,:,:)
  REAL(KIND=8),ALLOCATABLE,OPTIONAL :: work_array5(:,:,:),work_array6(:,:,:),work_array7(:,:,:),work_array8(:,:,:)

  INTEGER :: j,k,l
  REAL(KIND=8),ALLOCATABLE :: vertexx(:),vertexy(:),vertexz(:)
  REAL(KIND=8),ALLOCATABLE :: cellx(:),celly(:),cellz(:)
  REAL(KIND=8) :: radius,theta,x,y,z,mult,dx,dy,dz,sound_speed_squared,v,pressurebyenergy,pressurebyvolume,width
  REAL(KIND=8)  :: ugradx1,ugradx2,vgrady1,vgrady2,wgradz1,wgradz2
  REAL(KIND=8)  :: ugrady1,ugrady2,vgradx1,vgradx2,wgradx1,wgradx2
  REAL(KIND=8)  :: ugradz1,ugradz2,vgradz1,vgradz2,wgrady1,wgrady2
  REAL(KIND=8)  :: xx,yy,zz,xy,xz,yz
  REAL(KIND=8)  :: grad2,pgradx,pgrady,pgradz,pgradx2,pgrady2,pgradz2,grad     &
                  ,ygrad,pgrad,xgrad,zgrad,div,limiter

  ! Set the initial data

  dx=(10.0_8)/float(x_max-x_min+1)
  dy=(10.0_8)/float(y_max-y_min+1)
  dz=(10.0_8)/float(z_max-z_min+1)

  ALLOCATE(vertexx(x_min-2:x_max+3))
  ALLOCATE(vertexy(y_min-2:y_max+3))
  ALLOCATE(vertexz(z_min-2:z_max+3))
  ALLOCATE(cellx(x_min-2:x_max+2))
  ALLOCATE(celly(y_min-2:y_max+2))
  ALLOCATE(cellz(z_min-2:z_max+2))
  IF(PRESENT(vertexdx)) THEN
    ALLOCATE(vertexdx(x_min-2:x_max+3))
  ENDIF
  IF(PRESENT(vertexdy)) THEN
    ALLOCATE(vertexdy(y_min-2:y_max+3))
  ENDIF
  IF(PRESENT(vertexdz)) THEN
    ALLOCATE(vertexdz(z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(celldx)) THEN
    ALLOCATE(celldx(x_min-2:x_max+2))
  ENDIF
  IF(PRESENT(celldy)) THEN
    ALLOCATE(celldy(y_min-2:y_max+2))
  ENDIF
  IF(PRESENT(celldz)) THEN
    ALLOCATE(celldz(y_min-2:z_max+2))
  ENDIF
  IF(PRESENT(xarea)) THEN
    ALLOCATE(xarea(x_min-2:x_max+3 ,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(yarea)) THEN
    ALLOCATE(yarea(x_min-2:x_max+2 ,y_min-2:y_max+3,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(zarea)) THEN
    ALLOCATE(zarea(x_min-2:x_max+2 ,y_min-2:y_max+2,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(volume)) THEN
    ALLOCATE(volume(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(density0)) THEN
    ALLOCATE(density0(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(density1)) THEN
    ALLOCATE(density1(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(energy0)) THEN
    ALLOCATE(energy0(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(energy1)) THEN
    ALLOCATE(energy1(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(pressure)) THEN
    ALLOCATE(pressure(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(soundspeed)) THEN
    ALLOCATE(soundspeed(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(viscosity)) THEN
    ALLOCATE(viscosity(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(xvel0)) THEN
    ALLOCATE(xvel0(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(yvel0)) THEN
    ALLOCATE(yvel0(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(zvel0)) THEN
    ALLOCATE(zvel0(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(xvel1)) THEN
    ALLOCATE(xvel1(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(yvel1)) THEN
    ALLOCATE(yvel1(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(zvel1)) THEN
    ALLOCATE(zvel1(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array1)) THEN
    ALLOCATE(work_array1(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array2)) THEN
    ALLOCATE(work_array2(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array3)) THEN
    ALLOCATE(work_array3(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array4)) THEN
    ALLOCATE(work_array4(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array5)) THEN
    ALLOCATE(work_array5(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array6)) THEN
    ALLOCATE(work_array6(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array7)) THEN
    ALLOCATE(work_array7(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(work_array8)) THEN
    ALLOCATE(work_array8(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(vol_flux_x)) THEN
    ALLOCATE(vol_flux_x(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(vol_flux_y)) THEN
    ALLOCATE(vol_flux_y(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(vol_flux_z)) THEN
    ALLOCATE(vol_flux_z(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3))
  ENDIF
  IF(PRESENT(mass_flux_x)) THEN
    ALLOCATE(mass_flux_x(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(mass_flux_y)) THEN
    ALLOCATE(mass_flux_y(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2))
  ENDIF
  IF(PRESENT(mass_flux_z)) THEN
    ALLOCATE(mass_flux_z(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3))
  ENDIF

!$OMP PARALLEL

!$OMP DO
  DO j=x_min-2,x_max+3
     vertexx(j)=0.0_8+dx*float(j-x_min)
  ENDDO
!$OMP ENDDO

!$OMP DO
  DO k=y_min-2,y_max+3
     vertexy(k)=0.0_8+dy*float(k-y_min)
  ENDDO
!$OMP ENDDO

!$OMP DO
  DO l=z_min-2,z_max+3
     vertexz(l)=0.0_8+dz*float(l-z_min)
  ENDDO
!$OMP ENDDO

!$OMP DO
  DO j=x_min-2,x_max+2
    cellx(j)=0.5*(vertexx(j)+vertexx(j+1))
  ENDDO
!$OMP ENDDO

!$OMP DO
  DO k=y_min-2,y_max+2
    celly(k)=0.5*(vertexy(k)+vertexy(k+1))
  ENDDO
!$OMP ENDDO

!$OMP DO
  DO l=z_min-2,z_max+2
    cellz(l)=0.5*(vertexz(l)+vertexz(l+1))
  ENDDO
!$OMP ENDDO

  IF(PRESENT(vertexdx)) THEN
!$OMP DO
    DO j=x_min-2,x_max+3
      vertexdx(j)=dx
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vertexdy)) THEN
!$OMP DO
    DO k=y_min-2,y_max+3
      vertexdy(k)=dy
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vertexdz)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      vertexdz(l)=dz
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(celldx)) THEN
!$OMP DO
    DO j=x_min-2,x_max+2
      celldx(j)=dx
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(celldy)) THEN
!$OMP DO
    DO k=y_min-2,y_max+2
      celldy(k)=dy
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(celldz)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      celldz(l)=dz
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(xarea)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          xarea(j,k,l)=celldy(k)*celldz(l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(yarea)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          yarea(j,k,l)=celldx(j)*celldz(l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(zarea)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          zarea(j,k,l)=celldx(j)*celldy(k)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(volume)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          volume(j,k,l)=dx*dy*dz
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(density0)) THEN
!$OMP DO PRIVATE(radius,x,y,z)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          x=(cellx(j)-5.0_8)
          y=(celly(k)-5.0_8)
          z=(cellz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(radius.LE.2.5_8) THEN
            density0(j,k,l)=2.0_8-(radius*2.0_8/10.0_8)
          ELSE
            density0(j,k,l)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(density1)) THEN
!$OMP DO PRIVATE(radius,x,y,z)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          x=(cellx(j)-5.0_8)
          y=(celly(k)-5.0_8)
          z=(cellz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(radius.LE.2.5_8) THEN
            density1(j,k,l)=2.0_8-(radius*2.0_8/10.0_8)
          ELSE
            density1(j,k,l)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF
 
  IF(PRESENT(energy0)) THEN
!$OMP DO PRIVATE(radius,x,y,z)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          x=(cellx(j)-5.0_8)
          y=(celly(k)-5.0_8)
          z=(cellz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(radius.LE.2.5_8) THEN
            energy0(j,k,l)=2.0_8-(radius*2.0_8/10.0_8)
          ELSE
            energy0(j,k,l)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(energy1)) THEN
!$OMP DO PRIVATE(radius,x,y,z)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          x=(cellx(j)-5.0_8)
          y=(celly(k)-5.0_8)
          z=(cellz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(radius.LE.2.5_8) THEN
            energy1(j,k,l)=2.0_8-(radius*2.0_8/10.0_8)
          ELSE
            energy1(j,k,l)=1.0_8-((radius-5.0_8)*1.0_8/20.0_8)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(pressure)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          pressure(j,k,l)=(1.4_8-1.0_8)*density0(j,k,l)*energy0(j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(soundspeed)) THEN
!$OMP DO PRIVATE(v,pressurebyenergy,pressurebyvolume,sound_speed_squared)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          v=1.0_8/density0(j,k,l)
          pressurebyenergy=(1.4_8-1.0_8)*density0(j,k,l)
          pressurebyvolume=-density0(j,k,l)*pressure(j,k,l)
          sound_speed_squared=v*v*(pressure(j,k,l)*pressurebyenergy-pressurebyvolume)
          soundspeed(j,k,l)=SQRT(sound_speed_squared)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

! OMP THIS WITH A REDUCTION
  IF(PRESENT(dt)) THEN
    dt=0.0_8
    width=MIN(dx,dy,dz)
!$OMP DO REDUCTION(MAX : dt)
    DO l=z_min,z_max
      DO k=y_min,y_max
        DO j=x_min,x_max
          IF(soundspeed(j,k,l).GT.dt) dt=soundspeed(j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
    dt=width*0.7_8/dt
  ENDIF

  IF(PRESENT(xvel0)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,z,mult)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          x=(vertexx(j)-5.0_8)
          y=(vertexy(k)-5.0_8)
          z=(vertexz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.NE.0.0_8) THEN
            theta=atan(y/x)
          ELSE
            theta=atan(y/-0.000000001_8)
          ENDIF
          IF(radius.LE.2.5_8) THEN
            xvel0(j,k,l)=mult*(2.0_8-(radius*2.0_8/10.0_8))*sin(theta)
          ELSE
            xvel0(j,k,l)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*sin(theta)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(yvel0)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,z,mult)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          x=(vertexx(j)-5.0_8)
          y=(vertexy(k)-5.0_8)
          z=(vertexz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.NE.0.0_8) THEN
            theta=atan(y/x)
          ELSE
            theta=atan(y/-0.000000001_8)
          ENDIF
          IF(radius.LE.2.5_8) THEN
            yvel0(j,k,l)=mult*(2.0_8-(radius*2.0_8/10.0_8))*cos(theta)
          ELSE
            yvel0(j,k,l)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*cos(theta)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(zvel0)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,z,mult)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          x=(vertexx(j)-5.0_8)
          y=(vertexy(k)-5.0_8)
          z=(vertexz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.NE.0.0_8) THEN
            theta=atan(z/x)
          ELSE
            theta=atan(z/-0.000000001_8)
          ENDIF
          IF(radius.LE.2.5_8) THEN
            zvel0(j,k,l)=mult*(2.0_8-(radius*2.0_8/10.0_8))*cos(theta)
          ELSE
            zvel0(j,k,l)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*cos(theta)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(xvel1)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,z,mult)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          x=(vertexx(j)-5.0_8)
          y=(vertexy(k)-5.0_8)
          z=(vertexz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.NE.0.0_8) THEN
            theta=atan(y/x)
          ELSE
            theta=atan(y/-0.000000001_8)
          ENDIF
          IF(radius.LE.2.5_8) THEN
            xvel1(j,k,l)=mult*(2.0_8-(radius*2.0_8/10.0_8))*sin(theta)
          ELSE
            xvel1(j,k,l)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*sin(theta)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(yvel1)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,z,mult)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          x=(vertexx(j)-5.0_8)
          y=(vertexy(k)-5.0_8)
          z=(vertexz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.NE.0.0_8) THEN
            theta=atan(y/x)
          ELSE
            theta=atan(y/-0.000000001_8)
          ENDIF
          IF(radius.LE.2.5_8) THEN
            yvel1(j,k,l)=mult*(2.0_8-(radius*2.0_8/10.0_8))*cos(theta)
          ELSE
            yvel1(j,k,l)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*cos(theta)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(zvel1)) THEN
!$OMP DO PRIVATE(radius,theta,x,y,z,mult)
    DO l=z_min-2,z_max+2
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          x=(vertexx(j)-5.0_8)
          y=(vertexy(k)-5.0_8)
          z=(vertexz(l)-5.0_8)
          radius=sqrt(x**2.0_8+y**2.0_8+z**2.0)
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.LE.0.0_8) mult=1.0_8
          IF(x.LE.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.LE.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=-1.0_8
          IF(x.GT.0.0_8.AND.y.LE.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.GT.0.0_8.AND.y.GT.0.0_8.AND.z.GT.0.0_8) mult=1.0_8
          IF(x.NE.0.0_8) THEN
            theta=atan(z/x)
          ELSE
            theta=atan(z/-0.000000001_8)
          ENDIF
          IF(radius.LE.2.5_8) THEN
            zvel1(j,k,l)=mult*(2.0_8-(radius*2.0_8/10.0_8))*cos(theta)
          ELSE
            zvel1(j,k,l)=mult*(1.0_8-((radius-5.0_8)*1.0_8/20.0_8))*cos(theta)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(viscosity)) THEN
!$OMP DO PRIVATE(ugradx1,ugradx2,vgrady1,vgrady2,wgradz1,wgradz2,div,                     &
!$OMP            ugrady1,ugrady2,vgradx1,vgradx2,wgradx1,wgradx2,                         &
!$OMP            ugradz1,ugradz2,vgradz1,vgradz2,wgrady1,wgrady2,                         &
!$OMP            pgradx,pgrady,pgradz,pgradx2,pgrady2,pgradz2,limiter,                    &
!$OMP            pgrad,xgrad,ygrad,zgrad,grad,grad2,xx,yy,zz,xy,xz,yz)
  DO l=z_min,z_max
    DO k=y_min,y_max
      DO j=x_min,x_max
        ugradx1=xvel0(j  ,k  ,l  )+xvel0(j  ,k+1,l  )+xvel0(j  ,k  ,l+1)+xvel0(j  ,k+1,l+1)
        ugradx2=xvel0(j+1,k  ,l  )+xvel0(j+1,k+1,l  )+xvel0(j+1,k  ,l+1)+xvel0(j+1,k+1,l+1)
        ugrady1=xvel0(j  ,k  ,l  )+xvel0(j+1,k  ,l  )+xvel0(j  ,k  ,l+1)+xvel0(j+1,k  ,l+1)
        ugrady2=xvel0(j  ,k+1,l  )+xvel0(j+1,k+1,l  )+xvel0(j  ,k+1,l+1)+xvel0(j+1,k+1,l+1)
        ugradz1=xvel0(j  ,k  ,l  )+xvel0(j+1,k  ,l  )+xvel0(j  ,k+1,l  )+xvel0(j+1,k+1,l  )
        ugradz2=xvel0(j  ,k  ,l+1)+xvel0(j+1,k  ,l+1)+xvel0(j  ,k+1,l+1)+xvel0(j+1,k+1,l+1)

        vgradx1=yvel0(j  ,k  ,l  )+yvel0(j  ,k+1,l  )+yvel0(j  ,k  ,l+1)+yvel0(j  ,k+1,l+1)
        vgradx2=yvel0(j+1,k  ,l  )+yvel0(j+1,k+1,l  )+yvel0(j+1,k  ,l+1)+yvel0(j+1,k+1,l+1)
        vgrady1=yvel0(j  ,k  ,l  )+yvel0(j+1,k  ,l  )+yvel0(j  ,k  ,l+1)+yvel0(j+1,k  ,l+1)
        vgrady2=yvel0(j  ,k+1,l  )+yvel0(j+1,k+1,l  )+yvel0(j  ,k+1,l+1)+yvel0(j+1,k+1,l+1)
        vgradz1=yvel0(j  ,k  ,l  )+yvel0(j+1,k  ,l  )+yvel0(j  ,k+1,l  )+yvel0(j+1,k+1,l  )
        vgradz2=yvel0(j  ,k  ,l+1)+yvel0(j+1,k  ,l+1)+yvel0(j  ,k+1,l+1)+yvel0(j+1,k+1,l+1)

        wgradx1=zvel0(j  ,k  ,l  )+zvel0(j  ,k+1,l  )+zvel0(j  ,k  ,l+1)+zvel0(j  ,k+1,l+1)
        wgradx2=zvel0(j+1,k  ,l  )+zvel0(j+1,k+1,l  )+zvel0(j+1,k  ,l+1)+zvel0(j+1,k+1,l+1)
        wgrady1=zvel0(j  ,k  ,l  )+zvel0(j+1,k  ,l  )+zvel0(j  ,k  ,l+1)+zvel0(j+1,k  ,l+1)
        wgrady2=zvel0(j  ,k+1,l  )+zvel0(j+1,k+1,l  )+zvel0(j  ,k+1,l+1)+zvel0(j+1,k+1,l+1)
        wgradz1=zvel0(j  ,k  ,l  )+zvel0(j+1,k  ,l  )+zvel0(j  ,k+1,l  )+zvel0(j+1,k+1,l  )
        wgradz2=zvel0(j  ,k  ,l+1)+zvel0(j+1,k  ,l+1)+zvel0(j  ,k+1,l+1)+zvel0(j+1,k+1,l+1)

        div = (xarea(j,k,l)*(ugradx2-ugradx1)+  yarea(j,k,l)*(vgrady2-vgrady1))+ zarea(j,k,l)*(wgradz2-wgradz1)

        xx = 0.25_8*(ugradx2-ugradx1)/(celldx(j))
        yy = 0.25_8*(vgrady2-vgrady1)/(celldy(k))
        zz = 0.25_8*(wgradz2-wgradz1)/(celldz(l))
        xy = 0.25_8*(ugrady2-ugrady1)/(celldy(k))+0.25_8*(vgradx2-vgradx1)/(celldx(j))
        xz = 0.25_8*(ugradz2-ugradz1)/(celldz(l))+0.25_8*(wgradx2-wgradx1)/(celldx(j))
        yz = 0.25_8*(vgradz2-vgradz1)/(celldz(l))+0.25_8*(wgrady2-wgrady1)/(celldy(k))

        pgradx=(pressure(j+1,k,l)-pressure(j-1,k,l))/(celldx(j)+celldx(j+1))
        pgrady=(pressure(j,k+1,l)-pressure(j,k-1,l))/(celldy(k)+celldy(k+1))
        pgradz=(pressure(j,k,l+1)-pressure(j,k,l-1))/(celldz(l)+celldz(l+1))

        pgradx2 = pgradx*pgradx
        pgrady2 = pgrady*pgrady
        pgradz2 = pgradz*pgradz

        limiter = (xx*pgradx2+yy*pgrady2+zz*pgradz2                    &
                +  xy*pgradx*pgrady+xz*pgradx*pgradz+yz*pgrady*pgradz) &
                / MAX(pgradx2+pgrady2+pgradz2,1.0e-16_8)

        IF ((limiter.GT.0.0).OR.(div.GE.0.0))THEN
          viscosity(j,k,l) = 0.0
        ELSE
          pgradx = SIGN(MAX(1.0e-16_8,ABS(pgradx)),pgradx)
          pgrady = SIGN(MAX(1.0e-16_8,ABS(pgrady)),pgrady)
          pgradz = SIGN(MAX(1.0e-16_8,ABS(pgradz)),pgradz)
          pgrad = SQRT(pgradx*pgradx+pgrady*pgrady+pgradz*pgradz)
          xgrad = ABS(celldx(j)*pgrad/pgradx)
          ygrad = ABS(celldy(k)*pgrad/pgrady)
          zgrad = ABS(celldz(l)*pgrad/pgradz)
          grad  = MIN(xgrad,ygrad,zgrad)
          grad2 = grad*grad

          viscosity(j,k,l)=2.0_8*density0(j,k,l)*grad2*limiter*limiter
        ENDIF

      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF


  IF(PRESENT(work_array1)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array1(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array2)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array2(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array3)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array3(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array4)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array4(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array5)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array5(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array6)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array6(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array7)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array7(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(work_array8)) THEN
!$OMP DO
    DO l=z_min-2,z_max+3
      DO k=y_min-2,y_max+3
        DO j=x_min-2,x_max+3
          work_array8(j,k,l)=0.0_8
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vol_flux_x)) THEN
!$OMP DO
    DO l=z_min,z_max
      DO k=y_min,y_max
        DO j=x_min,x_max+1 
          vol_flux_x(j,k,l)=0.25_8*dt*xarea(j,k,l)                  &
                         *(xvel0(j,k,l)+xvel0(j,k+1,l)+xvel1(j,k,l)+xvel1(j,k+1,l))
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vol_flux_y)) THEN
!$OMP DO
    DO l=z_min,z_max
      DO k=y_min,y_max+1
        DO j=x_min,x_max
          vol_flux_y(j,k,l)=0.25_8*dt*yarea(j,k,l)                  &
                         *(yvel0(j,k,l)+yvel0(j+1,k,l)+yvel1(j,k,l)+yvel1(j+1,k,l))
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(vol_flux_z)) THEN
!$OMP DO
    DO l=z_min,z_max+1
      DO k=y_min,y_max
        DO j=x_min,x_max
          vol_flux_z(j,k,l)=0.25_8*dt*zarea(j,k,l)                  &
                         *(zvel0(j,k,l)+zvel0(j,k,l+1)+zvel1(j,k,l)+zvel1(j,k,l+1))
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(mass_flux_x)) THEN
!$OMP DO
    DO l=z_min-2,z_max+2
      DO k=y_min,y_max
        DO j=x_min,x_max+1 
          ! j-1 could be j, depending on the flow
          mass_flux_x(j,k,l)=vol_flux_x(j,k,l)*density1(j-1,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(mass_flux_y)) THEN
!$OMP DO
    DO l=z_min,z_max
      DO k=y_min,y_max+1
        DO j=x_min,x_max
          ! k-1 could be k, depending on the flow
          mass_flux_y(j,k,l)=vol_flux_y(j,k,l)*density1(j,k-1,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

  IF(PRESENT(mass_flux_z)) THEN
!$OMP DO
    DO l=z_min,z_max+1
      DO k=y_min,y_max
        DO j=x_min,x_max
          ! l-1 could be l, depending on the flow
          mass_flux_z(j,k,l)=vol_flux_z(j,k,l)*density1(j,k-1,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP ENDDO
  ENDIF

!$OMP END PARALLEL

END SUBROUTINE set_data

END MODULE set_data_module

