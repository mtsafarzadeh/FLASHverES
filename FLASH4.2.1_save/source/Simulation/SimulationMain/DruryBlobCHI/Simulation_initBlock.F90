!!****if* source/Simulation/SimulationMain/Sod/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID, 
!!                       integer(IN) :: myPE)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!!  sim_rhoambient      Density in the ambient medium
!!  sim_rhoblob         Density in the blob
!!  sim_tempambient
!!  sim_tempblob
!!  sim_blobradius      Radius of the blob
!!  sim_velambient     fluid velocity
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Eos_data, ONLY: eos_singleSpeciesA
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  Use Eos_interface, ONLY: Eos, Eos_wrapped


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  

  integer :: iMax, jMax, kMax
  
  
 !!$ Arguments -------------------
  integer, intent(in) :: blockID! , myPE
  !!$------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ,ii
  integer :: iq, jq, kq
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: del
  real :: xx, xxL, xxR, yy, zz, lposn,dx,dy,dz,u,t,height, width, taperwidth
  real, allocatable, dimension(:) :: xCoord, yCoord, yCoordL, yCoordR, zCoord, xCoordL, xCoordR
  real :: enerZone, ekinZone, eintZone, velxzone, velyzone,velzzone,magxzone,magyzone,magzzone,rhozone,crayzone
  real :: preszone,gamezone,gamczone,lp
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData
  real:: factor, err, sim_fieldLoopRadius, radius, sim_yCtr,x1,x2
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC+1,GRID_KHI_GC+1) :: Az,Ax,Ay
#else
  real, allocatable, dimension(:,:,:) :: Az !,Ax,Ay
#endif
  integer, dimension(MDIM) :: axis

  ! dump some output to stdout listing the paramters
!  if (myPE == MASTER_PE) then
!1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!  endif

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX), stat=istat)
  allocate(xCoordL(sizeX),stat=istat)
  allocate(xCoordR(sizeX),stat=istat)
  allocate(yCoord(sizeY), stat=istat)
  allocate(yCoordL(sizeY), stat=istat)
  allocate(yCoordR(sizeY), stat=istat)
!  allocate(zCoord(sizeZ), stat=istat)

  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER,sim_gCell, zCoord, sizeZ)
  call Grid_getCellCoords(JAXIS, blockID, CENTER,sim_gCell, yCoord, sizeY)
  call Grid_getCellCoords(JAXIS, blockID, LEFT_EDGE,sim_gCell, yCoordL, sizeY)
  call Grid_getCellCoords(JAXIS, blockID, RIGHT_EDGE,sim_gCell, yCoordR, sizeY)
  call Grid_getCellCoords(IAXIS, blockID, CENTER,     sim_gCell, xCoord,  sizeX)
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE,  sim_gCell, xCoordL, sizeX)
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, sim_gCell, xCoordR, sizeX)


#ifndef FIXEDBLOCKSIZE
  if (NDIM == 2) then
!     allocate(Ax(sizeX+1,sizeY+1,1),stat=istat)
!     allocate(Ay(sizeX+1,sizeY+1,1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,1),stat=istat)
  elseif (NDIM == 3) then
!     allocate(Ax(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
!     allocate(Ay(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
     allocate(Az(sizeX+1,sizeY+1,sizeZ+1),stat=istat)
  endif
#endif

    sim_yCtr =1.5e23
     sim_fieldLoopRadius = 0.2e23
      dx = xCoord(LOW+1)-xCoord(LOW)

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


     lposn = sim_posn
     height= sim_chi-1.
     width = sim_fieldLoopRadius
     taperwidth = sim_fieldLoopRadius/4.

     k = 1
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1

#if NFACE_VARS > 1
           ! x Coord at cell corner

           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoordL(i)
           else
              x1 = xCoordR(i-1)
           endif

           ! y Coord at cell corner

           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoordL(j)
           else
              x2 = yCoordR(j-1)
           endif
#else
           ! x Coord at cell center
           if (i <=blkLimitsGC(HIGH,IAXIS)) then
              x1 = xCoord(i)
           else

              x1 = xCoord(i-1) + dx
           endif

           ! y Coord at cell center

           if (j <=blkLimitsGC(HIGH,JAXIS)) then
              x2 = yCoord(j)
           else
              x2 = yCoord(j-1) + dx
           endif
#endif
           ! define radius of the field loop

           radius = sqrt((x1-lposn)**2 + (x2-sim_yCtr)**2)

           if (radius <= 1.3*sim_fieldLoopRadius) then
!              Ax(i,j,k) = 0.
!              Ay(i,j,k) = 0.

              Az(i,j,k) = 5.e-9*(1.3*sim_fieldLoopRadius - radius)
!              Az(i,j,k) = 5.e-6*(1.3*sim_fieldLoopRadius - radius)
!              Az(i,j,k) = 1.e-5*(1.3*sim_fieldLoopRadius - radius)

           else

!              Ax(i,j,k) = 0.
!              Ay(i,j,k) = 0.
                    Az(i,j,k) = 0.

           endif
        enddo
     enddo



 
  !ES  compute the magnitude of the b_field  and the x and y bfield strengths
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
!     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCoord(i)

           presZone = sim_pamb
           radius = sqrt((xx-lposn)**2 + (yy-sim_yCtr)**2)
           rhozone  = sim_rhoamb

!           if (radius .lt. width) rhoZone=sim_rhoamb*sim_chi

           lp = 0.8e23
           crayzone = 2.*sim_pamb-(2.*sim_pamb)*xx/100./sim_xmax

           if (xx .gt. lp) crayzone = 2.*sim_pamb-(2.*sim_pamb)*lp/100./sim_xmax -(2.*sim_pamb)/sim_fieldLoopRadius*(xx-lp)
           if (xx .gt. lp+sim_fieldLoopRadius) crayzone = 0.

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           velxZone = 0.
           velyZone = 0.
           velzZone = 0.

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           

           ! specific internal energy
           eintZone = presZone/(sim_gamma-1.)/rhoZone

           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
	   call Grid_putPointdata(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)
	   call Grid_putPointdata(blockId, CENTER, CRAY_MSCALAR, EXTERIOR, axis, crayZone)
#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif
!ES      add the initial MHD stuff
!           call Grid_putPointdata(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, 1.e-6)
!           call Grid_putPointdata(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, 1.e-6)
!           call Grid_putPointdata(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, 0.)
!           call Grid_putPointdata(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, 0.)
!           call Grid_putPointdata(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, 4.e-6)

!ES        Cell face-centered variables for StaggeredMesh scheme
!#if NFACE_VARS > 0
!           call Grid_putPointdata(blockId, FACEX, MAG_FACE_VAR,  EXTERIOR, axis, 1.e-6)
!           call Grid_putPointdata(blockId, FACEY, MAG_FACE_VAR,  EXTERIOR, axis, 1.e-6)
!           call Grid_putPointdata(blockId, FACEZ, MAG_FACE_VAR,  EXTERIOR, axis, 0.)
!#endif
        enddo
     enddo
  enddo

!! Cleanup!  Must deallocate arrays
  call Eos_wrapped(MODE_DENS_EI,blkLimitsGC,blockID)

  deallocate(xCoord)
  deallocate(xCoordL)
  deallocate(xCoordR)
  deallocate(yCoord)
  deallocate(yCoordL)
  deallocate(yCoordR)
!  deallocate(zCoord)
#ifndef FIXEDBLOCKSIZE
  deallocate(Az)
#endif

  return
end subroutine Simulation_initBlock



