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
  
  integer, intent(in) :: blockID

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  
  real :: xx, yy,  zz, xxL, xxR, dx
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone, eintZone, tempZone
  
  real :: A0,h0, radius
  logical :: gcell = .true.
  integer :: istat 

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData
  real :: smallx, ddelt, blobZone
  real :: bmag,bx,by,zero,magpZone 
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX),stat=istat)
  allocate(xRight(sizeX),stat=istat)
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0
  smallx = 1.0E-20
  ddelt = (sim_xmax-sim_xmin)/100.

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

!------------------------------------------------------------------------------

      dx = xCenter(LOW+1)-xCenter(LOW)
! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

      ! these are for setting L on either side of the boundary
      A0 = abs(sim_rhoambient-sim_rhoblob)/(sim_rhoambient+sim_rhoblob)
      
      h0 = 0.06*A0*0E-5*0E-5

     tempzone=sim_tempambient
     rhozone=sim_rhoambient
     velxzone=0. !sim_velambient
     velyzone= 0.
 
  !ES  compute the magnitude of the b_field  and the x and y bfield strengths
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)

           if (sim_square.ne.0.0) then
             radius=sqrt(xx*xx+zz*zz)
             radius=max(abs(yy),radius) 
           else
	     radius=sqrt(xx*xx+yy*yy+zz*zz)
           endif
          
           if (radius <= sim_blobradius) then
		
     	      rhoZone  = sim_rhoblob
              tempZone = sim_tempblob
              velxZone = 0.
              velyZone = 0.
              velzZone = 0.0E0 !
              blobZone = 1.
	   else
              rhoZone =  sim_rhoambient
              tempZone = sim_tempambient
              velyZone = sim_velambient 
              velxZone = 0.0E0 
              velzZone = 0.0E0 !
              blobZone = 0.00

           endif


           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           

           ! specific internal energy
           presZone= rhoZone*tempZone*sim_grv_boltz/sim_amu/eos_singleSpeciesA 
           eintZone = presZone/(sim_gamma-1.)/rhoZone

           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
	   call Grid_putPointdata(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)
	   call Grid_putPointdata(blockId, CENTER, BLOB_MSCALAR, EXTERIOR, axis, blobZone)

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
           call Grid_putPointdata(blockId, CENTER, MAGX_VAR, EXTERIOR, axis, sim_bx)
           call Grid_putPointdata(blockId, CENTER, MAGY_VAR, EXTERIOR, axis, sim_by)
           call Grid_putPointdata(blockId, CENTER, MAGZ_VAR, EXTERIOR, axis, zero)
           call Grid_putPointdata(blockId, CENTER, DIVB_VAR, EXTERIOR, axis, zero)
           call Grid_putPointdata(blockId, CENTER, MAGP_VAR, EXTERIOR, axis, sim_magp)

!ES        Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
           call Grid_putPointdata(blockId, FACEX, MAG_FACE_VAR,  EXTERIOR, axis, sim_bx)
           call Grid_putPointdata(blockId, FACEY, MAG_FACE_VAR,  EXTERIOR, axis, sim_by)
           call Grid_putPointdata(blockId, FACEZ, MAG_FACE_VAR,  EXTERIOR, axis, zero)
#endif
        enddo
     enddo
  enddo

!! Cleanup!  Must deallocate arrays
  call Eos_wrapped(MODE_DENS_EI,blkLimitsGC,blockID)

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

  return
end subroutine Simulation_initBlock



