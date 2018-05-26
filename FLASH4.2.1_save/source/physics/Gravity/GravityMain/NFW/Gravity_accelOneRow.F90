!!****if* source/physics/Gravity/GravityMain/PointMass/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2)(IN) :: pos,
!!                      integer(IN)    :: sweepDir,
!!                      integer(IN)    :: blockID,
!!                      integer(IN)    :: numCells,
!!                      real(:)(OUT)   :: grav,
!!                      integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y and SWEEP_X. These values are defined
!!              in constants.h
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in pointmass gravity
!! 
!!***

subroutine Gravity_accelOneRow (pos,sweepDir,blockID, numCells, grav, &
      potentialIndex)

!=======================================================================

  use Gravity_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex

!==========================================================================


#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
#endif
  real :: dr32, tmpdr32
  real :: r, rtilde

  integer :: sizeX,sizeY,sizez

  integer :: ii,j,k
  logical :: gcell = .true.

!==============================================================================

  if (.NOT.useGravity) return

  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
  endif
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
  

  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells

        dr32 = sqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
        rtilde = dr32/grv_rs
        if(rtilde .ge. grv_nfwc) rtilde = grv_nfwc
        dr32 = (dr32+grv_rs*0.001)**3

        grav(ii) = grv_factor*xCenter(ii)/dr32*(alog(1.+rtilde)-rtilde/(1.+rtilde))
!        print *,'dr',dr32**(1./3.),grav(ii)*dr32**(1./3.)
!        print *,'ggg',grv_factor,dr32**(1./3.)/3.08E21,rtilde
!        print *,'grav',ii,grav(ii)
     enddo


  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
        rtilde = dr32/grv_rs
        if(rtilde .ge. grv_nfwc) rtilde = grv_nfwc
        dr32 = (dr32+grv_rs*0.001)**3

        grav(ii) = grv_factor*yCenter(ii)/dr32*(alog(1.+rtilde)-rtilde/(1.+rtilde))
     enddo

  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

     do ii = 1, numCells

        dr32 = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)
        rtilde = dr32/grv_rs
        if(rtilde .ge. grv_nfwc) rtilde = grv_nfwc
        dr32 = (dr32+grv_rs*0.001)**3        
        
        grav(ii) = grv_factor*zCenter(ii)/dr32*(alog(1.+rtilde)-rtilde/(1.+rtilde))
     enddo

  endif

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif

  return

end subroutine Gravity_accelOneRow
