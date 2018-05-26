!!****if* source/physics/Gravity/GravityMain/PointMass/Gravity_accelListOfBlocks
!!
!! NAME
!!
!!  Gravity_accelListOfBlocks  
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelListOfBlocks(integer(IN) :: blockCount,
!!                         integer(IN)    :: blockList(blockCount),
!!                         integer(IN)    :: component,
!!                         integer(IN)    :: accelIndex,
!!                         integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!! ARGUMENTS
!!
!!   blockCount   - The number of blocks in the list
!!   blockList    - The list of blocks on which to calculate acceleration.
!!   component         - The component of the acceleration to compute.
!!                  Permitted values are IAXIS, JAXIS, KAXIS.  In
!!                     theory, ALLDIR will be implemented.  At the moment
!!                     this routine aborts in a messy way if called with ALLDIR.
!!   accelIndex - variable # to store the acceleration
!!   potentialIndex :   Variable # to take as potential if present, 
!!                     not applicable to pointmass
!!
!!***

subroutine Gravity_accelListOfBlocks (blockCount,blockList,component, &
     accelIndex, potentialIndex)

!==============================================================================

  use Gravity_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getCellCoords, Grid_releaseBlkPtr
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"

  integer,intent(IN)                      :: blockCount
  integer,dimension(MAXBLOCKS), intent(IN)     :: blockList
  integer, INTENT(in) ::  component
  integer, INTENT(in) ::  accelIndex
  integer,intent(IN),optional :: potentialIndex

!==============================================================================



  real          :: dr32, y2, z2
  real          :: rtilde,Frtilde
  integer       :: i, j, k, lb
  integer       :: csize
  logical       :: gcell = .true.
#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_IHI_GC) :: xCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_KHI_GC) :: zCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real, pointer :: solnVec(:,:,:,:)
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!==============================================================================
  if (.NOT.useGravity) return

  do lb = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     
#ifndef FIXEDBLOCKSIZE
     allocate(xCenter(blkLimitsGC(HIGH,IAXIS)))
     allocate(yCenter(blkLimitsGC(HIGH,JAXIS)))
     allocate(zCenter(blkLimitsGC(HIGH,KAXIS)))
#endif
     call Grid_getBlkPtr(blockList(lb), solnVec)
!------------------------------------------------------------------------------

     csize = blkLimitsGC(HIGH, IAXIS)
     call Grid_getCellCoords(IAXIS,blockList(lb),CENTER, gcell,xCenter,csize)
     yCenter = 0.
     zCenter = 0.
     if (NDIM >= 2) then
        csize = blkLimitsGC(HIGH, JAXIS)
        call Grid_getCellCoords(JAXIS,blockList(lb),CENTER,gcell,yCenter,csize)
     endif
     if (NDIM == 3) then 
        csize = blkLimitsGC(HIGH, KAXIS)
        call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, gcell, zCenter, csize)
     endif

!------------------------------------------------------------------------------
  
     if (component == IAXIS) then                    ! x-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           z2 = zCenter(k)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              y2 = yCenter(j)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                 dr32 = sqrt(xCenter(i)**2 + y2 + z2)
                 rtilde = dr32/grv_rs
                 if(rtilde .ge. grv_nfwc) rtilde = grv_nfwc
                 dr32 = dr32**3
                 solnVec(accelIndex, i, j, k) = grv_factor*xCenter(i)/dr32 &
                                               *(alog(1+rtilde)-rtilde/(1+rtilde))
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     Elseif (component == JAXIS) then                ! y-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           z2 = zCenter(k)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              y2 = yCenter(j)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                 dr32 = sqrt(xCenter(i)**2 + y2 + z2)
                 rtilde = dr32/grv_rs
                 if(rtilde .ge. grv_nfwc) rtilde = grv_nfwc
                 dr32 = (dr32+grv_rs*0.05)**3
                 solnVec(accelIndex, i, j, k) = grv_factor*yCenter(j)/dr32 &
                                               *(alog(1+rtilde)-rtilde/(1+rtilde))
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     elseif (component == KAXIS) then                ! z-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           z2 = zCenter(k)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              y2 = yCenter(j)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                 dr32    = sqrt(xCenter(i)**2 + y2 + z2)
                 rtilde  = dr32/grv_rs
                 Frtilde = (alog(1+rtilde)-rtilde/(1+rtilde))
                 if(rtilde .ge. grv_nfwc) Frtilde = grv_Fc
                 if(rtilde .le. 0.05)     Frtilde = rtilde**2/2.-2./3.*rtilde**3
                 dr32 = dr32**3
                 solnVec(accelIndex, i, j, k) = grv_factor*zCenter(k)/dr32 
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     else                                        ! ALLAXIS
        call Driver_abortFlash &
             ('[Gravity_accelListOfBlocks] ALLAXIS not supported!')
     endif

!------------------------------------------------------------------------------
     
#ifndef FIXEDBLOCKSIZE
    deallocate(xCenter)
    deallocate(yCenter)
    deallocate(zCenter)
#endif
    call Grid_releaseBlkPtr(blockList(lb),solnVec)
    
 enddo

!==============================================================================

 return

end subroutine Gravity_accelListOfBlocks
