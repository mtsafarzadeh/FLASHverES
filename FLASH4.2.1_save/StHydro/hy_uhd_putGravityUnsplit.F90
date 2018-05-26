!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_putGravityUnsplit
!!
!! NAME
!!
!!  hy_uhd_PutGravityUnsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_putGravityUnsplit( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(2,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real   (IN) :: dt,
!!                            real   (IN) :: dtOld,
!!                            real(OUT)   :: gravX(:,:,:),
!!                            real(OUT)   :: gravY(:,:,:),
!!                            real(OUT)   :: gravZ(:,:,:))
!!
!! ARGUMENTS
!!
!!  blockID     - a current block ID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dataSize    - dimensions for gravX, gravY and gravZ arrays
!!  dt          - timestep
!!  dtOld       - old timestep (needed for temporal extrapolations of gravity)
!!  gravX       - gravity components in x-direcition at time steps n
!!  gravY       - gravity components in y-direcition at time steps n
!!  gravZ       - gravity components in z-direcition at time steps n
!!
!! DESCRIPTION
!!
!!  This routine puts gravity components to arrays gravX, gravY and gravZ:
!!  gravX(:,:,:) includes gravity components at time step n
!!
!!*** 

Subroutine hy_uhd_putGravityUnsplit&
     (blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)

  use Gravity_interface, ONLY : Gravity_accelOneRow

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
  integer, dimension(MDIM), intent(IN) :: dataSize
  real,    intent(IN) :: dt, dtOld

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(OUT) :: &
       gravX,gravY,gravZ
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(OUT) :: &
       gravX,gravY,gravZ
#endif
  !! -----------------------------------------------------

  integer, dimension(2) :: gravPos
  integer :: ix, iy, iz

  do iz=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do iy=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        gravPos(1)=iy
        gravPos(2)=iz
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
        ! Gravity implementation defines FLASH_GRAVITY_TIMEDEP -> time-dependent gravity field
        ! gravity at time step n
        call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz),GPOT_VAR)
#else
        ! FLASH_GRAVITY_TIMEDEP not defined -> assume time-independent gravity field.
        ! Also if GPOT_VAR is defined -> use current accel without time
        ! interpolation, i.e., handle like time-independent gravity field - KW
        call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),gravX(:,iy,iz))
#endif
     enddo
  enddo


  if (NDIM >= 2) then
     do iz=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        do ix=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           gravPos(1)=ix
           gravPos(2)=iz
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
           ! gravity at time step n
           call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz),GPOT_VAR)
#else
           call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(ix,:,iz))
#endif
        enddo
     enddo


     if (NDIM == 3) then
        do iy=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do ix=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              gravPos(1)=ix
              gravPos(2)=iy
#if defined(GPOT_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
              ! gravity at time step n
              call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:),GPOT_VAR)
#else
              call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(ix,iy,:))
#endif
           enddo
        enddo
     endif
  endif
End Subroutine hy_uhd_putGravityUnsplit
