!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2),intent(in):: pos, 
!!                      integer, intent(in) :: sweepDir, 
!!                      integer, intent(in) :: blockID, 
!!                      integer, intent(in) :: numCells, 
!!                      real(numCells),intent(out) :: grav, 
!!                      integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  test against sweep_x,
!!                                 sweep_y, and sweep_z
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!
!! NOTES
!!
!!  If certain variables declared by the sink particles inplementation are
!!  declared, it is assumed that sink particles are in use
!!  The sets of variables to make this determination are
!!    o  {SGXO_VAR, SGYO_VAR, SGZO_VAR}  if potentialIndex is GPOL_VAR
!!    o  {SGAX_VAR, SGAY_VAR, SGAZ_VAR}  otherwise.
!!  If it is assumed that sink particles are in used, then the acceleration
!!  returned in the grav array will have the appropriated sink particle
!!  acceleration component added to the acceleration computed by differencing
!!  the potential variable given by potentialIndex.
!!***

!!REORDER(4): solnVec

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, potentialIndex)


  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits
!  use Driver_interface, ONLY : Driver_abortFlash
  use Gravity_data, ONLY: grv_defaultGpotVar

  implicit none

#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling per-block drift logging for this routine because it is called too much
#undef Grid_releaseBlkPtr
#endif

#include "constants.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer, intent(IN),optional        :: potentialIndex
  real            :: blockSize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer         :: ii, iimin, iimax
  real            :: gpot(numCells), delxinv
  real, parameter :: onesixth = 1.e0/6.e0
  integer         :: potVar, nxbBlock, nybBlock, nzbBlock
  integer         :: sink_ax_index, sink_ay_index, sink_az_index

  !==================================================
  
  
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  
  
  call Grid_getBlkPtr(blockID, solnVec)


  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
  nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1


!! IF a variable index is explicitly specified, assume that as the potential
!! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else if (grv_defaultGpotVar > 0) then
     potVar=grv_defaultGpotVar
  else
     potVar=GPOT_VAR
  end if

  ! select the current or the old sink acceleration containers
  if (potVar .eq. GPOT_VAR) then
#if defined(SGAX_VAR) && defined(SGAY_VAR) && defined(SGAZ_VAR)
     sink_ax_index = SGAX_VAR
     sink_ay_index = SGAY_VAR
     sink_az_index = SGAZ_VAR
#endif
  endif
  if (potVar .eq. GPOL_VAR) then
#if defined(SGXO_VAR) && defined(SGYO_VAR) && defined(SGZO_VAR)
     sink_ax_index = SGXO_VAR
     sink_ay_index = SGYO_VAR
     sink_az_index = SGZO_VAR
#endif
  endif
  if ((potVar .ne. GPOT_VAR) .and. (potVar .ne. GPOL_VAR)) then
     print *, "Gravity_accelOneRow called with neither GPOT_VAR nor GPOL_VAR."
     call Driver_abortFlash("Gravity_accelOneRow called with neither GPOT_VAR nor GPOL_VAR.")
  endif


  iimin   = 1
  iimax   = numCells

  !Get row of potential values and compute inverse of zone spacing  
  if (sweepDir == SWEEP_X) then                     ! x-direction

     delxinv = real(nxbBlock) / blockSize(IAXIS)
     
     gpot(:) = solnVec(potVar,:,pos(1),pos(2))
#ifdef SGAX_VAR
     ! acceleration due to sink particles
     grav(iimin:iimax) = solnVec(sink_ax_index,:,pos(1),pos(2))
#else
     grav(iimin:iimax) = 0.0
#endif

  elseif (sweepDir == SWEEP_Y) then                 ! y-direction
     
     delxinv = real(nybBlock) / blockSize(JAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),:,pos(2))
#ifdef SGAY_VAR
     ! acceleration due to sink particles
     grav(iimin:iimax) = solnVec(sink_ay_index,pos(1),:,pos(2))
#else
     grav(iimin:iimax) = 0.0
#endif

  else                                            ! z-direction
     
     delxinv = real(nzbBlock) / blockSize(KAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),pos(2),:)
#ifdef SGAZ_VAR
     ! acceleration due to sink particles
     grav(iimin:iimax) = solnVec(sink_az_index,pos(1),pos(2),:)
#else
     grav(iimin:iimax) = 0.0
#endif

  endif
  
  !-------------------------------------------------------------------------------
  
  !               Compute gravitational acceleration
  
  
  !**************** first-order differences
  !                 preserves conservation
  
  delxinv = 0.5e0 * delxinv
  
  do ii = iimin+1, iimax-1
     grav(ii) = grav(ii) + delxinv * (gpot(ii-1) - gpot(ii+1))
  enddo
  
  grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
  grav(iimax) = grav(iimax-1)
  
  
  call Grid_releaseBlkPtr(blockID, solnVec)
  
  return
   
end subroutine Gravity_accelOneRow


