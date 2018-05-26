!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/RadTrans_computeDt
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  RadTrans_computeDt(integer(IN) :: blockID,
!!                     integer(IN) :: blkLimits(2,MDIM)
!!                     integer(IN) :: blkLimitsGC(2,MDIM)
!!                     real(IN),pointer::  solnData(:,:,:,:),   
!!                     real(OUT)   :: dt_radtrans, 
!!                     real(OUT)   :: dt_minloc(5)) 
!!  DESCRIPTION 
!!    Compute radiative transfer time step.
!!    Does NOTHING for leakage.  This exists because the version of this
!!    routine in RadTransMain is MGD-specific.  This is in complete 
!!    violation of FLASH coding standards but, hey, I didn't write MGD!
!!
!!  ARGUMENTS
!!    blockID       --  local block ID
!!    blkLimits     --  the indices for the interior endpoints of the block
!!    blkLimitsGC   --  the indices for endpoints including the guardcells
!!    solnData      --  the physical, solution data from grid
!!    dt_radtrans   --  variable to hold timestep constraint
!!    dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!
!!***
subroutine RadTrans_computeDt(blockID,  blkLimits,blkLimitsGC, &
     solnData, dt_radtrans, dt_minloc)
  implicit none

#include "constants.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: blkLimits(2,MDIM)
  integer, intent(IN) :: blkLimitsGC(2,MDIM)
  real, pointer :: solnData(:,:,:,:) 
  real, intent(INOUT) :: dt_radtrans
  integer, intent(INOUT)  :: dt_minloc(5)

  ! Stub implementation
  return

end subroutine RadTrans_computeDt
