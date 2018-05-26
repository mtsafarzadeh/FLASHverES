!!****if* source/Grid/GridSolvers/BHTree/Wunsch/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!!
!! SYNOPSIS
!!
!!   Grid_solvePoisson(integer(IN) :: ipotvar,
!!           integer(IN) :: idensvar, 
!!           integer(6)(IN) :: bcTypes,
!!           real(2,6)(IN) :: bcValues,
!!           real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the tree
!!   summation method for isolated problems.  Periodic problems are
!!   not supported; 
!!
!!
!! ARGUMENTS
!!
!!  ipotvar -  index to variable containing potential
!!  idensvar - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!             only used in verifying that they are isolated
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact -  factor to be used in calculation, for gravity: 4*PI*G
!!
!!***



subroutine Grid_solvePoisson (ipotvar, idensvar, bcTypes, bcValues, poisfact)

  use tree, ONLY: grid_changed
  use gr_bhData, ONLY: gr_bhBndType, gr_bhGravFac, gr_bhDensVar, gr_bhGpotVar, &
    gr_bhPhysMACTW, gr_bhPhysMACComm, gr_bhPhysMACTW_step, gr_bhPhysMACComm_step, &
    gr_bhUseRelAccErr, gr_bhFirstCall
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_bhInterface, ONLY : gr_bhComBlkProperties, &
    gr_bhBuildTree, gr_bhExchangeTrees, gr_bhTreeWalk, gr_bhDestroyTree, &
    gr_bhInitFieldVar

  implicit none
#include "constants.h"

  integer, intent(in)    :: ipotvar, idensvar
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact

  
  !=======================================================================
  

  call Timers_start("gravity_tree")

  gr_bhBndType = bcTypes
  gr_bhGravFac = poisfact / (4*PI)
  gr_bhDensVar = idensvar
  gr_bhGpotVar = ipotvar

  ! set MAC for the actual step (*_step) to values,
  ! to ensure that their change in Gravity_init was valid only for the first time step
  if (gr_bhFirstCall) then
    gr_bhFirstCall = .false.
    if (gr_bhUseRelAccErr) then
      ! set MAC for the actual step (*_step) to .false., because
      ! acceleration from the previous time-step is not known
      gr_bhPhysMACTW_step = .false.
      gr_bhPhysMACComm_step = .false.
    endif

  else
    gr_bhPhysMACTW_step = gr_bhPhysMACTW
    gr_bhPhysMACComm_step = gr_bhPhysMACComm
  endif
  

  call gr_bhInitFieldVar(ipotvar)

  if (grid_changed .eq. 1) then
    call gr_bhComBlkProperties()
  endif

  call gr_bhBuildTree()
  call gr_bhExchangeTrees()
  call gr_bhTreeWalk()
  call gr_bhDestroyTree()

  call Timers_stop("gravity_tree")
  !=========================================================================
  
  return
end subroutine Grid_solvePoisson
