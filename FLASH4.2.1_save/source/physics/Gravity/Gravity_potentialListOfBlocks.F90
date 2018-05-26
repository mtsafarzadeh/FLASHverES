!!****f* source/physics/Gravity/Gravity_potentialListOfBlocks
!!
!! NAME
!! 
!!  Gravity_potentialListOfBlocks 
!!
!! SYNOPSIS
!!
!!  Gravity_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                integer(IN) :: blockList(blockCount))
!!
!! DESCRIPTION
!!
!!  Compute the zone-averaged gravitational potential on all blocks
!!  specified in the list.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  gravitational potential.  Invokes a solver (of the Poisson equation)
!!  if necessary. On return,
!!     GPOT_VAR:  contains potential for the current simulation time.
!!     GPOL_VAR (if defined): contains potential at the previous simulation time.
!!
!!  May affect other variables related to particle properties if particles
!!  are included in the simulation.  In particular,
!!     PDEN_VAR (if defined): may get updated to the current density from
!!                particles if particles have mass.
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Gravity implementation.
!!  The following information is subject to change without notice.
!!  For the Multigrid implementation:
!!     ISLS_VAR (residual)
!!     ICOR_VAR (correction)
!!     IMGM_VAR (image mass)
!!     IMGP_VAR (image potential)
!!  For the Multipole implementation:
!!     (none)
!!
!! NOTES
!!
!!  This call does nothing for those Gravity implementations provided with
!!  FLASH that just provide a given (time-independent) potential for other
!!  physics modules and don't compute the potential by invoking a solver:
!!  Constant, PointMass, and PlanePar.
!!***

subroutine Gravity_potentialListOfBlocks(blockCount,blockList)

!=============================================================================
  implicit none


  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList

  return
end subroutine Gravity_potentialListOfBlocks
