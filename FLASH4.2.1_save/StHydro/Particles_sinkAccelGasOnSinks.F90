!!****f* source/Particles/Particles_sinkAccelGasOnSinks
!!
!! NAME
!!
!!  Particles_sinkAccelGasOnSinks
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelGasOnSinks()
!!
!! DESCRIPTION
!!
!!  Computes gas -> sinks gravitational accelerations by direct summation
!!  over all sink particles and grid cells.
!!  For cosmology, will also want to get contribution from PDE
!!  (mapped DM delegate particle density).
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged and renamed to reflect symmetry with Particles_sinkAccelSinksOnGas (Christoph Federrath, 2013)
!!
!!***

subroutine Particles_sinkAccelGasOnSinks()
  implicit none
end subroutine Particles_sinkAccelGasOnSinks
