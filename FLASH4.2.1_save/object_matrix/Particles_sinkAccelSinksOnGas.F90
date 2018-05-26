!!****f* source/Particles/Particles_sinkAccelSinksOnGas
!!
!! NAME
!!
!!  Particles_sinkAccelSinksOnGas
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelSinksOnGas(integer,intent(IN) :: blockcount,
!!              integer,dimension(blockCount),intent(IN)  :: blocklist)
!!
!! DESCRIPTION
!!
!!  Computes SGAX,SGAY,SGAY unk vars from sink particles. (sinks -> gas accelerations).
!!  Computes sinks -> gas gravitational accelerations by direct summation
!!  over all sink particles and grid cells.
!!
!! ARGUMENTS
!!
!!   blockcount - the number of blocks on this processor
!!
!!   blocklist - the list of blocks held by this processor
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged by Christoph Federrath, 2013
!!
!!***

subroutine Particles_sinkAccelSinksOnGas(blockCount,blockList)
  implicit none
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
end subroutine Particles_sinkAccelSinksOnGas
  
