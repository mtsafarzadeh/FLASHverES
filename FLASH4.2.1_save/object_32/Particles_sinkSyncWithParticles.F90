!!****f* source/Particles/Particles_sinkSyncWithParticles
!!
!! NAME
!!
!!  Particles_sinkSyncWithParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSyncWithParticles(logical(in) :: sink_to_part)
!!
!! DESCRIPTION
!!
!!  Synchronizes global particle array with sink particle array.
!!
!! ARGUMENTS
!!
!!   sink_to_part -  logical flag indicating whether to sync with global
!!                   particle array or not
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   cleaned by Christoph Federrath, 2013
!!
!!***

subroutine Particles_sinkSyncWithParticles(sink_to_part)
  implicit none
  logical, intent(in) :: sink_to_part
end subroutine Particles_sinkSyncWithParticles
