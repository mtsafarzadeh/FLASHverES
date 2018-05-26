!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkSendOutputData
!!
!! NAME
!!
!!  Particles_sinkSendOutputData
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSendOutputData()
!!
!! DESCRIPTION
!!
!!  Calls pt_sinkGatherGlobal() and sets current globalNumSinkParticles.
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
!!
!!***

subroutine Particles_sinkSendOutputData()

  use Particles_sinkData, only: localnpf
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use IO_interface, only: IO_setScalar
  implicit none

  call pt_sinkGatherGlobal()
  call IO_setScalar("globalNumSinkParticles", localnpf)

end subroutine Particles_sinkSendOutputData
