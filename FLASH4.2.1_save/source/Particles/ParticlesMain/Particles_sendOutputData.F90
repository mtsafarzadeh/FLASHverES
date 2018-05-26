!!****if* source/Particles/ParticlesMain/Particles_sendOutputData
!!
!! NAME
!!  Particles_sendOutputData
!!
!! SYNOPSIS
!!
!!  call Particles_sendOutputData()
!!  
!! DESCRIPTION 
!!
!!  This routine allows the Particles unit to checkpoint any scalar data
!!  stored in the Particles_data Fortran modules and is called by the
!!  routine IO_updateScalars before checkpointing.  To send data to
!!  the IO unit this routine calls IO_setScalar.  In addition this
!!  routine may prepare any other data the IO unit needs for
!!  checkpointing which the Particles unit owns.
!!
!!
!!  ARGUMENTS  
!!
!!  SEE ALSO
!!   
!!   IO_setScalar, IO_updateScalars
!!
!!***

subroutine Particles_sendOutputData()

  use Particles_data, ONLY : pt_numLocal
  use Particles_interface, ONLY : Particles_getGlobalNum, Particles_sinkSendOutputData
  use IO_interface, ONLY : IO_setScalar
  implicit none

#include "Flash_mpi.h"

  integer :: globalNumParticles

  call Particles_getGlobalNum(globalNumParticles)
  call IO_setScalar("globalNumParticles", globalNumParticles)
  
  call Particles_sinkSendOutputData()
  
end subroutine Particles_sendOutputData
