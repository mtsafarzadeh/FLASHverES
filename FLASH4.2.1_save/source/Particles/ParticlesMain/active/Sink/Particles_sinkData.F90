!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkData
!!
!! NAME
!!
!!    Particles_sinkData
!!
!! SYNOPSIS
!!
!!    Particles_sinkData()
!!
!! DESCRIPTION
!!
!!    Module to hold local variables and data types for sink particle unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!***

module Particles_sinkData

  implicit none

#include "Particles.h"
#include "Flash.h"

  public ::  MAX_MSGS, maxsinks, n_empty
  logical, save :: RunningParticles = .false.
  integer, save :: n_empty
  integer, save :: sinks_maxSinks
  logical, save :: useSinkParticles

  integer, parameter :: maxsinks = 2048
  integer, save, pointer, dimension(:)  :: NumParticlesPerBlock
  integer, parameter :: MAX_MSGS = 12
  integer, parameter :: nrep_pbc = 2

  real,  dimension(MAX_MSGS) :: send_buff, recv_buff

  integer, save :: ipx, ipy, ipz, ipvx, ipvy, ipvz, ipm, iptag
  integer, save :: ipblk, iplx, iply, iplz, iplx_old, iply_old, iplz_old, ipmdot, ipt, ipsc
  integer, save :: ipdtold, ipcpu, iold_pmass, ipraccr, ipmgas
  integer, save :: ipbflx, ipbfly, ipbflz

  integer, parameter :: pt_sinkParticleProps = NPART_PROPS

  ! particles_local and particles_global refer to 
  ! sink particles - the local list and the global list
  real, save, allocatable, dimension(:,:) :: particles_local
  real, save, allocatable, dimension(:,:) :: particles_global
  
  integer, save :: local_tag_number
  integer, save :: localnp, localnpf

end module
