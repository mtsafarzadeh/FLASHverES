!!****if* source/Simulation/SimulationMain/Sod/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the Sod problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_rhoambient      Density in the ambient medium
!!  sim_rhoblob         Density in the blob
!!  sim_tempambient
!!  sim_tempblob
!!  sim_blobradius      Radius of the blob
!!  sim_velambient    fluid velocity 
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_rhoambient, sim_rhoblob, sim_tempambient, sim_tempblob,sim_blobradius,sim_velambient, sim_velframe
  real, save :: sim_gamma, sim_smallP, sim_smallX,sim_grv_boltz,sim_amu,sim_grv_const
  real, save :: sim_xmax, sim_xmin
  real, save :: sim_tcc,sim_shiftTime ! this the cloud crushing time and the time
                                      ! we last shifted the frame
  real, save :: sim_blobx, sim_bloby  ! this is average y and abs(x) of the blob
  real, save :: sim_square
  logical, save :: sim_gCell

!   real, save :: sim_xH,sim_xHP,sim_xHM,sim_xD,sim_xDP,sim_xDM
!    real, save :: sim_xHE,sim_xHEP,sim_xHEPP,sim_xH2,sim_xH2P
!    real, save :: sim_xHD,sim_xHDP,sim_xELEC

!  real, save :: sim_meta, sim_spin, sim_ang

!    real, save :: sim_cool_time

end module Simulation_data


