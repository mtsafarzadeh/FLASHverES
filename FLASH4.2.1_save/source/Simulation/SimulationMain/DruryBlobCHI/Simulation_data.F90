!!****if* source/Simulation/SimulationMain/magnetoHD/BrioWu/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: BrioWu
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!    
  real, save   :: sim_posn, sim_gamma, sim_diff_time, sim_boltz
  real, save   :: sim_uRight, sim_uLeft, sim_vRight, sim_vLeft, sim_wLeft, sim_wRight
  real, save   :: sim_rhoamb, sim_rhocloud, sim_pamb
  real, save   :: sim_byLeft, sim_byRight, sim_bzLeft, sim_bzRight
  real, save   :: sim_xmin,sim_xmax,sim_ymin,sim_ymax,sim_zmin,sim_zmax
  real, save    :: sim_beta, sim_chi, sim_amu
  integer, save :: sim_igrid, sim_jgrid, sim_kgrid

  !! Simulation variables
  real, save    :: sim_xangle, sim_yangle, sim_xcos, sim_ycos, sim_zcos
  real, save    :: sim_smallx, sim_smallP
  logical, save :: sim_gCell, sim_killdivb

integer, save :: sim_meshMe

  real, save, dimension (64,128) :: posx,posy,randomgrid

end module Simulation_data
