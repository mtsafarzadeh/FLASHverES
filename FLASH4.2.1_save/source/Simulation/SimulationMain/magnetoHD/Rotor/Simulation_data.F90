!!****if* source/Simulation/SimulationMain/magnetoHD/Rotor/Simulation_data
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
!!  Stores the local data for Simulation setup: Rotor
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma, sim_smallX, sim_smallRho, sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_xCtr, sim_yCtr, sim_zCtr, sim_Radius, sim_perturbZ

  logical, save :: sim_gCell, sim_killdivb

integer, save :: sim_meshMe
end module Simulation_data
