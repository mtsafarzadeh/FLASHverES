!!****if* source/Simulation/SimulationMain/Life/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***
module Simulation_data

  implicit none

#include "constants.h"

  ! The total mass of the target:
  real, save :: sim_targetMass

  !! *** Runtime Parameters *** !!  
  integer, save :: sim_meshComm
  integer, save :: sim_geometry
  real,    save :: sim_targetRadius
  real,    save :: sim_targetHeight
  real,    save :: sim_targetOffset
  real,    save :: sim_targetZOffset
  
  real,    save :: sim_rhoTarg  
  real,    save :: sim_teleTarg 
  real,    save :: sim_tionTarg 
  real,    save :: sim_tradTarg 

  real,    save :: sim_rhoCham  
  real,    save :: sim_teleCham 
  real,    save :: sim_tionCham 
  real,    save :: sim_tradCham 
  real,    save :: sim_velxCham
  real,    save :: sim_velyCham

  real,    save :: sim_smallX
  real,    save :: sim_pulseLength
  real,    save :: sim_inputEnergy
  integer, save :: sim_ndiv
  logical, save :: sim_oneSpec

  character(len=MAX_STRING_LENGTH), save :: sim_targetGeom
  character(len=MAX_STRING_LENGTH), save :: sim_driverType

end module Simulation_data


