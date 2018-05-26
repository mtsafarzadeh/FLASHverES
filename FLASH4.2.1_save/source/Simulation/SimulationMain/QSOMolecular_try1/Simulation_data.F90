!!****if* source/Simulation/SimulationMain/StirTurbScalar/Simulation_data
!!  SS : change the above directory to StirTurbScalar instead of StirTurb
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!   
!!   use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: StirTurb
!! 
!!  SS : introduced scalars and magnetic field 
!!
!!***


module Simulation_data

  implicit none
#include "Eos.h"
#include "Flash.h"
  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!
  real, save :: sim_smallX
  real, save :: sim_rhoAmbient, sim_cAmbient, sim_mach, sim_gamma, sim_tempAmbient
  real, save :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save :: sim_noise_amplitude
  real, save :: sim_amu,sim_boltz
  real, save :: sim_isothermal
  logical, save :: sim_gCell
  integer, save :: sim_vecLen, sim_mode


!!Molecular Stuff 
  real, save :: sim_xH, sim_xHP, sim_xHM, sim_xH2, sim_xH2P, sim_xH3P, sim_xHE, sim_xHEP, sim_xC, sim_xCP
  real, save :: sim_xCM, sim_xO,sim_xOP, sim_xOM, sim_xC2, sim_xO2, sim_xO2P, sim_xOH, sim_xOHP, sim_xCO
  real, save :: sim_xCOP, sim_xCH, sim_xCHP, sim_xCH2, sim_xCH2P, sim_xHCOP,sim_xHOCP, sim_xH2O, sim_xH2OP
  real, save :: sim_xH3OP, sim_xCH3P, sim_xELEC

!!CHEMCOOL
!!    real, save :: sim_shock_time

end module Simulation_data
