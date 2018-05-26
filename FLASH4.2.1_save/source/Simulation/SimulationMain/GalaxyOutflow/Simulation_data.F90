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
  real, save :: sim_cAmbient, sim_mach, sim_gamma, sim_tempAmbient
  real, save :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save :: sim_noise_amplitude
  real, save :: sim_amu,sim_boltz
  real, save :: sim_init_dust
  real, save :: sim_isothermal
  real, save :: sim_scalefactor
  real, save :: sim_dustinit
  logical, save :: sim_gCell
  integer, save :: sim_vecLen, sim_mode


end module Simulation_data
