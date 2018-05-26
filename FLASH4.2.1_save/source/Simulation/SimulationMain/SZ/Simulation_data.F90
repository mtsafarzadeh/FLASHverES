!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_data
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
!!  Store the simulation data for the SodSpherical problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_yangle     Angle made by diaphragm normal w/y-axis (deg)
!!  sim_shockpos   Point of intersection between the shock plane and the x-axis
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!
  !! Virial mass, radius, and temperature
  real, save :: sim_mvir, sim_rvir, sim_tvir
  !! concetration paramater, F(c) and 2 c/F(c)
  real, save :: sim_nfwc, sim_Fc, sim_Ac
  !!  central gas density
  real, save :: sim_rho0

  !! energy and radius of the QSO energy added to the sim 
  real, save :: sim_Eblast
  real, save :: sim_rblast
  real, save :: sim_mblast

integer, save :: sim_meshMe
end module Simulation_data


