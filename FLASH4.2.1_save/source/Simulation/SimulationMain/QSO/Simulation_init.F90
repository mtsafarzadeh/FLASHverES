!!****if* source/Simulation/SimulationMain/StirTurbScalar/Simulation_init
!!  SS : change the above directory to create a new StirTurbScalar 
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer sim_meshMe)
!!
!! ARGUMENTS
!!
!!    sim_meshMe      Current Processor Number
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Also initializes initial conditions for StirTurb problem
!!
!!  SS : added to get runtime parameters for the magnetic field
!!       Need to add a call to runtime parameters for the MHD solver
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use IO_interface, ONLY : IO_getScalar
  use Driver_data, ONLY : dr_restart
  implicit none
#include "Flash.h"
#include "Eos.h"
#include "constants.h"  

  call Driver_getMype(MESH_COMM, sim_meshMe)
! SS : added a call to the runtime xmin, xmax etc. 
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get('temp_ambient', sim_tempAmbient)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('isothermal_init',sim_isothermal)
  call RuntimeParameters_get('noise_amplitude', sim_noise_amplitude)

  call PhysicalConstants_get("proton mass",sim_amu)
  call PhysicalConstants_get("Boltzmann",sim_boltz)

  sim_vecLen = 1
  sim_mode = MODE_DENS_TEMP
  sim_gCell = .true.

end subroutine Simulation_init
