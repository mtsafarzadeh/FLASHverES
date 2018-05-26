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
  use Eos_interface, ONLY : Eos, Eos_wrapped
  use Driver_data, ONLY : dr_restart
  implicit none
#include "Flash.h"
#include "constants.h"  
#include "Multispecies.h"
#include "Eos.h"

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

!    Chemistry stuff is here

     call RunTimeParameters_get("sim_xH", sim_xH)
     call RunTimeParameters_get("sim_xHP", sim_xHP)
     call RunTimeParameters_get("sim_xHM", sim_xHM)
     call RunTimeParameters_get("sim_xH2", sim_xH2)
     call RunTimeParameters_get("sim_xH2P", sim_xH2P)
     call RunTimeParameters_get("sim_xH3P", sim_xH3P)
     call RunTimeParameters_get("sim_xHE", sim_xHE)
     call RunTimeParameters_get("sim_xHEP", sim_xHEP)
     call RunTimeParameters_get("sim_xC", sim_xC)
     call RunTimeParameters_get("sim_xCP", sim_xCP)
     call RunTimeParameters_get("sim_xCM", sim_xCM)
     call RunTimeParameters_get("sim_xO", sim_xO)
     call RunTimeParameters_get("sim_xOP", sim_xOP)
     call RunTimeParameters_get("sim_xOM", sim_xOM)
     call RunTimeParameters_get("sim_xC2", sim_xC2)
     call RunTimeParameters_get("sim_xO2", sim_xO2)
     call RunTimeParameters_get("sim_xO2P", sim_xO2P)
     call RunTimeParameters_get("sim_xOH", sim_xOH)
     call RunTimeParameters_get("sim_xOHP", sim_xOHP)
     call RunTimeParameters_get("sim_xCO", sim_xCO)
     call RunTimeParameters_get("sim_xCOP", sim_xCOP)
     call RunTimeParameters_get("sim_xCH", sim_xCH)
     call RunTimeParameters_get("sim_xCHP", sim_xCHP)
     call RunTimeParameters_get("sim_xCH2", sim_xCH2)
     call RunTimeParameters_get("sim_xCH2P", sim_xCH2P)
     call RunTimeParameters_get("sim_xHCOP", sim_xHCOP)
     call RunTimeParameters_get("sim_xHOCP", sim_xHOCP)
     call RunTimeParameters_get("sim_xH2O", sim_xH2O)
     call RunTimeParameters_get("sim_xH2OP", sim_xH2OP)
     call RunTimeParameters_get("sim_xH3OP", sim_xH3OP)
     call RunTimeParameters_get("sim_xCH3P", sim_xCH3P)
     call RunTimeParameters_get("sim_xELEC", sim_xELEC)
     !call RunTimeParameters_get("sim_shock_time", sim_shock_time)

end subroutine Simulation_init
