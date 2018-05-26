!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-jeans/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer myPE)
!!
!! ARGUMENTS
!!
!!    myPE      Current Processor Number
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Also initializes initial conditions for StirTurb problem
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype

  implicit none
#include "Flash.h"
#include "Eos.h"
#include "constants.h"
#include "Flash_mpi.h"


  integer :: i, istat
  real :: xmax, xmin, ymax, ymin, zmax, zmin

  call Driver_getMype(GLOBAL_COMM, sim_MyPE)
  call Driver_getComm(GLOBAL_COMM, sim_Comm)

  call RuntimeParameters_get( 'sim_rho0' ,  sim_rho0 )
  call RuntimeParameters_get( 'sim_T0'   ,  sim_T0 )
  call RuntimeParameters_get( 'sim_delta',  sim_delta )
  call RuntimeParameters_get( 'sim_hx',     sim_hx )
  call RuntimeParameters_get( 'sim_hy',     sim_hy )
  call RuntimeParameters_get( 'sim_hz',     sim_hz )

  call RuntimeParameters_get( 'abar_1',  sim_abar_1 )
  call RuntimeParameters_get( 'gamma_1', sim_gamma_1 )

  call RuntimeParameters_get( 'smlrho', smlrho)
  call RuntimeParameters_get( 'smallp', smallp)
  call RuntimeParameters_get( 'smallX', smallX)
  call RuntimeParameters_get( 'xmax', xmax)
  call RuntimeParameters_get( 'xmin', xmin)
  call RuntimeParameters_get( 'ymax', ymax)
  call RuntimeParameters_get( 'ymin', ymin)
  call RuntimeParameters_get( 'zmax', zmax)
  call RuntimeParameters_get( 'zmin', zmin)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance1', sim_solutionErrorTolerance1)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance2', sim_solutionErrorTolerance2)
  sim_Lx = xmax - xmin
  sim_Ly = ymax - ymin
  sim_Lz = zmax - zmin
  
  call PhysicalConstants_get( 'Boltzmann', sim_boltz)
  call PhysicalConstants_get( 'proton mass', sim_mH)
  call PhysicalConstants_get( 'pi', sim_pi)
  call PhysicalConstants_get("Newton", newton)


  sim_vecLen = 1
  sim_mode = MODE_DENS_PRES

  sim_p0 = sim_rho0*sim_boltz*sim_T0/sim_abar_1/sim_mH


end subroutine Simulation_init
