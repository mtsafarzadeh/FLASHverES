!!****if* source/Simulation/SimulationMain/StirTurb/Simulation_init
!!
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
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Simulation_interface, ONLY : st_getsold
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  real :: tinit
  real :: sold
#include "Flash.h"
#include "Eos.h"
#include "constants.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get( 'smallx', sim_smallX)
  call RuntimeParameters_get( 'rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get( 'c_ambient', sim_cAmbient)
  call RuntimeParameters_get( 'mach', sim_mach)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('writematrix',sim_writematrix)
  call RuntimeParameters_get('st_energy',sim_init_sten)
  call RuntimeParameters_get('mach',sim_rms_mach_target)

  sim_vecLen = 1
  sim_mode = MODE_DENS_PRES

!  Transition matrix stuff
   call RuntimeParameters_get('tinitial', tinit)
   call RuntimeParameters_get('npdfstart',npdfstart)
 
!ES  these are the # of timesteps between files 10, 20, 40, 80
!   steps_one   = 250
!   steps_two   = 500
!   steps_four  = 1000
!   steps_eight = 2000

! These are the values for the 256 runs 
!  steps_one   = 125
!   steps_two   = 250
!   steps_four  = 500
!   steps_eight = 1000

!   These are the short values
!    steps_one   = 15
!    steps_two   = 30
!    steps_four  = 60
!    steps_eight = 120

!   These are the shorter values
    steps_one   = 1
    steps_two   = 2
    steps_four  = 4
    steps_eight = 8



   n_one   = npdfstart/steps_one
   ntwo    = npdfstart/steps_two
   nfour = npdfstart/steps_four
   neight = npdfstart/steps_eight
   tone   = tinit
   ttwo    = tinit
   tfour = tinit
   teight = tinit
   call st_getsold(0.1,0.1,sold)
   print *,sold,sold*exp(-sold),'0.1'
   call st_getsold(5.1,5.1,sold)
   print *,sold,sold*exp(-sold),'5.1'
   call st_getsold(.32,5.1,sold)
   print *,sold,sold*exp(-sold),'.32'
   call st_getsold(.32,0.1,sold)
   print *,sold,sold*exp(-sold),'.32'


end subroutine Simulation_init
