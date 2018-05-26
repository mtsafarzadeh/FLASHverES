!!****if* source/Simulation/SimulationMain/Sod/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(in) :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Sod shock tube
!!  problem
!!
!! ARGUMENTS
!!
!!  myPE - the local processor number
!!
!! PARAMETERS
!!
!!  sim_rhoambient    	Density in the ambient medium
!!  sim_rhoblob   	Density in the blob
!!  sim_tempambient    
!!  sim_tempblob   
!!  sim_blobradius 	Radius of the blob
!!  sim_velambient    fluid velocity
!!
!!***

subroutine Simulation_init(myPE)
  
  use Simulation_data
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"

  integer, intent(in) :: myPE

  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 

  call RuntimeParameters_get('xmax', sim_xmax)
  call RuntimeParameters_get('xmin', sim_xmin)
  
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call RuntimeParameters_get('sim_rhoblob', sim_rhoblob)
  call RuntimeParameters_get('sim_rhoambient', sim_rhoambient)
  
  call RuntimeParameters_get('sim_tempblob', sim_tempblob)
  call RuntimeParameters_get('sim_tempambient', sim_tempambient)
  
  call RuntimeParameters_get('sim_velambient', sim_velambient)
  call RuntimeParameters_get('sim_blobradius', sim_blobradius)

  call RuntimeParameters_get('gconst', sim_grv_const)

  call PhysicalConstants_get("proton mass",sim_amu)
  call PhysicalConstants_get("Boltzmann",sim_grv_boltz)
 
  call Logfile_stamp(myPE, "initializing blob in wind tunnel with subgrid model",  &
       "[Simulation_init]")
     

  sim_grv_const= - sim_grv_const

!   call RunTimeParameters_get("sim_shock_time", sim_shock_time) 
!   call RunTimeParameters_get("sim_cool_time", sim_cool_time)


end subroutine Simulation_init







