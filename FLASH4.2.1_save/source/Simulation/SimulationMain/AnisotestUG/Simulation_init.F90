!!****if* source/Simulation/SimulationMain/Sod/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Sod shock tube
!!  problem
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

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_data, ONLY :  dr_tmax
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"


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
!  call RuntimeParameters_get('sim_velframe', sim_velframe)
  call RuntimeParameters_get('sim_bNormal', sim_bNormal)
  call RuntimeParameters_get('sim_byLeft', sim_byLeft)

  call RuntimeParameters_get('sim_square',sim_square)

!  call RuntimeParameters_get('gconst', sim_grv_const)

  call PhysicalConstants_get("proton mass",sim_amu)
  call PhysicalConstants_get("Boltzmann",sim_grv_boltz)

  sim_grv_const = - sim_grv_const

!  if (sim_velambient .gt. 0.) then
!  sim_tcc       = sqrt(sim_rhoblob/sim_rhoambient)*sim_blobradius/sim_velambient
!  else
!  sim_tcc = 1.
!  endif

!  sim_shiftTime = 0.
!  dr_tmax = dr_tmax*sim_tcc



end subroutine Simulation_init







