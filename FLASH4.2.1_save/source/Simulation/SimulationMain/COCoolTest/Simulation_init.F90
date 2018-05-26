!!****if* source/Simulation/SimulationMain/Marcus/Simulation_init
!!
!! NAME
!!
!!
!!  Simulation_init(integer(in) :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Marcus' cluster problem
!!  
!!
!! ARGUMENTS
!!
!!  myPE - the local processor number
!!
!! PARAMETERS
!!

!!
!!***

subroutine Simulation_init(myPE)


!!***used modules from FLASH.***
  
  use Simulation_data
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, RuntimeParameters_getPrev
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
#include "Flash.h"
#include "constants.h"
  integer, intent(in) :: myPE


!! These are the other variables we use
   real :: nblockx,nblocky,nblockz
   real :: sim_lrefine_max
 
!! These variables check to see if the Geometry is right 
  integer :: meshGeom
  logical :: validGeom

!!  This is a dummy variable
  real :: dummy


!!   ***make sure geometry is supported. I'm not sure it will run
!!      for 2D cartesian.***

     print *,'This is SimInit'
      call Grid_getGeometry(meshGeom)
     validGeom = (( meshGeom == CARTESIAN .AND. NDIM == 3) .OR.&
          &       ( meshGeom == CARTESIAN .AND. NDIM == 2))
     if (.NOT. validGeom) then
        print *, "ERROR: cluster only works for Cartesian 2/3-d  geometries"
        call Driver_abortFlash("ERROR: cluster only works for Cartesian 2/3-d geometries")
     endif


!! +++++Grab parameters from databases. Calculate some things.***
 

     print *,'Getting Physical Constants'
      
     print *,'Getting RuntimeParameters'
     call RuntimeParameters_get("xmin", sim_xMin)
     call RuntimeParameters_get("ymin", sim_yMin)
  
     call RuntimeParameters_get("xmax", sim_xMax)
     call RuntimeParameters_get("ymax", sim_yMax)
     call RunTimeParameters_get("sim_meta", sim_meta)
     call PhysicalConstants_get("ideal gas constant", sim_gasconst)


     print *,'Getting RuntimeParameters'
!    This is the number of blocks in the base grid
     call RuntimeParameters_get("sim_nblockx", nblockx)
     call RuntimeParameters_get("sim_nblocky", nblocky)

     if(NDIM == 3) then
        call RuntimeParameters_get("zmin", sim_zMin)
        call RuntimeParameters_get("zmax", sim_zMax)
        call RuntimeParameters_get("sim_nblockz", nblockz)
     endif

     print *,'Still getting RuntimeParameters'

     call RunTimeParameters_get("sim_c_temp", sim_c_temp)
     call RunTimeParameters_get("sim_c_den", sim_c_den)

	 
!!Here is where I will add my calls for my new shock parameters
    
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
 

     call RunTimeParameters_get("sim_shock_time", sim_shock_time) 




return

end subroutine Simulation_init






