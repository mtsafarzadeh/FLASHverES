!!****f* source/physics/TreeCol/TreeCol_init
!!
!! NAME
!!
!!  TreeCol_init
!!  
!! SYNOPSIS
!!
!!  TreeCol_init()
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables in the TreeCol unit, which are typically the 
!!  runtime parameters.  This routine must be called once by Driver_initFlash.F90. 
!!  Calling multiple times will not cause any harm but is unnecessary.
!!
!!  TreeCol unit will implement algorithm for estimating column densities developed
!!  by Clark et al. (2012, MNRAS, 420, 745). 
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!! 
!!  useTreeCol  BOOLEAN true  Controls turning on/off the compiled gravity unit
!!
!! NOTES
!!   
!!  Each implementation of TreeCol has its own runtime parameters.  Be sure to check
!!  the documentation or Config files to see them.
!!
!!***

subroutine TreeCol_init ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
   

  logical, save :: testUseTreeCol

!==============================================================================

  !! It is a failure to invoke the stub when useTreeCol is set TRUE.

  call RuntimeParameters_get ("useTreeCol", testUseTreeCol)
  if (testUseTreeCol) then
     call Driver_abortFlash("TreeCol unit seems not to be compiled in, and the TreeCol_init stub does not &
          &allow the value of useTreeCol to be TRUE.")
  end if

end subroutine TreeCol_init
