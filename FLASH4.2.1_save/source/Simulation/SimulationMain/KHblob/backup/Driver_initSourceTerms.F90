!!****if* source/Driver/DriverMain/Driver_initSourceTerms
!!
!! NAME
!!   
!!  Driver_initSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_initSourceTerms(integer(in) :: myPE,
!!                         logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes all source terms Units by 
!!   calling their respective initialization routines
!!   viz. Stir_init, Burn_init, Heat_init, Cool_init, etc.
!!  
!! ARGUMENTS
!!   myPE - current processor number
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_initSourceTerms(myPE, restart)

  use Burn_interface, ONLY:  Burn_init
  use Stir_interface, ONLY : Stir_init
  use Heat_interface, ONLY : Heat_init
  use Cool_interface, ONLY : Cool_init
  use Diffuse_interface, ONLY : Diffuse_init

  implicit none
  
  integer, intent(in) :: myPE
  logical, intent(in) :: restart

  call Stir_init(myPE, restart)
  call Cool_init(myPE)
  call Diffuse_init(myPE)
  call Heat_init(myPE)
  call Turbulence_init(myPE)


!!$  call Ioniz_init(myPE)
  call Burn_init(myPE)

end subroutine Driver_initSourceTerms
