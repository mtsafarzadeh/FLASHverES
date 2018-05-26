!!****if* source/Driver/DriverMain/Driver_finalizeSourceTerms
!!
!! NAME
!!   
!!  Driver_finalizeSourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_finalizeSourceTerms(integer(in) :: myPE,
!!                         logical(in) :: restart)
!!
!! DESCRIPTION
!!
!!   Finalizes all source terms Units by 
!!   calling their respective termination routines
!!   viz. Stir_finalize, Burn_finalize, Heat_finalize, Cool_finalize, etc.
!!  
!! ARGUMENTS
!!   myPE - current processor number
!!   restart - indicates if run is starting from scratch (.false.)
!!             or restarting from checkpoint (.true.)
!!
!!***

subroutine Driver_finalizeSourceTerms(myPE, restart)


  implicit none
  
  integer, intent(in) :: myPE
  logical, intent(in) :: restart


end subroutine Driver_finalizeSourceTerms
