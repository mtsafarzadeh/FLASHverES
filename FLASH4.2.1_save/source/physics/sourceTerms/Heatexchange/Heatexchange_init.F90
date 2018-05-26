!!****f* source/physics/sourceTerms/Heatexchange/Heatexchange_init
!!
!! NAME
!!
!!  Heatexchange_init
!!
!! SYNOPSIS
!!
!!  call Heatexchange_init(logical, intent(IN)  :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes the Heatexchange unit.
!!
!! ARGUMENTS
!!
!!   restart : indicates whether the run is restarting from a checkpoint
!!
!!
!!
!!***

subroutine Heatexchange_init(restart)
  implicit none
  
  logical, intent(IN) :: restart


end subroutine Heatexchange_init
