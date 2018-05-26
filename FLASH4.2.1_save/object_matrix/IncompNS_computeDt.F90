!!****f* source/physics/IncompNS/IncompNS_computeDt
!!
!! NAME
!!
!!  IncompNS_computeDt
!!
!!
!! SYNOPSIS
!!
!!  
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***


subroutine IncompNS_computeDt(ins_mindt,ins_minloc)

  implicit none
  real,    intent(INOUT) :: ins_mindt
  integer, intent(INOUT) :: ins_minloc(5)

  real, PARAMETER :: MAX_TSTEP = huge(1.0)

  ! Set to large value:
  ins_mindt  = MAX_TSTEP
  ins_minloc(:) = 0

  return

end subroutine IncompNS_computeDt
