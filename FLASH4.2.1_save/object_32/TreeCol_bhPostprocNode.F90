!!****f* source/physics/TreeCol/TreeCol_bhPostprocNode
!!
!! NAME
!!
!!  TreeCol_bhPostprocNode
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine TreeCol_bhPostprocNode(ndSize, node)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, intent(IN) :: ndSize
  real, dimension(:), intent(INOUT)  :: node

  return
end subroutine TreeCol_bhPostprocNode

