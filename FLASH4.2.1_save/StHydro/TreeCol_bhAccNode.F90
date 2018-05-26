!!****f* source/physics/TreeCol/TreeCol_bhAccNode
!!
!! NAME
!!
!!  TreeCol_bhAccNode
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

subroutine TreeCol_bhAccNode(subnode, accnode)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN)  :: subnode
  real, dimension(:), intent(INOUT) :: accnode
  return
end subroutine TreeCol_bhAccNode
