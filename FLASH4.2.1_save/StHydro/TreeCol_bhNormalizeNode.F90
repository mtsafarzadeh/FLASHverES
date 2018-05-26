!!****f* source/physics/TreeCol/TreeCol_bhNormalizeNode
!!
!! NAME
!!
!!  TreeCol_bhNormalizeNode
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

subroutine TreeCol_bhNormalizeNode(smr, node)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(MDIM), intent(IN) :: smr
  real, dimension(:), intent(INOUT) :: node

  return
end subroutine TreeCol_bhNormalizeNode

