!!****f* source/physics/TreeCol/TreeCol_bhFillBotNode
!!
!! NAME
!!
!!  TreeCol_bhFillBotNode
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

subroutine TreeCol_bhFillBotNode(blockno, point, solnData, botnode)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(INOUT) :: botnode
  return
end subroutine TreeCol_bhFillBotNode
