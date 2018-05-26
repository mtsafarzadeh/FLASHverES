!!****f* source/physics/TreeCol/TreeCol_bhAccBotNode
!!
!! NAME
!!
!!  TreeCol_bhAccBotNode
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

subroutine TreeCol_bhAccBotNode(blockno, point, blkLimits, solnData, botnode, accnode)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, dimension(:), intent(IN) :: botnode
  real, dimension(:), intent(INOUT) :: accnode
  return
end subroutine TreeCol_bhAccBotNode
