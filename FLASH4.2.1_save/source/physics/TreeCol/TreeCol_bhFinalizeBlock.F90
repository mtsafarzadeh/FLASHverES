!!****f* source/physics/TreeCol/TreeCol_bhFinalizeBlock
!!
!! NAME
!!
!!  TreeCol_bhFinalizeBlock
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

subroutine TreeCol_bhFinalizeBlock(blockno, blkLimits, solnData)
  implicit none
#include "constants.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

end subroutine TreeCol_bhFinalizeBlock
