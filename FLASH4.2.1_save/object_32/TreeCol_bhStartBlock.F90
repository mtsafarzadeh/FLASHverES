!!****f* source/physics/TreeCol/TreeCol_bhStartBlock
!!
!! NAME
!!
!!  TreeCol_bhStartBlock
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

subroutine TreeCol_bhStartBlock(blockno, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  return
end subroutine TreeCol_bhStartBlock

