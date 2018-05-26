!!****f* source/physics/TreeCol/TreeCol_bhSelfContrib
!!
!! NAME
!!
!!  TreeCol_bhSelfContrib
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

subroutine TreeCol_bhSelfContrib(cellsize, blockno, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(MDIM), intent(IN) :: cellsize
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  return
end subroutine TreeCol_bhSelfContrib

