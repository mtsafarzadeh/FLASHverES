!!****f* source/physics/TreeCol/TreeCol_bhNodeContrib
!!
!! NAME
!!
!!  TreeCol_bhNodeContrib
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

subroutine TreeCol_bhNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize  
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  return
end subroutine TreeCol_bhNodeContrib
