!!****f* source/physics/TreeCol/TreeCol_bhMAC
!!
!! NAME
!!
!!  TreeCol_bhMAC
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

logical function TreeCol_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize2
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, intent(IN) :: blockno
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  TreeCol_bhMAC = .true.
  return
end function TreeCol_bhMAC


