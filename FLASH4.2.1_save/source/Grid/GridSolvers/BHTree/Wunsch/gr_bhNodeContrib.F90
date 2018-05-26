!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhNodeContrib
!!
!! NAME
!!
!!  gr_bhNodeContrib
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhNodeContrib(
!!                          real(in)       :: node(:),
!!                          real(in)       :: ndSize,
!!                          real(in)       :: dr(MDIM+2),
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Calculates contribution of the node during the tree walk by calling
!!  corresponding *_bhNodeContrib subroutines of physical modules.
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
!!                contribution is added to the gravitational potential
!!  ndSize      : physical size of the node (the largest extent of the node)
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine gr_bhNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
  use gr_bhData, ONLY : gr_bhOAAvg, gr_bhOAMin, gr_bhOAMax, gr_bhOACnt
  use Gravity_interface, ONLY : Gravity_bhNodeContrib
  use TreeCol_interface, ONLY : TreeCol_bhNodeContrib
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize  
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, dimension(2,MDIM), intent(IN) :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  real oa

  ! call _bhNodeContrib subroutines of physical modules
  call Gravity_bhNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
  call TreeCol_bhNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)

  ! calculate opening angle - for monitoring
  oa = ndSize*dr(MDIM+2)
  gr_bhOAAvg = gr_bhOAAvg + oa
  if (oa < gr_bhOAMin) gr_bhOAMin = oa
  if (oa > gr_bhOAMax) gr_bhOAMax = oa
  gr_bhOACnt = gr_bhOACnt + 1

  return
end subroutine gr_bhNodeContrib

