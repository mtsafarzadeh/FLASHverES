!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBotNodeContrib
!!
!! NAME
!!
!!  gr_bhBotNodeContrib
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhBotNodeContrib(
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
!!  Calculates contribution of the bottom-most node. Calls *_bhBotNodeContrib 
!!  subroutines of physical modules (recently Gravity and TreeCol) for specific 
!!  calculation.
!!
!! ARGUMENTS
!!
!!  node        : array of the bottom-most (leaf) node of the tree, whose
!!                contribution is added to the gravitational potential
!!  ndSize      : physical size of the node (the largest extent of the grid cell)
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

subroutine gr_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
  use gr_bhData, ONLY : gr_bhOAAvg, gr_bhOAMin, gr_bhOAMax, gr_bhOACnt
  use Gravity_interface, ONLY : Gravity_bhBotNodeContrib
  use TreeCol_interface, ONLY : TreeCol_bhBotNodeContrib
  implicit none
#include "FortranLangFeatures.fh"
#include "constants.h"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize  
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  integer, dimension(2,MDIM), intent(IN) :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  real :: oa

  ! call _bhBotNodeContrib subroutines of physical modules
  call Gravity_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
  call TreeCol_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)

  ! calculate opening angle - for monitoring
  oa = ndSize*dr(MDIM+2)
  gr_bhOAAvg = gr_bhOAAvg + oa
  if (oa < gr_bhOAMin) gr_bhOAMin = oa
  if (oa > gr_bhOAMax) gr_bhOAMax = oa
  gr_bhOACnt = gr_bhOACnt + 1
  
  return
end subroutine gr_bhBotNodeContrib


