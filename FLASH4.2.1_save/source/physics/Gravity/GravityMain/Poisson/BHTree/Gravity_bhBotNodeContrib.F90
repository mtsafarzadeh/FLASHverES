!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhBotNodeContrib
!!
!! NAME
!!
!!  Gravity_bhBotNodeContrib
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhBotNodeContrib(
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
!!  Calculates contribution of the bottom-most node to the gravitational
!!  potential and adds it to the solnData array (with index grv_defaultGpotVar).
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

subroutine Gravity_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
  use Gravity_data, ONLY : useGravity, grv_bhNewton, grv_bhUseEwald, &
    grv_defaultGpotVar, grv_bhIM
  use Gravity_interface, ONLY : Gravity_bhEwald
  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"
  real, dimension(:), intent(IN) :: node
  real, intent(IN) :: ndSize  
  integer, intent(IN) :: blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData

  if (.not. useGravity) return

  ! add the contribution to the potential
  if (grv_bhUseEwald) then

!!     print*,'[Gravity_bhBotNodeContrib]: &
!!  &  gr_bhBndTypes',gr_bhBndType(1),gr_bhBndType(3), gr_bhBndType(5)
!!     print*,'Calling Ewald field!',dr

    solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & - grv_bhNewton*node(grv_bhIM) * Gravity_bhEwald(abs(dr(IAXIS)), abs(dr(JAXIS)), abs(dr(KAXIS)))
  else
    solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & - grv_bhNewton*node(grv_bhIM)*dr(MDIM+2)
  endif

  return
end subroutine Gravity_bhBotNodeContrib
