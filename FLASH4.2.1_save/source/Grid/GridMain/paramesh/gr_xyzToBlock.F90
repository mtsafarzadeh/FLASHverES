!!****if* source/Grid/GridMain/paramesh/gr_xyzToBlock
!!
!! NAME
!!  gr_xyzToBlock
!!
!! SYNOPSIS
!!
!!  gr_xyzToBlock( real, intent(IN): xyz,
!!                 integer, intent(OUT) : procID,
!!                 integer, intent(OUT) : blkID)
!!  
!! DESCRIPTION 
!!  
!!  This routine returns the identity of the block, on
!!  which the specified physical location falls. It also
!!  returns processor ID on which the block is residing
!!  
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_xyzToBlock(xyz, procID, blkID)
  use gr_interface, ONLY : gr_xyzToBlockLevel
  use Grid_data, ONLY : gr_meshNumProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lrefine_max
#ifdef BITTREE
  use bittree, only : amr_identify_block
#endif

  implicit none

  real, dimension(MDIM),intent(IN) :: xyz
  integer, intent(OUT) :: procID
  integer, intent(OUT) :: blkID

  integer, dimension(MDIM) :: ijk
  integer :: lev, proc, blk

  lev = lrefine_max

#ifdef BITTREE  
  call gr_xyzToBlockLevel(lev, xyz, ijk)
  call amr_identify_block(gr_meshNumProcs, lev, ijk, proc, blk)
#else
  call Driver_abortFlash("gr_xyzToBlock works only when bittree is enabled")
#endif

  procID=proc
  blkID=blk
End Subroutine gr_xyzToBlock
