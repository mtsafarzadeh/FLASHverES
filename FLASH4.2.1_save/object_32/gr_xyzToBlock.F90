!!****if* source/Grid/localAPI/gr_xyzToBlock
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

Subroutine gr_xyzToBlock(xyz, procID, blkID)

  implicit none

  real, dimension(MDIM),intent(IN) :: xyz
  integer, intent(OUT) :: procID
  integer, intent(OUT) :: blkID

  procID=0
  blkID=0  
End Subroutine gr_xyzToBlock
