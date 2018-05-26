!!****if* source/Grid/GridSolvers/HYPRE/UG/gr_hypreGridStatus
!!
!!  NAME 
!!
!!  gr_hypreGridStatus
!!
!!  SYNOPSIS
!!
!!  call gr_hypreGridStatus (integer, intent(IN):: blockCount,
!!                           integer, dimension(blockCount),intent(IN):: blockList)
!!
!!  DESCRIPTION 
!!      With AMR mesh verifies if the grid has been modified, if
!!      modified it resets the HYPRE grid object. If called first time 
!!      it sets up the HYPRE grid. 
!!
!! ARGUMENTS
!!
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated. 
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!   HYPRE grid is setup only once in UG. 
!!
!!***

!!REORDER(4): solnVec

#include "Flash.h"

subroutine gr_hypreGridStatus (blockCount, blockList)
  
  use gr_hypreData,     ONLY : gr_hypreGridIsSetUp
  
  implicit none
  
#include "constants.h" 
  
  integer, intent(IN):: blockCount  
  integer, dimension(blockCount),intent(IN):: blockList
  
  
  if (.not. gr_hypreGridIsSetUp) then 
     !! setup grid first time
     call gr_hypreSetupGrid (blockCount, blockList) 
  end if
  
  return
  
end subroutine gr_hypreGridStatus
