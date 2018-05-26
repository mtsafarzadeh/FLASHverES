!!****if* source/Grid/localAPI/gr_hypreGridStatus
!!
!!  NAME 
!!
!!  gr_hypreGridStatus
!!
!!  SYNOPSIS
!!
!!  call gr_hypreGridStatus(integer, intent(IN):: blockCount,
!!                          integer, dimension(blockCount),intent(IN):: blockList)
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
!!   Stub implementation.
!!
!!***


!!REORDER(4): solnVec

subroutine gr_hypreGridStatus (blockCount, blockList)
  
  implicit none  
  
  integer, intent(IN):: blockCount
  integer, dimension(blockCount),intent(IN):: blockList

  
  return
  
end subroutine gr_hypreGridStatus
