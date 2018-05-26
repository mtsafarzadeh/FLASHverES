!!****if* source/Grid/localAPI/gr_hypreSetupGrid
!!
!!  NAME 
!!
!! gr_hypreSetupGrid
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetupGrid(blockCount, blockList)
!!
!!
!!  DESCRIPTION 
!! This routine sets up the HYPRE Grid. 
!! Called only once in UG.
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  
!!
!!***

!!REORDER(4): solnVec

subroutine gr_hypreSetupGrid (blockCount, blockList)  
  
  implicit none 
  
  integer,                      intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  
  
  
end subroutine gr_hypreSetupGrid
