!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreGridStatus
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
  
  use gr_hypreData,      ONLY : gr_hypreGridIsSetUp, &
                                gr_hypreNStep
  use Driver_interface, ONLY : Driver_getNStep

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use tree, ONLY: grid_changed
  
  implicit none
  
#include "constants.h" 

  integer, intent(IN):: blockCount
  integer, dimension(blockCount),intent(IN):: blockList

  integer :: NStep

  call Timers_start ("gr_hypreGridStatus")

  
  if (.not. gr_hypreGridIsSetUp) then 
     call gr_hypreSetupGrid (blockCount, blockList)
     
     call Driver_getNStep(gr_hypreNStep)
  else
     if (grid_changed == 1) then
        
        call Driver_getNStep(NStep)
        
        if (Nstep /= gr_hypreNStep) then
           call gr_hypreDestroyGrid ()
           call gr_hypreSetupGrid (blockCount, blockList)                   
           gr_hypreNStep = NStep
        end if
        
     end if
  end if
 
  
  call Timers_stop ("gr_hypreGridStatus")

  
  return
  
end subroutine gr_hypreGridStatus
