!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_isDone
!!
!! NAME
!!  
!!  Pipeline_isDone
!!
!! SYNOPSIS
!! 
!!  call Pipeline_isDone (logical, intent (out) :: isDone)
!!
!! DESCRIPTION
!!
!!  The routine 
!!  smooth shutting down.
!!
!! ARGUMENTS
!!
!!  isDone : is set true, once the routine exits -> closing is in progress
!!
!!***

subroutine Pipeline_isDone (isDone)

  use Pipeline_data,     ONLY : pl_itemCount
  use pl_interface,      ONLY : pl_isCommDone

  include "Flash_mpi.h"

  implicit none

  logical, intent (out) :: isDone

  integer :: channel
!
!
!     ...Identify the channel on the pipeline corresponding to the item's procID.
!        Abort run if not found.
!
!
  if (pl_numChannels > 0) then

      do channel = 1, pl_numChannels
         if (pl_sendState (channel) == OPEN_STATE) then
             pl_sendState (channel) == PROMISE_TO_CLOSE_STATE
         end if
      end do

      call pl_progressSendComm ()

  end if

  isClosing = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_closeSendChannels
