!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_closeSendChannels
!!
!! NAME
!!  
!!  Pipeline_closeSendChannels
!!
!! SYNOPSIS
!! 
!!  call Pipeline_closeSendChannels (logical, intent (out) :: isClosing)
!!
!! DESCRIPTION
!!
!!  The routine invokes a closing of all sending channels in the pipeline for
!!  smooth shutting down.
!!
!! ARGUMENTS
!!
!!  isClosing : is set true, once the routine exits -> closing is in progress
!!
!!***

subroutine Pipeline_closeSendChannels (isClosing)

  use Pipeline_data,     ONLY : pl_numChannels,        &
                                pl_sendState,          &
                                OPEN_STATE,            &
                                PROMISE_TO_CLOSE_STATE

  use pl_interface,      ONLY : pl_progressSendComm

  include "Flash_mpi.h"

  implicit none

  logical, intent (out) :: isClosing

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
