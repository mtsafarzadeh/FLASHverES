!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_progressClosePromise
!!
!! NAME
!!  
!!  pl_progressClosePromise
!!
!! SYNOPSIS
!! 
!!  call pl_progressClosePromise ()
!!
!! DESCRIPTION
!!
!!  Fulfills the close promise on all those channels in such a state by sending a zero-byte
!!  notification message.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pl_progressClosePromise ()

  use Pipeline_data,     ONLY : pl_numChannels,         &
                                pl_sendCount,           &
                                pl_sendRequest,         &
                                pl_sendState,           &
                                pl_doLog,               &
                                pl_logUnit,             &
                                PROMISE_TO_CLOSE_STATE, &
                                WAITING_TO_CLOSE_STATE

  use pl_interface,      ONLY : pl_postSendMsg

  include "Flash_mpi.h"

  implicit none

  integer :: channel
!
!
!     ...Check all channels for closing promises.
!
!
  do channel = 1, pl_numChannels

     if ( pl_sendState   (channel) == PROMISE_TO_CLOSE_STATE .and. &
          pl_sendRequest (channel) == MPI_REQUEST_NULL             ) then

          pl_sendCount (channel) = 0        ! for a zero-byte message

          call pl_postSendMsg (channel)     ! send the zero-byte message

          pl_sendState (channel) = WAITING_TO_CLOSE_STATE

     end if

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pl_progressClosePromise
