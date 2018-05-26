!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_initComm
!!
!! NAME
!!  
!!  pl_initComm
!!
!! SYNOPSIS
!! 
!!  call pl_initComm ()
!!
!! DESCRIPTION
!!
!!  Initializes the communicattion environment inside the Pipeline unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  The Pipeline unit needs to be initialized before calling this routine.
!!
!!***

subroutine pl_initComm ()

  use Pipeline_data,     ONLY : pl_isCommDone,        &
                                pl_isCommInitialized, &
                                pl_isInitialized,     &
                                pl_isRecvCommDone,    &
                                pl_isSendCommDone,    &
                                pl_itemCount,         &
                                pl_numChannels,       &
                                pl_recvCount,         &
                                pl_sendCount,         &
                                pl_sendState

  use Driver_interface,  ONLY : Driver_abortFlash

  use pl_interface,      ONLY : pl_postRecvMsg

  implicit none

  integer :: channel
!
!
!     ...Abort, if pipeline unit has not been initialized.
!
!
  if (.not. pl_isInitialized) then
       call Driver_abortFlash ('[pl_initComm] ERROR: pipeline unit not initialized!')
  end if
!
!
!     ...Initialize the communication environment:
!
!         1) Sets all sending/receiving counts to 0
!         2) Sets the sending status of all channels to open
!            (waiting for the sending messages to be posted)
!         3) Posts all receiving messages on all channels
!            (waiting for the sending messages to arrive)
!
!
  pl_itemCount = 0

  if (pl_numChannels > 0) then

      pl_sendCount (1:pl_numChannels) = 0
      pl_recvCount (1:pl_numChannels) = 0

      do channel = 1, pl_numChannels
         pl_sendState (channel) = OPEN_STATE
         call pl_postRecvMsg (channel)
      end do

      pl_isCommDone     = .false.
      pl_isSendCommDone = .false.
      pl_isRecvCommDone = .false.
  else
      pl_isCommDone     = .true.
      pl_isSendCommDone = .true.
      pl_isRecvCommDone = .true.
  end if
!
!
!     ...Set communication environment initialization indicator.
!
!    
  pl_isCommInitialized = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine pl_initComm
