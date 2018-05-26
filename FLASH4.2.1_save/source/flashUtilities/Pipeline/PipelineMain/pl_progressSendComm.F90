!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_progressSendComm
!!
!! NAME
!!  
!!  pl_progressSendComm
!!
!! SYNOPSIS
!! 
!!  call pl_progressSendComm ()
!!
!! DESCRIPTION
!!
!!  Determines progress of send communications inside the Pipeline unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pl_progressSendComm ()

  use Pipeline_data,     ONLY : pl_isSendCommDone,     &
                                pl_numChannels,        &
                                pl_sendCount,          &
                                pl_sendIndex,          &
                                pl_sendRequest,        &
                                pl_sendStatus,         &
                                pl_sendState,          &
                                pl_procList,           &
                                pl_doLog,              &
                                pl_logUnit,            &
                                CLOSE_STATE,           &
                                WAITING_TO_CLOSE_STATE

  use Driver_interface,  ONLY : Driver_abortFlash,       &
                                Driver_checkMPIErrorCode

  use pl_interface,      ONLY : pl_progressClosePromise

  include "Flash_mpi.h"

  implicit none

  logical :: openSends

  integer :: nCompleted
  integer :: error
  integer :: channel
  integer :: i
  integer :: procID
!
!
!     ...Check for completed receives on all channels.
!
!
  if (pl_numChannels > 0 .and. .not.pl_isSendCommDone) then

      call pl_progressClosePromise ()  ! handle all those sending channels which
                                       ! have been promised to close. Their closing
                                       ! status will change to waiting.
!
!
!     ...Test all sending channels for completed sends.
!
!
      call MPI_Testsome (pl_numChannels, &    ! how many should be tested
                         pl_sendRequest, &    ! handles to be tested
                         nCompleted,     &    ! number of completed requests
                         pl_sendIndex,   &    ! 1d array of indices of operations that completed
                         pl_sendStatus,  &    ! 2d array of status objects for operations that completed
                         error           )

      call Driver_checkMPIErrorCode (error)

      do i = 1, nCompleted

         channel = pl_sendIndex  (i)              ! get the channel of the completed send
         procID  = pl_procList (channel)          ! the procID to where the send has been done

         pl_sendCount (channel) = 0               ! reset the channel item count to 0

         if (pl_sendState (channel) == WAITING_TO_CLOSE_STATE) then
             pl_sendState (channel) =  CLOSE_STATE
         end if

         if (pl_doLog) then
             write (pl_logUnit,'(a,i6)') ' Completed send to procID ',procID
         end if

      end do
!
!
!     ...Check for overall send completion and check for bad shutdown.
!
!
      pl_isSendCommDone = all (pl_sendState (1:pl_numChannels) == CLOSE_STATE)

      if (pl_isSendCommDone) then

          openSends = any ( pl_sendRequest (1:pl_numChannels) /= MPI_REQUEST_NULL .or. &
                            pl_sendCount   (1:pl_numChannels) /= 0                     )

          if (openSends) then
              call Driver_abortFlash ('[pl_progressSendComm] ERROR: bad shutdown!')
          end if
      end if
!
!
!    ...Ready!
!
!
  return
end subroutine pl_progressSendComm
