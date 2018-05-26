!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_postSendMsg
!!
!! NAME
!!  
!!  pl_postSendMsg
!!
!! SYNOPSIS
!! 
!!  call pl_postSendMsg (integer, intent (in) :: channel)
!!
!! DESCRIPTION
!!
!!  Posts a sending message for the specified channel.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to post the send
!!
!!***

subroutine pl_postSendMsg (channel)

  use Pipeline_data,     ONLY : pl_channelSize,    &
                                pl_comm,           &
                                pl_doLog,          &
                                pl_isSendCommDone, &
                                pl_itemSize,       &
                                pl_logUnit,        &
                                pl_procList,       &
                                pl_sendBuf,        &
                                pl_sendCount,      &
                                pl_sendRequest,    &
                                pl_tag
                                
  use Driver_interface,  ONLY : Driver_checkMPIErrorCode

  include "Flash_mpi.h"

  implicit none

  integer, intent (in) :: channel

  integer :: error
  integer :: msgSize
  integer :: procID
!
!
!     ...Post the send and write this action to the log file (if requested).
!
!
  msgSize = pl_sendCount (channel)

  if (msgSize >= 0) then

      procID = pl_procList (channel)

      call MPI_Isend (pl_sendBuf (1,1,channel),     &
                      pl_itemSize * msgSize,        &
                      FLASH_REAL,                   &
                      procID,                       &
                      pl_tag,                       &
                      pl_comm,                      &
                      pl_sendRequest (channel),     &
                      error                         )

      call Driver_checkMPIErrorCode (error)

      if (pl_doLog) then
          write (pl_logUnit,'(a,i6)') ' Posted sending message from proc ID ', procID
      end if

  end if
!
!
!     ...Set an indicator, if all sends were completed. At this stage, since
!        the sending has been posted but not checked for completion, this indicator
!        is false.
!
!    
  pl_isSendCommDone = .false.
!
!
!    ...Ready!
!
!
  return
end subroutine pl_postSendMsg
