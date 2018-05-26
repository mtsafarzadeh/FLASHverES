!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_postRecvMsg
!!
!! NAME
!!  
!!  pl_postRecvMsg
!!
!! SYNOPSIS
!! 
!!  call pl_postRecvMsg (integer, intent (in) :: channel)
!!
!! DESCRIPTION
!!
!!  Posts a receive message for the specified channel.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to post the receive
!!
!!***

subroutine pl_postRecvMsg (channel)

  use Pipeline_data,     ONLY : pl_channelSize,    &
                                pl_comm,           &
                                pl_doLog,          &
                                pl_isRecvCommDone, &
                                pl_itemSize,       &
                                pl_logUnit,        &
                                pl_procList,       &
                                pl_recvBuf,        &
                                pl_recvRequest,    &
                                pl_tag
                                
  use Driver_interface,  ONLY : Driver_checkMPIErrorCode

  include "Flash_mpi.h"

  implicit none

  integer, intent (in) :: channel

  integer :: error
  integer :: procID
!
!
!     ...Post the receive and write this action to the log file (if requested).
!
!
  procID = pl_procList (channel)

  call MPI_Irecv (pl_recvBuf (1,1,channel),     &
                  pl_itemSize * pl_channelSize, &
                  FLASH_REAL,                   &
                  procID,                       &
                  pl_tag,                       &
                  pl_comm,                      &
                  pl_recvRequest (channel),     &
                  error                         )

  call Driver_checkMPIErrorCode (error)

  if (pl_doLog) then
      write (pl_logUnit,'(a,i6)') ' Posted receive message from proc ID ', procID
  end if
!
!
!     ...Set an indicator, if all receives were completed. At this stage, since
!        the receive has been posted but not checked for completion, this indicator
!        is false.
!
!    
  pl_isRecvCommDone = .false.
!
!
!    ...Ready!
!
!
  return
end subroutine pl_postRecvMsg
