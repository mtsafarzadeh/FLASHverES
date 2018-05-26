!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_handleOldRecvMsg
!!
!! NAME
!!  
!!  pl_handleOldRecvMsg
!!
!! SYNOPSIS
!! 
!!  call pl_handleOldRecvMsg ()
!!
!! DESCRIPTION
!!
!!  Saves extra received leftover items from a previous receive. This can only
!!  happen, if the item buffer size is too small.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pl_handleOldRecvMsg ()

  use Pipeline_data,     ONLY : pl_numChannels,       &
                                pl_recvCount,         &
                                pl_recvRequest

  use pl_interface,      ONLY : pl_postRecvMsg,      &
                                pl_saveRecvItems

  include "Flash_mpi.h"

  implicit none

  logical :: isSaved
  integer :: channel
!
!
!     ...Saves all extra received items in case there were some left. This only
!        happens when there was previously insufficient space in pl_itemBuf.
!
!
  do channel = 1, pl_numChannels

     if (pl_recvRequest (channel) == MPI_REQUEST_NULL .and. &
         pl_recvCount   (channel)  > 0                      ) then

         call pl_saveRecvItems   (channel,  isSaved)

         if (isSaved) then
             call pl_postRecvMsg (channel)
         end if
     end if

  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pl_handleOldRecvMsg
