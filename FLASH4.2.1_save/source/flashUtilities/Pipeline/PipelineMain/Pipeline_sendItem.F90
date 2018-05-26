!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_sendItem
!!
!! NAME
!!  
!!  Pipeline_sendItem
!!
!! SYNOPSIS
!! 
!!  call Pipeline_sendItem (real,    intent (in)  :: item (:),
!!                          integer, intent (in)  :: procID,
!!                          logical, intent (out) :: isHandled)
!!
!! DESCRIPTION
!!
!!  The routine will try to add the item to the send buffer on the pipeline channel
!!  corresponding to the neighbor processor procID. This can be done, if there is
!!  no pending send on the send buffer and if there is enough space on the send buffer.
!!  This routine can be considered the feeder of the pipeline.
!!
!! ARGUMENTS
!!
!!  item          : the item (array of elements) to be added to the send buffer
!!  procID        : neighbor processor ID
!!  isHandled     : is true, if the item was successfully added to the buffer
!!
!! NOTES
!!
!!  If the procID corresponding to the item is not found within all the channels
!!  of the pipeline on the current MPI rank, then something must have gone wrong.
!!  Either the pipeline for the current MPI rank was set up incorrectly or procID
!!  does not correspond to a neighbor processor. In this case the only solution
!!  is to abort the program. 
!!
!!***

subroutine Pipeline_sendItem (item, procID,   isHandled)

  use Pipeline_data,     ONLY : pl_channelSize, &
                                pl_numChannels, &
                                pl_procList,    &
                                pl_sendBuf,     &
                                pl_sendCount,   &
                                pl_sendRequest, &
                                pl_sendState,   &
                                OPEN_STATE


  use Driver_interface,  ONLY : Driver_abortFlash

  use pl_interface,      ONLY : pl_progressSendComm, &
                                pl_postSendMsg

  include "Flash_mpi.h"

  implicit none

  real,    intent (in)  :: item (:)
  integer, intent (in)  :: procID
  logical, intent (out) :: isHandled

  integer :: channel, channelLocation
  integer :: i

  integer, parameter :: notFound = -1
!
!
!     ...Identify the channel on the pipeline corresponding to the item's procID.
!        Abort run if not found.
!
!
  channel = notFound

  do i = 1, pl_numChannels                        ! It may be necessary to change the
     if (pl_procList (i) == procID) then          ! pl_procList data structure to make
         channel = i                              ! the lookup faster (maybe pre-sorting
         exit                                     ! the list)
     end if
  end do

  if (channel == notFound) then
      call Driver_abortFlash ('[Pipeline_sendItem] ERROR: channel for item not found!')
  end if
!
!
!     ...If there is a pending send in our desired channel, we test all send channels.
!        Request values are reset to MPI_REQUEST_NULL when sends complete.
!
!
  if (pl_sendRequest (channel) /= MPI_REQUEST_NULL) then
      call pl_progressSendComm ()
  end if
!
!
!     ...We can safely add items to the send buffer if there is no pending send.
!        If the channel is full, post a sending message from that channel.
!
!
  if (pl_sendState   (channel) == OPEN_STATE       .and. &
      pl_sendRequest (channel) == MPI_REQUEST_NULL       ) then

      channelLocation = pl_sendCount (channel) + 1

      if (channelLocation > pl_channelSize) then
          call Driver_abortFlash ('[Pipeline_sendItem] ERROR: channel overflow!')
      end if

      pl_sendBuf   (:, channelLocation , channel) = item (:)
      pl_sendCount (                     channel) = channelLocation
          
      if (channelLocation == pl_channelSize) then
          call pl_postSendMsg (channel)
      end if

      isHandled = .true.
  else
      isHandled = .false.
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_sendItem
