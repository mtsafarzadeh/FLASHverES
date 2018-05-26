!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_progressRecvComm
!!
!! NAME
!!  
!!  pl_progressRecvComm
!!
!! SYNOPSIS
!! 
!!  call pl_progressRecvComm ()
!!
!! DESCRIPTION
!!
!!  Determines progress of receive communications inside the Pipeline unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine pl_progressRecvComm ()

  use Pipeline_data,     ONLY : pl_isRecvCommDone,    &
                                pl_numChannels,       &
                                pl_recvCount,         &
                                pl_recvIndex,         &
                                pl_recvRequest,       &
                                pl_recvStatus,        &
                                pl_procList,          &
                                pl_itemSize,          &
                                pl_doLog,             &
                                pl_logUnit

  use Driver_interface,  ONLY : Driver_abortFlash,       &
                                Driver_checkMPIErrorCode

  use pl_interface,      ONLY : pl_handleOldRecvMsg, &
                                pl_postRecvMsg,      &
                                pl_saveRecvItems

  include "Flash_mpi.h"

  implicit none

  logical :: isSaved

  integer :: nCompleted
  integer :: error
  integer :: i
  integer :: channel
  integer :: nElements, nItems, nItemsPerElement
  integer :: procID
!
!
!     ...Check for completed receives on all channels.
!
!
  if (pl_numChannels > 0 .and. .not.pl_isRecvCommDone) then

      call pl_handleOldRecvMsg ()  ! Copy old items from receive channels in which there
                                   ! are no pending MPI receives.  This only happens when
                                   ! there was previously insufficient space in pl_itemBuf.
!
!
!     ...Test all receive channels for new messages.  Save the corresponding items
!        and then post a new receive.
!
!
      call MPI_Testsome (pl_numChannels, &    ! how many should be tested
                         pl_recvRequest, &    ! handles to be tested
                         nCompleted,     &    ! number of completed requests
                         pl_recvIndex,   &    ! 1d array of indices of operations that completed
                         pl_recvStatus,  &    ! 2d array of status objects for operations that completed
                         error           )

      call Driver_checkMPIErrorCode (error)

      nItemsPerElement = pl_itemSize

      do i = 1, nCompleted

         channel = pl_recvIndex  (i)              ! get the channel of the completed receive
         procID  = pl_recvStatus (MPI_SOURCE,i)   ! get the procID from the source handle in status field

         if (procID /= pl_procList (channel)) then
             call Driver_abortFlash ('[pl_progressRecvComm] ERROR: procID mismatch!')
         end if

         call MPI_Get_count (pl_recvStatus (:,i), &   ! use the handles of the status field
                             FLASH_REAL,          &   ! to find out how many items (reals) were
                             nItems,              &   ! received from the completed receive
                             error                )   ! operation

         call Driver_checkMPIErrorCode (error)

         nElements = nItems / nItemsPerElement        ! each element can contain many items
         pl_recvCount (channel) = nElements           ! store number of elements received in channel

         if (pl_doLog) then
             write (pl_logUnit,'(2(a,i6))') ' Received ',nElements,' elements from procID ',procID
         end if

         if (nItems > 0) then
             call pl_saveRecvItems   (channel,  isSaved)
             if (isSaved) then
                 call pl_postRecvMsg (channel)
             end if
         end if

      end do
!
!
!     ...Check for overall completion.
!
!
      pl_isRecvCommDone = all ( pl_recvRequest (1:pl_numChannels) == MPI_REQUEST_NULL .and. &
                                pl_recvCount   (1:pl_numChannels) == 0                      )
!
!
!    ...Ready!
!
!
  return
end subroutine pl_progressRecvComm
