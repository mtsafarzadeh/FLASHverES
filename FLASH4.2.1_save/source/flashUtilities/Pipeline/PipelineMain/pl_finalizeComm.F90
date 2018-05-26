!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_finalizeComm
!!
!! NAME
!!  
!!  pl_finalizeComm
!!
!! SYNOPSIS
!! 
!!  call pl_finalizeComm (logical, optional, intent (in) :: doAsyncReturn)
!!
!! DESCRIPTION
!!
!!  Finalizes the communicattion environment inside the Pipeline unit.
!!
!! ARGUMENTS
!!
!!  doAsyncReturn : if set false, each processor has to wait for all others before returning
!!
!!***

subroutine pl_finalizeComm (doAsyncReturn)

  use Pipeline_data,     ONLY : pl_comm,              &
                                pl_isCommDone,        &
                                pl_isCommInitialized, &
                                pl_isRecvCommDone,    &
                                pl_isSendCommDone,    &
                                pl_numChannels,       &
                                pl_recvRequest,       &
                                pl_recvStatus,        &
                                pl_sendRequest,       &
                                pl_sendState,         &
                                pl_sendStatus,        &
                                pl_size,              &
                                CLOSE_STATE

  use Driver_interface,  ONLY : Driver_checkMPIErrorCode

  include "Flash_mpi.h"

  implicit none

  logical, optional, intent (in) :: doAsyncReturn

  logical :: doSyncReturn

  integer :: channel
  integer :: error
  integer :: handle
!
!
!     ...Finalize all channels. Wait for all sends and receives to finish.
!
!
  if (pl_numChannels > 0) then

      do channel = 1, pl_numChannels

         handle = pl_recvRequest (channel)

         if (handle /= MPI_REQUEST_NULL) then
             call MPI_Cancel       (handle, error)
             call Driver_checkMPIErrorCode (error)
         end if

      end do

      call MPI_Waitall (pl_numChannels, &
                        pl_recvRequest, &
                        pl_recvStatus,  &
                        error           )

      call Driver_checkMPIErrorCode (error)

      pl_isRecvCommDone = .true.

      do channel = 1, pl_numChannels

         handle = pl_sendRequest (channel)

         if (handle /= MPI_REQUEST_NULL) then
             call MPI_Cancel       (handle, error)
             call Driver_checkMPIErrorCode (error)
         end if

      end do

      call MPI_Waitall (pl_numChannels, &
                        pl_sendRequest, &
                        pl_sendStatus,  &
                        error           )

      call Driver_checkMPIErrorCode (error)

      pl_sendState (1:pl_numChannels) = CLOSE_STATE

      pl_isSendCommDone = .true.

  end if
!
!
!     ...Allow all proccesors to return at will (asynchronous) ?
!
!
  if (pl_size > 1) then

      if (present (doAsyncReturn)) then
          doSyncReturn = .not. doAsyncReturn
      else
          doSyncReturn = .true.
      end if

      if (doSyncReturn) then
          call MPI_Barrier     (pl_comm, error)
          call Driver_checkMPIErrorCode (error)
      end if

  end if
!
!
!     ...Set current state of communication environment.
!
!    
  pl_isCommDone        = .true.
  pl_isCommInitialized = .false.
!
!
!    ...Ready!
!
!
  return
end subroutine pl_finalizeComm
