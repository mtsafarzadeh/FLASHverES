!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_finalize
!!
!! NAME
!!  
!!  Pipeline_finalize
!!
!! SYNOPSIS
!! 
!!  call Pipeline_finalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the Pipeline unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Pipeline_finalize ()

  use Pipeline_data, ONLY : pl_isCommInitialized, &
                            pl_isInitialized,     &
                            pl_numChannels

  use pl_interface,  ONLY : pl_finalizeComm

  implicit none

  logical :: doAsyncReturn
!
!
!     ...Finalize first the communication environment.
!
!
  if (pl_isCommInitialized) then
      call pl_finalizeComm (doAsyncReturn = .true.)
  end if
!
!
!     ...Deallocate all arrays.
!
!
  if (pl_isInitialized .and. pl_numChannels > 0) then

      deallocate (pl_itemBuf)
      deallocate (pl_sendBuf)
      deallocate (pl_sendStatus)
      deallocate (pl_sendRequest)
      deallocate (pl_sendIndex)
      deallocate (pl_sendCount)
      deallocate (pl_sendState)
      deallocate (pl_recvBuf)
      deallocate (pl_recvStatus)
      deallocate (pl_recvRequest)
      deallocate (pl_recvIndex)
      deallocate (pl_recvCount)
      deallocate (pl_procList)

      pl_isInitialized = .false.

  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_finalize
