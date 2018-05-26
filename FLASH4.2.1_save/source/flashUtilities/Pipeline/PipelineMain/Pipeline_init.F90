!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_init
!!
!! NAME
!!  
!!  Pipeline_init
!!
!! SYNOPSIS
!! 
!!  call Pipeline_init (integer, intent (in)           :: itemSize,
!!                      integer, intent (in)           :: maxItems,
!!                      integer, intent (in)           :: channelSize,
!!                      integer, intent (in)           :: comm,
!!                      integer, intent (in)           :: numChannels,
!!                      integer, intent (in)           :: procList (:), 
!!                      integer, intent (in), optional :: logUnit)
!!
!! DESCRIPTION
!!
!!  Initializes the Pipeline unit. Needs to be called whenever the FLASH grid changes or whenever
!!  the pipeline code is used for another purpose.
!!
!! ARGUMENTS
!!
!!  itemSize      : number of items
!!  maxItems      : maximum number of items
!!  channelSize   : channel size (number of items each channel can hold at a time)
!!  comm          : communicator handle
!!  numChannels   : number of channels
!!  procList      : list of processor ID's
!!  logUnit       : unit ID number for creating a log file
!!
!!***

subroutine Pipeline_init (itemSize, maxItems, channelSize, comm, numChannels, procList, logUnit)

  use Pipeline_data

  use Driver_interface,  ONLY : Driver_checkMPIErrorCode
  use pl_interface,      ONLY : pl_initComm

  include "Flash_mpi.h"

  implicit none

  integer, intent (in)           :: itemSize
  integer, intent (in)           :: maxItems
  integer, intent (in)           :: channelSize
  integer, intent (in)           :: comm
  integer, intent (in)           :: numChannels
  integer, intent (in)           :: procList (:)
  integer, intent (in), optional :: logUnit

  integer :: error
!
!
!     ...Copy the arguments to the internal variables.
!
!
  pl_itemSize    = itemSize
  pl_maxItems    = maxItems
  pl_channelSize = channelSize
  pl_comm        = comm
  pl_numChannels = numChannels
!
!
!     ...Prepare the log file data (if needed).
!
!    
  pl_doLog   = .false.
  pl_logUnit = -1

  if (present (logUnit)) then
      pl_doLog = (logUnit >= 0)
      if (pl_doLog) pl_logUnit = logUnit
  end if
!
!
!     ...Allocate the necessary arrays.
!
!    
  if (numChannels > 0) then

      allocate (pl_sendRequest (1:numChannels))
      allocate (pl_sendIndex   (1:numChannels))
      allocate (pl_sendCount   (1:numChannels))
      allocate (pl_sendState   (1:numChannels))
      allocate (pl_recvRequest (1:numChannels))
      allocate (pl_recvIndex   (1:numChannels))
      allocate (pl_recvCount   (1:numChannels))
      allocate (pl_procList    (1:numChannels))
      allocate (pl_itemBuf     (1:itemSize , 1:maxItems))
      allocate (pl_sendStatus  (1:MPI_STATUS_SIZE , 1:numChannels))
      allocate (pl_recvStatus  (1:MPI_STATUS_SIZE , 1:numChannels))
      allocate (pl_sendBuf     (1:itemSize , 1:channelSize , 1:numChannels))
      allocate (pl_recvBuf     (1:itemSize , 1:channelSize , 1:numChannels))

      pl_procList    (1:numChannels) = procList (1:numChannels)
      pl_sendRequest (1:numChannels) = MPI_REQUEST_NULL
      pl_recvRequest (1:numChannels) = MPI_REQUEST_NULL

  end if
!
!
!     ...Get rank and size of communicator.
!
!    
  call MPI_Comm_rank            (pl_comm, pl_rank, error)
  call Driver_checkMPIErrorCode (error)
  call MPI_Comm_size            (pl_comm, pl_size, error)
  call Driver_checkMPIErrorCode (error)
!
!
!     ...Set initialization indicator.
!
!    
  pl_isInitialized = .true.
!
!
!     ...Now initialize the communicator environment.
!
!    
  call pl_initComm ()
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_init
