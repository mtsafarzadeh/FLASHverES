!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_data
!!
!! NAME
!!
!!  Pipeline_data
!!
!! SYNOPSIS
!!
!!  use Pipeline_data
!!  
!! DESCRIPTION
!!
!!  Data module for a Pipeline
!!  --------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (I) means data that is supplied as input at pipeline initialization (in argument list)
!!           (C) means data that is (C)alculated internally by the pipeline code
!!
!!   pl_channelSize        (R) : number of items each pipeline channel can hold
!!   pl_comm               (G) : the communicator within which the pipeline is defined
!!   pl_doLog              (C) : is set true, if a pipeline log file is requested
!!   pl_isCommDone
!!   pl_isCommInitialized
!!   pl_isInitialized
!!   pl_isRecvCommDone
!!   pl_isSendCommDone
!!   pl_itemBuf
!!   pl_itemCount          (C) : counts the total items processed by the pipeline at each stage
!!   pl_itemSize           (I) : the size (number of reals) defining each item
!!   pl_logUnit            (I) : the unit number for the pipeline log file
!!   pl_maxItems           (I) : the maximum number of items that the pipeline can handle
!!   pl_numChannels        (I) : the number of channels of the pipeline for the current rank 
!!   pl_procList           (I) : list of processors that will communicate with the current rank 
!!   pl_rank
!!   pl_recvBuf
!!   pl_recvCount
!!   pl_recvIndex
!!   pl_recvRequest
!!   pl_recvStatus
!!   pl_sendBuf
!!   pl_sendCount
!!   pl_sendIndex
!!   pl_sendRequest
!!   pl_sendState
!!   pl_sendStatus
!!   pl_size
!!   pl_tag                (P) : tag number
!!
!!***

Module Pipeline_data

  implicit none

  logical, save :: pl_doLog
  logical, save :: pl_isCommDone        = .false.
  logical, save :: pl_isCommInitialized = .false.
  logical, save :: pl_isInitialized     = .false.
  logical, save :: pl_isRecvCommDone
  logical, save :: pl_isSendCommDone

  integer, save :: pl_channelSize
  integer, save :: pl_comm
  integer, save :: pl_itemCount
  integer, save :: pl_itemSize
  integer, save :: pl_logUnit
  integer, save :: pl_maxItems
  integer, save :: pl_numChannels
  integer, save :: pl_rank
  integer, save :: pl_size

  integer, parameter :: pl_tag                 =  1235
  integer, parameter :: OPEN_STATE             = -3000
  integer, parameter :: PROMISE_TO_CLOSE_STATE = -4000
  integer, parameter :: WAITING_TO_CLOSE_STATE = -5000
  integer, parameter :: CLOSE_STATE            = -6000

  integer, save, allocatable :: pl_procList    (:)
  integer, save, allocatable :: pl_recvCount   (:)
  integer, save, allocatable :: pl_recvIndex   (:)
  integer, save, allocatable :: pl_recvRequest (:)
  integer, save, allocatable :: pl_sendCount   (:)
  integer, save, allocatable :: pl_sendIndex   (:)
  integer, save, allocatable :: pl_sendRequest (:)
  integer, save, allocatable :: pl_sendState   (:)

  real,    save, allocatable :: pl_itemBuf     (:,:)
  integer, save, allocatable :: pl_recvStatus  (:,:)
  integer, save, allocatable :: pl_sendStatus  (:,:)

  real,    save, allocatable :: pl_recvBuf     (:,:,:)
  real,    save, allocatable :: pl_sendBuf     (:,:,:)

end Module Pipeline_data
