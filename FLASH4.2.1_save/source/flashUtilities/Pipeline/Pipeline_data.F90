!!****f* source/flashUtilities/Pipeline/Pipeline_data
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
!!           (C) means data that is (C)alculated internally by the pipeline code
!!
!!   pl_channelSize
!!   pl_comm                : pipeline communicator handle
!!   pl_doLog
!!   pl_isCommDone
!!   pl_isCommInitialized
!!   pl_isInitialized
!!   pl_isRecvCommDone
!!   pl_isSendCommDone
!!   pl_itemBuf
!!   pl_itemCount
!!   pl_itemSize
!!   pl_logUnit
!!   pl_maxItems
!!   pl_numChannels
!!   pl_procList
!!   pl_rank                : rank of processor in pipeline communicator
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
!!   pl_size                : size of pipeline communicator
!!   pl_tag             (P) : tag number
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
