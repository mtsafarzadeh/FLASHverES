
module UTPipeline


  implicit none


contains



  subroutine UTPipeline_sendFullestChannel
    implicit none
    integer, parameter :: notFound = -1 !This value must be negative
    integer :: fullestChannel, bufSize, i

    if (utpipe_numChannels > 0) then
       if (any(utpipe_sendCount(:) > 0)) then
          fullestChannel = notFound
          bufSize = notFound
          do i = 1, utpipe_numChannels
             !Test for data that is not currently being sent
             if ( utpipe_sendState(i) == OPEN_STATE .and. &
                  utpipe_sendRequest(i) == MPI_REQUEST_NULL .and. &
                  utpipe_sendCount(i) > 0 ) then
                if (utpipe_sendCount(i) > bufSize) then
                   fullestChannel = i
                   bufSize = utpipe_sendCount(i)
                end if
             end if
          end do
          !It is possible that there are no sends meeting the above criteria
          if (fullestChannel >= 1 .and. fullestChannel <= utpipe_numChannels) then
             call utpipe_postSendMsg(fullestChannel)
          end if
       end if
    end if
  end subroutine UTPipeline_sendFullestChannel


  subroutine UTPipeline_progressComm(doFlush)
    implicit none
    logical, optional, intent(IN) :: doFlush
    call UTPipeline_progressRecvComm()
    call UTPipeline_progressSendComm()

    !The following code guarantees global progress.  It will normally
    !be executed when we are processing the last few items
    if (present(doFlush)) then
       if (doFlush .and. utpipe_itemCount == 0) then
          call UTPipeline_sendFullestChannel()
       end if
    end if
  end subroutine UTPipeline_progressComm


  !Rename UTPipeline_progressSendComm
  subroutine UTPipeline_progressSendComm()
    implicit none
    integer :: outcount, index, ierr, i

    if (utpipe_numChannels > 0 .and. .not.utpipe_isSendCommDone) then

       call utpipe_progressClosePromise()

       call MPI_Testsome(utpipe_NumChannels, utpipe_sendRequest, &
            outcount, utpipe_sendIndex, utpipe_sendStatus, ierr)
       ASSERT_MPI_SUCCESS(ierr)

       !Note that status objects are only meaningful for receive messages.
       do i = 1, outcount
          index = utpipe_sendIndex(i)

          utpipe_sendCount(index) = 0
          if (utpipe_sendState(index) == WAITING_TO_CLOSE_STATE) then
             utpipe_sendState(index) = CLOSE_STATE
          end if
          if (utpipe_doLog) then
             write(utpipe_logUnit,'(a,i6)') 'Completed send msg to ', &
                  utpipe_procList(index)
          end if
       end do

       !Check for completion
       utpipe_isSendCommDone = all(utpipe_sendState == CLOSE_STATE)
       if (utpipe_isSendCommDone .and. &
            any(utpipe_sendRequest /= MPI_REQUEST_NULL .or. &
            utpipe_sendCount /= 0)) then
          call Driver_abortFlash('Bad shutdown')
       end if
    end if
  end subroutine UTPipeline_progressSendComm





  !Promise to close the send channels.
  subroutine UTPipeline_closeSendChannels(isClosing)
    implicit none
    logical, intent(OUT) :: isClosing
    integer :: i
    if (utpipe_numChannels > 0) then
       do i = 1, utpipe_numChannels
          if (utpipe_sendState(i) == OPEN_STATE) then
             utpipe_sendState(i) = PROMISE_TO_CLOSE_STATE
          end if
       end do
       call UTPipeline_progressSendComm()
    end if
    isClosing = .true.
  end subroutine UTPipeline_closeSendChannels





  !Fulfill the close promise by sending a zero-byte notification message
  subroutine utpipe_progressClosePromise()
    implicit none
    integer :: i
    do i = 1, utpipe_numChannels
       if ( utpipe_sendState(i) == PROMISE_TO_CLOSE_STATE .and. &
            utpipe_sendRequest(i) == MPI_REQUEST_NULL ) then
          utpipe_sendCount(i) = 0 !For a zero-byte message
          call utpipe_postSendMsg(i)
          utpipe_sendState(i) = WAITING_TO_CLOSE_STATE
       end if
    end do
  end subroutine utpipe_progressClosePromise


  !Call this after UTPipeline_closeSendChannels
  subroutine UTPipeline_isCommDone(isCommDone)
    implicit none
    logical, intent(OUT) :: isCommDone

    if (.not.utpipe_isCommDone) then
       if (utpipe_isSendCommDone .and. utpipe_isRecvCommDone) then
          call UTPipeline_finalizeComm(doAsyncReturn=.true.)
       else
          call UTPipeline_progressComm()
       end if
    end if
    isCommDone = utpipe_isCommDone
  end subroutine UTPipeline_isCommDone


  subroutine UTPipeline_isDone(isDone)
    implicit none
    logical, intent(OUT) :: isDone
    logical :: isCommDone

    call UTPipeline_isCommDone(isCommDone)
    isDone = isCommDone .and. utpipe_itemCount == 0
  end subroutine UTPipeline_isDone


  subroutine UTPipeline_numItems(numItems)
    implicit none
    integer, intent(OUT) :: numItems
    numItems = utpipe_itemCount
  end subroutine UTPipeline_numItems



  subroutine UTPipeline_getItems(userArray, userMaxCount, userCount)
    implicit none
    real, dimension(:,:), intent(INOUT) :: userArray
    integer, intent(IN) :: userMaxCount
    integer, intent(INOUT) :: userCount
    integer :: freeSpace, itemsToCopy, firstItem

    freeSpace = userMaxCount - userCount
    itemsToCopy = min(utpipe_itemCount, freeSpace)
    if (itemsToCopy > 0) then
       !Copy from the end of utpipe_itemBuf to allow for a fast memcpy
       firstItem = utpipe_itemCount - itemsToCopy + 1

       if (utpipe_doLog) then
          write(utpipe_logUnit,'(2(a,2(i6)))') 'Copy from buf slice ', &
               firstItem, firstItem+itemsToCopy-1, ' to user slice ', &
               userCount+1, userCount+itemsToCopy
       end if

       userArray(:,userCount+1:userCount+itemsToCopy) = &
            utpipe_itemBuf(:,firstItem:firstItem+itemsToCopy-1)

       utpipe_itemCount = utpipe_itemCount - itemsToCopy
       userCount = userCount + itemsToCopy

       if (utpipe_doLog) then
          write(utpipe_logUnit,'(2(a,i6))') 'Elements in user array ', &
               userCount, ' Elements in buf ', utpipe_itemCount
       end if
    end if
  end subroutine UTPipeline_getItems





  !Caller should probably add the following: if (isHandled) item = NONEXISTENT
  subroutine UTPipeline_sendItem(item, procID, isHandled)
    implicit none
    real, dimension(:), intent(IN) :: item
    integer, intent(IN) :: procID
    logical, intent(OUT) :: isHandled
    integer :: channel, ptr, i
    integer, parameter :: notFound = -1

    !It may be necessary to change the utpipe_procList data structure
    !to make the lookup faster.
    channel = notFound
    do i = 1, utpipe_numChannels
       if (utpipe_procList(i) == procID) then
          channel = i
          exit
       end if
    end do
    if (channel == notFound) call Driver_abortFlash("Msg channel not found")

    !If there is a pending send in our desired channel we test all
    !send channels.  Request values are reset to MPI_REQUEST_NULL when
    !sends complete.
    if (utpipe_sendRequest(channel) /= MPI_REQUEST_NULL) then
       call UTPipeline_progressSendComm()
    end if

    !We can safetly add items to the send buffer if there is no pending send.
    if ( utpipe_sendState(channel) == OPEN_STATE .and. &
         utpipe_sendRequest(channel) == MPI_REQUEST_NULL ) then
       ptr = utpipe_sendCount(channel) + 1
       if (ptr > utpipe_channelSize) call Driver_abortFlash("Counting error")
       utpipe_sendBuf(:,ptr,channel) = item(:)
       utpipe_sendCount(channel) = ptr !Array is needed in utpipe_postSendMsg
          
       if (utpipe_sendCount(channel) == utpipe_channelSize) then
          call utpipe_postSendMsg(channel)
       end if
       isHandled = .true.
    else
       isHandled = .false.
    end if
  end subroutine UTPipeline_sendItem






  !We allow the caller to see the internal state of the pipeline
  !message exchange.
  subroutine UTPipeline_iterateItems(readOnlyFn)
    implicit none
    interface
       subroutine readOnlyFn(item, itemDescription)
         implicit none
         real, dimension(:), intent(IN) :: item
         character(len=*), intent(IN) :: itemDescription
       end subroutine readOnlyFn
    end interface
    integer :: i, n
    character(len=100) :: itemDescription

    do i = 1, utpipe_itemCount
       write (itemDescription,'(a,i10)') 'itemBuf: item ', i
       call readOnlyFn(utpipe_itemBuf(:,i), trim(itemDescription))
    end do

    !"The sender should not modify any part of the send buffer after a
    !nonblocking send operation is called, until the send completes."
    ![MPI-3 3.7.2].  "(the send operation itself leaves the content of
    !the send buffer unchanged)" [MPI-3 3.7.3]
    !... it should therefore always be OK to read what is there.
    do n = 1, utpipe_numChannels
       if (utpipe_sendCount(n) > 0) then
          do i = 1, utpipe_sendCount(n)
             write (itemDescription,'(2(a,i10))') 'sendBuf: channel ', n, &
                  & ', item ', i
             call readOnlyFn(utpipe_sendBuf(:,i,n), trim(itemDescription))
          end do
       end if
    end do

    !"The receiver should not access any part of the receive buffer
    !after a nonblocking receive operation is called, until the
    !receive completes." [MPI-3 3.7.2].
    !... it is not OK to read what is there.
    do n = 1, utpipe_numChannels
       if ( utpipe_recvCount(n) > 0 .and. &
            utpipe_recvRequest(n) == MPI_REQUEST_NULL ) then
          do i = 1, utpipe_recvCount(n)
             write (itemDescription,'(2(a,i10))') 'recvBuf: channel ', n, &
                  & ', item ', i
             call readOnlyFn(utpipe_recvBuf(:,i,n), trim(itemDescription))
          end do
       end if
    end do
  end subroutine UTPipeline_iterateItems

end module UTPipeline
