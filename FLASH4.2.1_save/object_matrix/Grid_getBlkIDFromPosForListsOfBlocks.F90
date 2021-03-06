!!****if* source/Grid/GridMain/Grid_getBlkIDFromPosForListsOfBlocks
!!
!! NAME
!!  Grid_getBlkIDFromPosForListsOfBlocks
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkIDFromPosForListsOfBlocks(real(IN) :: pos(:), 
!!                            integer(IN) :: blkList(:),
!!                            integer(IN) :: blkCount,
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID,
!!                   optional,integer(IN)  :: comm)
!!  
!!  call Grid_getBlkIDFromPos(real(IN)    :: pos(:), 
!!                            integer(IN) :: blkList(:),
!!                            integer(IN) :: blkCount,
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID,
!!                   optional,integer(IN)  :: comm)
!!
!! DESCRIPTION 
!! 
!!  Returns the processor and block ID
!!  containing the cell that overlaps with the 
!!  specified position co-ordinate
!! 
!! 
!! ARGUMENTS 
!!
!!  pos        :: co-ordinates of the point
!!  blkList    :: the list of blocks to search
!!  blkCount   :: the count of blocks in the list
!!  ansBlockID    :: the local block ID of the block that contains the point
!!  ansProcID     :: the processor ID that contains the point
!!  comm       :: if a communicator other than the default mesh communicator is
!!                desired, it should be specified here
!!
!! EXAMPLE
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkIDFromPosForListsOfBlocks(pos,blkList, blkCount,ansBlockID, ansProcID,comm)

  use Grid_data, ONLY : gr_minCellSizes, gr_globalDomain, gr_meshMe, gr_meshComm
  use Grid_interface, ONLY : Grid_getBlkCornerID
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  include "Flash_mpi.h"
  real, dimension(1:MDIM), intent(IN) :: pos
  integer,intent(IN)  :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  integer, intent(OUT) :: ansBlockID, ansProcID
  integer, optional, intent(IN) :: comm

  integer :: i,j,blk, ierr, mycomm
  integer,parameter :: bufsize=2
  integer,dimension(bufsize) :: sbuf, rbuf
  integer, dimension(MDIM) :: cornerID, stride, cornerIDHigh, posInd
  logical,dimension(MDIM) :: inBlk
  logical :: onUpperBoundary(NDIM), found
  real,dimension(MDIM) :: temp

  
  ansBlockID=NONEXISTENT
  ansProcID=NONEXISTENT
  temp(1:NDIM)=(gr_globalDomain(HIGH,1:NDIM) - pos(1:NDIM))
  onUpperBoundary(:) = (temp(1:NDIM)==0)
  if(minval(temp).ge.0) then
     temp(1:NDIM)=(pos(1:NDIM)-gr_globalDomain(LOW,1:NDIM))/gr_minCellSizes(1:NDIM)
     if(minval(temp).ge.0) then
        posInd(1:NDIM)=int(temp(1:NDIM))+1
        where (onUpperBoundary(:)) posInd(1:NDIM)=posInd(1:NDIM)-1
        sbuf=0
        do i = 1,blkCount
           ansBlockID=blkList(i)
           call Grid_getBlkCornerID(ansBlockID,cornerID,stride,cornerIDHigh)
           inBlk=.true.
           inBlk(1:NDIM)=(posInd(1:NDIM).ge.cornerID(1:NDIM)).and.&
                (posInd(1:NDIM).le.cornerIDHigh(1:NDIM))
           found=inBlk(IAXIS).and.inBlk(JAXIS).and.inBlk(KAXIS)
           if(found) then
              sbuf(1)=gr_meshMe+1 !! to compensate for the fact that 0 is a valid procID
              sbuf(2)=ansBlockID
           end if
        end do
        if(present(comm)) then
           mycomm=comm
        else
           mycomm=gr_meshComm
        end if
        call MPI_AllReduce(sbuf,rbuf, 1, MPI_2INTEGER, MPI_MAXLOC,mycomm, ierr)
        ansProcID=rbuf(1)-1 !! remove the added 1
        ansBlockID=rbuf(2)
     end if
  end if
end subroutine Grid_getBlkIDFromPosForListsOfBlocks
