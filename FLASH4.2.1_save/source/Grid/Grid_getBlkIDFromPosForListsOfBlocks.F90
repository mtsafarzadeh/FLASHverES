!!****f* source/Grid/Grid_getBlkIDFromPosForListsOfBlocks
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

subroutine Grid_getBlkIDFromPosForListsOfBlocks(pos,blkList, blkCount,ansBlockID, ansProcID,comm)

  implicit none

  real, dimension(1:MDIM), intent(IN) :: pos
  integer,intent(IN)  :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  integer, intent(OUT) :: ansBlockID, ansProcID
  integer,optional,intent(IN) :: comm
  
  ansProcID=NONEXISTENT
  ansBlockID=NONEXISTENT
  
end subroutine Grid_getBlkIDFromPosForListsOfBlocks
