!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhTreeWalk
!!
!! NAME
!!
!!  gr_bhTreeWalk
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhTreeWalk()
!!
!! DESCRIPTION
!!
!!   Executes the tree walk. Calls gr_treeTreeWalkBlock
!!   for all LEAF blocks.
!!
!! ARGUMENTS
!!
!!  idensvar : number of grid varible with density (recently not used)
!!  ipotvar  : number of grid varible with grav. potential (for compatibility with other solvers)
!!
!!***



subroutine gr_bhTreeWalk()

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement
  use gr_bhData, ONLY : gr_bhTreeDcount, gr_bhTreeZones, &
    gr_bhTreeMyPE, gr_bhComm, gr_bhTreeBS, &
    gr_bhTreeLimangle, gr_bhPhysMACTW_step, gr_bhUseUnifiedTW, &
    gr_bhTWType, GR_TREE_TWSTD, GR_TREE_TWUNI, GR_TREE_TWPQ, &
    gr_bhOAMin, gr_bhOAAvg, gr_bhOAMax, gr_bhOACnt !, gr_bhTotMass
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_bhLocalInterface, ONLY : gr_bhTreeWalkBlock, &
    gr_bhTreeWalkBlockUnified, gr_bhTreeWalkBlockPQ
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  
  integer :: tot_zones, ierr
  real    :: tot_dcount(1:3), tot_OAAvg, tot_OAMin, tot_OAMax, tot_OACnt
  character(len=256) :: strBuff

  call Timers_start("treewalk")

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! for calculating numebrs of interactions and average limit angle
  gr_bhTreeDcount = 0.0
  gr_bhTreeZones = 0
  gr_bhOAAvg = 0.0
  gr_bhOAMin = 10.0
  gr_bhOAMax = 0.0
  gr_bhOAcnt = 0
  !gr_bhTotMass = 0.0
  do blockID = 1, blockCount
    if (gr_bhTWType == GR_TREE_TWSTD) then
      call gr_bhTreeWalkBlock(blockList(blockID))
    else if (gr_bhTWType == GR_TREE_TWUNI) then
      call gr_bhTreeWalkBlockUnified(blockList(blockID))
    else if (gr_bhTWType == GR_TREE_TWPQ) then
      call gr_bhTreeWalkBlockPQ(blockList(blockID))
    else
    endif
    gr_bhTreeZones = gr_bhTreeZones + gr_bhTreeBS*gr_bhTreeBS*gr_bhTreeBS
  enddo
  !print *, "TW: ", gr_bhTreeMyPE, gr_bhTreeZones, gr_bhTreeDcount
  !print *, "TW totmass: ", gr_bhTreeMyPE, gr_bhTotMass

  call MPI_Reduce(gr_bhTreeDcount,tot_dcount,3,FLASH_REAL,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhTreeZones,tot_zones,1,FLASH_INTEGER,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOAAvg,tot_OAAvg,1,FLASH_REAL,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOAMin,tot_OAMin,1,FLASH_REAL,MPI_MIN,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOAMax,tot_OAMax,1,FLASH_REAL,MPI_MAX,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhOACnt,tot_OACnt,1,FLASH_REAL,FLASH_SUM,MASTER_PE,gr_bhComm, ierr)
  tot_OAAvg = tot_OAAvg / tot_OACnt

  if (gr_bhTreeMyPE == MASTER_PE) then
     write (strBuff, '("cell-cell distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(1), tot_dcount(1)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("cell-node distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(2), tot_dcount(2)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("cell-block distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(3), tot_dcount(3)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("opening angle (min,avg,max): ", e10.3, e10.3, e10.3)') &
     & tot_OAMin, tot_OAAvg, tot_OAMax
     call Logfile_stamp( strBuff, "[BHTree]")
     !write (strBuff, '("count check: (dcount vs. OACnt): ", e20.10, e20.10)') &
     !& tot_dcount(1)+tot_dcount(2)+tot_dcount(3), tot_OACnt
     !call Logfile_stamp( strBuff, "[BHTree]")
  end if

  call Timers_stop("treewalk")

  return
end subroutine gr_bhTreeWalk


