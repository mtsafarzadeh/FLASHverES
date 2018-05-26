!!****h* source/physics/sourceTerms/Chemistry/Chemistry_interface
!!
!! This is the header file for the heat module that defines its
!! public interfaces.
!!***
Module Chemistry_interface
#include "constants.h"
#include "Flash.h"

  interface Chemistry
     subroutine Chemistry (blockCount,blockList,dt)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN)::blockList
       real,intent(IN) :: dt
     end subroutine Chemistry
  end interface

  interface Chemistry_computeDt
     subroutine Chemistry_computeDt (block_no, myPE, &
          x, dx, uxgrid, &
          y, dy, uygrid, &
          z, dz, uzgrid, &
          blkLimits,blkLimitsGC,  &
          solnData,   &
          dt_check, dt_minloc )

       integer, intent(IN) :: block_no, myPE
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x
       real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y
       real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z
       real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: dx, uxgrid
       real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: dy, uygrid
       real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: dz, uzgrid
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x
       real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y
       real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: dx, uxgrid
       real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: dy, uygrid
       real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: dz, uzgrid
#endif
       real,INTENT(OUT)    :: dt_check
       integer,INTENT(OUT)    :: dt_minloc(5)
       real, pointer, dimension(:,:,:,:) :: solnData
     end subroutine Chemistry_computeDt

  end interface

  interface Chemistry_init
     subroutine Chemistry_init()
     end subroutine Chemistry_init
  end interface


  interface Chemistry_finalize
     subroutine Chemistry_finalize ()
     end subroutine Chemistry_finalize
  end interface

  interface Chemistry_sendOutputData
     subroutine Chemistry_sendOutputData()
     end subroutine Chemistry_sendOutputData
  end interface


end Module Chemistry_interface
