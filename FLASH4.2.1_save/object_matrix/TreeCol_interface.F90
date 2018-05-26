!!****h* source/physics/TreeCol/TreeCol_interface
!!
!! This is the header file for the gravity module that defines its
!! public interfaces.
!!***

Module TreeCol_interface

  interface TreeCol_finalize
     subroutine TreeCol_finalize()
     end subroutine TreeCol_finalize
  end interface

  interface TreeCol_init
     subroutine TreeCol_init()
       
     end subroutine TreeCol_init
  end interface

  interface TreeCol_potentialListOfBlocks
     subroutine TreeCol_potentialListOfBlocks(blockCount,blockList)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN) :: blockList
     end subroutine TreeCol_potentialListOfBlocks
  end interface

  interface TreeCol_unitTest
     subroutine TreeCol_unitTest( fileUnit, perfect)
       implicit none
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine TreeCol_unitTest
  end interface

! The following subroutines and functions are needed by the tree solver
#include "constants.h"
#include "FortranLangFeatures.fh"
  interface TreeCol_bhAccBotNode
    subroutine TreeCol_bhAccBotNode(blockno, point, blkLimits, locCoords, solnData, botnode, accnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, dimension(:,:,:) :: locCoords
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(IN) :: botnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine TreeCol_bhAccBotNode
  end interface

  interface TreeCol_bhAccNode
    subroutine TreeCol_bhAccNode(subnode, accnode)
      real, dimension(:), intent(IN)  :: subnode
      real, dimension(:), intent(INOUT) :: accnode
    end subroutine TreeCol_bhAccNode
  end interface

  interface TreeCol_bhBotNodeContrib
    subroutine TreeCol_bhBotNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
      real, dimension(:), intent(IN) :: node
      real, intent(IN) :: ndSize
      real, dimension(MDIM+2), intent(IN) :: dr
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeCol_bhBotNodeContrib
  end interface


  interface TreeCol_bhFillBotNode
    subroutine TreeCol_bhFillBotNode(blockno, point, solnData, botnode)
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      real, DIMENSION(:,:,:,:), POINTER :: solnData
      real, dimension(:), intent(INOUT) :: botnode
    end subroutine TreeCol_bhFillBotNode
  end interface

  interface TreeCol_bhNodeContrib
    subroutine TreeCol_bhNodeContrib(node, ndSize, dr, blockno, point, blkLimits, solnData)
      real, dimension(:), intent(IN) :: node
      real, intent(IN) :: ndSize
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      real, dimension(MDIM+2), intent(IN) :: dr
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeCol_bhNodeContrib
  end interface

  interface TreeCol_bhNormalizeNode
    subroutine TreeCol_bhNormalizeNode(smr, node)
      real, dimension(MDIM), intent(IN) :: smr
      real, dimension(:), intent(INOUT) :: node
    end subroutine TreeCol_bhNormalizeNode
  end interface

  interface TreeCol_bhPostprocNode
    subroutine TreeCol_bhPostprocNode(ndSize, node)
      real, intent(IN) :: ndSize
      real, dimension(:), intent(INOUT) :: node
    end subroutine TreeCol_bhPostprocNode
  end interface

  interface TreeCol_bhStartBlock
    subroutine TreeCol_bhStartBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeCol_bhStartBlock
  end interface

  interface TreeCol_bhFinalizeBlock
    subroutine TreeCol_bhFinalizeBlock(blockno, blkLimits, solnData)
      integer, intent(IN) :: blockno
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER :: solnData
    end subroutine TreeCol_bhFinalizeBlock
  end interface


  interface TreeCol_bhInitFieldVar
    subroutine TreeCol_bhInitFieldVar(gpotVar)
      implicit none
      integer, intent(IN) :: gpotVar
    end subroutine TreeCol_bhInitFieldVar
  end interface

  interface TreeCol_bhMAC
    logical function TreeCol_bhMAC(node, ndSize2, dr, blockno, point, blkLimits, solnData)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize2
    integer, dimension(MDIM), intent(IN) :: point
    real, dimension(MDIM+2), intent(IN) :: dr
    integer, intent(IN) :: blockno
    integer, dimension(2,MDIM)   :: blkLimits
    real, DIMENSION(:,:,:,:), POINTER :: solnData
    end function TreeCol_bhMAC
  end interface
  
  interface TreeCol_bhPartErr
    subroutine TreeCol_bhPartErr(node, ndSize, dr, perr)
    implicit none
    real, dimension(:), intent(IN) :: node
    real, intent(IN) :: ndSize
    real, dimension(MDIM+2), intent(IN) :: dr
    real :: perr
    end subroutine TreeCol_bhPartErr
  end interface
  
  interface TreeCol_bhSelfContrib
    subroutine TreeCol_bhSelfContrib(node, cellsize, blockno, point, blkLimits, solnData)
      implicit none
      real, dimension(:), intent(IN) :: node
      real, dimension(MDIM), intent(IN) :: cellsize
      integer, intent(IN) :: blockno
      integer, dimension(MDIM), intent(IN) :: point
      integer, dimension(2,MDIM)   :: blkLimits
      real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
    end subroutine TreeCol_bhSelfContrib
  end interface

  interface TreeCol_bhGetNodeStruct
    subroutine TreeCol_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
      implicit none
      integer, intent(IN) :: im, ix, iy, iz
      integer, intent(INOUT) :: bnsize, nsize
    end subroutine TreeCol_bhGetNodeStruct
  end interface




! end of tree solver subroutines/functions



end Module TreeCol_interface
