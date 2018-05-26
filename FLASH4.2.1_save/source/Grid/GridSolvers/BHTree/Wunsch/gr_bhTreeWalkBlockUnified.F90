!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhTreeWalkBlockUnified
!!
!! NAME
!!
!!  gr_bhTreeWalkBlock
!!
!!
!! SYNOPSIS
!!
!!   gr_bhTreeWalkBlock(
!!          integer,intent(in) :: block
!!          )
!!
!! DESCRIPTION
!!
!!   Traverses the tree in a specific block. In this, "unified tree walk",
!!   the MAC for all interactions is evaluated for each cell, as in the original
!!   Barnes&Hut algorithm.
!!
!! ARGUMENTS
!!
!!  block    - ID of a block where the contribution is calculated
!!
!!***



subroutine gr_bhTreeWalkBlockUnified(block)

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_bhLocalInterface, ONLY : gr_bhLeafContrib, gr_bhParentContrib &
    , gr_bhStartBlock, gr_bhFinalizeBlock, gr_bhMAC, gr_bhBotNodeContrib &
    , gr_bhNodeContrib, gr_bhSelfContrib, gr_bhPartErr, gr_bhPeriodicDr
  use Driver_interface, ONLY : Driver_abortFlash

  use gr_bhData, ONLY : gr_bhTreeLevels, &
    gr_bhTreeBS, gr_bhLocCoords, gr_bhTreeMyPE, &
    gr_bhTreeLrefine, gr_bhTreeNodeSize2, gr_bhPhysMACTW_step, &
    gr_bhTreeDcount, gr_bhTreeLimangle2i, &
    gr_bhTreeNFirstLev, gr_bhTreeChild, gr_bhTreeNodetype, &
    gr_bhTreeFirstLevBlocks, gr_bhTreeParentTree, &
    gr_bhTreeCellSize, gr_bhTWMaxQueueSize, gr_bhTreeArray, &
    GR_TREE_NSIZE, GR_TREE_BNSIZE, GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, &
    GR_TREE_BTP_POS ,GR_TREE_BTP_LEV, GR_TREE_BTP_I, GR_TREE_BTP_J, GR_TREE_BTP_K, &
    GR_TREE_BTP_C1, GR_TREE_BTP_C8, gr_bhTreeBCen, gr_bhLocRecvTreeLevels, &
    gr_bhBlockTreePos, gr_bhTreeNodeSize, gr_bhBlockTreeNodeCen
  use tree, ONLY : nodetype, child, mchild

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer, intent(IN) :: block
  integer       :: i, j, k, l, cpu, tr, sp, istat, level
  integer       :: ii, jj, kk, ii_mbLp1, jj_mbLp1, kk_mbLp1
  integer       :: node_btp, pos
  real          :: dcount(1:3)
  real, POINTER, DIMENSION(:,:,:,:) :: solnData
  integer :: stack(3, gr_bhTWMaxQueueSize)
  real    :: dist2
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real, dimension(MDIM+2) :: dr ! two additional elements dr^2 and 1/|dr|
  integer, dimension(MDIM) :: point
  real, dimension(MDIM) :: drDiff, drGC
  real, allocatable  :: node(:), botnode(:)
  logical       :: node_accepted, gcell = .false.
  !real :: ss, perr(2)

  ! THERE IS PROBABLY A BUG IN THIS ROUTINE
  ! RESULTS SLIGHTLY DIFFER AT THE VERY CENTRE FROM RESULTS
  ! OF TreeWalkBlock and TreeWalkBlockPQ


  ! allocate the node array
  allocate(node(GR_TREE_NSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhTreeWalkBlockUnified.F90")
  allocate(botnode(GR_TREE_BNSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate botnode in gr_bhTreeWalkBlockUnified.F90")

  ! get information about the block
  ! indeces necessary to write values into solnData
  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC) 
  ! pointer to the density and gp field
  call Grid_getBlkPtr(block,solnData,CENTER) 

  ! call physics routines to do block initialization
  call gr_bhStartBlock(block, blkLimits, solnData)
  
  ! reset distance counter
  dcount = 0.0

  do kk = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS) ! this block
    point(KAXIS) = kk
    kk_mbLp1 = kk - blkLimits(LOW, KAXIS) + 1
    do jj = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
      point(JAXIS) = jj
      jj_mbLp1 = jj - blkLimits(LOW, JAXIS) + 1
      do ii = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        point(IAXIS) = ii
        ii_mbLp1 = ii - blkLimits(LOW, IAXIS) + 1
        ! put 1st level blocks on stack
        stack(1:2,1:gr_bhTreeNFirstLev) = gr_bhTreeFirstLevBlocks(1:2,1:gr_bhTreeNFirstLev)
        stack(3,:) = 1
        sp = gr_bhTreeNFirstLev
        !ss = 0.0
   
        ! walk the tree - start with the amr tree, continue into block-trees
        do
   
          ! pop a block from the stack
          tr       = stack(1,sp)
          cpu      = stack(2,sp)
          node_btp = stack(3,sp)
          pos      = -1
          level    = 0
          sp = sp - 1

          ! determine distance between the cell and the tree node
          ! 1. subtract cell coordinates
          dr(IAXIS) = -gr_bhLocCoords(ii_mbLp1, IAXIS, block)
          dr(JAXIS) = -gr_bhLocCoords(jj_mbLp1, JAXIS, block)
          dr(KAXIS) = -gr_bhLocCoords(kk_mbLp1, KAXIS, block)
          if (gr_bhTreeNodetype(tr, cpu) .ne. 1) then ! 1 = LEAF block
            ! 2. parent block, add position of its mass centre
            dr(IAXIS) = dr(IAXIS) + gr_bhTreeParentTree(GR_TREE_IX, tr, cpu)
            dr(JAXIS) = dr(JAXIS) + gr_bhTreeParentTree(GR_TREE_IY, tr, cpu)
            dr(KAXIS) = dr(KAXIS) + gr_bhTreeParentTree(GR_TREE_IZ, tr, cpu)
            ! fill the node
            node(1:GR_TREE_NSIZE) = gr_bhTreeParentTree(1:GR_TREE_NSIZE, tr, cpu)
          else
            ! leaf block, get the position in gr_bhTreeArray and the node level
            pos = gr_bhBlockTreePos(GR_TREE_BTP_POS, node_btp)
            level = gr_bhBlockTreePos(GR_TREE_BTP_LEV, node_btp)

            if (level == gr_bhTreeLevels) then
              ! leaf level of block-tree: cell-cell interaction
              botnode = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_BNSIZE-1)

              ! get indeces of the interacting cell in the block
              i = gr_bhBlockTreePos(GR_TREE_BTP_I, node_btp)
              j = gr_bhBlockTreePos(GR_TREE_BTP_J, node_btp)
              k = gr_bhBlockTreePos(GR_TREE_BTP_K, node_btp)
            
              if ((gr_bhTreeMyPE .eq. cpu) .and. (block .eq. tr) &
              .and. (i .eq. (ii-blkLimits(LOW,IAXIS)+1)) &
              .and. (j .eq. (jj-blkLimits(LOW,JAXIS)+1)) &
              .and. (k .eq. (kk-blkLimits(LOW,KAXIS)+1))) then
                ! add self-contribution of the cell
                call gr_bhSelfContrib(botnode, gr_bhTreeCellSize(gr_bhTreeLrefine(tr,cpu),:) &
                & , block, point, blkLimits, solnData)
              else
                ! 2. add position of a leaf node (cell) of the block-tree
                dr(IAXIS) = dr(IAXIS) + gr_bhTreeBCen(IAXIS, tr, cpu) &
                & - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS)*(gr_bhTreeBS-2*i+1)
                dr(JAXIS) = dr(JAXIS) + gr_bhTreeBCen(JAXIS, tr, cpu) & 
                & - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS)*(gr_bhTreeBS-2*j+1)
                dr(KAXIS) = dr(KAXIS) + gr_bhTreeBCen(KAXIS, tr, cpu) &
                & - 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS)*(gr_bhTreeBS-2*k+1)
              endif ! test for self contribution
            else
              ! non-leaf nodes of the block tree
              ! 2. add position of the mass centre of the node
              dr(IAXIS) = dr(IAXIS) + gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IX-1)
              dr(JAXIS) = dr(JAXIS) + gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IY-1)
              dr(KAXIS) = dr(KAXIS) + gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IZ-1)
              ! fill the node
              node = gr_bhTreeArray(cpu, tr)%p(pos:pos+GR_TREE_NSIZE-1)
            endif ! test for the leaf level of the block-tree

          endif
   
          ! 3. correct for periodic boundaries
          call gr_bhPeriodicDr(dr)

          ! calculate square of distance and inverted distance (stored also in dr)
          dist2  = dr(IAXIS)*dr(IAXIS) + dr(JAXIS)*dr(JAXIS) + dr(KAXIS)*dr(KAXIS)
          dr(MDIM+1) = dist2
          dr(MDIM+2) = sqrt(1.0 / (dist2 + 1D-99))
   
          if (level == gr_bhTreeLevels) then
            ! at leaf level, just call BotNodeContrib and thats it
            call gr_bhBotNodeContrib(botnode, gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) &
            & , dr, block, point, blkLimits, solnData)
            dcount(1) = dcount(1) + 1
          else
            ! position vector to the node geometric centre
            ! drDiff is a vector from the node mass centre to its geometrical centre
            drDiff(:) = gr_bhTreeBCen(:, tr, cpu) &
            & + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), :) &
            & * gr_bhBlockTreeNodeCen(:, node_btp) &
            & - node(GR_TREE_IX:GR_TREE_IZ)
            !drDiff(IAXIS) = gr_bhTreeBCen(IAXIS, tr, cpu) &
            !& + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), IAXIS) &
            !& * gr_bhBlockTreeNodeCen(IAXIS, node_btp) &
            !& - node(GR_TREE_IX)
            !drDiff(JAXIS) = gr_bhTreeBCen(JAXIS, tr, cpu) &
            !& + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), JAXIS) &
            !& * gr_bhBlockTreeNodeCen(JAXIS, node_btp) &
            !& - node(GR_TREE_IY)
            !drDiff(KAXIS) = gr_bhTreeBCen(KAXIS, tr, cpu) &
            !& + gr_bhTreeCellSize(gr_bhTreeLrefine(tr, cpu), KAXIS) &
            !& * gr_bhBlockTreeNodeCen(KAXIS, node_btp) &
            !& - node(GR_TREE_IZ)
            drGC = dr(1:MDIM) + drDiff

            ! call the multipole acceptance criterion (MAC)
            node_accepted = gr_bhMAC(gr_bhPhysMACTW_step, node &
            & , gr_bhTreeNodeSize2(level+gr_bhTreeLrefine(tr,cpu)) &
            & , drGC, dr, block, point, blkLimits, solnData)
            if (node_accepted) then 
              ! node accepted, call contribution subroutine
              call gr_bhNodeContrib(node, gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) &
              & , dr, block, point, blkLimits, solnData)
              ! increase the distance counters
              if (level > 0) then
                ! cell-node interaction
                dcount(2) = dcount(2) + 1
              else
                ! cell-block interaction
                dcount(3) = dcount(3) + 1
              endif
              !call gr_bhPartErr(node, gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)), dr, perr)
              !ss = ss + perr(1)**2
            else
              ! node not accepted, put all its children on the stack
              if (gr_bhTreeNodetype(tr, cpu) == 1) then
                ! LEAF block, at first test if the requeted level is present on my cpu
                if (cpu /= gr_bhTreeMyPE) then
                  if ((level+1) > gr_bhLocRecvTreeLevels(tr, cpu)) then
                    print *, "Asking for higher level than present on a given cpu: "
                    print *, "TWB2: ", gr_bhTreeMyPE, block, tr, cpu, sqrt(dr(MDIM+1)) &
                    , gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) &
                    , gr_bhTreeNodeSize(level+gr_bhTreeLrefine(tr,cpu)) / sqrt(dr(MDIM+1)) &
                    , "|dr:", dr &
                    , "|Cell:", &
                      gr_bhLocCoords(ii_mbLp1, IAXIS, block), &
                      gr_bhLocCoords(jj_mbLp1, JAXIS, block), &
                      gr_bhLocCoords(kk_mbLp1, KAXIS, block) &
                    , "|TA:", gr_bhTreeArray(cpu, tr)%p(pos+GR_TREE_IX-1:pos+GR_TREE_IZ-1)
                    call Driver_abortFlash("TREE_ARRAY ACCESS FAILED!!!")
                  endif
                endif
                ! put all children on the stack
                do l = GR_TREE_BTP_C1,GR_TREE_BTP_C8
                  sp = sp + 1
                  stack(1,sp) = tr
                  stack(2,sp) = cpu
                  stack(3,sp) = gr_bhBlockTreePos(l, node_btp)
                enddo
              else
                ! 
                do l = 1,mchild
                  sp = sp + 1
                  stack(1,sp) = gr_bhTreeChild(1,l,tr,cpu)
                  stack(2,sp) = gr_bhTreeChild(2,l,tr,cpu)
                  stack(3,sp) = 1
                enddo
              endif
   
            endif ! testing MAC
            
          endif ! LEAF block
   
   
          ! stack is empty - exiting
          if (sp == 0) exit
        enddo ! end tree walk
        !print *, "SS = ", sqrt(ss)
      enddo
    enddo
  enddo

  ! call physics routines to do block finalization
  call gr_bhFinalizeBlock(block, blkLimits, solnData)

  ! release the block pointer
  call Grid_releaseBlkPtr(block,solnData, CENTER)
  gr_bhTreeDcount = gr_bhTreeDcount + dcount

  deallocate(node)
  deallocate(botnode)

  return
end subroutine gr_bhTreeWalkBlockUnified

