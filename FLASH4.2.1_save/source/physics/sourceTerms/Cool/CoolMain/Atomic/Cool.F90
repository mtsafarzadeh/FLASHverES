!!****f* source/physics/sourceterms/Cool/CoolMain.F90
!!
!! NAME
!!
!! Mix
!!
!!
!! SYNOPSIS
!!
!!   call Cool ( integer, intent(IN)    :: blockCount,
!!                integer(:), intent(IN) :: blockList,
!!                real, intent(IN)       ::  dt  )
!!
!! DESCRIPTION
!!
!!  Apply mix module to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks which should receive burning
!!   dt  --       passed to the internal bn_burner module
!!
!! PARAMETERS
!!
!! NOTES
!!
!!
!!***

subroutine Cool ( blockCount, blockList, dt,time  )
  use Cool_data
  implicit none
  integer, intent(in) :: blockCount
  integer, intent(in),dimension(blockCount) :: blockList
  real, intent(in)    :: dt,time

  integer :: thisBlock,blockID
  external Cool_block

!  print *,'running cool.f90 blockCount:',blockCOunt
  ! loop over list of blocks passed in
if (nocool > 0.) then 
  do thisBlock = 1, blockCount
    blockID = blockList(thisBlock)
        call Cool_block(blockID)
  enddo
endif

return
end subroutine Cool
