!!****f* source/physics/sourceTerms/Heatexchange/Heatexchange
!!
!! NAME
!!
!!  Heatexchange
!!
!!
!! SYNOPSIS
!!
!!   call Heatexchange ( integer(IN) :: blockCount, 
!!                       integer(IN) :: blockList(blockCount), 
!!                       real(IN)    ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply thermal head exchange among temperature components
!!  to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks where componenets should exchange
!!                heat
!!   dt  --       passed to the internal hx_burner module  
!!
!! PARAMETERS
!!
!!  useHeatexchange -- Boolean, True.  Turns on Heatexchange unit
!!
!!
!!***


subroutine Heatexchange ( blockCount, blockList, dt )
  
  implicit none

  !args
  integer, INTENT(in)                        :: blockCount
  integer, INTENT(in), DIMENSION(blockCount)  :: blockList
  real,    INTENT(in)                        :: dt

  return
  
end subroutine Heatexchange
