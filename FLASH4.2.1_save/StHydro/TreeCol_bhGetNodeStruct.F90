!!****f* source/physics/TreeCol/TreeCol_bhGetNodeStruct
!!
!! NAME
!!
!!  TreeCol_bhGetNodeStruct
!!
!!
!! SYNOPSIS
!!
!!   call TreeCol_bhGetNodeStruct(botnodesize, nodesize)
!!
!! DESCRIPTION
!!
!!   Calculates structure (size, indeces of variables) of the 
!!   treecol section of the tree node
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine TreeCol_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
  implicit none
  integer, intent(IN) :: im, ix, iy, iz
  integer, intent(INOUT) :: bnsize, nsize
  
  return
end subroutine TreeCol_bhGetNodeStruct
