! See fl_fsInterface for a description subroutines
!
! Dean Townsley 2008
!

subroutine fl_fsGcMask(fl_gcMask,fl_gcDoEos)

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  ! this flame speed doesn't need anything in the guardcells
  ! so we leave both alone, as they default to false

  return

end subroutine
