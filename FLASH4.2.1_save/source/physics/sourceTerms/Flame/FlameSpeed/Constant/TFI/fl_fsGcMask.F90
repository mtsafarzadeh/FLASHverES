! See fl_fsInterface for a description subroutines
!
! Dean Townsley 2008
!

subroutine fl_fsGcMask(fl_gcMask,fl_gcDoEos)

  use fl_fsTFIInterface, ONLY : fl_fsTFIGcMask

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  ! this flame speed doesn't need anything in the guardcells
  ! so we leave both alone, as they default to false

  call fl_fsTFIGcMask(fl_gcMask, fl_gcDoEos)

  return

end subroutine
