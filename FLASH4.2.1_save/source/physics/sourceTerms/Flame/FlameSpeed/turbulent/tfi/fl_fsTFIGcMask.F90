! See fl_fsTFIInterface for a description of subroutines
!
! Aaron Jackson 2010
!

#include "Flash.h"
subroutine fl_fsTFIGcMask(fl_gcMask,fl_gcDoEos)

  implicit none

  logical, dimension(:), intent(inout) ::  fl_gcMask
  logical, intent(inout) :: fl_gcDoEos

  ! turbulent strength is needed to calculate the TFI-enhanced flame speed
  fl_gcMask(TURB_VAR) = .true.

  return

end subroutine
