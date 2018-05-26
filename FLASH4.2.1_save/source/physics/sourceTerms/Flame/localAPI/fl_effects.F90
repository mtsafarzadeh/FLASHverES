! Dean Townsley 2008
!
! this is a stub for flame effects that don't do anything locally

#include "FortranLangFeatures.fh"

subroutine fl_effects( solnData, flamdot, dt, blockID)

  implicit none

  real, dimension(:,:,:,:),POINTER_INTENT_IN  :: solnData
  real,dimension(:,:,:), intent(in)     :: flamdot
  real,intent(in)                       :: dt
  integer, intent(in)                   :: blockID

  return
end subroutine
