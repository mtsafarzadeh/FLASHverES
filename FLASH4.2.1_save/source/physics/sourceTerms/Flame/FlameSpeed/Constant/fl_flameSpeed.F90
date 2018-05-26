! Dean Townsley 2008

subroutine fl_flameSpeed( solnData, flamespeed, blockID, nlayers)

#include "Flash.h"
#include "constants.h"
#include "FortranLangFeatures.fh"

  use fl_fsData, only : fl_fsConstFlameSpeed, fl_fsConstFlameWidth
  implicit none
  real, dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  real, dimension(:,:,:),intent(out) :: flamespeed
  integer, intent(in) :: blockID, nlayers

  flamespeed(:,:,:) = fl_fsConstFlameSpeed
#ifdef FSPD_VAR
  solndata(FSPD_VAR,:,:,:) = fl_fsConstFlameSpeed
#endif

end subroutine
