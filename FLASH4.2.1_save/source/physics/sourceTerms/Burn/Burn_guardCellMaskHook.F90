#include "Flash.h"

subroutine Burn_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

end subroutine Burn_guardCellMaskHook

