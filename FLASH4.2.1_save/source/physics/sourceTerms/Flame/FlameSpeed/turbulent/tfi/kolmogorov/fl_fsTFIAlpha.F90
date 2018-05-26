! See fl_fsTFIInterface for subroutine functional description
!
! Aaron Jackson 2010
!
! This subroutine calculates the model constant alpha in the
! Colin et al. (2000) TFI perscription

subroutine fl_fsTFIAlpha(alpha, up, de, dl0)

#include "Flash.h"

  use fl_fsTFIData, ONLY : fl_fsTFIBeta, fl_fsTFIViscosity

  implicit none

  real, intent(out) :: alpha
  real, intent(in) :: up, de, dl0

  real :: Re

  Re = up * de / fl_fsTFIViscosity

  ! fl_fsBeta already includes 2 ln2 / 3 c_ms pre-factor
  ! this assumes Kolmogorov limit, depends on Re
  alpha = fl_fsTFIBeta / sqrt(Re)

  return

end subroutine fl_fsTFIAlpha
