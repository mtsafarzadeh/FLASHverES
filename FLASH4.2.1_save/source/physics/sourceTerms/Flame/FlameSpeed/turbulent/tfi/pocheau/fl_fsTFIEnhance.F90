! Aaron Jackson 2010
!
! This subroutine calculates the enhancement factor to the flame speed.

subroutine fl_fsTFIEnhance(E, up, s, de, dl0, de_over_dl1, E_lim)

  use fl_fsTFIData, only : fl_fsTFICt

  implicit none

  real, intent(out) :: E
  real, intent(in) :: up, s, de, dl0, de_over_dl1
  real, intent(out), optional :: E_lim

  real :: up_over_s

  ! use negative to represent no quenching
  if (present(E_lim)) E_lim = -1.0e0

  up_over_s = up / s

  ! calculate enhancement
  E = sqrt( 1.0e0 + fl_fsTFICt * up_over_s * up_over_s )

  return
end subroutine fl_fsTFIEnhance
