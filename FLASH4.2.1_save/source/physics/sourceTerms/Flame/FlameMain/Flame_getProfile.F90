! see Flame_interface.F90 at top level for function description
!
! Dean Townsley 2008
!
subroutine Flame_getProfile(x, f)

  use Flame_data, ONLY : fl_width, fl_initProfileAdjustWidth

  implicit none
  real, intent(in)  :: x
  real, intent(out) :: f

  ! This is an approximate profile form based on widths determined in
  ! Vladimirova et al.  Here the "width" is approximately the
  ! distance between where phi=0.12 and 0.88.  Over twice width, phi
  ! goes from 0.02 to 0.98.
  ! The fl_initProfileAdjustmentWidth is to allow compatibility with
  ! slight variations if necessary

  ! tanh is slow, but no need to be super efficient here since this
  ! should only be called at init
  f = 0.5 * (1.0 - tanh(x/fl_width/0.5/fl_initProfileAdjustWidth))

  return

end subroutine Flame_getProfile
