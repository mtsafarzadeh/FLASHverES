!
! Dean Townsley 008
!

subroutine Flame_laminarSpeed(dens, s, ds, info)

  use fl_fsData

  implicit none

  real, intent(in)   :: dens
  real, intent(out)  :: s
  real, optional, intent(out) :: ds
  real, dimension(:), optional, intent(in) :: info

  s = fl_fsConstFlameSpeed

  if (present(ds)) ds = fl_fsConstFlameWidth

end subroutine Flame_laminarSpeed
