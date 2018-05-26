! see Flame_interface.F90 at top level for function description
!
! Dean Townsley 2008
!
subroutine Flame_getWidth(laminarWidth)

  use Flame_data, ONLY : fl_width

  implicit none
  real, intent(OUT) :: laminarWidth

  laminarWidth=fl_width
  return

end subroutine Flame_getWidth
