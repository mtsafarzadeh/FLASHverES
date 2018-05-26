!
! Aaron Jackson 2009
!
! this correspondence is good for dens >~ 1.5e5 (degenerate)
! 0.92 is adjusted for fit
! constants for non-rel and relativistic degenerate e^- gas from 
!   Hansen & Kawaler

subroutine fl_fsUnburnDens(pres, ye, dens)

  implicit none

  real, intent(IN)  :: pres, ye
  real, intent(OUT) :: dens

  dens = 0.92/ye*sqrt( (pres/1.243e15)**(6.0/4.0) + &
                     (pres/1.004e13)**(6.0/5.0) )

  return
end subroutine fl_fsUnburnDens
