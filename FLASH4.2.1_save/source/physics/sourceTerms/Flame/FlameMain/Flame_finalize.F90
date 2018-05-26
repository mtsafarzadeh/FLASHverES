!  see Flame_interface.F90 at top level for function description
!
! Dean Townsley 2008
!
subroutine Flame_finalize()

  use Flame_data, ONLY: fl_useFlame

  implicit none

  if(fl_useFlame) then 
     call fl_fsFinalize
     call fl_effFinalize
  end if

  return

end subroutine Flame_finalize
