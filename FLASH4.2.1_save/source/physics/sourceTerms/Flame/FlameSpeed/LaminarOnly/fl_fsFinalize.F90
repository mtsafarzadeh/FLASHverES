!
! Dean Townsley 2008
!

subroutine fl_fsFinalize()

   implicit none

   call fl_fsLaminarFinalize
   call fl_fsTFIFinalize

   return

end subroutine fl_fsFinalize
