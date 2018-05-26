! Deallocates Atwood Table
! Aaron Jackson 2008

subroutine fl_fsAtwoodEndTable()

  use fl_fsAtwoodData, ONLY : fl_fsAtwoodTabA
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer :: istat

  deallocate(fl_fsAtwoodTabA,stat=istat)
     if (istat /= 0) call Driver_abortFlash("Cannot deallocate AtwoodTabA in fl_fsAtwoodEndTable")

  return
end subroutine
