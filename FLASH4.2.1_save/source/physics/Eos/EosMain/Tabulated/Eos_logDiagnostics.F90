subroutine Eos_logDiagnostics(force)
  use eos_tabInterface,ONLY: eos_tabLogOutsideCounts, eos_tabZeroOutsideCounts
  implicit none
  logical, intent(IN) :: force

  call eos_tabLogOutsideCounts(force)
  call eos_tabZeroOutsideCounts()

end subroutine Eos_logDiagnostics
