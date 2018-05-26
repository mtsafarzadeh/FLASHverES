subroutine Grid_guardCellMaskHook(ccMask, needEos)
  use Burn_interface, ONLY: Burn_guardCellMaskHook
  use Grid_data, ONLY: gr_enableMaskedGCFill
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

  if (gr_enableMaskedGCFill) then

     call Burn_guardCellMaskHook(ccMask,needEos)
     call Simulation_guardCellMaskHook(ccMask,needEos)

     ! ... add calls to other units if necessary ...

  end if

end subroutine Grid_guardCellMaskHook
