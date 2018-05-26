#include "Flash.h"

subroutine Simulation_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

  !!  Additional logic necessary due to the hydrostatic boundary conditions.
  ! For constant isothermal BC we need temperature and material info (ye, sumy) in interior.
  ! This temperature setting does not trigger the "needEos" flag.
  if ( ccMask(DENS_VAR) .or. ccMask(EINT_VAR) .or. ccMask(TEMP_VAR) &
       .or. ccMask(PRES_VAR) .or. ccMask(ENER_VAR) &
       .or. ccMask(GAMC_VAR) .or. ccMask(GAME_VAR) ) then
     ccMask(DENS_VAR) = .true.
     ccMask(TEMP_VAR) = .true.
     ! material information (ye, yi) stored in flame variable
     ccMask(FLAM_MSCALAR) = .true.
  endif

end subroutine Simulation_guardCellMaskHook

