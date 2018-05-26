! Aaron Jackson 2010
!
! The default behavior is to do nothing

subroutine fl_fsTFIFlameSpeedBlock(solnData, s, ds, dx, compLimits, &
                                   quench_limit)

#include "Flash.h"
#include "constants.h"

  implicit none

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(:,:,:), intent(inout) :: s
  real, dimension(:,:,:), intent(in) :: ds
  real, intent(in) :: dx
  integer, dimension(LOW:HIGH,MDIM), intent(in) :: compLimits
  real, dimension(:,:,:), intent(out), optional :: quench_limit

  if (present(quench_limit)) quench_limit(:,:,:) = -1.0e0

end subroutine fl_fsTFIFlameSpeedBlock
