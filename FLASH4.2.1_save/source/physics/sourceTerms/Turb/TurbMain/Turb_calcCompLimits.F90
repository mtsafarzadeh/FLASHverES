! See Turb_interface for subroutine description.
!
! Aaron Jackson 2010
!

#include "constants.h"
#include "Flash.h"
subroutine Turb_calcCompLimits(blkLimits, compLimits, local_stepSize)

  use Turb_data, only : turb_stepSize, turb_fillGC

  implicit none
  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: blkLimits
  integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: compLimits
  integer, intent(IN) :: local_stepSize

  ! check if we are not at flame refinement level
  ! or if we are performing guard cell fill between operators
  if ( turb_fillGC .or. (local_stepSize .ne. turb_stepSize) ) then

     ! we only need to compute for interior cells
     compLimits(:,:) = blkLimits(:,:)

  else

     ! we both operators can be performed without guard cell filling
     ! need to compute first operator into guard cells for second
     ! operator
     compLimits(LOW,IAXIS) = blkLimits(LOW,IAXIS)-2*turb_stepSize
     compLimits(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+2*turb_stepSize
     compLimits(LOW,JAXIS) = blkLimits(LOW,JAXIS)-K2D*2*turb_stepSize
     compLimits(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+K2D*2*turb_stepSize
     compLimits(LOW,KAXIS) = blkLimits(LOW,KAXIS)-K3D*2*turb_stepSize
     compLimits(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+K3D*2*turb_stepSize

  endif

  return
end subroutine Turb_calcCompLimits
