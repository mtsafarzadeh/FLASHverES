!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_getRenormGroup(mscalar,group)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(out) ::group 
   integer, intent(in) :: mscalar
   integer :: ctr

   !! values is an array of possible input strings (with possibly extra spaces)
   !! groups is an array of corresponding answers

   integer, dimension(1:1), parameter :: values = (/&
   & 0 &
   &/)

   integer, dimension(1:1), parameter :: groups = (/&
   & 0 &
   &/)

   !! The group numbers correspond to the following names used in the Config files
   !! 
   !! This is for informational purposes only

   group = 0  ! By default group 0 doesn't exist
   do ctr = LBOUND(values,1), UBOUND(values,1)
      if ( values(ctr) .eq. mscalar ) group = groups(ctr)
   end do

end subroutine Simulation_getRenormGroup

