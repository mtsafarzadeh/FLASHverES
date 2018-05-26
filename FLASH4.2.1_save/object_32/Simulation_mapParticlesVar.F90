!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!

subroutine Simulation_mapParticlesVar(part_key, var_key, var_type)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(in)  :: part_key
   integer, intent(out) :: var_key, var_type


   integer :: ctr, tmp_var_key, tmp_var_type

   integer, dimension(1:1), parameter :: part_keys = (/&
   & NONEXISTENT&
   &/)

   integer, dimension(1:1), parameter :: var_keys = (/&
   & NONEXISTENT&
   &/)

   integer, dimension(1:1), parameter :: var_types = (/&
   & NONEXISTENT&
   &/)

   tmp_var_key = NONEXISTENT
   tmp_var_type = NONEXISTENT

   do ctr = LBOUND(part_keys,1), UBOUND(part_keys,1)
      if (part_keys(ctr) .eq. part_key) then
         tmp_var_key = var_keys(ctr)
         tmp_var_type = var_types(ctr)
      end if
   end do
   var_key = tmp_var_key
   var_type = tmp_var_type

end subroutine Simulation_mapParticlesVar


