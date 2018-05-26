!!****f* source/physics/Gravity/Gravity_bhGenEwaldField
!!
!! NAME
!!
!!  Gravity_bhGenEwaldField
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhGenEwaldField()
!!
!! DESCRIPTION
!!
!!  Generates the Ewald field. It is done in parallel, the whole field is
!!  distributed among all cpus. The field is eventually saved in file 
!!  grv_bhEwaldFName. In case grv_bhEwaldAlwaysGenerate is FALSE and
!!  the file with the Ewald field exists, the Ewald field is read from it.
!!  Ewald field is stored on a nested grid (with the highest density at small
!!  distances), one level of the EF grid is calculated by
!!  grv_bhGenEwaldFieldLevel which is called by this subroutine.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***




subroutine Gravity_bhGenEwaldField()
  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  return
end subroutine Gravity_bhGenEwaldField

