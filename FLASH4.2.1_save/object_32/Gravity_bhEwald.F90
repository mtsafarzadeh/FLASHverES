!!****f* source/physics/Gravity/Gravity_bhEwald
!!
!! NAME
!!
!!  Gravity_bhEwald
!!
!!
!! SYNOPSIS
!!
!!   real field = Gravity_bhEwald(
!!                           real(in)    :: x,
!!                           real(in)    :: y,
!!                           real(in)    :: z
!!        )
!!
!! DESCRIPTION
!!
!!   Interpolates in the Ewald field and returns its value for the point x,y,z.
!!
!! ARGUMENTS
!!
!!  x   : x-coordinate of the point where the Ewald field is determined
!!  y   : y-coordinate of the point where the Ewald field is determined
!!  z   : z-coordinate of the point where the Ewald field is determined
!!
!! RESULT
!!
!!  Value of the Ewald field in a given point.
!!
!! NOTES
!!
!!***



real function Gravity_bhEwald(x, y, z)
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real,intent(in) :: x, y, z

  Gravity_bhEwald = 0.0
  return
end

