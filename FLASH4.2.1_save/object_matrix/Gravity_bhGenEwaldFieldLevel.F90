!!****f* source/physics/Gravity/Gravity_bhGenEwaldFieldLevel
!!
!! NAME
!!
!!  Gravity_bhGenEwaldFieldLevel
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhGenEwaldFieldLevel(
!!                      integer(in) :: ewald_periodicity,
!!                      integer(in) :: nx,
!!                      integer(in) :: ny,
!!                      integer(in) :: nz,
!!                      real(in)    :: ewald_xmax,
!!                      real(in)    :: ewald_ymax,
!!                      real(in)    :: ewald_zmax,
!!                      real(inout) :: field(-1:,-1:,-1:)
!!        )
!!
!! DESCRIPTION
!!
!!   Generates one level of the nested grid of the Ewald field. Called by
!!   grv_bhGenerateEwaldField.
!!
!! ARGUMENTS
!!
!!   ewald_periodicity : bit array denoting which direction has periodic (1)
!!                       or isolated (0) boundaries
!!   nx : number of points in ewald field in direction x
!!   ny : number of points in ewald field in direction y
!!   nz : number of points in ewald field in direction z
!!   ewald_xmax : size of ewald field in direction x
!!   ewald_ymax : size of ewald field in direction y
!!   ewald_zmax : size of ewald field in direction z
!!   field - ewald field array
!!
!!***

subroutine Gravity_bhGenEwaldFieldLevel(nx, ny, nz, ewald_xmax, ewald_ymax, ewald_zmax, &
           & ewald_periodicity, field)
  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer, intent(IN) :: ewald_periodicity, nx, ny, nz
  real, intent(IN) ::  ewald_xmax, ewald_ymax, ewald_zmax
  real, intent(INOUT) :: field(-1:,-1:,-1:)

  return
end subroutine Gravity_bhGenEwaldFieldLevel

