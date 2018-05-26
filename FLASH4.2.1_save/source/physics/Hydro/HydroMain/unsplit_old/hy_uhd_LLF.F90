!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_LLF
!!
!! NAME
!!
!!  hy_uhd_LLF
!!
!! SYNOPSIS
!!
!!  hy_uhd_LLF( integer(IN) :: dir,
!!             real(IN)    :: Vm(HY_VARINUMMAX),
!!             real(IN)    :: Vp(HY_VARINUMMAX),
!!             real(OUT)   :: Fstar(HY_VARINUM))
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  Vm     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,(MAGX,MAGY,MAGZ) + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,(MAGX,MAGY,MAGZ) + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!
!! DESCRIPTION
!! 
!!   This routine computes local Lax-Friedrichs fluxes that are conservative and consistent at intercells.
!!   This flux has the maximum amount of dissipation allowed by the stability condition and
!!   exhibits odd-even decoupling. (See e.g., Laney, Computational Gasdynamics, pp.315)
!!
!!   This implementation is a stub implementation at this level.
!!
!! REFERENCE
!!
!!  * Laney, Computational Gasdynamics, Cambridge University Press
!!
!!***

Subroutine hy_uhd_LLF(dir,Vm,Vp,Fstar)

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM), intent(OUT):: Fstar
  !! --------------------------------------


End Subroutine hy_uhd_LLF
