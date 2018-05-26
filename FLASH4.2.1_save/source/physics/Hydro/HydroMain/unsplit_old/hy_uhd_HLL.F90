!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_HLL
!!
!! NAME
!!
!!  hy_uhd_HLL
!!
!! SYNOPSIS
!!
!!  hy_uhd_HLL( integer(IN) :: dir,
!!              real(IN)    :: Vm(HY_VARINUMMAX),
!!              real(IN)    :: Vp(HY_VARINUMMAX),
!!              real(OUT)   :: Fstar(HY_VARINUM))
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
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLL Riemann fan:
!!
!!            SL                  SR
!!             \                 /
!!              \               /
!!               \      U*     /
!!                \           /
!!                 \         /
!!                  \       /
!!           UL      \     /       UR
!!                    \   /
!!                     \ /
!!   --------------------------------------
!!
!!  This implementation is a stub implementation at this level.
!!
!! REFERENCES
!!
!!  * Harten, Lax and van Leer, SIAM  Review, 25(1):35--61, 1983
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***

Subroutine hy_uhd_HLL(dir,Vm,Vp,Fstar)
  
  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
  !! --------------------------------------

End Subroutine hy_uhd_HLL
