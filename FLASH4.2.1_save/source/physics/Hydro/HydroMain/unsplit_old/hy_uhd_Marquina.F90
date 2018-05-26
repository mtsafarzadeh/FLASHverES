!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_Marquina
!!
!! NAME
!!
!!  hy_uhd_Marquina
!!
!! SYNOPSIS
!!
!!  hy_uhd_Marquina( integer(IN) :: dir,
!!                   real(IN)    :: Vm(HY_VARINUMMAX),
!!                   real(IN)    :: Vp(HY_VARINUMMAX),
!!                   real(OUT)   :: Fstar(HY_VARINUM))
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  Vm     - primitive variable for left state  
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state 
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!
!! DESCRIPTION
!! 
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!! REFERENCES
!!
!!  * Donat and Marquina, JCP, 125:42-58, 1996
!!  * Stiriba and Donat, Computers & Mathematics with Applications, 46:719-739, 2003
!!  * Leveque, Mihalas, Dorfi, and Muller, Computational Methods for Astrophysical Flows, Springer, 1997
!!  * Cunningham et al., AstroBear
!!
!!***

Subroutine hy_uhd_Marquina(dir,Vm,Vp,Fstar)

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
  real, intent(OUT) :: vint
  !! --------------------------------------

End Subroutine hy_uhd_Marquina
