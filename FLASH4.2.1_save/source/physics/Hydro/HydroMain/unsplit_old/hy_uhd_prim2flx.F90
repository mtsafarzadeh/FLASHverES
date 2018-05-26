!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_prim2flx
!!
!! NAME
!!
!!  hy_uhd_prim2flx
!!
!! SYNOPSIS
!!
!!  hy_uhd_prim2flx( integer(IN) :: dir,
!!                   real(IN)    :: V(HY_VARINUM2),
!!                   real(OUT)   :: F(HY_VARINUM))
!!
!! ARGUMENTS
!!
!! dir- directional index
!! V  - primitive variables  + GAMC,GAME
!! F  - flux
!!
!! DESCRIPTION
!!
!!  This routine calculates conversion from primitive variables to fluxes.
!!
!!***

Subroutine hy_uhd_prim2flx(dir,V,F)

  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN)  :: dir
  real, dimension(HY_VARINUM2), intent(IN) :: V
  real, dimension(HY_VARINUM),  intent(OUT) :: F
  !! --------------------------------------

  real  :: u2,E
  real  :: B2,UB,ptot

  B2 = 0.
  UB = 0.
  ptot = 0.

  u2 = dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))
  E   = 0.5*V(HY_DENS)*u2+V(HY_PRES)/(V(HY_GAME)-1.0)

#ifdef FLASH_USM_MHD
  B2 = dot_product(V(HY_MAGX:HY_MAGZ),V(HY_MAGX:HY_MAGZ))
  UB = dot_product(V(HY_VELX:HY_VELZ),V(HY_MAGX:HY_MAGZ))
  ptot= V(HY_PRES) + 0.5 *B2
  E   = E + 0.5*B2
#endif



  if (dir==DIR_X) then
     F(F01DENS_FLUX)= V(HY_DENS)*V(HY_VELX)
#ifndef FLASH_USM_MHD
     F(F01DENS_FLUX)= V(HY_DENS)*V(HY_VELX)
     F(F02XMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELX)+V(HY_PRES)
     F(F03YMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELY)
     F(F04ZMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELZ)
     F(F05ENER_FLUX)= (E+V(HY_PRES))*V(HY_VELX)
#else
     F(F02XMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELX)-V(HY_MAGX)*V(HY_MAGX)+ptot
     F(F03YMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELY)-V(HY_MAGX)*V(HY_MAGY)
     F(F04ZMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELZ)-V(HY_MAGX)*V(HY_MAGZ)
     F(F05ENER_FLUX)= (E+ptot)*V(HY_VELX)-V(HY_MAGX)*UB
     F(F06MAGX_FLUX)= 0.
     F(F07MAGY_FLUX)= V(HY_VELX)*V(HY_MAGY)-V(HY_VELY)*V(HY_MAGX)
     F(F08MAGZ_FLUX)= V(HY_VELX)*V(HY_MAGZ)-V(HY_VELZ)*V(HY_MAGX)
#endif

  elseif (dir==DIR_Y) then
     F(F01DENS_FLUX)= V(HY_DENS)*V(HY_VELY)
#ifndef FLASH_USM_MHD
     F(F01DENS_FLUX)= V(HY_DENS)*V(HY_VELY)
     F(F02XMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELX)
     F(F03YMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELY)+V(HY_PRES)
     F(F04ZMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELZ)
     F(F05ENER_FLUX)= (E+V(HY_PRES))*V(HY_VELY)
#else
     F(F02XMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELX)-V(HY_MAGY)*V(HY_MAGX)
     F(F03YMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELY)-V(HY_MAGY)*V(HY_MAGY)+ptot
     F(F04ZMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELZ)-V(HY_MAGY)*V(HY_MAGZ)
     F(F05ENER_FLUX)= (E+ptot)*V(HY_VELY)-V(HY_MAGY)*UB
     F(F06MAGX_FLUX)= V(HY_VELY)*V(HY_MAGX)-V(HY_VELX)*V(HY_MAGY)
     F(F07MAGY_FLUX)= 0.
     F(F08MAGZ_FLUX)= V(HY_VELY)*V(HY_MAGZ)-V(HY_VELZ)*V(HY_MAGY)
#endif

  elseif (dir==DIR_Z) then
     F(F01DENS_FLUX)= V(HY_DENS)*V(HY_VELZ)

#ifndef FLASH_USM_MHD
     F(F01DENS_FLUX)= V(HY_DENS)*V(HY_VELZ)
     F(F02XMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELX)
     F(F03YMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELY)
     F(F04ZMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELZ)+V(HY_PRES)
     F(F05ENER_FLUX)= (E+V(HY_PRES))*V(HY_VELZ)
#else
     F(F02XMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELX)-V(HY_MAGZ)*V(HY_MAGX)
     F(F03YMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELY)-V(HY_MAGZ)*V(HY_MAGY)
     F(F04ZMOM_FLUX)= F(F01DENS_FLUX)*V(HY_VELZ)-V(HY_MAGZ)*V(HY_MAGZ)+ptot
     F(F05ENER_FLUX)= (E+ptot)*V(HY_VELZ)-V(HY_MAGZ)*UB
     F(F06MAGX_FLUX)= V(HY_VELZ)*V(HY_MAGX)-V(HY_VELX)*V(HY_MAGZ)
     F(F07MAGY_FLUX)= V(HY_VELZ)*V(HY_MAGY)-V(HY_VELY)*V(HY_MAGZ)
     F(F08MAGZ_FLUX)= 0.0
#endif

  endif

End Subroutine hy_uhd_prim2flx
