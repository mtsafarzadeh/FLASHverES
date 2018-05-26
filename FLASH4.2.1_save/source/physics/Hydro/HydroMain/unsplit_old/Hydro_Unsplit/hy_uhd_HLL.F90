!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/hy_uhd_HLL
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
!!  Vm     - primitive variable for left state (DENS,VELX,VELY,VELZ,PRES + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state(DENS,VELX,VELY,VELZ,PRES + GAMC,GAME,EINT,TEMP)
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
!! REFERENCES
!!
!!  * Harten, Lax and van Leer, SIAM  Review, 25(1):35--61, 1983
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***

Subroutine hy_uhd_HLL(dir,Vm,Vp,Fstar,ierr)
  
  use hy_uhd_interface, ONLY : hy_uhd_prim2con,hy_uhd_prim2flx
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM),   intent(OUT):: Fstar
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real :: densL,presL,velxL,velyL,velzL
  real :: densR,presR,velxR,velyR,velzR
  real :: SL,SR,cfL,cfR,aL2,aR2,gamcL,gamcR,velNL,velNR
  real, dimension(HY_VARINUM) :: UL,UR,FL,FR

  ierr = 0 ! Set no error

  ! Begin HLL
  densL = Vm(HY_DENS)
  velxL = Vm(HY_VELX)
  velyL = Vm(HY_VELY)
  velzL = Vm(HY_VELZ)
  presL = Vm(HY_PRES)
  gamcL = Vm(HY_GAMC)

  densR = Vp(HY_DENS)
  velxR = Vp(HY_VELX)
  velyR = Vp(HY_VELY)
  velzR = Vp(HY_VELZ)
  presR = Vp(HY_PRES)
  gamcR = Vp(HY_GAMC)

  ! Get normal velocity component
  select case (dir)
  case (DIR_X)
     velNL = velxL
     velNR = velxR
  case (DIR_Y)
     velNL = velyL
     velNR = velyR
  case (DIR_Z)
     velNL = velzL
     velNR = velzR
  end select


  ! Some parameters
  aL2   = gamcL*presL/densL
  aR2   = gamcR*presR/densR

  if (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  else
     cfL   = sqrt(aL2)
     cfR   = sqrt(aR2)
  endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:HY_ENER))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:HY_ENER))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F05ENER_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F05ENER_FLUX))


  if (SL > 0.) then
     !Ustar(HY_DENS:HY_PRES) = UL(HY_DENS:HY_PRES)
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = FL(F01DENS_FLUX:F05ENER_FLUX)
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
     !Ustar(HY_DENS:HY_PRES) = (SR*UR - SL*UL - FR + FL)/(SR - SL)
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = (SR*FL(F01DENS_FLUX:F05ENER_FLUX) &
                                       - SL*FR(F01DENS_FLUX:F05ENER_FLUX) &
                                       + SR*SL*(UR(HY_DENS:HY_ENER) - UL(HY_DENS:HY_ENER)))/(SR - SL)
  else
     !Ustar(HY_DENS:HY_PRES) = UR(HY_DENS:HY_PRES)
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = FR(F01DENS_FLUX:F05ENER_FLUX)
  endif

End Subroutine hy_uhd_HLL
