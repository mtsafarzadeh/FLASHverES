!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/hy_uhd_HLLC
!!
!! NAME
!!
!!  hy_uhd_HLLC
!!
!! SYNOPSIS
!!
!!  hy_uhd_HLLC( integer(IN) :: dir,
!!               real(IN)    :: Vm(HY_VARINUMMAX),
!!               real(IN)    :: Vp(HY_VARINUMMAX),
!!               real(OUT)   :: Fstar(HY_VARINUM))
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  Vm     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!
!! DESCRIPTION
!! 
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLLC Riemann fan:
!!
!!            SL      SM=qStar    SR
!!             \        |        /
!!              \       |       /
!!               \  UL* | UR*  /
!!                \     |     /
!!                 \    |    /
!!                  \   |   /
!!           UL      \  |  /       UR
!!                    \ | /
!!                     \|/
!!   --------------------------------------
!!
!! REFERENCES
!!
!!  * S. Li, JCP, 203:344-357, 2005
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!
!!***

Subroutine hy_uhd_HLLC(dir,Vm,Vp,Fstar,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_prim2con,hy_uhd_prim2flx
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
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
  real :: dStarL,dStarR,totalPresL,totalPresR
  real :: pStar,qStar
  real :: denomL,denomR,numerL,numerR
  real, dimension(HY_VARINUM) :: UL,UR,FL,FR,Uhll,UCstarR,UCstarL

  ierr = 0 ! Set no error

  ! Begin HLLC
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
     cfR   = sqrt(aR2)
     cfL   = sqrt(aL2)
  endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Total pressure
  totalPresL = presL
  totalPresR = presR

  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:HY_ENER))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:HY_ENER))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F05ENER_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F05ENER_FLUX))


  ! Get HLL states for later use
  if (SL > 0.) then
     Uhll(HY_DENS:HY_PRES) = UL(HY_DENS:HY_PRES)
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
     Uhll(HY_DENS:HY_PRES) = (SR*UR - SL*UL - FR + FL)/(SR - SL)
  else
     Uhll(HY_DENS:HY_PRES) = UR(HY_DENS:HY_PRES)
  endif

  ! (1) Normal velocity component
  ! qStarL = qStarR = qStar
  qStar =( densR*velNR*(SR-velNR) - densL*velNL*(SL-velNL)  &
          +totalPresL - totalPresR  )/&
         (densR*(SR-velNR) - densL*(SL-velNL))

  ! Convenient parameters
  numerL = SL-velNL
  denomL = SL-qStar
  numerR = SR-velNR
  denomR = SR-qStar

  ! (2) Total pressure
  ! pStarL = pStarR = pStar
  pStar = densL*numerL*(qStar-velNL)+totalPresL

  ! (3) Density
  dStarL = UL(HY_DENS)*numerL/denomL
  dStarR = UR(HY_DENS)*numerR/denomR

  ! (4) Conserved variables in the two-state (left & right) star regions
  UCstarL(HY_DENS)  = dStarL
  UCstarL(HY_ENER)  = UL(HY_ENER)*numerL/denomL + &
               ((pStar*qStar - totalPresL*velNL))/denomL

  UCstarR(HY_DENS)  = dStarR
  UCstarR(HY_ENER)  = UR(HY_ENER)*numerR/denomR + &
               ((pStar*qStar - totalPresR*velNR))/denomR


  select case (dir)
  case (DIR_X)
     UCstarL(HY_XMOM) = dStarL*qStar
     UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
     UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

     UCstarR(HY_XMOM) = dStarR*qStar
     UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
     UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

  case (DIR_Y)
     UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
     UCstarL(HY_YMOM) = dStarL*qStar
     UCstarL(HY_ZMOM) = UL(HY_ZMOM)*numerL/denomL

     UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
     UCstarR(HY_YMOM) = dStarR*qStar
     UCstarR(HY_ZMOM) = UR(HY_ZMOM)*numerR/denomR

  case (DIR_Z)
     UCstarL(HY_XMOM) = UL(HY_XMOM)*numerL/denomL
     UCstarL(HY_YMOM) = UL(HY_YMOM)*numerL/denomL
     UCstarL(HY_ZMOM) = dStarL*qStar

     UCstarR(HY_XMOM) = UR(HY_XMOM)*numerR/denomR
     UCstarR(HY_YMOM) = UR(HY_YMOM)*numerR/denomR
     UCstarR(HY_ZMOM) = dStarR*qStar
  end select
  ! End of calculating HLLC intermediate states ---------------------------


  ! (5) Finally, calculate HLLC fluxes
  if (SL >= 0.) then
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = FL(F01DENS_FLUX:F05ENER_FLUX)

  elseif ((SL < 0.).and. (qStar >= 0.)) then
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = FL(F01DENS_FLUX:F05ENER_FLUX) &
          + SL*(UCstarL(HY_DENS:HY_ENER) - UL(HY_DENS:HY_ENER))

  elseif ((qStar <0.) .and. (SR >= 0.)) then
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = FR(F01DENS_FLUX:F05ENER_FLUX) &
          + SR*(UCstarR(HY_DENS:HY_ENER) - UR(HY_DENS:HY_ENER))

  else
     Fstar(F01DENS_FLUX:F05ENER_FLUX) = FR(F01DENS_FLUX:F05ENER_FLUX)

  endif

End Subroutine hy_uhd_HLLC
