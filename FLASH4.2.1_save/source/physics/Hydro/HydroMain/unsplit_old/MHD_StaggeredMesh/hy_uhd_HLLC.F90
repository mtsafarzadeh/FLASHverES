!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_HLLC
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
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state 
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
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
  real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real :: densL,presL,velxL,velyL,velzL,magxL,magyL,magzL
  real :: densR,presR,velxR,velyR,velzR,magxR,magyR,magzR
  real :: SL,SR,cfL,cfR,aL2,aR2,gamcL,gamcR,magNL,magNR,magBL2,magBR2,velNL,velNR
  real :: dStarL,dStarR,totalPresL,totalPresR
  real :: BxStar,ByStar,BzStar,pStar,qStar,Bn_hll
  real :: denomL,denomR,numerL,numerR
  real, dimension(HY_VARINUM) :: UL,UR,FL,FR,Uhll,UCstarR,UCstarL

!!$  real :: avgPres,avgMagz,avgVelz

  ierr = 0 ! Set no error

  ! Begin HLLC
  densL = Vm(HY_DENS)
  velxL = Vm(HY_VELX)
  velyL = Vm(HY_VELY)
  velzL = Vm(HY_VELZ)
  magxL = Vm(HY_MAGX)
  magyL = Vm(HY_MAGY)
  magzL = Vm(HY_MAGZ)
  presL = Vm(HY_PRES)
  gamcL = Vm(HY_GAMC)

  densR = Vp(HY_DENS)
  velxR = Vp(HY_VELX)
  velyR = Vp(HY_VELY)
  velzR = Vp(HY_VELZ)
  magxR = Vp(HY_MAGX)
  magyR = Vp(HY_MAGY)
  magzR = Vp(HY_MAGZ)
  presR = Vp(HY_PRES)
  gamcR = Vp(HY_GAMC)

  ! Get normal velocity component
  select case (dir)
  case (DIR_X)
     magNL = magxL
     magNR = magxR
     velNL = velxL
     velNR = velxR
  case (DIR_Y)
     magNL = magyL
     magNR = magyR
     velNL = velyL
     velNR = velyR
  case (DIR_Z)
     magNL = magzL
     magNR = magzR
     velNL = velzL
     velNR = velzR
  end select


  ! Some parameters
  magBL2= dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))/densL
  aL2   = gamcL*presL/densL

  magBR2= dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))/densR
  aR2   = gamcR*presR/densR

  if (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  else
     cfL   = sqrt(0.5*( aL2 + magBL2 + sqrt((aL2 + magBL2 )**2 - 4.*aL2*magNL*magNL/densL)))
     cfR   = sqrt(0.5*( aR2 + magBR2 + sqrt((aR2 + magBR2 )**2 - 4.*aR2*magNR*magNR/densR)))
  endif

  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velNL - cfL, velNR - cfR)
  SR = max(velNL + cfL, velNR + cfR)

  ! Total pressure
  totalPresL = presL + 0.5*dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))
  totalPresR = presR + 0.5*dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))

  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F08MAGZ_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F08MAGZ_FLUX))


!!$if (NDIM == 2) then
!!$   if ( Vm(HY_MAGZ) == 0. .and. Vp(HY_MAGZ) == 0. .and. Vm(HY_VELZ) == 0. .and. Vp(HY_VELZ) == 0.) then
!!$      FL(8) = 0.
!!$      FR(8) = 0.
!!$      FL(4) = 0.
!!$      FR(4) = 0.
!!$
!!$      UL(8) = 0.
!!$      UR(8) = 0.
!!$      UL(4) = 0.
!!$      UR(4) = 0.
!!$
!!$
!!$   endif
!!$endif


!!$print*,'------------'
!!$print*,FL(F08MAGZ_FLUX)
!!$print*,FR(F08MAGZ_FLUX)
!!$print*,'------------'
!print*,'------------------'
!!$if (abs(FL(8)) > 0.) print*,'FL=',FL(8) 
!!$if (abs(FR(8)) > 0.) print*,'FR=',FR(8) 
!!$if (abs(UL(8)) > 0.) print*,'UL=',UL(8) 
!!$if (abs(UR(8)) > 0.) print*,'UR=',UR(8)
!print*,'------------------'
!!$FL(8) = 0.
!!$FR(8) = 0.
!!$UL(8) = 0.
!!$UR(8) = 0.



  ! Get HLL states for later use
  if (SL > 0.) then
     Uhll(HY_DENS:HY_MAGZ) = UL(HY_DENS:HY_MAGZ)
  elseif ((SL <= 0.) .and. (SR >= 0.)) then
     Uhll(HY_DENS:HY_MAGZ) = (SR*UR(HY_DENS:HY_MAGZ) &
                            - SL*UL(HY_DENS:HY_MAGZ) &
                            - FR(F01DENS_FLUX:F08MAGZ_FLUX) &
                            + FL(F01DENS_FLUX:F08MAGZ_FLUX))/(SR - SL)
  else
     Uhll(HY_DENS:HY_MAGZ) = UR(HY_DENS:HY_MAGZ)
  endif


  ! Calculate intermediate states ---------------------------------------------
  Bn_hll = Uhll(HY_PRES+dir) !=(SR*magNR-SL*magNL)/(SR-SL)
  BxStar = Uhll(HY_MAGX)     !BxStarL = BxStarR = BxHLL
  ByStar = Uhll(HY_MAGY)     !ByStarL = ByStarR = ByHLL
  BzStar = Uhll(HY_MAGZ)     !BzStarL = BzStarR = BzHLL

  ! (1) Normal velocity component
  ! qStarL = qStarR = qStar
  qStar =( densR*velNR*(SR-velNR) - densL*velNL*(SL-velNL)  &
          +totalPresL - totalPresR - magNL**2 + magNR**2  )/&
         (densR*(SR-velNR) - densL*(SL-velNL))

  ! Convenient parameters
  numerL = SL-velNL
  denomL = SL-qStar
  numerR = SR-velNR
  denomR = SR-qStar

  ! (2) Total pressure
  ! pStarL = pStarR = pStar
  pStar = densL*numerL*(qStar-velNL)+totalPresL-magNL**2+Bn_hll**2

  ! (3) Density
  dStarL = UL(HY_DENS)*numerL/denomL
  dStarR = UR(HY_DENS)*numerR/denomR

  ! (4) Conserved variables in the two-state (left & right) star regions
  UCstarL(HY_DENS)  = dStarL
  UCstarL(HY_MAGX:HY_MAGZ)= Uhll(HY_MAGX:HY_MAGZ)
  UCstarL(HY_ENER)  = UL(HY_ENER)*numerL/denomL + &
               ((pStar*qStar - totalPresL*velNL) - &
                (Bn_hll*dot_product(Uhll(HY_MAGX:HY_MAGZ),Uhll(HY_VELX:HY_VELZ))/Uhll(HY_DENS)-&
                  magNL*dot_product(Vm(HY_MAGX:HY_MAGZ),  Vm(HY_VELX:HY_VELZ)) ))/denomL

  UCstarR(HY_DENS)  = dStarR
  UCstarR(HY_MAGX:HY_MAGZ)= Uhll(HY_MAGX:HY_MAGZ)
  UCstarR(HY_ENER)  = UR(HY_ENER)*numerR/denomR + &
               ((pStar*qStar - totalPresR*velNR) - &
                (Bn_hll*dot_product(Uhll(HY_MAGX:HY_MAGZ),Uhll(HY_VELX:HY_VELZ))/Uhll(HY_DENS)-&
                  magNR*dot_product(Vp(HY_MAGX:HY_MAGZ),  Vp(HY_VELX:HY_VELZ)) ))/denomR


  select case (dir)
  case (DIR_X)
     UCstarL(HY_VELX) = dStarL*qStar
     UCstarL(HY_VELY) = UL(HY_VELY)*numerL/denomL - (BxStar*ByStar-magxL*magyL)/denomL
     UCstarL(HY_VELZ) = UL(HY_VELZ)*numerL/denomL - (BxStar*BzStar-magxL*magzL)/denomL

     UCstarR(HY_VELX) = dStarR*qStar
     UCstarR(HY_VELY) = UR(HY_VELY)*numerR/denomR - (BxStar*ByStar-magxR*magyR)/denomR
     UCstarR(HY_VELZ) = UR(HY_VELZ)*numerR/denomR - (BxStar*BzStar-magxR*magzR)/denomR

  case (DIR_Y)
     UCstarL(HY_VELX) = UL(HY_VELX)*numerL/denomL - (ByStar*BxStar-magyL*magxL)/denomL
     UCstarL(HY_VELY) = dStarL*qStar
     UCstarL(HY_VELZ) = UL(HY_VELZ)*numerL/denomL - (ByStar*BzStar-magyL*magzL)/denomL

     UCstarR(HY_VELX) = UR(HY_VELX)*numerR/denomR - (ByStar*BxStar-magyR*magxR)/denomR
     UCstarR(HY_VELY) = dStarR*qStar
     UCstarR(HY_VELZ) = UR(HY_VELZ)*numerR/denomR - (ByStar*BzStar-magyR*magzR)/denomR

  case (DIR_Z)
     UCstarL(HY_VELX) = UL(HY_VELX)*numerL/denomL - (BzStar*BxStar-magzL*magxL)/denomL
     UCstarL(HY_VELY) = UL(HY_VELY)*numerL/denomL - (BzStar*ByStar-magzL*magyL)/denomL
     UCstarL(HY_VELZ) = dStarL*qStar

     UCstarR(HY_VELX) = UR(HY_VELX)*numerR/denomR - (BzStar*BxStar-magzR*magxR)/denomR
     UCstarR(HY_VELY) = UR(HY_VELY)*numerR/denomR - (BzStar*ByStar-magzR*magyR)/denomR
     UCstarR(HY_VELZ) = dStarR*qStar
  end select
  ! End of calculating HLLC intermediate states ---------------------------
!UL(8) = 0.;
!UR(8) = 0.;
!UCstarL(8) = 0.;
!UCstarR(8) = 0.;
  ! (5) Finally, calculate HLLC fluxes
  if (SL >= 0.) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX)

  elseif ((SL < 0.).and. (qStar >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SL*(UCstarL(HY_DENS:HY_MAGZ) - UL(HY_DENS:HY_MAGZ))

  elseif ((qStar <0.) .and. (SR >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SR*(UCstarR(HY_DENS:HY_MAGZ) - UR(HY_DENS:HY_MAGZ))

  else
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX)

  endif
!Fstar(8) = 0.

!!$avgPres=(Vp(HY_PRES)+Vm(HY_PRES))*0.5
!!$avgMagz=(Vp(HY_MAGZ)+Vm(HY_MAGZ))**2*0.25
!!$avgVelz=((Vp(HY_VELZ)+Vm(HY_VELZ))**2)*(Vp(HY_DENS)+Vm(HY_DENS))*0.125
!!$
!!$
!!$if (NDIM == 2) then
!!$   if ( avgMagz**2 < 1.e-16*avgPres .and. &
!!$        avgVelz**2 < 1.e-16*avgPres) then
!!$      Fstar(4) = 0.
!!$      Fstar(8) = 0.
!!$   endif
!!$endif


End Subroutine hy_uhd_HLLC
