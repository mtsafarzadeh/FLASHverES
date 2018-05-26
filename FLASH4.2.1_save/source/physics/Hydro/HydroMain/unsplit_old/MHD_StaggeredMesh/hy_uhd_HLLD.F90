!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_HLLD
!!
!! NAME
!!
!!  hy_uhd_HLLD
!!
!! SYNOPSIS
!!
!!  hy_uhd_HLLD( integer(IN) :: dir,
!!               real(IN)    :: Vm(HY_VARINUM4),
!!               real(IN)    :: Vp(HY_VARINUM4),
!!               real(OUT)   :: Fstar(HY_VARINUM),
!!               real(OUT)   :: vint  )
!!
!! ARGUMENTS
!!
!!  dir    - a spatial direction for which the flux is being considered and computed
!!  Vm     - primitive variable for left state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES,MAGX,MAGY,MAGZ + GAMC,GAME,EINT,TEMP)
!!  Fstar  - computed flux data
!!  vint   - interface velocity
!!
!! DESCRIPTION
!! 
!!   This routine computes high-order Godunov fluxes based on the left and right Riemann states.
!!
!!   The HLLD Riemann fan:
!!
!!            SL*       SM       SR*
!!   SL        \        |        /        SR      
!!     *        \       |       /        *
!!       *   UL* \ UL** | UR** / UR*   *
!!         *      \     |     /      *
!!           *     \    |    /     *
!!             *    \   |   /    *
!!           UL  *   \  |  /   *   UR
!!                 *  \ | /  *
!!                   * \|/ *
!!   --------------------------------------
!!
!! REFERENCE
!!
!!  * Miyoshi & Kusano, JCP, 208:315-344, 2005
!!
!!***

Subroutine hy_uhd_HLLD(dir,Vm,Vp,Fstar,ierr)
  
  use hy_uhd_interface, ONLY : hy_uhd_prim2con,hy_uhd_prim2flx,hy_uhd_HLLC
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM), intent(OUT):: Fstar
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real, parameter :: epsilon = 1.e-4
  real :: densL,presL,velxL,velyL,velzL,magxL,magyL,magzL
  real :: densR,presR,velxR,velyR,velzR,magxR,magyR,magzR
  real :: cfL,cfR,aL2,aR2,gamcL,gamcR,magBL2,magBR2
  real :: magnL,magnR,magtL,magtR,velnL,velnR
  real :: SM,SL,SR,SL2,SR2
  real :: dStarL,dStarR,prestL,prestR,pres,Bn_hll
  real :: scrch1L,scrch1R,scrch2L,scrch2R
  real :: scrch3L,scrch3R,scrch4L,scrch4R,scrch5L,scrch5R
  real :: velxStarL,velyStarL,velzStarL,velxStarR,velyStarR,velzStarR
  real :: magxStarL,magyStarL,magzStarL,magxStarR,magyStarR,magzStarR
  real :: velxStar2,velyStar2,velzStar2,magxStar2,magyStar2,magzStar2
  real :: signumBn
  real, dimension(HY_VARINUM) :: UL,UR,FL,FR,Uhll,UCstarL,UCstarR,UCstar2L,UCstar2R
  real :: BxStar,ByStar,BzStar
  logical :: degeneracyHLLD

  ierr = 0 ! Set no error
  degeneracyHLLD=.false.


  ! Begin HLLD
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
     magnL = magxL
     magnR = magxR
     magtL = magyL*magyL+magzL*magzL
     magtR = magyR*magyR+magzR*magzR
     velnL = velxL
     velnR = velxR
  case (DIR_Y)
     magnL = magyL
     magnR = magyR
     magtL = magxL*magxL+magzL*magzL
     magtR = magxR*magxR+magzR*magzR
     velnL = velyL
     velnR = velyR
  case (DIR_Z)
     magnL = magzL
     magnR = magzR
     magtL = magxL*magxL+magyL*magyL
     magtR = magxR*magxR+magyR*magyR
     velnL = velzL
     velnR = velzR
  end select

  ! Some parameters
  magBL2= dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))/densL
  aL2   = gamcL*presL/densL
  cfL   = 0.5*(aL2+magBL2+sqrt((aL2+magBL2)**2-4.*aL2*magnL*magnL/densL))
  cfL   = sqrt(cfL)

  magBR2= dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))/densR
  aR2   = gamcR*presR/densR
  cfR   = 0.5*(aR2+magBR2+sqrt((aR2+magBR2)**2-4.*aR2*magnR*magnR/densR))
  cfR   = sqrt(cfR)

  if (aL2 < 0. .or. aR2 < 0.) then
     call Driver_abortFlash&
          ("[hy_uhd_HLLD]: Imaginary sound speed has obtained! "//&
           "Please try other (more diffusive) slope limiter, flux, order, cfl, etc. "//&
           "in order to increase numerical stability.")
  endif


  ! Get left/right going fastest wave speeds SL & SR for the left and right states
  ! by S. F. Davis, SIAM J. Sci. Stat, Comput., 9(1988) 445.
  ! Also see Miyoshi, Kusano, JCP, 208 (2005)
  SL = min(velnL - cfL, velnR - cfR)
  SR = max(velnL + cfL, velnR + cfR)


  ! Total pressure
  prestL = presL + 0.5*dot_product(Vm(HY_MAGX:HY_MAGZ),Vm(HY_MAGX:HY_MAGZ))
  prestR = presR + 0.5*dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))

  ! Convert primitive variables to conservative variables
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),UL(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),UR(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F08MAGZ_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F08MAGZ_FLUX))


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
  Bn_hll = Uhll(HY_PRES+dir) !=(SR*magnR-SL*magnL)/(SR-SL)
  BxStar = Uhll(HY_MAGX)     !BxStarL = BxStarR = BxHLL
  ByStar = Uhll(HY_MAGY)     !ByStarL = ByStarR = ByHLL
  BzStar = Uhll(HY_MAGZ)     !BzStarL = BzStarR = BzHLL


  !!***************************************
  !! (I)    UL* and UR* regions           *
  !!***************************************
  ! Normal velocity component and the middle wave SM
  ! SM = u*L = u**L = u*R = u**R
  SM =(densR*velnR*(SR-velnR)-densL*velnL*(SL-velnL)&
       -prestR+prestL-magnL*magnL+magnR*magnR)/&
      (densR*(SR-velnR)-densL*(SL-velnL))

  ! Convenient parameters
  scrch1L = SL-velnL
  scrch2L = SL-SM
  scrch3L = SM-velnL

  scrch1R = SR-velnR
  scrch2R = SR-SM
  scrch3R = SM-velnR

  ! Total pressure in the whole Riemann fan
  ! pres*L = pres*R = pres**L = pres**R = pres
  pres = densL*scrch1L*(SM-velnL)+prestL-magnL*magnL+Bn_hll*Bn_hll

  ! Densities in UL* and UR*
  dStarL = UL(HY_DENS)*scrch1L/scrch2L
  dStarR = UR(HY_DENS)*scrch1R/scrch2R

  SL2 = SM - abs(magnL)/sqrt(dStarL) ! = SL*
  SR2 = SM + abs(magnR)/sqrt(dStarR) ! = SR*

  ! Check if degeneracy happens: 
  ! This is the case of cf=ca, ca>a
  ! (1) transverse B=0 (note: in general, when trans B=0, cf=ca when ca>a; cs=a when ca<a), and
  ! (2) normal B .ne. 0, and strong to give ca>a.
  ! If this happens, we use HLLC (see PLUTO implementation as well)
  if ( (SL2-SL) < epsilon*(SM-SL) ) then
     degeneracyHLLD=.true.
  endif
  if ( (SR-SR2) < epsilon*(SR-SM) ) then
     degeneracyHLLD=.true.
  endif
  if (degeneracyHLLD) then
     call hy_uhd_HLLC(dir,Vm,Vp,Fstar,ierr)
     return
  endif

  ! Proceed to calculate left star regions if there is no degeneracy
  scrch4L = scrch3L/(densL*scrch1L*scrch2L-magnL*magnL)
  scrch5L = (densL*scrch1L*scrch1L-magnL*magnL)/(densL*scrch1L*scrch2L-magnL*magnL)

  select case (dir)
  case (DIR_X)
     ! Left primitive variables
     velxStarL = SM
     velyStarL = velyL - magxL*magyL*scrch4L
     velzStarL = velzL - magxL*magzL*scrch4L

     magxStarL = Bn_hll
     magyStarL = magyL*scrch5L
     magzStarL = magzL*scrch5L
  case (DIR_Y)
     ! Left
     velxStarL = velxL - magyL*magxL*scrch4L
     velyStarL = SM
     velzStarL = velzL - magyL*magzL*scrch4L

     magxStarL = magxL*scrch5L
     magyStarL = Bn_hll
     magzStarL = magzL*scrch5L
  case (DIR_Z)
     ! Left
     velxStarL = velxL - magzL*magxL*scrch4L
     velyStarL = velyL - magzL*magyL*scrch4L
     velzStarL = SM

     magxStarL = magxL*scrch5L
     magyStarL = magyL*scrch5L
     magzStarL = Bn_hll
  end select

  ! Left conserved variables
  UCstarL(HY_DENS)= dStarL
  UCstarL(HY_VELX)= UCstarL(HY_DENS)*velxStarL
  UCstarL(HY_VELY)= UCstarL(HY_DENS)*velyStarL
  UCstarL(HY_VELZ)= UCstarL(HY_DENS)*velzStarL

  UCstarL(HY_MAGX)= magxStarL
  UCstarL(HY_MAGY)= magyStarL
  UCstarL(HY_MAGZ)= magzStarL
  UCstarL(HY_ENER)= (scrch1L*UL(HY_ENER)-prestL*velnL+pres*SM&
                    +magnL*( Vm(HY_VELX)*Vm(HY_MAGX)&
                            +Vm(HY_VELY)*Vm(HY_MAGY)&
                            +Vm(HY_VELZ)*Vm(HY_MAGZ))&
                    -magnL*( UCstarL(HY_VELX)*UCstarL(HY_MAGX)&
                            +UCstarL(HY_VELY)*UCstarL(HY_MAGY)&
                            +UCstarL(HY_VELZ)*UCstarL(HY_MAGZ))&
                            /UCstarL(HY_DENS))/scrch2L


  ! Proceed to calculate right star regions if there is no degeneracy
  scrch4R = scrch3R/(densR*scrch1R*scrch2R-magnR*magnR)
  scrch5R = (densR*scrch1R*scrch1R-magnR*magnR)/(densR*scrch1R*scrch2R-magnR*magnR)

  select case (dir)
  case (DIR_X)
     ! Right primitive variables
     velxStarR = SM
     velyStarR = velyR - magxR*magyR*scrch4R
     velzStarR = velzR - magxR*magzR*scrch4R

     magxStarR = Bn_hll
     magyStarR = magyR*scrch5R
     magzStarR = magzR*scrch5R
  case (DIR_Y)
     ! Right
     velxStarR = velxR - magyR*magxR*scrch4R
     velyStarR = SM
     velzStarR = velzR - magyR*magzR*scrch4R

     magxStarR = magxR*scrch5R
     magyStarR = Bn_hll
     magzStarR = magzR*scrch5R
  case (DIR_Z)
     ! Right
     velxStarR = velxR - magzR*magxR*scrch4R
     velyStarR = velyR - magzR*magyR*scrch4R
     velzStarR = SM

     magxStarR = magxR*scrch5R
     magyStarR = magyR*scrch5R
     magzStarR = Bn_hll
  end select

  ! Right conserved variables
  UCstarR(HY_DENS)= dStarR
  UCstarR(HY_VELX)= UCstarR(HY_DENS)*velxStarR
  UCstarR(HY_VELY)= UCstarR(HY_DENS)*velyStarR
  UCstarR(HY_VELZ)= UCstarR(HY_DENS)*velzStarR

  UCstarR(HY_MAGX)= magxStarR
  UCstarR(HY_MAGY)= magyStarR
  UCstarR(HY_MAGZ)= magzStarR
  UCstarR(HY_ENER)= (scrch1R*UR(HY_ENER)-prestR*velnR+pres*SM&
                    +magnR*( Vp(HY_VELX)*Vp(HY_MAGX)&
                            +Vp(HY_VELY)*Vp(HY_MAGY)&
                            +Vp(HY_VELZ)*Vp(HY_MAGZ))&
                    -magnR*( UCstarR(HY_VELX)*UCstarR(HY_MAGX)&
                            +UCstarR(HY_VELY)*UCstarR(HY_MAGY)&
                            +UCstarR(HY_VELZ)*UCstarR(HY_MAGZ))&
                            /UCstarR(HY_DENS))/scrch2R
  !! Done with calculating UL* and UR* regions !!




  !!***************************************
  !! (II)    UL** and UR** regions        *
  !!***************************************
  ! Densities
  UCstar2L(HY_DENS) = UCstarL(HY_DENS)
  UCstar2R(HY_DENS) = UCstarR(HY_DENS)

  scrch1L = sqrt(UCstarL(HY_DENS))
  scrch1R = sqrt(UCstarR(HY_DENS))
  scrch2L = 1./(scrch1L + scrch1R)
  scrch2R = scrch2L

  signumBn = sign(1.,Bn_hll)

  select case (dir)
  case (DIR_X)
     ! Left primitive variables
     velxStar2 = SM
     velyStar2 = (scrch1L*velyStarL+scrch1R*velyStarR&
                 +(UCstarR(HY_MAGY)-UCstarL(HY_MAGY))*signumBn)*scrch2L
     velzStar2 = (scrch1L*velzStarL+scrch1R*velzStarR&
                 +(UCstarR(HY_MAGZ)-UCstarL(HY_MAGZ))*signumBn)*scrch2L

     magxStar2 = Bn_hll
     magyStar2 = (scrch1L*magyStarR+scrch1R*magyStarL&
                 +scrch1L*scrch1R*(velyStarR-velyStarL)*signumBn)&
                 *scrch2L
     magzStar2 = (scrch1L*magzStarR+scrch1R*magzStarL&
                 +scrch1L*scrch1R*(velzStarR-velzStarL)*signumBn)&
                 *scrch2L

  case (DIR_Y)
     ! Left primitive variables
     velxStar2 = (scrch1L*velxStarL+scrch1R*velxStarR&
                 +(UCstarR(HY_MAGX)-UCstarL(HY_MAGX))*signumBn)*scrch2L
     velyStar2 = SM
     velzStar2 = (scrch1L*velzStarL+scrch1R*velzStarR&
                 +(UCstarR(HY_MAGZ)-UCstarL(HY_MAGZ))*signumBn)*scrch2L

     magxStar2 = (scrch1L*magxStarR+scrch1R*magxStarL&
                 +scrch1L*scrch1R*(velxStarR-velxStarL)*signumBn)&
                 *scrch2L
     magyStar2 = Bn_hll
     magzStar2 = (scrch1L*magzStarR+scrch1R*magzStarL&
                 +scrch1L*scrch1R*(velzStarR-velzStarL)*signumBn)&
                 *scrch2L

  case (DIR_Z)
     ! Left primitive variables
     velxStar2 = (scrch1L*velxStarL+scrch1R*velxStarR&
                 +(UCstarR(HY_MAGX)-UCstarL(HY_MAGX))*signumBn)*scrch2L
     velyStar2 = (scrch1L*velyStarL+scrch1R*velyStarR&
                 +(UCstarR(HY_MAGY)-UCstarL(HY_MAGY))*signumBn)*scrch2L
     velzStar2 = SM

     magxStar2 = (scrch1L*magxStarR+scrch1R*magxStarL&
                 +scrch1L*scrch1R*(velxStarR-velxStarL)*signumBn)&
                 *scrch2L
     magyStar2 = (scrch1L*magyStarR+scrch1R*magyStarL&
                 +scrch1L*scrch1R*(velyStarR-velyStarL)*signumBn)&
                 *scrch2L
     magzStar2 = Bn_hll

  end select

  ! Left conservative variables
  UCstar2L(HY_VELX) = UCstar2L(HY_DENS)*velxStar2
  UCstar2L(HY_VELY) = UCstar2L(HY_DENS)*velyStar2
  UCstar2L(HY_VELZ) = UCstar2L(HY_DENS)*velzStar2

  UCstar2L(HY_MAGX) = magxStar2
  UCstar2L(HY_MAGY) = magyStar2
  UCstar2L(HY_MAGZ) = magzStar2
  UCstar2L(HY_ENER) = UCstarL(HY_ENER)-sqrt(UCstarL(HY_DENS))*signumBn*&
                (( UCstarL(HY_VELX)*UCstarL(HY_MAGX)&
                  +UCstarL(HY_VELY)*UCstarL(HY_MAGY)&
                  +UCstarL(HY_VELZ)*UCstarL(HY_MAGZ))/UCstarL(HY_DENS)&
                -( UCstar2L(HY_VELX)*UCstar2L(HY_MAGX)&
                  +UCstar2L(HY_VELY)*UCstar2L(HY_MAGY)&
                  +UCstar2L(HY_VELZ)*UCstar2L(HY_MAGZ))/UCstar2L(HY_DENS))

  ! Right conservative variables
  UCstar2R(HY_VELX) = UCstar2R(HY_DENS)*velxStar2
  UCstar2R(HY_VELY) = UCstar2R(HY_DENS)*velyStar2
  UCstar2R(HY_VELZ) = UCstar2R(HY_DENS)*velzStar2

  UCstar2R(HY_MAGX) = magxStar2
  UCstar2R(HY_MAGY) = magyStar2
  UCstar2R(HY_MAGZ) = magzStar2
  UCstar2R(HY_ENER) = UCstarR(HY_ENER)+sqrt(UCstarR(HY_DENS))*signumBn*&
                (( UCstarR(HY_VELX)*UCstarR(HY_MAGX)&
                  +UCstarR(HY_VELY)*UCstarR(HY_MAGY)&
                  +UCstarR(HY_VELZ)*UCstarR(HY_MAGZ))/UCstarR(HY_DENS)&
                -( UCstar2R(HY_VELX)*UCstar2R(HY_MAGX)&
                  +UCstar2R(HY_VELY)*UCstar2R(HY_MAGY)&
                  +UCstar2R(HY_VELZ)*UCstar2R(HY_MAGZ))/UCstar2R(HY_DENS))

  !! END of calculating all HLLD states !!

  !!***************************************
  !! (III) HLLD fluxes                    *
  !!***************************************
  if (SL >= 0.) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX)

  elseif ((SL < 0.) .and. (SL2 >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SL*(UCstarL(HY_DENS:HY_MAGZ) - UL(HY_DENS:HY_MAGZ))

  elseif ((SL2 < 0.) .and. (SM >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FL(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SL2*(UCstar2L(HY_DENS:HY_MAGZ) - UCstarL(HY_DENS:HY_MAGZ))&
          +  SL*(UCstarL(HY_DENS:HY_MAGZ) - UL(HY_DENS:HY_MAGZ))

  elseif ((SM < 0.) .and. (SR2 >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SR2*(UCstar2R(HY_DENS:HY_MAGZ) - UCstarR(HY_DENS:HY_MAGZ))&
          +  SR*(UCstarR(HY_DENS:HY_MAGZ) - UR(HY_DENS:HY_MAGZ))

  elseif ((SR2 < 0.) .and. (SR >= 0.)) then
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX) &
          + SR*(UCstarR(HY_DENS:HY_MAGZ) - UR(HY_DENS:HY_MAGZ))

  else
     Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FR(F01DENS_FLUX:F08MAGZ_FLUX)

  endif

End Subroutine hy_uhd_HLLD
