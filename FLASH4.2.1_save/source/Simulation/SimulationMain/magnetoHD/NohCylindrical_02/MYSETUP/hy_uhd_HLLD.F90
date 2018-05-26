!!****if* source/Simulation/SimulationMain/magnetoHD/NohCylindrical_02/MYSETUP/hy_uhd_HLLD
!!
!! NAME
!!
!!  hy_uhd_HLLD
!!
!! SYNOPSIS
!!
!!  hy_uhd_HLLD( integer(IN) :: dir,
!!               real(IN)    :: Vm(HY_VARINUMMAX),
!!               real(IN)    :: Vp(HY_VARINUMMAX),
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
  use hy_uhd_slopeLimiters, ONLY : signum

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

  real :: densL,presL,velxL,velyL,velzL,magxL,magyL,magzL
  real :: densR,presR,velxR,velyR,velzR,magxR,magyR,magzR
  real :: cfL,cfR,aL2,aR2,gamcL,gamcR,magBL2,magBR2
  real :: magnL,magnR,magtL,magtR,velnL,velnR
  real :: SM,SL,SR,SL2,SR2
  real :: dStarL,dStarR,prestL,prestR,pres,Bn_hll
  real :: scrch1L,scrch1R,scrch2L,scrch2R
  real :: scrch3L,scrch3R,scrch4L,scrch4R,scrch5L,scrch5R,denom
  real :: velxStarL,velyStarL,velzStarL,velxStarR,velyStarR,velzStarR
  real :: magxStarL,magyStarL,magzStarL,magxStarR,magyStarR,magzStarR
  real :: velxStar2,velyStar2,velzStar2,magxStar2,magyStar2,magzStar2
  real :: signumBn
  real, dimension(HY_VARINUM) :: UL,UR,FL,FR,Uhll,UCstarL,UCstarR,UCstar2L,UCstar2R
  logical :: degeneracyHLLD=.false.
  real :: epsilon=1.e-16
real :: avgPres,avgMagz,avgVelz

  ! Set no error
  ierr = 0

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
 

  magBR2= dot_product(Vp(HY_MAGX:HY_MAGZ),Vp(HY_MAGX:HY_MAGZ))/densR
  aR2   = gamcR*presR/densR
  cfR   = 0.5*(aR2+magBR2+sqrt((aR2+magBR2)**2-4.*aR2*magnR*magnR/densR))


  if (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  else
     cfL   = sqrt(cfL)
     cfR   = sqrt(cfR)
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


  !!***************************************
  !! (I)    UL* and UR* regions           *
  !!***************************************
  ! Normal velocity component and the middle wave SM
  ! SM = u*L = u**L = u*R = u**R
  SM =(densR*velnR*(SR-velnR)-densL*velnL*(SL-velnL)&
       -prestR+prestL-magnL*magnL+magnR*magnR)/&
      (densR*(SR-velnR)-densL*(SL-velnL))


  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! (Ia) First degeneracy check                  !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! This is the case when any of fields become zero.
  !! Compare B^2 fields with epsilon*pressure to 
  !! check numerical zero.
  scrch1L = epsilon*min(presL,presR)
  if ( magnL**2 .le. scrch1L .or. &
       magtL    .le. scrch1L .or. &
       magnR**2 .le. scrch1L .or. &
       magtR    .le. scrch1L) then
     degeneracyHLLD = .true.
  endif


  ! Convenient parameters
  scrch1L = SL-velnL
  scrch2L = SL-SM
  scrch3L = SM-velnL

  scrch1R = SR-velnR
  scrch2R = SR-SM
  scrch3R = SM-velnR


  ! Total pressure in the whole Riemann fan
  ! pres*L = pres*R = pres**L = pres**R = pres
  pres = scrch1R*densR*prestL &
        -scrch1L*densL*prestR &
        +densL*densR*scrch1R*scrch1L*(velnR-velnL)

  pres = pres/(scrch1R*densR-scrch1L*densL)

  ! Densities in UL* and UR*
  dStarL = UL(HY_DENS)*scrch1L/scrch2L
  dStarR = UR(HY_DENS)*scrch1R/scrch2R

  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! (Ib) Second degeneracy check                 !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! Check degeneracy cases on the left and right regions:
  !! This is when  (i) all transverse fields are zero, and
  !!              (ii) C_alfven > C_sound.
  !! Then we get a degeracy in which we have
  !!               (a) C_fast = C_alfven, and 
  !!               (b) C_slow = C_sound  

  !! *** Left star region ***
  denom=densL*scrch1L*scrch2L-magnL*magnL
  if (abs(denom) < epsilon*aL2*densL) then
     if (.not. degeneracyHLLD) degeneracyHLLD = .true.
  else
     scrch4L = scrch3L/denom
     scrch5L = (densL*scrch1L*scrch1L-magnL*magnL)/denom
  endif

  !! *** Right star region ***
  denom=densR*scrch1R*scrch2R-magnR*magnR
  if (abs(denom) < epsilon*aR2*densR) then
     if (.not. degeneracyHLLD) degeneracyHLLD = .true.
  else
     scrch4R = scrch3R/denom
     scrch5R = (densR*scrch1R*scrch1R-magnR*magnR)/denom
  endif


  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! (Ic) Proceed calculating left star regions   !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  ! Left primitive variables
  select case (dir)
  case (DIR_X)
     magxStarL = Bn_hll
     velxStarL = SM
     if (.not. degeneracyHLLD) then
        magyStarL = magyL*scrch5L
        magzStarL = magzL*scrch5L
        velyStarL = velyL - magxL*magyL*scrch4L
        velzStarL = velzL - magxL*magzL*scrch4L
     else
        ! if degeneracy happens, we use a more diffusive HLLC approach
        if (magtL > 0.) then
           magyStarL = Uhll(HY_MAGY)
           magzStarL = Uhll(HY_MAGZ)
           velyStarL = velyL + Bn_hll*(magyL-magyStarL)/scrch1L/densL
           velzStarL = velzL + Bn_hll*(magzL-magzStarL)/scrch1L/densL
        else
           magyStarL = 0.
           magzStarL = 0.
           velyStarL = velyL
           velzStarL = velzL
        endif
     endif
  case (DIR_Y)
     magyStarL = Bn_hll
     velyStarL = SM
     if (.not. degeneracyHLLD) then
        magxStarL = magxL*scrch5L
        magzStarL = magzL*scrch5L
        velxStarL = velxL - magyL*magxL*scrch4L
        velzStarL = velzL - magyL*magzL*scrch4L
     else
        if (magtL > 0.) then
           magxStarL = Uhll(HY_MAGX)
           magzStarL = Uhll(HY_MAGZ)
           velxStarL = velxL + Bn_hll*(magxL-magxStarL)/scrch1L/densL
           velzStarL = velzL + Bn_hll*(magzL-magzStarL)/scrch1L/densL
        else
           magxStarL = 0.
           magzStarL = 0.
           velxStarL = velxL
           velzStarL = velzL
        endif
     endif
  case (DIR_Z)
     magzStarL = Bn_hll
     velzStarL = SM
     if (.not. degeneracyHLLD) then
        magxStarL = magxL*scrch5L
        magyStarL = magyL*scrch5L
        velxStarL = velxL - magzL*magxL*scrch4L
        velyStarL = velyL - magzL*magyL*scrch4L
     else
        if (magtL > 0.) then
           magxStarL = Uhll(HY_MAGX)
           magyStarL = Uhll(HY_MAGY)
           velxStarL = velxL + Bn_hll*(magxL-magxStarL)/scrch1L/densL
           velyStarL = velyL + Bn_hll*(magyL-magyStarL)/scrch1L/densL
        else
           magxStarL = 0.
           magyStarL = 0.
           velxStarL = velxL
           velyStarL = velyL
        endif
     endif
  end select

  ! Left conserved variables
  UCstarL(HY_DENS)= dStarL
  UCstarL(HY_VELX)= UCstarL(HY_DENS)*velxStarL
  UCstarL(HY_VELY)= UCstarL(HY_DENS)*velyStarL
  UCstarL(HY_VELZ)= UCstarL(HY_DENS)*velzStarL

  UCstarL(HY_MAGX)= magxStarL
  UCstarL(HY_MAGY)= magyStarL
  UCstarL(HY_MAGZ)= magzStarL
  UCstarL(HY_ENER)= scrch1L*UL(HY_ENER)-prestL*velnL+pres*SM+&
                    Bn_hll*( Vm(HY_VELX)*Vm(HY_MAGX)&
                            +Vm(HY_VELY)*Vm(HY_MAGY)&
                            +Vm(HY_VELZ)*Vm(HY_MAGZ)&
                            -velxStarL*UCstarL(HY_MAGX)&
                            -velyStarL*UCstarL(HY_MAGY)&
                            -velzStarL*UCstarL(HY_MAGZ))
  UCstarL(HY_ENER) =  UCstarL(HY_ENER)/scrch2L



  !! +++++++++++++++++++++++++++++++++++++++++++++!
  !! (Id) Proceed calculating right star regions  !
  !! +++++++++++++++++++++++++++++++++++++++++++++!
  ! Right primitive variables
  select case (dir)
  case (DIR_X)
     magxStarR = Bn_hll
     velxStarR = SM
     if (.not. degeneracyHLLD) then
        magyStarR = magyR*scrch5R
        magzStarR = magzR*scrch5R
        velyStarR = velyR - magxR*magyR*scrch4R
        velzStarR = velzR - magxR*magzR*scrch4R
     else
        if (magtR > 0.) then
           ! if degeneracy happens, we use a more diffusive HLLC approach
           magyStarR = Uhll(HY_MAGY)
           magzStarR = Uhll(HY_MAGZ)
           velyStarR = velyR + Bn_hll*(magyR-magyStarR)/scrch1R/densR
           velzStarR = velzR + Bn_hll*(magzR-magzStarR)/scrch1R/densR
        else
           magyStarR = 0.
           magzStarR = 0.
           velyStarR = velyR
           velzStarR = velzR
        endif
     endif
  case (DIR_Y)
     magyStarR = Bn_hll
     velyStarR = SM
     if (.not. degeneracyHLLD) then
        magxStarR = magxR*scrch5R
        magzStarR = magzR*scrch5R
        velxStarR = velxR - magyR*magxR*scrch4R
        velzStarR = velzR - magyR*magzR*scrch4R
     else
        if (magtR > 0.) then
           magxStarR = Uhll(HY_MAGX)
           magzStarR = Uhll(HY_MAGZ)
           velxStarR = velxR + Bn_hll*(magxR-magxStarR)/scrch1R/densR
           velzStarR = velzR + Bn_hll*(magzR-magzStarR)/scrch1R/densR
        else
           magxStarR = 0.
           magzStarR = 0.
           velxStarR = velxR
           velzStarR = velzR
        endif
     endif
  case (DIR_Z)
     magzStarR = Bn_hll
     velzStarR = SM
     if (.not. degeneracyHLLD) then
        magxStarR = magxR*scrch5R
        magyStarR = magyR*scrch5R
        velxStarR = velxR - magzR*magxR*scrch4R
        velyStarR = velyR - magzR*magyR*scrch4R
     else
        if (magtR > 0.) then
           magxStarR = Uhll(HY_MAGX)
           magyStarR = Uhll(HY_MAGY)
           velxStarR = velxR + Bn_hll*(magxR-magxStarR)/scrch1R/densR
           velyStarR = velyR + Bn_hll*(magyR-magyStarR)/scrch1R/densR
        else
           magxStarR = 0.
           magyStarR = 0.
           velxStarR = velxR
           velyStarR = velyR
        endif
     endif
  end select

  ! Right conserved variables
  UCstarR(HY_DENS)= dStarR
  UCstarR(HY_VELX)= UCstarR(HY_DENS)*velxStarR
  UCstarR(HY_VELY)= UCstarR(HY_DENS)*velyStarR
  UCstarR(HY_VELZ)= UCstarR(HY_DENS)*velzStarR

  UCstarR(HY_MAGX)= magxStarR
  UCstarR(HY_MAGY)= magyStarR
  UCstarR(HY_MAGZ)= magzStarR
  UCstarR(HY_ENER)= scrch1R*UR(HY_ENER)-prestR*velnR+pres*SM+&
                    Bn_hll*(Vp(HY_VELX)*Vp(HY_MAGX)&
                           +Vp(HY_VELY)*Vp(HY_MAGY)&
                           +Vp(HY_VELZ)*Vp(HY_MAGZ)&
                           -velxStarR*UCstarR(HY_MAGX)&
                           -velyStarR*UCstarR(HY_MAGY)&
                           -velzStarR*UCstarR(HY_MAGZ))
  UCstarR(HY_ENER) = UCstarR(HY_ENER)/scrch2R
  !! Done with calculating UL* and UR* regions !!



  !!***************************************
  !! (II)    UL** and UR** regions        *
  !!***************************************
  ! First calculate SL* and SR*
  SL2 = SM - abs(Bn_hll)/sqrt(UCstarL(HY_DENS)) ! = SL*
  SR2 = SM + abs(Bn_hll)/sqrt(UCstarR(HY_DENS)) ! = SR*

  !! Last degeneracy check:
  !! This is when normal field becomes zero (i.e., C_alfven = 0).
  !! As a result, we get
  !!       (a) C_fast^2 = C_sound^2 + (B_trans/sqrt(dens))^2, and
  !!       (b) C_slow = C_alfven = 0.
  !! When this happens, we convert to HLLC as HLLD reduces to HLLC 
  !! which resolves the intermediate state when missing Alfven wave.
  if ( Bn_hll**2 .le. epsilon*pres) then
     call hy_uhd_HLLC(dir,Vm,Vp,Fstar,ierr)
     return
  else

     ! Densities
     UCstar2L(HY_DENS) = UCstarL(HY_DENS)
     UCstar2R(HY_DENS) = UCstarR(HY_DENS)

     scrch1L = sqrt(UCstarL(HY_DENS))
     scrch1R = sqrt(UCstarR(HY_DENS))
     scrch2L = 1./(scrch1L + scrch1R)
     scrch2R = scrch2L

     signumBn = signum(Bn_hll)

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
                         (velxStarL*UCstarL (HY_MAGX)&
                         +velyStarL*UCstarL (HY_MAGY)&
                         +velzStarL*UCstarL (HY_MAGZ)&
                         -velxStar2*UCstar2L(HY_MAGX)&
                         -velyStar2*UCstar2L(HY_MAGY)&
                         -velzStar2*UCstar2L(HY_MAGZ))

     ! Right conservative variables
     UCstar2R(HY_VELX) = UCstar2R(HY_DENS)*velxStar2
     UCstar2R(HY_VELY) = UCstar2R(HY_DENS)*velyStar2
     UCstar2R(HY_VELZ) = UCstar2R(HY_DENS)*velzStar2

     UCstar2R(HY_MAGX) = magxStar2
     UCstar2R(HY_MAGY) = magyStar2
     UCstar2R(HY_MAGZ) = magzStar2
     UCstar2R(HY_ENER) = UCstarR(HY_ENER)+sqrt(UCstarR(HY_DENS))*signumBn*&
                         (velxStarR*UCstarR (HY_MAGX)&
                         +velyStarR*UCstarR (HY_MAGY)&
                         +velzStarR*UCstarR (HY_MAGZ)&
                         -velxStar2*UCstar2R(HY_MAGX)&
                         -velyStar2*UCstar2R(HY_MAGY)&
                         -velzStar2*UCstar2R(HY_MAGZ))
     !! END of calculating all HLLD states !!
  endif

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

  !!!***************************************
  !! (IV) Enforce zero for corresponding   *
  !!      magnetic field components        *
  !!!***************************************
  Fstar(F06MAGX_FLUX+dir-1) = 0.


avgPres=(Vp(HY_PRES)+Vm(HY_PRES))*0.5
avgMagz=(Vp(HY_MAGZ)+Vm(HY_MAGZ))**2*0.25
avgVelz=((Vp(HY_VELZ)+Vm(HY_VELZ))**2)*(Vp(HY_DENS)+Vm(HY_DENS))*0.125


if (NDIM == 2) then
   if ( avgMagz**2 < 1.e-16*avgPres .and. &
        avgVelz**2 < 1.e-16*avgPres) then
      Fstar(4) = 0.
      Fstar(8) = 0.
   endif
endif

End Subroutine hy_uhd_HLLD
