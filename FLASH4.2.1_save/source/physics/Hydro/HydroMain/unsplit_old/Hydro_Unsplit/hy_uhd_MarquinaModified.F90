!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/hy_uhd_MarquinaModified
!!
!! NAME
!!
!!  hy_uhd_MarquinaModified
!!
!! SYNOPSIS
!!
!!  hy_uhd_MarquinaModified( integer(IN) :: dir,
!!                           real(IN)    :: Vm(HY_VARINUM+4),
!!                           real(IN)    :: Vp(HY_VARINUM+4),
!!                           real(OUT)   :: Fstar(HY_VARINUM))
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
#define CUSTOM
Subroutine hy_uhd_MarquinaModified(dir,Vm,Vp,Fstar) !,ix,iy,iz,blkLimitsGC,LambdaFix)

  use hy_uhd_interface, ONLY : hy_uhd_avgState,       &
                               hy_uhd_eigenParameters,&
                               hy_uhd_eigenValue,     &
                               hy_uhd_eigenVector,    &
                               hy_uhd_prim2con,       &
                               hy_uhd_prim2flx

  use Hydro_data,       ONLY : hy_entropy,hy_shockInstabilityFix
!!$  use hy_uhd_slopeLimiters, ONLY : minmod,mc,vanLeer
  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
!!$  integer, intent(IN) :: ix,iy,iz
!!$  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
!!$#ifdef FIXEDBLOCKSIZE
!!$  real, dimension(NDIM,2, &
!!$                  GRID_ILO_GC:GRID_IHI_GC, &
!!$                  GRID_JLO_GC:GRID_JHI_GC, &
!!$                  GRID_KLO_GC:GRID_KHI_GC),&
!!$                  intent(IN) :: LambdaFix
!!$#else
!!$  real, dimension(NDIM,2,&
!!$                  blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!$                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!$                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
!!$                  intent(IN) :: LambdaFix
!!$#endif
  !! --------------------------------------

  real,dimension(HY_VARINUM2) :: Vavg,Uavg!DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM)   :: Um,Up,sigF,vec,FL,FR,FLtemp,FRtemp,vecL,vecR
  real,dimension(HY_WAVENUM)   :: lambda,lambdaL,lambdaR,lambL,lambR
  real,dimension(MDIM)         :: beta
  real,dimension(HY_WAVENUM,HY_VARINUM):: leig,leigL,leigR
  real,dimension(HY_VARINUM,HY_WAVENUM):: reig,reigL,reigR
  integer :: k
  logical :: cons
  real :: cs,ca,cf,as,af,bbN,uN


  ! Second order Godunov flux
  cons=.true.

!!$  ! Avg state for MHD
!!$  call hy_uhd_avgState(Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))
!!$  call hy_uhd_eigenParameters(Vavg(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
!!$  call hy_uhd_eigenValue(lambda,uN,cf,C_alfn=ca,C_slow=cs)
!!$  call hy_uhd_eigenVector(leig,reig,Vavg(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
!!$
!!$  vint = Vavg(dir+1)

  ! Left state
  call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,cons,uN,cf)
  call hy_uhd_eigenValue(lambdaL,uN,cf)
  call hy_uhd_eigenVector(leigL,reigL,Vm(HY_DENS:HY_GAME),dir,cons,cf)

  ! Right state
  call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,cons,uN,cf)
  call hy_uhd_eigenValue(lambdaR,uN,cf)
  call hy_uhd_eigenVector(leigR,reigR,Vp(HY_DENS:HY_GAME),dir,cons,cf)

!!$  vint = 0.5*(Vm(dir+1)+Vp(dir+1))

!!$  call hy_uhd_prim2con(Vavg(HY_DENS:HY_GAME),Uavg(HY_DENS:HY_ENER))

  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:HY_ENER))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:HY_ENER))

  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F05ENER_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F05ENER_FLUX))


!!$  ! averaging for more diffusion for MHD
!!$  leigL = leig
!!$  leigR = leig
!!$  reigL = reig
!!$  reigR = reig
!!$  lambdaL = lambda
!!$  lambdaR = lambda

#ifndef CUSTOM
  FLtemp = 0.
  FRtemp = 0.
  sigF = 0.

  do k=1,HY_WAVENUM
     if (lambdaL(k)*lambdaR(k) > 0.) then
        if (lambdaL(k) > 0.) then
           vecL(k)=dot_product(leigL(k,HY_DENS:HY_ENER),FL(F01DENS_FLUX:F05ENER_FLUX))
           vecR(k)=0.
        else
           vecL(k)=0.
           vecR(k)=dot_product(leigR(k,HY_DENS:HY_ENER),FR(F01DENS_FLUX:F05ENER_FLUX))
        endif
     else
        vecL(k)=0.5*(max(abs(lambdaL(k)),abs(lambdaR(k)))*dot_product(leigL(k,HY_DENS:HY_ENER),Um(HY_DENS:HY_ENER)))
        vecR(k)=0.5*(max(abs(lambdaL(k)),abs(lambdaR(k)))*dot_product(leigR(k,HY_DENS:HY_ENER),Up(HY_DENS:HY_ENER)))
!!$        vecL(HY_DENS:HY_ENER)=0.5*( dot_product(leigL(k,HY_DENS:HY_ENER),FL(F01DENS_FLUX:F05ENER_FLUX))&
!!$                                   +max(abs(lambdaL(k)),abs(lambdaR(k)))&
!!$                                   *dot_product(leigL(k,HY_DENS:HY_ENER),Um(HY_DENS:HY_ENER)))
!!$
!!$        vecR(HY_DENS:HY_ENER)=0.5*( dot_product(leigR(k,HY_DENS:HY_ENER),FR(F01DENS_FLUX:F05ENER_FLUX))&
!!$                                   -max(abs(lambdaL(k)),abs(lambdaR(k)))&
!!$                                   *dot_product(leigR(k,HY_DENS:HY_ENER),Up(HY_DENS:HY_ENER)))
     endif

!!$     FLtemp(F01DENS_FLUX:F05ENER_FLUX) = FLtemp(F01DENS_FLUX:F05ENER_FLUX)+vecL(HY_DENS:HY_ENER)*reigL(HY_DENS:HY_ENER,k)
!!$     FRtemp(F01DENS_FLUX:F05ENER_FLUX) = FRtemp(F01DENS_FLUX:F05ENER_FLUX)+vecR(HY_DENS:HY_ENER)*reigR(HY_DENS:HY_ENER,k)

     FLtemp(F01DENS_FLUX:F05ENER_FLUX) = &
           vecL(HY_FASTLEFT)*reigL(HY_DENS:HY_ENER,HY_FASTLEFT)&
          +vecL(HY_SLOWLEFT)*reigL(HY_DENS:HY_ENER,HY_SLOWLEFT)&
          +vecL(HY_ENTROPY )*reigL(HY_DENS:HY_ENER,HY_ENTROPY) &
          +vecL(HY_SLOWRGHT)*reigL(HY_DENS:HY_ENER,HY_SLOWRGHT)&
          +vecL(HY_FASTRGHT)*reigL(HY_DENS:HY_ENER,HY_FASTRGHT)

     FRtemp(F01DENS_FLUX:F05ENER_FLUX) = &
           vecR(HY_FASTLEFT)*reigR(HY_DENS:HY_ENER,HY_FASTLEFT)&
          +vecR(HY_SLOWLEFT)*reigR(HY_DENS:HY_ENER,HY_SLOWLEFT)&
          +vecR(HY_ENTROPY )*reigR(HY_DENS:HY_ENER,HY_ENTROPY) &
          +vecR(HY_SLOWRGHT)*reigR(HY_DENS:HY_ENER,HY_SLOWRGHT)&
          +vecR(HY_FASTRGHT)*reigR(HY_DENS:HY_ENER,HY_FASTRGHT)


!$     vec(HY_DENS:HY_ENER) = 0.5*(abs(lambda(k))*(1.-dt/del(1)*abs(lambda(k)))*&
!$          reig(HY_DENS:HY_ENER,k)*dot_product(leig(k,HY_DENS:HY_ENER),Uavg(HY_DENS:HY_ENER)))
!$
!$     sigF(F01DENS_FLUX:F05ENER_FLUX) = sigF(F01DENS_FLUX:F05ENER_FLUX) + vec(HY_DENS:HY_ENER)


  enddo

  Fstar(F01DENS_FLUX:F05ENER_FLUX) = FLtemp(F01DENS_FLUX:F05ENER_FLUX)+FRtemp(F01DENS_FLUX:F05ENER_FLUX)


  ! in conjunction with high-resolution method using flux-difference splitting
!$  Fstar(F01DENS_FLUX:F05ENER_FLUX) = Fstar(F01DENS_FLUX:F05ENER_FLUX) &
!$                                   - FL(F01DENS_FLUX:F05ENER_FLUX)&
!$                                   -sigF(F01DENS_FLUX:F05ENER_FLUX)

!$  Fstar(F01DENS_FLUX:F05ENER_FLUX) = FR(F01DENS_FLUX:F05ENER_FLUX)-Fstar(F01DENS_FLUX:F05ENER_FLUX) &
!$                                   -sigF(F01DENS_FLUX:F05ENER_FLUX)

#endif


#ifdef CUSTOM
  sigF =0.
  vec  =0.

  do k=1,HY_WAVENUM
     vec(HY_DENS:HY_ENER)&
!!$          = max(abs(lambdaL(k)),abs(lambdaR(k)))*& ! this produces a spottie region
          = 0.5*(abs(lambdaL(k))+abs(lambdaR(k)))*& ! this give circular symmetries
!!$          = max(abs(lambdaL(k)),abs(lambdaR(k)))*& ! this give circular symmetries
!!$          = mc(abs(lambdaL(k)),abs(lambdaR(k)))*& ! this is similar to the arithmetic avg
          ( reigR(HY_DENS:HY_ENER,k)*dot_product(leigR(k,HY_DENS:HY_ENER),Up(HY_DENS:HY_ENER)) &
           -reigL(HY_DENS:HY_ENER,k)*dot_product(leigL(k,HY_DENS:HY_ENER),Um(HY_DENS:HY_ENER)) )

     sigF(F01DENS_FLUX:F05ENER_FLUX) = sigF(F01DENS_FLUX:F05ENER_FLUX) + vec(HY_DENS:HY_ENER)
  enddo

  ! Godunov Flux
  Fstar(F01DENS_FLUX:F05ENER_FLUX) = .5*(   FR(F01DENS_FLUX:F05ENER_FLUX) &
                                         +  FL(F01DENS_FLUX:F05ENER_FLUX) &
                                         -sigF(F01DENS_FLUX:F05ENER_FLUX) )
#endif

End Subroutine hy_uhd_MarquinaModified
