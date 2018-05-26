!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_Marquina
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

Subroutine hy_uhd_Marquina(dir,Vm,Vp,Fstar,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_avgState,       &
                               hy_uhd_eigenParameters,&
                               hy_uhd_eigenValue,     &
                               hy_uhd_eigenVector,    &
                               hy_uhd_prim2con,       &
                               hy_uhd_prim2flx

  use Hydro_data,       ONLY : hy_entropy

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

  real,dimension(HY_VARINUM2) :: Vavg,Uavg!DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM)   :: Um,Up,sigF,vec,FL,FR,FLtemp,FRtemp,vecL,vecR
  real,dimension(HY_WAVENUM)   :: lambdaL,lambdaR
  real,dimension(MDIM)         :: beta
  real,dimension(HY_WAVENUM,HY_VARINUM):: leig
  real,dimension(HY_VARINUM,HY_WAVENUM):: reig
  integer :: k
  logical :: cons
  real    :: cs,ca,cf,as,af,bbN,uN
  real    :: aL2,aR2

  ierr = 0 ! Set no error

  aL2 = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2 = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)

  if (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  endif

  ! Godunov flux
  cons=.true.

  ! Avg state for MHD
  call hy_uhd_avgState(dir,Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))
  call hy_uhd_eigenParameters(Vavg(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenVector(leig,reig,Vavg(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)

  ! Left state
  call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambdaL,uN,cf,C_alfn=ca,C_slow=cs)

  ! Right state
  call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambdaR,uN,cf,C_alfn=ca,C_slow=cs)

  ! Fluxes
  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F08MAGZ_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F08MAGZ_FLUX))

  FLtemp = 0.
  FRtemp = 0.


  ! Note: eigenvectors are not "sided (i.e., left and right)" local system but averaged
  !       at the interfaces. This averaging is used to introduce more diffusion for
  !       MHD system. The averaging is not required for gas dynamics.
  do k=1,HY_WAVENUM
     if (lambdaL(k)*lambdaR(k) > 0.) then
        if (lambdaL(k) > 0.) then
           vecL(k)=dot_product(leig(k,HY_DENS:HY_MAGZ),FL(F01DENS_FLUX:F08MAGZ_FLUX))
           vecR(k)=0.
        else
           vecL(k)=0.
           vecR(k)=dot_product(leig(k,HY_DENS:HY_MAGZ),FR(F01DENS_FLUX:F08MAGZ_FLUX))
        endif
     else

        vecL(k)=0.5*( dot_product(leig(k,HY_DENS:HY_MAGZ),FL(F01DENS_FLUX:F08MAGZ_FLUX))&
                                   +max(maxval(abs(lambdaL)),maxval(abs(lambdaR)))&
                                   *dot_product(leig(k,HY_DENS:HY_MAGZ),Um(HY_DENS:HY_MAGZ)))

        vecR(k)=0.5*( dot_product(leig(k,HY_DENS:HY_MAGZ),FR(F01DENS_FLUX:F08MAGZ_FLUX))&
                                   -max(maxval(abs(lambdaL)),maxval(abs(lambdaR)))&
                                   *dot_product(leig(k,HY_DENS:HY_MAGZ),Up(HY_DENS:HY_MAGZ)))
     endif
  enddo

  FLtemp(F01DENS_FLUX:F08MAGZ_FLUX) = &
        vecL(HY_FASTLEFT)*reig(HY_DENS:HY_MAGZ,HY_FASTLEFT)&
       +vecL(HY_ALFNLEFT)*reig(HY_DENS:HY_MAGZ,HY_ALFNLEFT)&
       +vecL(HY_SLOWLEFT)*reig(HY_DENS:HY_MAGZ,HY_SLOWLEFT)&
       +vecL(HY_ENTROPY )*reig(HY_DENS:HY_MAGZ,HY_ENTROPY) &
       +vecL(HY_SLOWRGHT)*reig(HY_DENS:HY_MAGZ,HY_SLOWRGHT)&
       +vecL(HY_ALFNRGHT)*reig(HY_DENS:HY_MAGZ,HY_ALFNRGHT)&
       +vecL(HY_FASTRGHT)*reig(HY_DENS:HY_MAGZ,HY_FASTRGHT)

  FRtemp(F01DENS_FLUX:F08MAGZ_FLUX) = &
        vecR(HY_FASTLEFT)*reig(HY_DENS:HY_MAGZ,HY_FASTLEFT)&
       +vecR(HY_ALFNLEFT)*reig(HY_DENS:HY_MAGZ,HY_ALFNLEFT)&
       +vecR(HY_SLOWLEFT)*reig(HY_DENS:HY_MAGZ,HY_SLOWLEFT)&
       +vecR(HY_ENTROPY )*reig(HY_DENS:HY_MAGZ,HY_ENTROPY) &
       +vecR(HY_SLOWRGHT)*reig(HY_DENS:HY_MAGZ,HY_SLOWRGHT)&
       +vecR(HY_ALFNRGHT)*reig(HY_DENS:HY_MAGZ,HY_ALFNRGHT)&
       +vecR(HY_FASTRGHT)*reig(HY_DENS:HY_MAGZ,HY_FASTRGHT)

  Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = FLtemp(F01DENS_FLUX:F08MAGZ_FLUX)+FRtemp(F01DENS_FLUX:F08MAGZ_FLUX)


End Subroutine hy_uhd_Marquina
