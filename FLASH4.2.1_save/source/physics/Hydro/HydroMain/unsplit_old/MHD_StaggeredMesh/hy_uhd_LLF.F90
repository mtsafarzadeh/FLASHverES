!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_LLF
!!
!! NAME
!!
!!  hy_uhd_LLF
!!
!! SYNOPSIS
!!
!!  hy_uhd_LLF(integer(IN) :: dir,
!!             real(IN)    :: Vm(HY_VARINUMMAX),
!!             real(IN)    :: Vp(HY_VARINUMMAX),
!!             real(OUT)   :: Fstar(HY_VARINUM))
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
!!   This routine computes local Lax-Friedrichs fluxes that are conservative and consistent at intercells.
!!   This flux has the maximum amount of dissipation allowed by the stability condition and
!!   exhibits odd-even decoupling. (See e.g., Laney, Computational Gasdynamics, pp.315)
!!
!! REFERENCE
!!
!!  * Laney, Computational Gasdynamics, Cambridge University Press
!!
!!***

Subroutine hy_uhd_LLF(dir,Vm,Vp,Fstar,ierr)

  use hy_uhd_interface, ONLY : hy_uhd_avgState,       &
                               hy_uhd_prim2con,       &
                               hy_uhd_prim2flx,       &
                               hy_uhd_eigenParameters,&
                               hy_uhd_eigenValue

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration -----------
  integer, intent(IN) :: dir
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM), intent(OUT):: Fstar
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real,dimension(HY_VARINUM2) :: Vavg !DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM)   :: Um,Up,FL,FR
  real,dimension(HY_WAVENUM)   :: lambda,lambdaL,lambdaR
  real,dimension(3)  :: beta
  real :: cs,ca,cf,as,af,bbN,uN
  logical :: cons
  real :: aL2,aR2

  ierr = 0 ! Set no error

  aL2 = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2 = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)

  if (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  endif


  cons = .true.
  call hy_uhd_avgState(dir,Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))
  call hy_uhd_eigenParameters(Vavg(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambda,uN,cf,C_alfn=ca,C_slow=cs)

  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F08MAGZ_FLUX))
  call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambdaL,uN,cf,C_alfn=ca,C_slow=cs)

  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F08MAGZ_FLUX))
  call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambdaR,uN,cf,C_alfn=ca,C_slow=cs)


  !! local Lax-Friedrichs Flux
  Fstar(F01DENS_FLUX:F08MAGZ_FLUX)  = &
       .5*(FR(F01DENS_FLUX:F08MAGZ_FLUX)  + FL(F01DENS_FLUX:F08MAGZ_FLUX) &
       -max(maxval(abs(lambda)),maxval(abs(lambdaL)),maxval(abs(lambdaR)))&
       *(Up(HY_DENS:HY_MAGZ)-Um(HY_DENS:HY_MAGZ)))

End Subroutine hy_uhd_LLF
