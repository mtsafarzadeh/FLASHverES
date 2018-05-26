!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_Roe
!!
!! NAME
!!
!!  hy_uhd_Roe
!!
!! SYNOPSIS
!!
!!  hy_uhd_Roe( integer(IN) :: dir,
!!              real(IN)    :: Vm(HY_VARINUMMAX),
!!              real(IN)    :: Vp(HY_VARINUMMAX),
!!              real(OUT)   :: Fstar(HY_VARINUM))
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
!!  * Roe, JCP, 43:357-372, 1981
!!
!!***

Subroutine hy_uhd_Roe(dir,Vm,Vp,Fstar,ierr)

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
  real, dimension(HY_VARINUMMAX), intent(IN)  :: Vm, Vp
  real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real,dimension(HY_VARINUM2)  :: Vavg !DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM)   :: Um,Up,sigF,vec,FL,FR
  real,dimension(HY_WAVENUM)   :: lambda,lambdaL,lambdaR
  real,dimension(MDIM)         :: beta
  real,dimension(HY_WAVENUM,HY_VARINUM):: leig
  real,dimension(HY_VARINUM,HY_WAVENUM):: reig
  integer :: k
  logical :: cons
  real :: cs,ca,cf,as,af,bbN,uN
  real :: nu
  real :: aL2,aR2

!!$real :: avgPres,avgMagz,avgVelz

  ierr = 0 ! Set no error

  aL2 = Vm(HY_GAMC)*Vm(HY_PRES)/Vm(HY_DENS)
  aR2 = Vp(HY_GAMC)*Vp(HY_PRES)/Vp(HY_DENS)

  if (aL2 < 0. .or. aR2 < 0.) then
     ierr = 1
     return
  endif

  ! Godunov flux
  cons=.true.
  call hy_uhd_avgState(dir,Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))
  call hy_uhd_eigenParameters(Vavg(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambda,uN,cf,C_alfn=ca,C_slow=cs)
  call hy_uhd_eigenVector(leig,reig,Vavg(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)


  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:HY_MAGZ))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F08MAGZ_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F08MAGZ_FLUX))

  if (hy_entropy) then
     cons=.false.
     call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
     call hy_uhd_eigenValue(lambdaL,uN,cf,C_alfn=ca,C_slow=cs)
     call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
     call hy_uhd_eigenValue(lambdaR,uN,cf,C_alfn=ca,C_slow=cs)
     call hy_uhd_entropyFix(lambda,lambdaL,lambdaR)
  endif

  sigF =0.
  vec  =0.

  do k=1,HY_WAVENUM
     vec(HY_DENS:HY_MAGZ)&
          = abs(lambda(k))*reig(HY_DENS:HY_MAGZ,k)&
          *dot_product(leig(k,HY_DENS:HY_MAGZ),Up(HY_DENS:HY_MAGZ)-Um(HY_DENS:HY_MAGZ))

     sigF(F01DENS_FLUX:F08MAGZ_FLUX) = sigF(F01DENS_FLUX:F08MAGZ_FLUX) + vec(HY_DENS:HY_MAGZ)
  enddo

  ! Godunov Flux
  Fstar(F01DENS_FLUX:F08MAGZ_FLUX) = .5*(FR(F01DENS_FLUX:F08MAGZ_FLUX) &
                                       + FL(F01DENS_FLUX:F08MAGZ_FLUX) &
                                     - sigF(F01DENS_FLUX:F08MAGZ_FLUX))


!!$avgPres=(Vp(HY_PRES)+Vm(HY_PRES))*0.5
!!$avgMagz=(Vp(HY_MAGZ)+Vm(HY_MAGZ))**2*0.25
!!$avgVelz=((Vp(HY_VELZ)+Vm(HY_VELZ))**2)*(Vp(HY_DENS)+Vm(HY_DENS))*0.125


!!$if (NDIM == 2) then
!!$   if ( avgMagz**2 < 1.e-16*avgPres .and. &
!!$        avgVelz**2 < 1.e-16*avgPres) then
!!$      Fstar(4) = 0.
!!$      Fstar(8) = 0.
!!$   endif
!!$endif
End Subroutine hy_uhd_Roe
