!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/hy_uhd_Roe
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
!!            (DENS,VELX,VELY,VELZ,PRES + GAMC,GAME,EINT,TEMP)
!!  Vp     - primitive variable for right state
!!            (DENS,VELX,VELY,VELZ,PRES + GAMC,GAME,EINT,TEMP)
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
  real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
  real, dimension(HY_VARINUM),   intent(OUT):: Fstar
  integer, intent(OUT) :: ierr
  !! --------------------------------------

  real,dimension(HY_VARINUM2) :: Vavg !DENS,VELX,VELY,VELZ,PRES + GAMC,GAME
  real,dimension(HY_VARINUM)   :: Um,Up,sigF,vec,FL,FR
  real,dimension(HY_WAVENUM)   :: lambda,lambL,lambR
  real,dimension(HY_WAVENUM,HY_VARINUM) :: leig
  real,dimension(HY_VARINUM,HY_WAVENUM) :: reig
  integer :: k
  logical :: cons
  real    :: cf,uN
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
  call hy_uhd_avgState(dir,Vm(HY_DENS:HY_EINT),Vp(HY_DENS:HY_EINT),Vavg(HY_DENS:HY_GAME))
  call hy_uhd_eigenParameters(Vavg(HY_DENS:HY_GAME),dir,cons,uN,cf)
  call hy_uhd_eigenValue(lambda,uN,cf)
  call hy_uhd_eigenVector(leig,reig,Vavg(HY_DENS:HY_GAME),dir,cons,cf)

  call hy_uhd_prim2con(Vm(HY_DENS:HY_GAME),Um(HY_DENS:HY_ENER))
  call hy_uhd_prim2con(Vp(HY_DENS:HY_GAME),Up(HY_DENS:HY_ENER))
  call hy_uhd_prim2flx(dir,Vm(HY_DENS:HY_GAME),FL(F01DENS_FLUX:F05ENER_FLUX))
  call hy_uhd_prim2flx(dir,Vp(HY_DENS:HY_GAME),FR(F01DENS_FLUX:F05ENER_FLUX))

  if (hy_entropy) then
     cons=.false.
     ! Entropy fix for low density region
     call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,cons,uN,cf)
     call hy_uhd_eigenValue(lambL,uN,cf)
     call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,cons,uN,cf)
     call hy_uhd_eigenValue(lambR,uN,cf)
     call hy_uhd_entropyFix(lambda,lambL,lambR)
  endif

  sigF =0.
  vec  =0.

  do k=1,HY_WAVENUM
     vec(HY_DENS:HY_ENER)&
          = abs(lambda(k))*reig(HY_DENS:HY_PRES,k)&
           *dot_product(leig(k,HY_DENS:HY_PRES),Up(HY_DENS:HY_ENER)-Um(HY_DENS:HY_ENER))

     sigF(F01DENS_FLUX:F05ENER_FLUX) = sigF(F01DENS_FLUX:F05ENER_FLUX) + vec(HY_DENS:HY_ENER)
  enddo

  ! Godunov Flux
  Fstar(F01DENS_FLUX:F05ENER_FLUX) = .5*(FR(F01DENS_FLUX:F05ENER_FLUX) &
                                       + FL(F01DENS_FLUX:F05ENER_FLUX) &
                                     - sigF(F01DENS_FLUX:F05ENER_FLUX))
End Subroutine hy_uhd_Roe
