!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_upwindTransverseFlux
!!
!! NAME
!!
!!  hy_uhd_upwindTransverseFlux
!!
!! SYNOPSIS
!!
!!  hy_uhd_upwindTransverseFlux( integer(IN) :: dir,
!!                               integer(IN) :: order,
!!                               real,pointer:: vm1,
!!                               real,pointer:: vc0,
!!                               real,pointer:: vp1,
!!                               real(IN)    :: lambda(HY_WAVENUM),
!!                               real(IN)    :: leig(HY_VARINUM,HY_WAVENUM),
!!                               real(IN)    :: reig(HY_VARINUM,HY_WAVENUM),
!!                               integer(IN) :: sigSize,
!!                               real(OUT)   :: sig,
!!                               logical(IN),optional :: speciesScalar)
!!
!!
!! ARGUMENTS
!!
!!  dir      - direction of transverse flux vector
!!  order    - order for spatial discretization
!!  vm1      - pointer for the data values at i-1 cell
!!  vc0      - pointer for the data values at i cell
!!  vp1      - pointer for the data values at i+1 cell
!!  lambda   - eigen values
!!  leig     - left eigen vectors
!!  reig     - right eigen vectors
!!  sigSize  - size of the transverse flux vector
!!  sig      - transverse flux vector
!!  speciesScalar - this is present when updating species and mass scalars
!!  
!!
!! DESCRIPTION
!!
!!  This routine advances species by locally an Eulerian algorithm in an unspit way.
!!
!!*** 

Subroutine hy_uhd_upwindTransverseFlux&
     (dir,order,vm1,vc0,vp1,lambda,leig,reig,sigSize,sig,speciesScalar)

  use Hydro_data,           ONLY : hy_transOrder, hy_useAuxEintEqn
  use hy_uhd_slopeLimiters, ONLY : minmod, mc, vanLeer, signum
  use Driver_interface,     ONLY : Driver_abortFlash

  implicit none

#include "UHD.h"
#include "Flash.h"

  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: dir,order
  real,pointer,dimension(:)  :: vm1,vc0,vp1
  real,intent(IN),dimension(HY_WAVENUM) :: lambda
  real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
  real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
  integer, intent(IN) :: sigSize
  real,intent(OUT), dimension(sigSize) :: sig
  logical, intent(IN), optional :: speciesScalar
  !!------------------------------------------------------------------------

  integer :: n,hyVarEnd,hyBeg,hyConEnd
  real, dimension(HY_END_VARS) :: TransFlux
  real, dimension(HY_WAVENUM)     :: delbar
  real :: sigSpd
  logical :: useSpeciesScalarLocal


  if (present(speciesScalar)) then
     useSpeciesScalarLocal = speciesScalar
  else
     useSpeciesScalarLocal = .FALSE.
  endif


  !! First set ranges for indices
  if ( .not. useSpeciesScalarLocal ) then
     ! (1) hydro/MHD updates (species updates are treated separately)
     hyBeg    = HY_DENS
     hyConEnd = HY_ENER
#ifdef FLASH_USM_MHD /* for USM-MHD */
     hyConEnd = HY_MAGZ
#endif
     hyVarEnd = hyConEnd

     ! (2) For gamc, game, eint, and 3T variables;
#ifndef FLASH_EOS_GAMMA
     !! Constant gammas don't need to be calculated anyway.
     hyVarEnd = HY_GAME
#endif
     ! Adjust hyEnd depending on updates of gamc, game, eint, grav & 3T vars
     if (hy_useAuxEintEqn) hyVarEnd = HY_EINT
#ifdef FLASH_UHD_3T
     hyVarEnd = HY_END_VARS
#endif

  else
     ! (2)for species and mass scalar advection
     !allocate(TransFlux(HY_NSPEC))
     hyBeg    = 1        !HY_SPEC_BEG
     hyConEnd = HY_NSPEC !HY_SPEC_END
     hyVarEnd = hyConEnd
  endif


  if (order > 1) then
     call Driver_abortFlash&
          ("[hy_uhd_upwindTransverseFlux]: No high-order transFlux support! Please set transOrder=1.")
  endif


  ! Initialize sig
  sig = 0.

  !! Perform transverse flux calculation now.
  if ( .not. useSpeciesScalarLocal ) then

     !! (1) For hydro/MHD variables:
     !!    (1a) Calculate upwind transverse fluxes for conservative variables
     do n=1,HY_WAVENUM
        ! Upwinding
        ! (NOTE: Using this If-else-endif is much faster than using the signum approach as before)
        if (lambda(n) > 0.) then
           TransFlux(hyBeg:hyConEnd) = vc0(hyBeg:hyConEnd)-vm1(hyBeg:hyConEnd)
        else
           TransFlux(hyBeg:hyConEnd) = vp1(hyBeg:hyConEnd)-vc0(hyBeg:hyConEnd)
        endif

        ! Make sigma sums (or the transverse fluxes) for primitive variables 
        ! except for gamc, game, eint, gravity, and 3T variables.
        ! gamc, game, eint, grav, and 3T are treated separately in the below, (1b).
        delbar(n) = dot_product(leig(hyBeg:hyConEnd,n),TransFlux(hyBeg:hyConEnd))
        TransFlux(hyBeg:hyConEnd) = lambda(n)*reig(hyBeg:hyConEnd,n)*delbar(n)
        sig(hyBeg:hyConEnd) = sig(hyBeg:hyConEnd) + TransFlux(hyBeg:hyConEnd)
     enddo ! End of do n=1,HY_WAVENUM

        !! (1) For hydro/MHD variables: 
        !!     (1b) Calculate transverse fluxes for gamc, game, eint, and 3T variables;
        !!          they are simply advected with flow velocity
#ifndef FLASH_EOS_GAMMA
        hyBeg = HY_GAMC
#else
        if (hy_useAuxEintEqn) then
           hyBeg = HY_EINT
        else
           hyBeg = HY_EINT+1
        endif
#endif

     if (hyVarEnd > hyConEnd .and. hyBeg .ge. hyVarEnd) then
        if (lambda(HY_ENTROPY) > 0.) then
           TransFlux(hyBeg:hyVarEnd) = vc0(hyBeg:hyVarEnd)-vm1(hyBeg:hyVarEnd)
        else
           TransFlux(hyBeg:hyVarEnd) = vp1(hyBeg:hyVarEnd)-vc0(hyBeg:hyVarEnd)
        endif
        sig(hyBeg:hyVarEnd) = lambda(HY_ENTROPY)*TransFlux(hyBeg:hyVarEnd)

        !! Gravity components should not have ANY transverse fluxes. They are source terms!
#ifdef GRAVITY
        sig(HY_GRAV) = 0.
#endif
     endif ! endif of if (hyVarEnd > hyConEnd)


  else !if ( .not. useSpeciesScalarLocal ) then
     !! (2) For species advection:
     !!     Transverse fluxes for species; they are simply advected with fluid velocity

     if (lambda(HY_ENTROPY) > 0.) then
        sig(hyBeg:hyConEnd) = vc0(hyBeg:hyConEnd)-vm1(hyBeg:hyConEnd)
     else
        sig(hyBeg:hyConEnd) = vp1(hyBeg:hyConEnd)-vc0(hyBeg:hyConEnd)

     endif
     sig(hyBeg:hyConEnd) = lambda(HY_ENTROPY)*sig(hyBeg:hyConEnd)
  endif



End Subroutine hy_uhd_upwindTransverseFlux
