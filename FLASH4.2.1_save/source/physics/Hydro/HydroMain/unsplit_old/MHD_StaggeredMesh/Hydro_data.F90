!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/Hydro_data
!!
!! NAME
!!
!!  Hydro_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE Hydro_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data that are specific to the Unsplit MHD_StaggeredMesh unit.
!!
!!***
 
 
Module Hydro_data

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  logical, save :: hy_killdivb,          hy_E_modification,        &
                   hy_fluxCorrect,       hy_charLimiting,          &
                   hy_entropy,           hy_forceHydroLimit,       &
                   hy_shockDetectOn,     hy_flattening,            &
                   hy_facevar2ndOrder,   hy_useResistivity,        &
                   hy_useDiffuse,        hy_useViscosity,          &
                   hy_useConductivity,   hy_useMagneticResistivity,&
                   hy_ContactSteepening, hy_EOSforRiemann,         &
                   hy_upwindTVD,         hy_use_avisc,             &
                   hy_energyFixSwitch,   hy_hybridRiemannOnly,     &
                   hy_gravConsv,         hy_useGravHalfUpdate,     &
                   hy_useGravPotUpdate,  hy_useGravity,            &
                   hy_updateHydroFluxes, hy_useHydro,              &
                   hy_use3dFullCTU,      hy_addThermalFlux,        &
                   hy_E_upwind,          hy_useHybridOrder,        &
                   hy_useVaryingCFL,     hy_conserveAngMom,        &
                   hy_useAuxEintEqn
                   

  integer, save :: hy_geometry, hy_limiter, hy_irenorm
  integer, save :: hy_unsplitEosMode = MODE_DENS_EI
  integer, save :: hy_EosModeAfter = MODE_DENS_EI
  integer, save :: hy_prol_method, hy_RiemannSolver,hy_RiemannSolverLoc, hy_entropyFixMethod
  integer, parameter :: hy_numXN = NSPECIES+NMASS_SCALARS

  character(len=MAX_STRING_LENGTH), save :: hy_limiter_str,      &
                                            hy_prol_method_str,  &
                                            hy_RiemannSolver_str,&
                                            hy_entropyFixMethod_str

  real, save :: hy_cfl
  real, save :: hy_cfl_original
  real, save :: hy_tiny=1.e-32
  real, save :: hy_Rconst
  real, save :: hy_eswitch
  real, save :: hy_smalldens
  real, save :: hy_smallpres
  real, save :: hy_smallE
  real, save :: hy_maxMagDiff
  real, save :: hy_LimitedSlopeBeta
  real, save :: hy_cvisc
  real, save :: hy_hybridOrderKappa

  ! System of units used
  character(4), save :: hy_units

  ! Everybody should know these!
  integer, save :: hy_meshNumProcs, hy_meshMe

  ! Constants for non-dimensionalization
  real, save :: hy_xref
  real, save :: hy_tref
  real, save :: hy_dref
  real, save :: hy_vref
  real, save :: hy_pref
  real, save :: hy_eref
  real, save :: hy_qref
  real, save :: hy_bref
  real, save :: hy_gref
  real, save :: hy_mref
  real, save :: hy_nref
  real, save :: hy_kref

  ! Order of Accuracy
  ! hy_order=1 ==> First order Godunov
  ! hy_order=2 ==> Muscl-Hancock
  ! hy_order=3 ==> PPM
  ! hy_order=5 ==> WENO: Note when using WENO, users should use 6 guardcells 
  !                      along with increased sizes of nxb, nyb, and nzb that are
  !                      larger than 2*NGUARD=12 at least.
  integer, save :: hy_order, hy_transOrder, hy_3Torder
  integer,save  :: hy_gcMaskSize=NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(NUNK_VARS+NDIM*NFACE_VARS),save :: hy_gcMask
  integer, dimension(NFLUXES), save :: hy_fluxCorVars

!#ifdef FLASH_UHD_3T
  ! For Biermann Battery Terms
  real,    save :: hy_avogadro, hy_qele
  logical, save :: hy_useBiermann1T
  logical, save :: hy_useBiermann
  logical, save :: hy_biermannSource
  real,    save :: hy_bier1TZ
  real,    save :: hy_bier1TA
  real,    save :: hy_biermannCoef
  real,    save :: hy_speedOfLight
  logical, save :: hy_hallVelocity
  logical, save :: hy_conserveAngField 
!#endif

  logical, save :: hy_threadBlockList = .false.
  logical, save :: hy_threadWithinBlock = .false.
  integer, save :: hy_3TMode = HY3T_NONE

End Module Hydro_data
