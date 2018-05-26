!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  Hydro_init()
!!  
!!
!! DESCRIPTION
!! 
!!  This routine initializes unit scope variables which are typically the runtime parameters.
!!  The routine must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!
!!***

Subroutine Hydro_init()

  use Hydro_data
  use Driver_interface,            ONLY : Driver_abortFlash, Driver_getMype, &
                                          Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stampMessage, &
                                          Logfile_stampVarMask, &
                                          Logfile_stamp
  use Grid_interface,              ONLY : Grid_setFluxHandling, &
                                          Grid_getBlkIndexLimits
  use IO_interface,                ONLY : io_getScalar

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  character(len=MAX_STRING_LENGTH) :: str_geometry
  integer :: i
  logical :: threadBlockListBuild, threadWithinBlockBuild
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  
  

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,hy_meshMe)
  call Driver_getNumProcs(MESH_COMM,hy_meshNumProcs)

  call RuntimeParameters_get("useHydro",hy_useHydro) !!DEV: not yet systematically honored

  call RuntimeParameters_get("cfl", hy_cfl)
  hy_cfl_original = hy_cfl

  call RuntimeParameters_get("restart", hy_restart)
  if (hy_restart) then
     call IO_getScalar("dt", hy_dt)
     hy_dtmin = hy_dt
  else
     hy_dtmin = HUGE(1.0)
  endif


  call RuntimeParameters_get("UnitSystem",          hy_units)
  call RuntimeParameters_get("order",               hy_order)
  call RuntimeParameters_get("hy_3Torder",          hy_3Torder)
  if (hy_3Torder == -1) hy_3Torder = hy_order
  call RuntimeParameters_get("transOrder",          hy_transOrder)
  if ((NGUARD <= 4) .and. (hy_order > 3)) then
     if (hy_order .ne. 6) call Driver_abortFlash&
          ("[Hydro_init]: Hydro requires more guardcells for the given hy_order method.")
  endif
  hy_useVaryingCFL = .false.
  call RuntimeParameters_get("use_hybridOrder",     hy_useHybridOrder)
  if (hy_useHybridOrder) hy_useVaryingCFL = .true.
  call RuntimeParameters_get("hybridOrderKappa",    hy_hybridOrderKappa)

#ifdef BDRY_VAR
  if (.not. hy_useVaryingCFL) hy_useVaryingCFL = .true.
#endif

#ifdef FLASH_GRID_PARAMESH
  if ((NGUARD > 4) .and. (NXB < 2*NGUARD)) then
     call Driver_abortFlash&
          ("[Hydro_init]: Hydro requires larger NXB, etc. for the given number of guardcells.")
  endif
#endif

  call RuntimeParameters_get("entropy",             hy_entropy)
  call RuntimeParameters_get('entropyFixMethod',    hy_entropyFixMethod_str)
  call RuntimeParameters_get("charLimiting",        hy_charLimiting)
  call RuntimeParameters_get("eintSwitch",          hy_eswitch)
  call RuntimeParameters_get("slopeLimiter",        hy_limiter_str)
  call RuntimeParameters_get("smlrho",              hy_smalldens)
  call RuntimeParameters_get("smallp",              hy_smallpres)
  call RuntimeParameters_get("smallE",              hy_smallE)
  call RuntimeParameters_get('irenorm',             hy_irenorm)
  call RuntimeParameters_get('RiemannSolver',       hy_RiemannSolver_str)
  call RuntimeParameters_get('shockDetect',         hy_shockDetectOn)
  if (hy_shockDetectOn) hy_useVaryingCFL = .true.

  call RuntimeParameters_get("use_steepening",      hy_ContactSteepening)
  call RuntimeParameters_get("use_flattening",      hy_flattening)
  call RuntimeParameters_get("use_upwindTVD",       hy_upwindTVD)
  if (NGUARD <= 4) hy_upwindTVD = .false.
  call RuntimeParameters_get("use_auxEintEqn",      hy_useAuxEintEqn)
  call RuntimeParameters_get("hydroComputeDtOption",hy_hydroComputeDtOption)
  call RuntimeParameters_get("use_avisc",           hy_use_avisc)
  call RuntimeParameters_get("cvisc",               hy_cvisc)
  call RuntimeParameters_get("updateHydroFluxes",   hy_updateHydrofluxes)
  call RuntimeParameters_get("addThermalFlux",      hy_addThermalFlux)
  call RuntimeParameters_get("conserveAngMom",      hy_conserveAngMom) 
  if (NDIM == 3) then
     call RuntimeParameters_get("use_3dFullCTU",    hy_use3dFullCTU)
  else
     hy_use3dFullCTU = .false.
  endif

  !! Extra definitions for MHD
#ifdef FLASH_USM_MHD
  call RuntimeParameters_get("killdivb",            hy_killdivb)
  call RuntimeParameters_get("E_modification",      hy_E_modification)
  call RuntimeParameters_get("E_upwind",            hy_E_upwind)
  call RuntimeParameters_get("energyFix",           hy_energyFixSwitch)
  call RuntimeParameters_get('prolMethod',          hy_prol_method_str)
  call RuntimeParameters_get('ForceHydroLimit',     hy_forceHydroLimit)

  !! Hydro limit B=0 -----------------------------------------------------------
  if (hy_forceHydroLimit) then
     if (hy_killdivb) then
        hy_killdivb = .false.
        if (hy_meshMe .EQ. MASTER_PE) print*, &
             "[Hydro_init] Pure hydro mode is chosen and hy_killdivb is turned off."
        call Logfile_stampMessage &
             ('[Hydro_init] Pure hydro mode: hy_killdivb is turned off.')
     endif
  endif

  !! DivB=0 constraint ---------------------------------------------------------
  if ((NDIM == 1) .OR. (hy_order == 1)) then
     hy_killdivb = .false.
     hy_energyFixSwitch = .false.
     if (hy_meshMe .EQ. MASTER_PE) print*, &
          "[Hydro_init] Parameters killdivb, and energyFix are not "// &
          "required and now they are turned off."
     call Logfile_stampMessage &
          ('[Hydro_init] Parameters killdivb, and energyFix are now turned off')
  endif


  !! Prolongation Method for divergence-free face-centered B fields ------------
  if(trim(hy_prol_method_str) == "injection_prol" .or. &
     trim(hy_prol_method_str) == "INJECTION_PROL" ) then
     hy_prol_method = INJECTION_PROL
  elseif (trim(hy_prol_method_str) == "balsara_prol" .or. &
          trim(hy_prol_method_str) == "BALSARA_PROL" ) then
     hy_prol_method = BALSARA_PROL
  else
     call Driver_abortFlash&
          ("[Hydro_init]: The prolongation method is of unknown type: " // &
           "Please choose one of 'injection_prol(2D & 3D)' or 'balsara_prol(3D)'.")
  endif
  
  ! In cylindrical, injection  prolongation is not correct use modified Balsara instead
  if ((hy_geometry == CYLINDRICAL) .and. (hy_prol_method .ne. BALSARA_PROL)) then
     hy_prol_method = BALSARA_PROL
     call Logfile_stampMessage &
          ('[Hydro_init] Balsara prolongation method is required for cylindrical AMR.')
  endif

#endif /* For USM-MHD */

  !! Gravity -------------------------------------------------------------------
  hy_useGravity = .false.
#ifdef GRAVITY
  call RuntimeParameters_get("useGravity", hy_useGravity)
  if (hy_useGravity) then
     call RuntimeParameters_get("use_gravHalfUpdate", hy_useGravHalfUpdate)
  endif
#endif


  !! Non-ideal diffusions ------------------------------------------------------
  call RuntimeParameters_get("useViscosity",    hy_useViscosity)
  call RuntimeParameters_get("useConductivity", hy_useConductivity)
  hy_useDiffuse = .false.
  if (hy_useViscosity .or. hy_useConductivity) then
     hy_useDiffuse = .true.
  endif
#if defined(FLASH_USM_MHD)
  call RuntimeParameters_get("useMagneticResistivity", hy_useMagneticResistivity)
#ifdef FLASH_USM_MHD
  if (hy_useMagneticResistivity) then
     hy_useDiffuse = .true.
     hy_E_upwind   = .false.
     if (hy_meshMe .EQ. MASTER_PE) print*, &
          "[Hydro_init] The upwind electric field construction is NOT allowed for diffusion."
  endif
#endif
#endif


  !! Entropy Fix ---------------------------------------------------------------
  if (hy_entropy) then
     if(trim(hy_entropyFixMethod_str) == "HartenHyman" .or. &
        trim(hy_limiter_str) == "hartenhyman") then
        hy_entropyFixMethod = HARTENHYMAN
     elseif(trim(hy_entropyFixMethod_str) == "Harten" .or. &
            trim(hy_limiter_str) == "harten") then
        hy_entropyFixMethod = HARTEN
     endif
  endif


  !! Slope Limiter -------------------------------------------------------------
  if(trim(hy_limiter_str) == "minmod" .or. &
     trim(hy_limiter_str) == "MINMOD") then
     hy_limiter = MINMOD
  else if(trim(hy_limiter_str) == "mc" .or. &
          trim(hy_limiter_str) == "MC" ) then
     hy_limiter = MC
  else if (trim(hy_limiter_str) == "hybrid" .or. &
           trim(hy_limiter_str) == "HYBRID" ) then
     hy_limiter = HYBRID
  else if (trim(hy_limiter_str) == "vanLeer" .or. &
           trim(hy_limiter_str) == "VANLEER") then
     hy_limiter = VANLEER
  else if (trim(hy_limiter_str) == "limited" .or. &
           trim(hy_limiter_str) == "LIMITED") then
     hy_limiter = LIMITED
     call RuntimeParameters_get('LimitedSlopeBeta', hy_LimitedSlopeBeta)
  else
     call Driver_abortFlash&
          ("[Hydro_init]: The hy_limter of unknown type! It should be one of" // &
           "'minmod','mc', 'vanLeer', 'hybrid' or 'limited'.")
  end if


  !! Riemann Solver ------------------------------------------------------------
  hy_hybridRiemannOnly = .false.

  if(trim(hy_RiemannSolver_str) == "Roe" .or. &
     trim(hy_RiemannSolver_str) == "roe" .or. &
     trim(hy_RiemannSolver_str) == "ROE" ) then
     hy_RiemannSolver = ROE
  elseif (trim(hy_RiemannSolver_str) == "hll" .or. &
          trim(hy_RiemannSolver_str) == "HLL" ) then
     hy_RiemannSolver = HLL
  elseif (trim(hy_RiemannSolver_str) == "hllc" .or. &
          trim(hy_RiemannSolver_str) == "HLLC" ) then
     hy_RiemannSolver = HLLC
#if defined(FLASH_USM_MHD)
  elseif (trim(hy_RiemannSolver_str) == "hlld" .or. &
          trim(hy_RiemannSolver_str) == "HLLD" ) then
     hy_RiemannSolver = HLLD
#endif
  elseif (trim(hy_RiemannSolver_str) == "marquina" .or. &
          trim(hy_RiemannSolver_str) == "Marquina" ) then
     hy_RiemannSolver = MARQ
  elseif (trim(hy_RiemannSolver_str) == "marquinaModified" .or. &
          trim(hy_RiemannSolver_str) == "marquinamodified" .or. &
          trim(hy_RiemannSolver_str) == "MarquinaModified" ) then
     hy_RiemannSolver = MARM
  elseif (trim(hy_RiemannSolver_str) == "LocalLaxFriedrichs" .or. &
          trim(hy_RiemannSolver_str) == "llf"  .or. &
          trim(hy_RiemannSolver_str) == "LLF" ) then
     hy_RiemannSolver = LLF
  elseif (trim(hy_RiemannSolver_str) == "HYBRID" .or. &
          trim(hy_RiemannSolver_str) == "hybrid"  .or. &
          trim(hy_RiemannSolver_str) == "Hybrid" ) then
     hy_RiemannSolver = HYBR
     if (.not.hy_shockDetectOn) then
        hy_hybridRiemannOnly = .true.
        hy_shockDetectOn = .true.
     endif
  else
     call Driver_abortFlash&
          ("[Hydro_init]: The Riemann Solver is of unknown type: " // &
           "Options are LLF, HLL, HLLC, HLLD (only for MHD), Roe, "//& 
           "hybrid, Marquina, or MarquinaModified.")
  endif

  !! Geometry ------------------------------------------------------------------
  call RuntimeParameters_get("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, hy_geometry)
  if ((hy_geometry .NE. CARTESIAN) .AND. (hy_geometry .NE. CYLINDRICAL) )  then
     call Driver_abortFlash&
          ("[Hydro_init]: Spherical coordinates are not supported in the unsplit solvers")
  endif
!!$  if (hy_geometry .NE. CARTESIAN )  then
!!$     if (hy_meshME == MASTER_PE) then
!!$        print *, "[Hydro_init]: Using non-Cartesian Geometry!"
!!$!$        if (.not. hy_useAuxEintEqn) then
!!$!$           print *, "[Hydro_init]: Using non-Cartesian Geometry requires to turn on hy_useAuxEintEqn!"
!!$!$        endif
!!$     endif
!!$!$     hy_useAuxEintEqn = .true.
!!$  endif
  call RuntimeParameters_get("flux_correct", hy_fluxCorrect)
  if (NDIM > 1) then
     if (hy_fluxCorrect) then
        if (hy_geometry == CARTESIAN) then
           call Grid_setFluxHandling('consv_flux_densities')
        else
           call Grid_setFluxHandling('consv_fluxes')
        endif
     end if
  end if

  call PhysicalConstants_get("ideal gas constant", hy_Rconst)

  !! For correct flux correction in non-Cartesian geometry----------------------
  do i = 1, NFLUXES
     hy_fluxCorVars(i) = i
  enddo
     

  !! System units --------------------------------------------------------------
  hy_xref = 1.0
  hy_vref = 1.0
  hy_dref = 1.0

  hy_mref = hy_xref*hy_vref
  hy_tref = hy_xref/hy_vref
  hy_eref = hy_vref*hy_vref
  hy_nref = hy_dref*hy_vref*hy_xref
  hy_pref = hy_dref*hy_vref*hy_vref
  hy_gref = hy_vref*hy_vref/hy_xref

  hy_qref = hy_vref*hy_vref/hy_Rconst
  hy_kref = hy_dref*hy_vref*hy_xref*hy_Rconst


#if defined(FLASH_USM_MHD)
  if ( hy_units == "SI" .or. hy_units == "si" ) then
     hy_bref = hy_vref*sqrt(4.0*PI*hy_dref*1.e-7)
  else if ( hy_units == "CGS" .or. hy_units == "cgs" ) then
     hy_bref = hy_vref*sqrt(4.0*PI*hy_dref)
  else
     hy_bref = hy_vref*sqrt(hy_dref)
  end if
#endif

  !! Allow selective guardcell fill calls ---------------------------------------
  hy_gcMaskSize = NUNK_VARS
  hy_gcMask = .TRUE.

  ! Begin extra definitions for MHD
#ifdef FLASH_USM_MHD /* Extra setups for MHD */
  hy_gcMaskSize = NUNK_VARS+NDIM*NFACE_VARS

#ifdef DIVB_VAR
  hy_gcMask(DIVB_VAR) = .FALSE.
#endif
#ifdef MAGP_VAR
  hy_gcMask(MAGP_VAR) = .FALSE.
#endif
#ifdef TOTP_VAR
  hy_gcMask(TOTP_VAR) = .FALSE.
#endif
#ifdef BETA_VAR
  hy_gcMask(BETA_VAR) = .FALSE.
#endif
#ifdef VECZ_VAR
  hy_gcMask(VECZ_VAR) = .FALSE.
#endif
#ifdef CURX_VAR
  hy_gcMask(CURX_VAR) = .FALSE.
#endif
#ifdef CURY_VAR
  hy_gcMask(CURY_VAR) = .FALSE.
#endif
#ifdef CURZ_VAR
  hy_gcMask(CURZ_VAR) = .FALSE.
#endif
#endif /* for MHD */

#if NSPECIES == 1
#ifdef SPECIES_BEGIN
  hy_gcMask(SPECIES_BEGIN) = .FALSE.
#endif
#endif

!!$!! Species and mass scalars need to use full scratch arrays
!!$#if (NSPECIES+NMASS_SCALARS) > 0
!!$#ifndef FLASH_UHD_NEED_SCRATCHVARS
!!$  call Driver_abortFlash&
!!$       ("[Hydro_init]: Species and mass scalar update for UHD requires '+uhdWithSpeciesMassScalar' in setup"//&
!!$       " instead of '+uhd'.")
!!$#endif
!!$#endif

#ifdef FLASH_UHD_3T
  call hy_uhd_MultiTempInit()
#endif

  call Logfile_stampVarMask(hy_gcMask, .FALSE., '[Hydro_init]', 'gcNeed')


!#ifdef FLASH_UHD_3T
  call PhysicalConstants_get('Avogadro',        hy_avogadro)
  call PhysicalConstants_get('electron charge', hy_qele)
  call PhysicalConstants_get('speed of light',  hy_speedOfLight)
  hy_qele = hy_qele * sqrt(4*PI)/hy_speedOfLight !this is the correct charge

#if defined(FLASH_USM_MHD)
  call RuntimeParameters_get('use_Biermann',      hy_useBiermann)
  call RuntimeParameters_get('use_Biermann1T',    hy_useBiermann1T)
  call RuntimeParameters_get('hy_biermannSource', hy_biermannSource)
  call RuntimeParameters_get('hy_bier1TZ',        hy_bier1TZ)
  call RuntimeParameters_get('hy_bier1TA',        hy_bier1TA)
  call RuntimeParameters_get('hy_biermannCoef',   hy_biermannCoef)
#endif
!#endif


  call RuntimeParameters_get("threadBlockListBuild", threadBlockListBuild)
  call RuntimeParameters_get("threadHydroBlockList", hy_threadBlockList)

  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadHydroWithinBlock", hy_threadWithinBlock)

  if (hy_threadBlockList .and. .not. threadBlockListBuild) then
     call Logfile_stamp('WARNING! Turning off block list threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadBlockList = .false.
  end if
  if (hy_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Hydro_init]')
     hy_threadWithinBlock = .false.
  end if


  hy_needScrchVars = .false.
#ifdef FLASH_UHD_NEED_SCRATCHVARS
  hy_needScrchVars = .true.
#endif


  !! GP initialization
!!  if (hy_order == 6) then

!!     if (NDIM < 2) then
!!        call Driver_abortFlash("[Hydro_init]: GP is designed to support 2D and 3D.")
!!     endif


    ! get runtime parameters
!!    call RuntimeParameters_get("radiusGP", hy_radiusGP)
!!    call RuntimeParameters_get("sigmaGP",  hy_sigmaGP)
!!    call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)
    
    ! get counter size
!!    call hy_uhd_counterGP_init(hy_radiusGP, hy_counterGP, blkLimitsGC)

    ! allocate the arrays and matrices
!!    allocate(hy_RinvGP(hy_counterGP,hy_counterGP))
!!    allocate(hy_WpGP(NDIM,hy_counterGP))
!!    allocate(hy_WmGP(NDIM,hy_counterGP))
    
    ! get weights and Rinv     
!!    call hy_uhd_initGP(hy_RinvGP, hy_WpGP, hy_WmGP, blkLimitsGC)
!! endif

End Subroutine Hydro_init
