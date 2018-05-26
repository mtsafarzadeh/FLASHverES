!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/threadBlockList/hy_uhd_unsplit
!!
!! NAME
!!
!!  hy_uhd_unsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_unsplit( integer (IN) :: blockCount,
!!                  integer (IN) :: blockList(blockCount),
!!                  real    (IN) :: dt,
!!                  real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!
!!  Performs MHD update in a directionally unsplit fashion over a set
!!  of blocks. Divergence cleaning control for the magnetic fields is
!!  handled by using the staggered mesh algorithm with a call to
!!  hy_uhd_staggeredDivb. Before calling this routine, electric fields are
!!  to be calculated with a call to hy_uhd_getElectricFields.
!!  blockList is an integer array of size blockCount that contains 
!!  the blocks over which to update.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - applies an eos to the guard cells; 
!!   - computes fluxes using a call to hy_uhd_getFaceFlux
!!   - if we're not doing flux correction (as controlled by the flux_correct
!!     runtime parameter), then we update all the cell values from the fluxes 
!!     (with a call to hy_uhd_unsplitUpdate), otherwise, we update just cells 
!!     not on the boundaries, and save fluxes for cells on the boundary;
!!   - and finally, we apply an eos to the block.
!! 
!!  After the main block loop, if doing flux correction, we have
!!  the Grid correct boundary fluxes for all blocks where approriate,
!!  and do another loop over blocks, updating the cell values for
!!  cells on the block boundaries using the corrected fluxes, and
!!  apply an eos on the block. 
!!  The same is true for the electric fields correction where
!!  there are block boundaries that are sharing different levels of
!!  refinements.
!!
!!
!! REFERENCES
!!
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (needed for temporal extrapolations of gravity)
!!
!!***

!!REORDER(4): U, scrch_Ctr, fl[xyz]

#ifdef DEBUG_ALL
#define DEBUG_MHD
#endif
#define DEBUG_GRID_GCMASK

Subroutine hy_uhd_unsplit ( blockCount, blockList, dt, dtOld )

  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_killdivb,         &
                         hy_gref,             &
                         hy_useGravity,       &
                         hy_order,            &
                         hy_units,            &
                         hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_forceHydroLimit,  &
                         hy_facevar2ndOrder,  &
                         hy_eosModeAfter,     &
                         hy_useGravHalfUpdate,&
                         hy_useGravPotUpdate, &
                         hy_gravConsv,        &
                         hy_updateHydroFluxes,&
                         hy_geometry,         &
                         hy_fluxCorVars,      &
                         hy_cfl,              &
                         hy_cfl_original,     &
                         hy_useBiermann,      &
                         hy_useBiermann1T,    &
                         hy_biermannSource,   &
                         hy_threadBlockList,  &
                         hy_shockDetectOn

  use Driver_interface, ONLY : Driver_abortFlash

  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_getElectricFields,&
                               hy_uhd_staggeredDivb,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_eigenParameters,  &
                               hy_uhd_eigenValue,       &
                               hy_uhd_putGravityUnsplit,&
                               hy_uhd_addGravityUnsplit,&
                               hy_uhd_biermannSource,   &
                               hy_uhd_shockDetect

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_getBlkData,        &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  use Eos_interface, ONLY : Eos_wrapped

  use Logfile_interface, ONLY : Logfile_stampVarMask, Logfile_stamp

!!$  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Gravity_interface, ONLY : Gravity_potentialListOfBlocks

  use hy_uhd_slopeLimiters, ONLY : minmod, get_upwind

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "UHD.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  real,    INTENT(IN) :: dt, dtOld
  !! -----------------------------------------------------

  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: i,blockID,level=0
  integer :: ix,iy,iz
  real, dimension(MDIM) :: del
  logical :: gcMask(hy_gcMaskSize)
  logical :: halfTimeAdvance
  integer, dimension(2,MDIM) :: eosRange

#ifdef FIXEDBLOCKSIZE
  real :: flx(NFLUXES,&
              GRID_ILO_GC:GRID_IHI_GC,     &
              GRID_JLO_GC:GRID_JHI_GC,     &
              GRID_KLO_GC:GRID_KHI_GC)
  real :: fly(NFLUXES,&
              GRID_ILO_GC:GRID_IHI_GC,     &
              GRID_JLO_GC:GRID_JHI_GC,     &
              GRID_KLO_GC:GRID_KHI_GC)
  real :: flz(NFLUXES,&
              GRID_ILO_GC:GRID_IHI_GC,     &
              GRID_JLO_GC:GRID_JHI_GC,     &
              GRID_KLO_GC:GRID_KHI_GC)


  real, dimension(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: gravX,gravY,gravZ
  real :: faceAreas(GRID_ILO_GC:GRID_IHI_GC,     &
                    GRID_JLO_GC:GRID_JHI_GC,     &
                    GRID_KLO_GC:GRID_KHI_GC)

#else
  real, allocatable, dimension(:,:,:,:) :: flx,fly,flz
  real, allocatable, dimension(:,:,:,:) :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)
#endif

!#ifdef HY_NEED_EXTRA_GCFILL
  real, pointer, dimension(:,:,:,:) :: U, scrch_Ctr
!#endif

  real :: csx,cax,cfx,csy,cay,cfy,as,af,bbN,uN 
  real,dimension(MDIM)         :: beta 
  real,dimension(HY_WAVENUM)   :: lambdaLoc 
  real,dimension(HY_VARINUM2)  :: Vloc
  integer :: updateMode=UPDATE_ALL
  integer :: iSize, jSize, kSize
  real    :: gravDtFactor
  logical :: halfTimeGravUpdate

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif
  
  !! End of data declaration ***********************************************


#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The staggeredMesh MHD solver only works with PARAMESH 3 or 4!")
#endif

#ifdef FLASH_GRID_UG
  hy_fluxCorrect = .false.
#endif



  !! ***************************************************************************
  !! Shock detection before MHD                                                *
  !! ***************************************************************************
  !! Call shock detect algorithm to determine tagging shocks before hydro begins:
  if (hy_shockDetectOn) then

     !! Call guardcell filling to properly detect shocks
     gcMask = .false.
     gcMask(DENS_VAR) = .true.   
     gcMask(PRES_VAR) = .true.   
     gcMask(GAMC_VAR) = .true.   
     gcMask(VELX_VAR:VELZ_VAR) = .true.

#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant[Detect]')
     end if
#endif

     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,&
          maskSize=NUNK_VARS, mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)

     !! Detect shocks
     do i=1,blockCount
        blockID = blockList(i)
        call hy_uhd_shockDetect(blockID)
     enddo
  endif


  !! ***************************************************************************
  !! Call guardcell filling with Eos before MHD                                *
  !! ***************************************************************************

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcNeed')
  end if
#endif

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=hy_gcMaskSize,mask=hy_gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)


  !! ***************************************************************************
  !! First part of advancement                                                 *
  !! ***************************************************************************
  !! Loop over the blocks

  !! Retain the original cfl that may have been changed in some leaf blocks.
  if (hy_updateHydroFluxes) then
     hy_cfl = hy_cfl_original
  endif

  
  !$omp parallel if (hy_threadBlockList) &
  !$omp default(none) &
  !$omp firstprivate(level,updateMode) &
  !$omp private(i,blockID,del,blkLimits,blkLimitsGC,eosRange,datasize,&
  !$omp flx,fly,flz,gravX,gravY,gravZ,U,scrch_Ctr,faceAreas, halftimeadvance,& 
  !$omp csx, cax, cfx, csy, cay, cfy, as, af, bbN, uN, Vloc, lambdaLoc, beta) &
  !$omp shared(blockCount,blockList,gcMask,gcMaskLogged,dt,dtOld,hy_units,&
  !$omp hy_unsplitEosMode,hy_useGravity,hy_gref,hy_fluxCorrect,hy_order,&
  !$omp hy_updateHydroFluxes,hy_eosModeAfter,hy_useGravHalfUpdate,&
  !$omp hy_useGravPotUpdate,hy_geometry,hy_forcehydrolimit,&
  !$omp hy_facevar2ndorder, hy_killdivb, hy_gcmasksize,&
  !$omp hy_usebiermann, hy_usebiermann1t, hy_biermannsource, hy_fluxCorVars)

    
  
  !$omp do schedule(static)
  do i=1,blockCount

     blockID = blockList(i)

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        call hy_uhd_unitConvert(blockID,FWDCONVERT)
     endif

     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1     

#ifndef FIXEDBLOCKSIZE
     allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravX(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! ************************************************************************
     !! Get gravity
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity .and. .not. hy_useGravPotUpdate) then
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref
     endif


     if (hy_updateHydroFluxes) then
        !! ************************************************************************
        !! Calculate Riemann (interface) states
        !! Note: gravX(1,:,:,:) - gravity at n
        !!       gravX(2,:,:,:) - gravity at n+1/2 (extrapolated)
        !!       gravX(3,:,:,:) - gravity at n+1 (exprapolated)

        call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                    gravX(1,:,:,:),gravY(1,:,:,:),gravZ(1,:,:,:),&
                                    gravX(2,:,:,:),gravY(2,:,:,:),gravZ(2,:,:,:),&
                                    normalFieldUpdate=.false.)

        ! Note: Two different ways of handling gravity:
        ! 1. With gravity calculated potential at n+1/2, Riemann states do not include gravity
        !    source terms at this point, and will include them in hy_uhd_addGravityUnsplit later
        !    to the primitive Riemann variables (not available for conservative formulation).
        ! 2. With gravity extrapolated from n-1 & n states, gravity source terms have been
        !    included to Riemann states in conservative formulation in hy_uhd_getRiemannState.

     endif !! End of if (hy_updateHydroFluxes) then

     !! Don't do this part when pure hydro (NFACE_VARS = 0) is considered
#if NFACE_VARS > 0
#if NDIM > 1
     !! Note: Do this intermediate part for multidimensional problem only.
     !!       One dimensional problem doesn't need to do this part as 1D setup
     !!       doesn't need to use magnetic fields (facevars) at cell face centers.

     !! ************************************************************************
     !! Evolve facevars n+1/2 time step to get 2nd order accuracy
     if ((.not. hy_forceHydroLimit) .and. hy_facevar2ndOrder)then
        !! Intermediate Godunov fluxes
        call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del,flx,fly,flz)
        if (hy_killdivb) then
           !! DivB=0 calls
           halfTimeAdvance = .true.
           !! Calculate Godunov fluxes one more time to calculate electric fields
           call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
        endif
     endif
     !! End of evolution of face-centered magnetic fields by n+1/2 time step
     !! ************************************************************************

#ifndef FIXEDBLOCKSIZE
     deallocate(flx)
     deallocate(fly)
     deallocate(flz)
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
#endif

  enddo
  !$omp end do                                                                                          


  !! End of the first loop of evolution by n+1/2 time step
  !! ***************************************************************************



#ifdef GPOT_VAR
  if (hy_useGravity .and. hy_useGravPotUpdate) then
     !! ***************************************************************************
     !! Gravity calculation at n+1/2 by calling Poisson solver                    *
     !! ***************************************************************************
     ! At this stage, U(DENS_VAR,:,:,:) is averaged extrapolated (or bilinear) density 
     ! at n+1/2 values to cell center (already done in getRiemann).
     ! This density is used to compute gravity components at n+1/2.
     ! Note: In case hy_useGravPotUpdate = .false., gravity components will be extrapolated to n+1/2.

     !$omp single     
     call Gravity_potentialListOfBlocks(blockCount, blockList)
     
     !! Fill guardcells for only the gpot variable
     gcMask  =.false.
     gcMask(GPOT_VAR) = .true.
#if defined(GPOL_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
     if (hy_useGravPotUpdate) gcMask(GPOL_VAR) = .true.
#endif
     gcMask(ENER_VAR) = .true.
#ifdef EINT_VAR
     gcMask(EINT_VAR) = .true.
#endif

#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcWant[Pot]')
     end if
#endif

     !! Guardcell filling for gpot
     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&
          maskSize=NUNK_VARS,mask=gcMask,makeMaskConsistent=.true.,&
          doLogMask=.NOT.gcMaskLogged)
     !$omp end single

  endif
#endif




  !! ***************************************************************************
  !! Exchange guardcell data for face-centered magnetic fields that will be
  !! used in the data reconstruction-evolution calculation.
  if ( hy_facevar2ndOrder .and. (.not. hy_forceHydroLimit)) then

     if (hy_fluxCorrect .and. hy_killdivb) then
        !! Call the electric field correction routine
        !$omp single     
        call Grid_conserveField()
        !$omp end single
     endif
  
     !! Calculate divergence-free magnetic fields at interfaces
     !$omp do schedule(static)
     do i=1,blockCount
        blockID = blockList(i)
        call Grid_getDeltas(blockID,del)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        halfTimeAdvance = .true.
        !! Evolve face-centered magnetic fields n+1/2 time step
        call hy_uhd_staggeredDivb(blockID,0.5*dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
     enddo
     !$omp end do

     !! Fill guardcells for only the face-centered field variable that are
     !! evolved to n+1/2 time step for date reconstruction

     !$omp single     
     gcMask  =.false.
#ifdef MAGI_FACE_VAR
     gcMask(MAGI_FACE_VAR:NFACE_VARS*NDIM:NFACE_VARS) = .true.
#endif
#ifdef MAG_FACE_VAR
     gcMask(MAG_FACE_VAR:NFACE_VARS*NDIM:NFACE_VARS) = .true.
#endif

#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stamp('gcNeed(MAGI_FACE_VAR,MAG_FACE_VAR) - FACES', '[hy_uhd_unsplit]')
     end if
#endif
     !! Guardcell filling for facevars
     call Grid_fillGuardCells(FACES,ALLDIR,&
          maskSize=(hy_gcMaskSize-NUNK_VARS),mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)
     !$omp end single
  endif


  !! ***************************************************************************
  !! Second part of advancement related to gravity using potential             *
  !! ***************************************************************************
  !! Loop over the blocks
  !$omp do schedule(static)
  do i=1,blockCount
     blockID = blockList(i)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

#ifndef FIXEDBLOCKSIZE
     allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravX(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! ************************************************************************
     !! Get gravity
     gravX = 0.
     gravY = 0.
     gravZ = 0.


     !! This part of the gravity update is only getting called when hy_useGravPotUpdate
     !! is chosen instead of extrapolating gravity components.
     if (hy_useGravity .and. hy_useGravHalfUpdate .and. hy_useGravPotUpdate) then
        call Grid_getBlkPtr(blockID,U,CENTER)
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
        
        ! Reset density in unk to n state values as we are done with
        ! computing gravities at n+1/2 using density at n+1/2.
        U(DENS_VAR,:,:,:) = scrch_Ctr(VAR2_SCRATCH_CENTER_VAR,:,:,:)

        call Grid_releaseBlkPtr(blockID,U,CENTER)
        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)

        
        ! Now update cell-center grav accelerations (using new potential)
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref

        halfTimeGravUpdate = .true.
        ! and add gravity source terms to Riemann states at n+1/2.
        call hy_uhd_addGravityUnsplit(blockID,blkLimits,dataSize,dt,&
             gravX(2,:,:,:),gravY(2,:,:,:),gravZ(2,:,:,:),halfTimeGravUpdate)

     endif

     ! We need to call gravity again here for cases (AMR or UG) WITH scratch vars when
     ! hy_useGravHalfUpdate=.false. and hy_useGravPotUpdate=.false.
     ! This gravity is for the final update in the corrector step.
     if (hy_useGravity .and. .not. hy_useGravPotUpdate) then
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref
     endif


     if (hy_updateHydroFluxes) then
        if ((hy_order >= 2) .and. hy_facevar2ndOrder) then
           !! *********************************************************************
           !! Calculate Riemann (interface) states with divergence-free fields at n+1/2
           !! Note: gravX(1,:,:,:) - gravity at n
           !!       gravX(2,:,:,:) - gravity at n+1/2 (extrapolated)
           !!       gravX(3,:,:,:) - gravity at n+1 (exprapolated)
           call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                       gravX(1,:,:,:),gravY(1,:,:,:),gravZ(1,:,:,:),&
                                       gravX(2,:,:,:),gravY(2,:,:,:),gravZ(2,:,:,:),&
                                       normalFieldUpdate=.true.)
        endif
     endif !end of if (hy_updateHydroFluxes) then
#endif
! end if of #if NDIM > 1
#endif 
! end if of #if NFACE_VARS > 0

     if (NDIM == 1) then
        if (hy_useGravity) then
           if (.not. hy_useGravPotUpdate) then
              ! We need to call gravity again here for cases (AMR or UG) WITH scratch vars when
              ! hy_useGravHalfUpdate=.false. and hy_useGravPotUpdate=.false.
              ! This gravity is for the final update in the corrector step.
              call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
              gravX = gravX/hy_gref
              gravY = gravY/hy_gref
              gravZ = gravZ/hy_gref
           else
              call Driver_abortFlash&
                   ("[hy_uhd_unsplit]: 1D gravity potential update is not YET available!")
              !! DEV-Dongwok: we need to put gravity for potential update in 1D
           endif
        endif
     endif
        
     

     !! ************************************************************************
     !! Calculate high order Godunov fluxes
     !! Initialize arrays with zero
     flx = 0.
     fly = 0.
     flz = 0.

     call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del,flx,fly,flz,&
                             lastCall=.true.)

     if (.not. hy_fluxCorrect) then
        !! ************************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,dtOld,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ)

#if NFACE_VARS > 0
        !! ************************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if (.not. hy_forceHydroLimit) then
           if (NDIM >= 2) then
              !! DivB=0 calls
              if (hy_killdivb) then
                 halfTimeAdvance = .false.
                 call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
                 call hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
              endif
           endif
        endif ! End of if (.not. hy_forceHydroLimit) then
#endif

#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,del,hy_unsplitEosMode)

#ifdef FLASH_UHD_3T
        call hy_uhd_unsplitUpdateMultiTemp&
             (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz)
#endif

        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert unit
           call hy_uhd_unitConvert(blockID,BWDCONVERT)
        endif


        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
#endif /* ifndef GRAVITY */

     else !! if Flux correction is used.
        !! Flux conservation calls on AMR:
        !! Correct fluxes at each block boundary where coarse and fine
        !! blocks are neighboring each other.

#if NFACE_VARS > 0
        !! ************************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if (.not. hy_forceHydroLimit) then
           if (NDIM >= 2) then
              !! DivB=0 calls
              if (hy_killdivb) then
                 halfTimeAdvance = .false.
                 call hy_uhd_getElectricFields(blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
              endif
           endif
        endif ! End of if (.not. hy_forceHydroLimit) then
#endif


        if (hy_geometry /= CARTESIAN) then
           ! we are using consv_fluxes and need to divide by face areas
#ifndef FIXEDBLOCKSIZE
           allocate(faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif
           faceAreas = 0.
           call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_putFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)

           if (NDIM >= 2) then
              faceAreas = 0.
              call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
                                   (/1,1,1/), faceAreas, datasize)
              call Grid_putFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)

              if (NDIM == 3) then
                 faceAreas = 0.
                 call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                      (/1,1,1/), faceAreas, datasize)
                 call Grid_putFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
              endif
           endif
#ifndef FIXEDBLOCKSIZE
           deallocate(faceAreas)
#endif
        else ! Cartesian geometry
           call Grid_putFluxData(blockID,IAXIS,flx,datasize)
           if (NDIM >= 2) then
              call Grid_putFluxData(blockID,JAXIS,fly,datasize)
              if (NDIM == 3) then
                 call Grid_putFluxData(blockID,KAXIS,flz,datasize)
              endif
           endif
        endif

     endif

     
#ifndef FIXEDBLOCKSIZE
     deallocate(flx)
     deallocate(fly)
     deallocate(flz)
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
#endif
  end do  
  !$omp end do


  !! End of leaf block do-loop before flux conserve call

  
  !! ***************************************************************************
  !! Third part of advancement                                                 *
  !! ***************************************************************************
  !! Do this part only if refining and flux correcting
  if (hy_fluxCorrect) then
     !$omp single
  
     !! ************************************************************************
     !! Conservation of Fluxes at each block boundary
     call Grid_conserveFluxes(ALLDIR,level)
     
#if NFACE_VARS > 0
     !! ************************************************************************
     !! Correct electric fields at fine-coarse boundaries
     if ((NDIM >= 2) .and. hy_killdivb .and. (.not. hy_forceHydroLimit)) then
        !! Call the electric field correction routine

        call Grid_conserveField()

     endif
#endif
     !$omp end single
     
     !! ************************************************************************
     !! Perform updates for every blocks
     !$omp do schedule(static)
     do i = 1,blockCount

        blockID = blockList(i)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getDeltas(blockID,del)

        datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
#ifndef FIXEDBLOCKSIZE
        allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(gravX(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(gravY(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(gravZ(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

        !! Initialize arrays with zero
        flx = 0.
        fly = 0.
        flz = 0.

        !! *********************************************************************
        !! Get gravity
        gravX = 0.
        gravY = 0.
        gravZ = 0.

        if (hy_useGravity) then
           call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
           gravX = gravX/hy_gref
           gravY = gravY/hy_gref
           gravZ = gravZ/hy_gref
        endif

        !! Get modified conserved flux values at block interfaces:
        !! This is important especially at the block interface at 
        !! different refinement levels of neighboring blocks
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)

        flx(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
             =scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)

        if (hy_geometry /= CARTESIAN) then
           faceAreas = 0.
           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,IAXIS,flx,datasize)
        endif

#if NDIM > 1
        fly(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
             =scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)

        if (hy_geometry /= CARTESIAN) then
           faceAreas = 0.
           call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,JAXIS,fly,datasize)
        endif
        
#if NDIM > 2
        flz(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
             =scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)

        if (hy_geometry /= CARTESIAN) then
           faceAreas = 0.
           call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,KAXIS,flz,datasize)
        endif
#endif
#endif

        !! *********************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,dtOld,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ)

#if NFACE_VARS > 0
        !! *********************************************************************
        !! Update face-centered magnetic fields from n to n+1 time step
        if (NDIM >=2) then
           if (hy_killdivb) then
              halfTimeAdvance = .false.
              call hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
           endif
        endif
#endif

#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! *********************************************************************
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,del,hy_unsplitEosMode)

#ifdef FLASH_UHD_3T
        call hy_uhd_unsplitUpdateMultiTemp&
             (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz)
#endif

        !! *********************************************************************
        !! Convert unit if necessary
        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert unit
           call hy_uhd_unitConvert(blockID,BWDCONVERT)
        endif

        
        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)

        
        !! *********************************************************************
        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
#endif /* if gravity is included, we calculate gravity for n+1 in the below */


#ifndef FIXEDBLOCKSIZE
        deallocate(flx)
        deallocate(fly)
        deallocate(flz)
        deallocate(gravX)
        deallocate(gravY)
        deallocate(gravZ)
        deallocate(faceAreas)
#endif
        
     end do !! End of the loop over blocks
     !$omp end do
  end if !! End of the flux conservation routine
  !$omp end parallel  

#ifdef GRAVITY /* Perform this only when gravity is used */
  !! ***************************************************************************
  !! Fourth part of advancement to compute gravity at n+1 state                *
  !! ***************************************************************************
#ifdef GPOT_VAR
  if (hy_useGravity) then

     !! Gravity calculation at n+1 by calling Poisson solver
     call Gravity_potentialListOfBlocks(blockCount, blockList)

     !! Fill guardcells for only the gpot variable
     gcMask  =.false.
     gcMask(GPOT_VAR) = .true.

#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant[Pot]')
     end if
#endif

     !! Guardcell filling for gpot
     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.false.,&
          maskSize=NUNK_VARS,mask=gcMask,makeMaskConsistent=.false.,&
          doLogMask=.NOT.gcMaskLogged)
  endif
#endif


  !! Proceed to couple the updated gravitational accelerations 
  !! to energy and momenta on each block (see hy_uhd_addGravityUnsplit)
  do i = 1,blockCount

     blockID = blockList(i)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getDeltas(blockID,del)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
#ifndef FIXEDBLOCKSIZE
     allocate(gravX(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! *********************************************************************
     !! Get and add gravity to hydro variables (momenta and energy)         *
     !! *********************************************************************
     if (hy_useGravity) then
        if (hy_order > 1 .and. hy_useGravHalfUpdate) then
           gravDtFactor = 1.
        else
           gravDtFactor = 2.
        endif

        gravX = 0.
        gravY = 0.
        gravZ = 0.
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref

        halfTimeGravUpdate=.false.
        call hy_uhd_addGravityUnsplit(blockID,blkLimits,dataSize,gravDtFactor*dt,&
             gravX(1,:,:,:),gravY(1,:,:,:),gravZ(1,:,:,:),halfTimeGravUpdate)
     endif

     !! *********************************************************************
     !! Correct energy if necessary
     call hy_uhd_energyFix(blockID,blkLimits,dt,del,hy_unsplitEosMode)

     !! *********************************************************************
     !! Convert unit if necessary
     if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
        !! Convert unit
        call hy_uhd_unitConvert(blockID,BWDCONVERT)
     endif



     !! *********************************************************************
     !! Call to Eos
     call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)

#ifndef FIXEDBLOCKSIZE
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
#endif

  end do !! End of the loop over blocks
#endif /* End of n+1 gravity coupling */






!!$#ifdef FAST_VAR 
!!$  do i=1,blockCount 
!!$      
!!$     blockID = blockList(i) 
!!$     call Grid_getBlkPtr(blockID,U,CENTER) 
!!$     call Grid_getDeltas(blockID,del) 
!!$     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
!!$      
!!$     do iz=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS) 
!!$        do iy=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) 
!!$           do ix=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS) 
!!$               
!!$              Vloc(HY_DENS)=U(DENS_VAR,ix,iy,iz) 
!!$              Vloc(HY_VELX:HY_VELZ)=U(VELX_VAR:VELZ_VAR,ix,iy,iz) 
!!$              Vloc(HY_MAGX:HY_MAGZ)=U(MAGX_VAR:MAGZ_VAR,ix,iy,iz) 
!!$              Vloc(HY_PRES)=U(PRES_VAR,ix,iy,iz) 
!!$              Vloc(HY_GAMC:HY_GAME)=U(GAMC_VAR:GAME_VAR,ix,iy,iz) 
!!$               
!!$              call hy_uhd_eigenParameters(Vloc(HY_DENS:HY_GAME),DIR_X,.true.,& 
!!$                   uN,cfx,C_alfn=cax,C_slow=csx,A_f=af,A_s=as,B_beta=beta) 
!!$              call hy_uhd_eigenValue(lambdaLoc,uN,cfx,C_alfn=cax,C_slow=csx) 
!!$               
!!$              call hy_uhd_eigenParameters(Vloc(HY_DENS:HY_GAME),DIR_Y,.true.,& 
!!$                   uN,cfy,C_alfn=cay,C_slow=csy,A_f=af,A_s=as,B_beta=beta) 
!!$              call hy_uhd_eigenValue(lambdaLoc,uN,cfy,C_alfn=cay,C_slow=csy) 
!!$               
!!$              U(FAST_VAR,ix,iy,iz) = max(cfx,cfy) 
!!$              U(SLOW_VAR,ix,iy,iz) = max(csx,csy) 
!!$
!!$           enddo 
!!$        enddo 
!!$     enddo 
!!$  enddo 
!!$#endif 


#if NFACE_VARS > 0
  !! ***********************************************************************
  !! Note: Need guardcell fills to compute current densities and divv.
  !!       Use this call only if needed as this is expensive after all!
#ifdef HY_NEED_EXTRA_GCFILL


#ifdef FLASH_GRID_PARAMESH
  gcMask  =.FALSE.
  gcMask(MAGX_VAR)=.true.
  gcMask(MAGY_VAR)=.true.
  gcMask(MAGZ_VAR)=.true.
  gcMask(VELX_VAR)=.true.
  gcMask(VELY_VAR)=.true.
  gcMask(VELZ_VAR)=.true.
#ifdef MAG_FACE_VAR
  gcMask(NUNK_VARS+MAG_FACE_VAR:NUNK_VARS+NFACE_VARS*NDIM:NFACE_VARS) = .true.
#endif

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(gcMask, .FALSE., '[hy_uhd_unsplit]', 'gcWant')
  end if
#endif

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,minLayers=1,doEos=.false.,&
       maskSize=hy_gcMaskSize,mask=gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)
#endif

#ifdef FLASH_GRID_UG
  !! Note that, on UG, guardcell masking is not applicable
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR)
#endif
  do i=1,blockCount

     blockID = blockList(i)
     call Grid_getBlkPtr(blockID,U,CENTER)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do iz=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do iy=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do ix=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              !! current density jx
#ifdef CURX_VAR
#if NDIM > 1
              U(CURX_VAR,ix,iy,iz)=&
                   0.5*(U(MAGZ_VAR,ix,iy+1,iz)-U(MAGZ_VAR,ix,iy-1,iz))/del(DIR_Y)
#if NDIM > 2
              U(CURX_VAR,ix,iy,iz)=U(CURX_VAR,ix,iy,iz)&
                   -0.5*(U(MAGY_VAR,ix,iy,iz+1)-U(MAGY_VAR,ix,iy,iz-1))/del(DIR_Z)
#endif
#endif
#endif

              !! current density jy
#ifdef CURY_VAR
              U(CURY_VAR,ix,iy,iz)=&
                   -0.5*(U(MAGZ_VAR,ix+1,iy,iz)-U(MAGZ_VAR,ix-1,iy,iz))/del(DIR_X)
#if NDIM > 2
              U(CURY_VAR,ix,iy,iz)=U(CURY_VAR,ix,iy,iz)&
                   +0.5*(U(MAGX_VAR,ix,iy,iz+1)-U(MAGX_VAR,ix,iy,iz-1))/del(DIR_Z)
#endif
#endif

              !! current density jz
#ifdef CURZ_VAR
              U(CURZ_VAR,ix,iy,iz)=&
                   0.5*(U(MAGY_VAR,ix+1,iy,iz)-U(MAGY_VAR,ix-1,iy,iz))/del(DIR_X)

#if NDIM > 1
              U(CURZ_VAR,ix,iy,iz)=U(CURZ_VAR,ix,iy,iz)&
                   -0.5*(U(MAGX_VAR,ix,iy+1,iz)-U(MAGX_VAR,ix,iy-1,iz))/del(DIR_Y)
#endif
#endif

#ifdef DIVV_VAR

              U(DIVV_VAR,ix,iy,iz) = (U(VELX_VAR,ix+1,iy,iz)-U(VELX_VAR,ix-1,iy,iz))/del(DIR_X)
              if (NDIM > 1) then
                 U(DIVV_VAR,ix,iy,iz) = U(DIVV_VAR,ix,iy,iz) + &
                                   +(U(VELY_VAR,ix,iy+1,iz)-U(VELY_VAR,ix,iy-1,iz))/del(DIR_Y)
              if (NDIM == 3) then
                 U(DIVV_VAR,ix,iy,iz) = U(DIVV_VAR,ix,iy,iz) + &
                                   +(U(VELZ_VAR,ix,iy,iz+1)-U(VELZ_VAR,ix,iy,iz-1))/del(DIR_Z)
              endif
              endif
#endif
           enddo
        enddo
     enddo
     
     call Grid_releaseBlkPtr(blockID,U,CENTER)
  enddo

#endif 
#endif


  ! ***************************************
  ! *                                     *
  ! *     ADD BIERMANN BATTERY SOURCE     *
  ! *                                     *
  ! ***************************************
  !
  ! This subroutine call modifies the magnetic field directly using an
  ! external Biermann Battery source term. It will only preserve
  ! div(B)=0 in xy and rz geometries ...
 
  if ( (hy_useBiermann .or. hy_useBiermann1T) .and. hy_biermannSource ) then

 
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stamp('gcNeed - CENTER unmasked', '[hy_uhd_unsplit]')
     end if
#endif

     call Grid_fillGuardCells(CENTER, ALLDIR)
     call hy_uhd_biermannSource(blockCount, blockList, dt)
     
  end if

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if  
#endif
  
End Subroutine hy_uhd_unsplit
