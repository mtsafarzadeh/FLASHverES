!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/threadWithinBlock/hy_uhd_unsplit
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
!!  Performs Hydro update in a directionally unsplit fashion over a set
!!  of blocks.
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
#define DEBUG_UHD
#endif
#define DEBUG_GRID_GCMASK

Subroutine hy_uhd_unsplit ( blockCount, blockList, dt, dtOld )

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkData

#include "Flash.h"
  use Hydro_data, ONLY : hy_fluxCorrect,      &
                         hy_gref,             &
                         hy_useGravity,       &
                         hy_order,            &
                         hy_units,            &
                         hy_gcMaskSize,       &
                         hy_gcMask,           &
                         hy_unsplitEosMode,   &
                         hy_eosModeAfter,     &
                         hy_useGravHalfUpdate,&
                         hy_useGravPotUpdate, &
                         hy_gravConsv,        &
                         hy_updateHydroFluxes,&
                         hy_geometry,         &
                         hy_fluxCorVars,      &
                         hy_cfl,              &
                         hy_cfl_original,     &
                         hy_needScrchVars,    &
                         hy_3TMode,           &
                         hy_threadWithinBlock,&
                         hy_shockDetectOn
#ifndef FLASH_UHD_NEED_SCRATCHVARS
  use Hydro_data, ONLY : scrch_Ctr => hy_scrchCtr
#endif

  use Driver_interface, ONLY : Driver_abortFlash

  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_putGravityUnsplit,&
                               hy_uhd_addGravityUnsplit,&
                               hy_uhd_shockDetect

  use Eos_interface, ONLY : Eos_wrapped

  use Logfile_interface, ONLY : Logfile_stampVarMask

!!$  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Gravity_interface, ONLY : Gravity_potentialListOfBlocks

  implicit none

#include "constants.h"
#include "Eos.h"
#include "UHD.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) ::  blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList
  real,    INTENT(IN) :: dt, dtOld
  !! -----------------------------------------------------

  integer, dimension(MDIM) :: datasize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: i,blockID,level
  integer :: ix,iy,iz
  real, dimension(MDIM) :: del
  logical :: gcMask(hy_gcMaskSize)
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
  real, allocatable, dimension(:,:,:,:)   :: flx,fly,flz
  real, allocatable, dimension(:,:,:,:)   :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)
#endif

  real, pointer, dimension(:,:,:,:) :: U
#ifdef FLASH_UHD_NEED_SCRATCHVARS
  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
#endif

#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif
  integer :: updateMode
  real    :: gravDtFactor
  logical :: halfTimeGravUpdate
  !! End of data declaration ***********************************************
  level = 0

#ifdef FLASH_GRID_PARAMESH2
  call Driver_abortFlash("The unsplit Hydro solver only works with PARAMESH 3 or 4!")
#endif


#ifdef FLASH_GRID_PARAMESH3OR4
  if (hy_fluxCorrect) then
     updateMode = UPDATE_INTERIOR
  else
     updateMode = UPDATE_ALL
  end if
#endif


#ifdef FLASH_GRID_UG
  updateMode = UPDATE_ALL
  hy_fluxCorrect = .false.
#endif


  !! ***************************************************************************
  !! Shock detection before hydro                                              *
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
  !! Call guardcell filling with Eos before hydro                              *
  !! ***************************************************************************

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     call Logfile_stampVarMask(hy_gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcNeed')
  end if
#endif

  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&
       maskSize=hy_gcMaskSize, mask=hy_gcMask,makeMaskConsistent=.true.,&
       doLogMask=.NOT.gcMaskLogged)


  !! ***************************************************************************
  !! First part of advancement                                                 *
  !! ***************************************************************************
  !! Loop over the blocks

  !! Retain the original cfl that may have been changed in some leaf blocks.
  if (hy_updateHydroFluxes) then
     hy_cfl = hy_cfl_original
  endif



  do i=1,blockCount

     blockID = blockList(i)

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        call hy_uhd_unitConvert(blockID,FWDCONVERT)
     endif

     call Grid_getDeltas(blockID,del)

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

#ifndef FLASH_UHD_NEED_SCRATCHVARS
        allocate(scrch_Ctr(HY_NSCRATCH_VARS,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        scrch_Ctr = 0.
#endif

#ifndef FIXEDBLOCKSIZE
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
                                    gravX(2,:,:,:),gravY(2,:,:,:),gravZ(2,:,:,:))

        ! Note: Two different ways of handling gravity:
        ! 1. With gravity calculated potential at n+1/2, Riemann states do not include gravity
        !    source terms at this point, and will include them in hy_uhd_addGravityUnsplit later
        !    to the primitive Riemann variables (not available for conservative formulation).
        ! 2. With gravity extrapolated from n-1 & n states, gravity source terms have been
        !    included to Riemann states in conservative formulation in hy_uhd_getRiemannState.

     endif !! End of if (hy_updateHydroFluxes) then


!!\ BEGIN EXITING OUT THIS BLOCK DO-LOOP FOR GRAVITATIONAL POTENTIAL UPDATE =============================
!!\ Perform Gravity potential update by exiting out this block do-loop and calling                      /
!!\ Poisson solver for the gravitational potential                                                      /
!!\ NOTE: IF GRAVITATIONAL POTENTIAL UPDATE IS NOT CHOSEN AT SETUP, THE FOLLOWING PART OF GRAVITATIONAL /
!!\       UPDATE WON'T BE EXECUTED.                                                                     /
#ifdef FLASH_UHD_NEED_SCRATCHVARS
!!\ BEGIN EXITING OUT THIS BLOCK DO-LOOP FOR GRAVITATIONAL POTENTIAL UPDATE                             !
!!\                                                                                                     !
!!                                                                                                      !
#ifndef FIXEDBLOCKSIZE
     deallocate(gravX)                                                                                  !
     deallocate(gravY)                                                                                  !
     deallocate(gravZ)                                                                                  !
#endif
                                                                                                        !
  enddo                                                                                                 !
                                                                                                        !
                                                                                                        !
                                                                                                        !
#ifdef GPOT_VAR
  if (hy_useGravity .and. hy_useGravPotUpdate) then                                                     !
     !! ***************************************************************************                     !
     !! Gravity calculation at n+1/2 by calling Poisson solver                    *                     !
     !! ***************************************************************************                     !
     ! At this stage, U(DENS_VAR,:,:,:) is averaged extrapolated (or bilinear) density                  !
     ! at n+1/2 values to cell center (already done in getRiemann).                                     !
     ! This density is used to compute gravity components at n+1/2.                                     !
     ! Note: In case hy_useGravPotUpdate = .false., gravity components will be                          !
     !       extrapolated to n+1/2.                                                                     !
     call Gravity_potentialListOfBlocks(blockCount, blockList)                                          !
                                                                                                        !
     !! Fill guardcells for only the gpot variable                                                      !
     gcMask  =.false.                                                                                   !
     gcMask(GPOT_VAR) = .true.                                                                          !
#if defined(GPOL_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
     if (hy_useGravPotUpdate) gcMask(GPOL_VAR) = .true.                                                 !
#endif
     gcMask(ENER_VAR) = .true.                                                                          !
#ifdef EINT_VAR
     gcMask(EINT_VAR) = .true.                                                                          !
#endif
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then                                                                        !
        call Logfile_stampVarMask(gcMask, .TRUE., '[hy_uhd_unsplit]', 'gcWant[Pot]')                    !
     end if                                                                                             !
#endif
                                                                                                        !
     !! Guardcell filling for gpot                                                                      !
     call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.,&                                              !
          maskSize=NUNK_VARS,mask=gcMask,makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskLogged)         !
                                                                                                        !
  endif                                                                                                 !
#endif
                                                                                                        !
                                                                                                        !
  !! ***************************************************************************                        !
  !! Second part of advancement related to gravity using potential             *                        !
  !! ***************************************************************************                        !
  !! Loop over the blocks                                                                               !
  do i=1,blockCount                                                                                     !
                                                                                                        !
     blockID = blockList(i)                                                                             !
     call Grid_getDeltas(blockID,del)                                                                   !
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)                                         !
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1                                !
#ifndef FIXEDBLOCKSIZE
     allocate(gravX(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))                                 !
     allocate(gravY(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))                                 !
     allocate(gravZ(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))                                 !
#endif
     !! ************************************************************************                        !
     !! Get gravity                                                                                     !
     gravX = 0.                                                                                         !
     gravY = 0.                                                                                         !
     gravZ = 0.                                                                                         !
                                                                                                        !
     !! This part of the gravity update is only getting called when hy_useGravPotUpdate                 !
     !! is chosen instead of extrapolating gravity components.                                          !
     if (hy_useGravity .and. hy_useGravHalfUpdate .and. hy_useGravPotUpdate) then                       !
        call Grid_getBlkPtr(blockID,U,CENTER)                                                           !
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)                                              !
                                                                                                        !
        ! Reset density in unk to n state values as we are done with                                    !
        ! computing gravities at n+1/2 using density at n+1/2.                                          !
        U(DENS_VAR,:,:,:) = scrch_Ctr(VAR2_SCRATCH_CENTER_VAR,:,:,:)                                    !
                                                                                                        !
        call Grid_releaseBlkPtr(blockID,U,CENTER)                                                       !
        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)                                          !
                                                                                                        !
        ! Now update cell-center grav accelerations (using new potential)                               !
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)          !
        gravX = gravX/hy_gref                                                                           !
        gravY = gravY/hy_gref                                                                           !
        gravZ = gravZ/hy_gref                                                                           !

        !half time update                                                                               !
        halfTimeGravUpdate = .true.                                                                     !
        ! and add gravity source terms to Riemann states at n+1/2.                                      !
        call hy_uhd_addGravityUnsplit(blockID,blkLimits,dataSize,dt,&                                   !
             gravX(2,:,:,:),gravY(2,:,:,:),gravZ(2,:,:,:),halfTimeGravUpdate)                           !
                                                                                                        !
     endif                                                                                              !
                                                                                                        !
     ! We need to call gravity again here for cases (AMR or UG) WITH scratch vars when 
     ! hy_useGravHalfUpdate=.false. and hy_useGravPotUpdate=.false.
     ! This gravity is for the final update in the corrector step.
     if (hy_useGravity .and. .not. hy_useGravPotUpdate) then
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)          !
        gravX = gravX/hy_gref                                                                           !
        gravY = gravY/hy_gref                                                                           !
        gravZ = gravZ/hy_gref                                                                           !
     endif


#endif
!!\ ENDING HERE:                                                                                        /
!!\ NOTE: IF GRAVITATIONAL POTENTIAL UPDATE IS NOT CHOSEN AT SETUP, THE ABOVE PART OF GRAVITATIONAL     /
!!\       UPDATE WON'T BE EXECUTED.                                                                     /
!!\ RETURNING TO CONTINUE THE FIRST BLOCK DO-LOOP FOR GRAVITATIONAL POTENTIAL UPDATE ===================/

#ifndef FIXEDBLOCKSIZE
     allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif


     !! ************************************************************************
     !! Calculate high order Godunov fluxes
     !! Initialize arrays with zero
     flx = 0.
     fly = 0.
     flz = 0.

     call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del,flx,fly,flz)


     !! ************************************************************************
     !! Unsplit update for conservative variables from n to n+1 time step
#ifndef FLASH_GRID_UG
     if ((.not. hy_needScrchVars) .or. (hy_needScrchVars .and. .not. hy_fluxCorrect)) then
#endif
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,dtOld,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ)

#ifdef FLASH_UHD_3T
        call hy_uhd_unsplitUpdateMultiTemp&
             (blockID,updateMode,blkLimits, dataSize, dt, del, flx, fly, flz)
#endif
#ifndef FLASH_GRID_UG
     endif
#endif

     if (.not. hy_fluxCorrect) then
#ifndef GRAVITY !! if gravity is included we delay energy fix until we update gravity at n+1 state 
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,del,hy_unsplitEosMode)

        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert unit
           call hy_uhd_unitConvert(blockID,BWDCONVERT)
        endif

        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
#endif !! ifndef GRAVITY 

     else !! if Flux correction is used.
          !! Flux conservation calls on AMR:
          !! Correct fluxes at each block boundary where coarse and fine
          !! blocks are neighboring each other.

        if (hy_geometry /= CARTESIAN) then
           ! we are using consv_fluxes and need to divide by face areas
           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
                                (/1,1,1/), &
                                faceAreas, datasize)
           call Grid_putFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
           if (NDIM > 1) then
              call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
                                   (/1,1,1/), &
                                   faceAreas, datasize)
              call Grid_putFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
              if (NDIM > 2) then
                 call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                      (/1,1,1/), &
                                      faceAreas, datasize)
                 call Grid_putFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
              endif
           endif
        else ! Cartesian geometry
           call Grid_putFluxData(blockID,IAXIS,flx,datasize)

           if (NDIM > 1) then
              call Grid_putFluxData(blockID,JAXIS,fly,datasize)
              if (NDIM > 2) then
                 call Grid_putFluxData(blockID,KAXIS,flz,datasize)
              endif
           endif
        endif
     end if

     
#ifndef FIXEDBLOCKSIZE
     deallocate(flx)
     deallocate(fly)
     deallocate(flz)
     deallocate(gravX)
     deallocate(gravY)
     deallocate(gravZ)
     deallocate(faceAreas)
#endif

#ifndef FLASH_UHD_NEED_SCRATCHVARS
     deallocate(scrch_Ctr)
#endif

  end do
  !! End of leaf block do-loop before flux conserve call


  !! ***************************************************************************
  !! Third part of advancement                                                 *
  !! ***************************************************************************
  !! Do this part only if refining and flux correcting
  if (hy_fluxCorrect) then

     !! ************************************************************************
     !! Conservation of Fluxes at each block boundary
     call Grid_conserveFluxes(ALLDIR,level)

     !! ************************************************************************
     !! Perform updates for every blocks
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

#ifdef FLASH_UHD_NEED_SCRATCHVARS
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
#endif

        !$omp parallel if (hy_threadWithinBlock) &
        !$omp default(none) &
        !$omp shared(flx,fly,flz,gravX,gravY,gravZ,hy_useGravity,hy_gref,&
        !$omp scrch_Ctr,blockID,blkLimitsGC,dataSize,dt,dtOld,hy_needScrchVars)

        !! Initialize arrays with zero
        !$omp workshare
        flx = 0.
        fly = 0.
        flz = 0.

        !! *********************************************************************
        !! Get gravity
        gravX = 0.
        gravY = 0.
        gravZ = 0.
        !$omp end workshare

        if (hy_useGravity) then
           !$omp single
           call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
           !$omp end single
           !$omp workshare
           gravX = gravX/hy_gref
           gravY = gravY/hy_gref
           gravZ = gravZ/hy_gref
           !$omp end workshare nowait
        endif


        !! Get modified conserved flux values at block interfaces:
        !! This is important especially at the block interface at 
        !! different refinement levels of neighboring blocks

        if (hy_needScrchVars) then
           !CD: The memory copies from scrch_Ctr into flx,fly,flz are the most
           !expensive parts of hy_uhd_unsplit.F90.  We group them together so
           !that they can be put into a single workshare.
           !$omp workshare
           flx(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
                =scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)
#if NDIM > 1
           fly(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
                =scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)
#if NDIM > 2
           flz(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
                =scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)
#endif
#endif
           !$omp end workshare nowait
        endif
        !$omp end parallel


        if (hy_geometry /= CARTESIAN) then
           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,IAXIS,flx,datasize)
        endif

#if NDIM > 1
        if (hy_geometry /= CARTESIAN) then
           call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,JAXIS,fly,datasize)
        endif

#if NDIM > 2
        if (hy_geometry /= CARTESIAN) then
           call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_getFluxData(blockID,KAXIS,flz,datasize,hy_fluxCorVars,faceAreas)
        else
           call Grid_getFluxData(blockID,KAXIS,flz,datasize)
        endif
#endif
#endif

#ifdef FLASH_UHD_NEED_SCRATCHVARS
        call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
#endif

        !! *********************************************************************
        !! Unsplit update for conservative variables from n to n+1 time step
        if (.not. hy_needScrchVars) then
           updateMode=UPDATE_BOUND
        else
           updateMode=UPDATE_ALL
        endif

        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,dtOld,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ)


#ifndef GRAVITY !! if gravity is included we delay energy fix until we update gravity at n+1 state 
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


        !! *********************************************************************
        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
#endif  !! if gravity is included, we calculate gravity for n+1 in the below 


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
  end if !! End of the flux conservation routine


#ifdef GRAVITY !! Perform this only when gravity is used 
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

        !$omp parallel if (hy_threadWithinBlock) &
        !$omp default(none) &
        !$omp shared(gravX,gravY,gravZ,hy_gref,&
        !$omp blockID,blkLimitsGC,dataSize,dt,dtOld)

        !! Initialize arrays with zero
        !$omp workshare
        gravX = 0.
        gravY = 0.
        gravZ = 0.
        !$omp end workshare
        !$omp single
        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        !$omp end single
        !$omp workshare
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref
        !$omp end workshare
        !$omp end parallel

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
#endif !! End of n+1 gravity coupling 




  ! Do 3T update...
#ifdef FLASH_UHD_3T
  call hy_uhd_multiTempAfter(blockCount, blockList, dt)
#endif

#ifdef DEBUG_GRID_GCMASK
  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
#endif

End Subroutine hy_uhd_unsplit
