!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_Unsplit/threadBlockList/hy_uhd_unsplit
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
!!  * Lee, D., "A solution accurate, efficient and stable unsplit staggered mesh scheme 
!!              for three dimensional magnetohydrodynamics", 243 (2013), 269-292, JCP
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

#include "Flash.h"
Subroutine hy_uhd_unsplit ( blockCount, blockList, dt, dtOld )

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
                         hy_updateHydroFluxes,&
                         hy_geometry,         &
                         hy_fluxCorVars,      &
                         hy_cfl,              &
                         hy_cfl_original,     &
                         hy_needScrchVars,    &
                         hy_dtmin,            &
                         hy_3TMode,           &
                         hy_shockDetectOn,    &
                         hy_threadBlockList

  use Driver_interface, ONLY : Driver_abortFlash

  use hy_uhd_interface, ONLY : hy_uhd_getRiemannState,  &
                               hy_uhd_getFaceFlux,      &
                               hy_uhd_unsplitUpdate,    &
                               hy_uhd_unitConvert,      &
                               hy_uhd_energyFix,        &
                               hy_uhd_putGravityUnsplit,&
                               hy_uhd_addGravityUnsplit,&
                               hy_uhd_shockDetect

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkData

  use Eos_interface, ONLY : Eos_wrapped

  use Logfile_interface, ONLY : Logfile_stampVarMask

  use Timers_interface, ONLY : Timers_start, Timers_stop

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

  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: &
       gravX,gravY,gravZ
  real :: faceAreas(GRID_ILO_GC:GRID_IHI_GC,     &
                    GRID_JLO_GC:GRID_JHI_GC,     &
                    GRID_KLO_GC:GRID_KHI_GC)
#else
  real, allocatable, dimension(:,:,:,:)   :: flx,fly,flz
  real, allocatable, dimension(:,:,:)   :: gravX, gravY, gravZ
  real, allocatable :: faceAreas(:,:,:)
#endif

  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
  real, pointer, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,save :: gcMaskLogged =.TRUE.
#endif

  real, pointer, dimension(:,:,:,:) :: U
  integer :: updateMode
  real    :: gravDtFactor
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
     if (1.2*hy_cfl < hy_cfl_original) then
        !! Slow recover (of factor of 1.2) to the original CFL once it gets to
        !! reduced to a smaller one in the presence of strong shocks.
        !! This variable CFL takes place in the following three cases using:
        !! (1) use_hybridOrder = .true.,
        !! (2) use_hybridOrder = .true., or
        !! (3) BDRY_VAR is defined and used for stationary objects.
        hy_cfl = 1.2*hy_cfl
     else
        hy_cfl = hy_cfl_original
     endif
     hy_dtmin = huge(1.0)
  endif


  !$omp parallel if (hy_threadBlockList) &
  !$omp default(none) &
  !$omp firstprivate(level,updateMode) &
  !$omp private(i,blockID,del,blkLimits,blkLimitsGC,datasize,&
  !$omp flx,fly,flz,gravX,gravY,gravZ,U,scrch_Ctr,faceAreas,hy_SpcR,hy_SpcL,hy_SpcSig) &
  !$omp shared(blockCount,blockList,dt,dtOld,gcMask,gcMaskLogged,hy_units,&
  !$omp hy_unsplitEosMode,hy_useGravity,hy_gref,hy_fluxCorrect,&
  !$omp hy_updateHydroFluxes,hy_eosModeAfter,hy_useGravHalfUpdate,&
  !$omp hy_geometry,hy_fluxCorVars,hy_needScrchVars,&
  !$omp hy_3TMode)

  !$omp do schedule(static)  
  do i=1,blockCount

     blockID = blockList(i)

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        call hy_uhd_unitConvert(blockID,FWDCONVERT)
     endif

     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1

#ifdef FLASH_UHD_NEED_SCRATCHVARS
     call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)      
#else
     allocate(scrch_Ctr(HY_NSCRATCH_VARS,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#if (NSPECIES+NMASS_SCALARS) > 0
     allocate(  hy_SpcR(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
     allocate(  hy_SpcL(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
     allocate(hy_SpcSig(HY_NSPEC,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS),NDIM))
#endif     
#endif

#ifndef FIXEDBLOCKSIZE
     allocate(gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! ************************************************************************
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


     if (hy_updateHydroFluxes) then
        !! ************************************************************************
        !! Calculate Riemann (interface) states
        !! Note: gravX(:,:,:) - gravity at n

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
        hy_SpcL=0.
        hy_SpcR=0.
        hy_SpcSig=0.
#endif
#endif
        call Timers_start("RiemannState")
        call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del, &
                                    gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:),scrch_Ctr,&
                                    hy_SpcR,hy_SpcL,hy_SpcSig)
        call Timers_stop("RiemannState")
     endif !! End of if (hy_updateHydroFluxes) then

#ifndef FIXEDBLOCKSIZE
     allocate(flx(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(fly(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(flz(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif


     !! ************************************************************************
     !! Calculate high order Godunov fluxes
     !! Initialize arrays with zero
     flx = 0.
     fly = 0.
     flz = 0.
     call Timers_start("getFaceFlux")
     call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del,flx,fly,flz,scrch_Ctr,&
                             hy_SpcR,hy_SpcL,hy_SpcSig)
     call Timers_stop("getFaceFlux")
     !! ************************************************************************
     !! Unsplit update for conservative variables from n to n+1 time step
#ifndef FLASH_GRID_UG
     if ((.not. hy_needScrchVars) .or. (hy_needScrchVars .and. .not. hy_fluxCorrect)) then
#endif
        call Timers_start("unsplitUpdate")
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ)
        call Timers_stop("unsplitUpdate")
#ifdef FLASH_UHD_3T
        call hy_uhd_unsplitUpdateMultiTemp&
             (blockID,updateMode,blkLimits,dataSize,dt,del,flx,fly,flz)
#endif
#ifndef FLASH_GRID_UG
     endif
#endif

     if (.not. hy_fluxCorrect) then
#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)

        if ( hy_units .NE. "none" .and. hy_units .NE. "NONE" ) then
           !! Convert unit
           call hy_uhd_unitConvert(blockID,BWDCONVERT)
        endif

!#ifndef FLASH_EOS_GAMMA
        !! Call to Eos
        call Eos_wrapped(hy_eosModeAfter, blkLimits, blockID)
!#endif
#endif /* ifndef GRAVITY */  

     else !! if Flux correction is used.
          !! Flux conservation calls on AMR:
          !! Correct fluxes at each block boundary where coarse and fine
          !! blocks are neighboring each other.

        if (hy_geometry /= CARTESIAN) then
           ! we are using consv_fluxes and need to divide by face areas
           call Grid_getBlkData(blockID,CELL_FACEAREA,ILO_FACE, EXTERIOR, &
                                (/1,1,1/), faceAreas, datasize)
           call Grid_putFluxData(blockID,IAXIS,flx,datasize,hy_fluxCorVars,faceAreas)
           if (NDIM > 1) then
              call Grid_getBlkData(blockID,CELL_FACEAREA,JLO_FACE, EXTERIOR, &
                                   (/1,1,1/), faceAreas, datasize)
              call Grid_putFluxData(blockID,JAXIS,fly,datasize,hy_fluxCorVars,faceAreas)
              if (NDIM > 2) then
                 call Grid_getBlkData(blockID,CELL_FACEAREA,KLO_FACE, EXTERIOR, &
                                      (/1,1,1/), faceAreas, datasize)
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

#ifdef FLASH_UHD_NEED_SCRATCHVARS
     call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)      
#else
     deallocate(scrch_Ctr)
#if (NSPECIES+NMASS_SCALARS) > 0
     deallocate(hy_SpcR)
     deallocate(hy_SpcL)
     deallocate(hy_SpcSig)
#endif     
#endif

  end do
  !$omp end do
  !! End of leaf block do-loop before flux conserve call


  !! ***************************************************************************
  !! Third part of advancement                                                 *
  !! ***************************************************************************
  !! Do this part only if refining and flux correcting
  if (hy_fluxCorrect) then

     !! ************************************************************************
     !! Conservation of Fluxes at each block boundary
     !$omp single
     call Grid_conserveFluxes(ALLDIR,level)
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
        allocate(    gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(    gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(    gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
        allocate(  faceAreas(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

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
        flx = 0.
        fly = 0.
        flz = 0.

#ifdef FLASH_UHD_NEED_SCRATCHVARS
        call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
#endif

        if (hy_needScrchVars) then
           flx(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
                =scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)
#if NDIM > 1
           fly(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
                =scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)
#if NDIM == 3
           flz(HY_DENS_FLUX:HY_END_FLUX,:,:,:)&
                =scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-2,:,:,:)
#endif
#endif
        endif

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
        call Timers_start("unsplitUpdate")
        call hy_uhd_unsplitUpdate(blockID,updateMode,dt,del,datasize,blkLimits,&
                                  blkLimitsGC,flx,fly,flz,gravX,gravY,gravZ)
        call Timers_stop("unsplitUpdate")

#ifndef GRAVITY /* if gravity is included we delay energy fix until we update gravity at n+1 state */
        !! *********************************************************************
        !! Correct energy if necessary
        call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)
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
     allocate(gravX(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravY(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(gravZ(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
#endif

     !! *********************************************************************
     !! Get and add gravity to hydro variables (momenta and energy)         *
     !! *********************************************************************
     gravX = 0.
     gravY = 0.
     gravZ = 0.
     if (hy_useGravity) then
        if (hy_order > 1 .and. hy_useGravHalfUpdate) then
           gravDtFactor = 1.
        else
           gravDtFactor = 2.
        endif

        call hy_uhd_putGravityUnsplit(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
        gravX = gravX/hy_gref
        gravY = gravY/hy_gref
        gravZ = gravZ/hy_gref

        call hy_uhd_addGravityUnsplit(blockID,blkLimits,dataSize,gravDtFactor*dt,&
             gravX(:,:,:),gravY(:,:,:),gravZ(:,:,:))
     endif

     !! *********************************************************************
     !! Correct energy if necessary
     call hy_uhd_energyFix(blockID,blkLimits,dt,dtOld,del,hy_unsplitEosMode)

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
