!!****if* source/physics/Hydro/HydroMain/unsplit_old/threadWithinBlock/hy_uhd_getRiemannState
!!
!! NAME
!!
!!  hy_uhd_getRiemannState
!!
!! SYNOPSIS
!!
!!  hy_uhd_getRiemannState( integer(IN) :: blockID,
!!                          integer(IN) :: blkLimits,
!!                          integer(IN) :: blkLimitsGC,
!!                          real(IN)    :: dt,
!!                          integer(IN) :: del(MDIM),
!!                          real(IN)    :: gravX(:,:,:),
!!                          real(IN)    :: gravY(:,:,:),
!!                          real(IN)    :: gravZ(:,:,:),
!!                          logical(IN), optional :: normalFieldUpdate)
!!
!! DESCRIPTION
!!
!!  This routine computes the Riemann state values at cell interfaces using 
!!  the cell centered variables and store them in the scratch arrays.
!!
!!  A 2D Cartesian configuration of a single cell is shown:
!!
!!
!!           ---------------------
!!          |          yp         |
!!          |                     |
!!          |                     |
!!          |                     |
!!          |                     |
!!          |xm      (i,j)      xp|
!!          |                     |
!!          |                     |
!!          |                     |
!!  y       |                     |
!!  |       |          ym         |
!!  |        ---------------------
!!  |______x
!!
!!
!! ARGUMENTS
!!
!!  blockID     - local block ID
!!  blkLimits   - an array that holds the lower and upper indices of the section
!!                of block without the guard cells
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dt          - a current time step
!!  del         - deltas in each direction 
!!  gravX       - gravity component in x-direction at n+1/2 step
!!  gravY       - gravity component in y-direction at n+1/2 step
!!  gravZ       - gravity component in z-direction at n+1/2 step
!!  normalFieldUpdate - a logical switch to choose normal magnetic fields updates only
!!                      (needed for MHD only)
!!   
!!***

!!REORDER(4):U, V0, scrch_Ctr, scrch_gpot, B[xyz]
Subroutine hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del,&
                                  ogravX,ogravY,ogravZ,&
                                  hgravX,hgravY,hgravZ,&
                                  normalFieldUpdate)

  use Grid_interface,       ONLY : Grid_getBlkPtr,     &
                                   Grid_releaseBlkPtr, &
                                   Grid_getCellCoords
#include "Flash.h"
  use Hydro_data,           ONLY : hy_shockDetectOn,     &
                                   hy_order,             &
                                   hy_forceHydroLimit,   &
                                   hy_flattening,        &
                                   hy_useGravity,        &
                                   hy_useGravHalfUpdate, &
                                   hy_gravConsv,         &
                                   hy_useGravPotUpdate,  &
                                   hy_transOrder,        &
                                   hy_use3dFullCTU,      &
                                   hy_geometry,          &  
                                   hy_useHybridOrder,    &
                                   hy_eswitch,           &
                                   hy_threadWithinBlock, &
                                   hy_killdivB

#ifndef FLASH_USM_MHD
#ifndef FLASH_UHD_NEED_SCRATCHVARS
  use Hydro_data,           ONLY : scrch_Ctr => hy_scrchCtr
#endif
#endif
  use hy_uhd_slopeLimiters, ONLY : mc
  use hy_uhd_interface,     ONLY : hy_uhd_dataReconstOnestep, &
                                   hy_uhd_shockDetect,        &
                                   hy_uhd_prim2con,           &
                                   hy_uhd_con2prim,           &
                                   hy_uhd_eigenParameters,    &
                                   hy_uhd_eigenValue,         &
                                   hy_uhd_eigenVector

  implicit none

#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------------------------------------
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
  real,    intent(IN)   :: dt
  real,    intent(IN),dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
       ogravX,ogravY,ogravZ,hgravX,hgravY,hgravZ
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: ogravX,ogravY,ogravZ,hgravX,hgravY,hgravZ
#endif
  logical, intent(IN), optional :: normalFieldUpdate
  !! ---------------------------------------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax, i, j, k
  integer,dimension(MDIM) :: dataSize
  real, dimension(HY_VARINUMMAX) :: Wxp,  Wxn,  Wyp,  Wyn,  Wzp,  Wzn
  real, dimension(HY_VARINUMMAX) :: Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn, &
                                    Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn,&
                                    Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn
  real, pointer, dimension(:,:,:,:) :: U, scrch_gpot
  real, allocatable, dimension(:,:,:,:) :: V0
#if defined(FLASH_UHD_NEED_SCRATCHVARS) || defined(FLASH_USM_MHD)
  real, pointer, dimension(:,:,:,:) :: scrch_Ctr ! Pointer array for Riemann States
#endif
  integer :: dir

! MHD only-------------------------------------------------------------------------------
!CD: Remove conditional compilation so that we only need a single version of
!the OpenMP parallel region.
!#ifdef FLASH_USM_MHD
!#if NFACE_VARS > 0
!#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz
!#endif
!#endif
!#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

!CD: Remove conditional compilation so that we only need a single version of
!the OpenMP parallel region.
!#ifdef GRAVITY
  real,allocatable,dimension(:,:,:) :: gravXP,gravXN,gravYP,gravYN,gravZP,gravZN
  real, dimension(HY_VARINUM) :: consVar
!#endif
#ifdef FIXEDBLOCKSIZE
  real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: FlatCoeff,FlatTilde
#else
  real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)) :: &
       FlatCoeff, FlatTilde
#endif
  real :: Sp, dv1, dp1, dp2, presL,presR,hdt

  real :: nOneHalfDens

#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_IHI_GC) :: xCenter  
#else  
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter  
#endif   
   real :: Rinv, geoFac, eta, enth   
   integer :: velPhi, velTht, magPhi, magZ
   real :: sGeo_dens, sGeo_velx, sGeo_velp, sGeo_pres, sGeo_trans, sGeo_eint, sGeo_magz,sGeo_magp
   real, allocatable,dimension(:,:,:) :: DivU,soundSpeed
   integer :: k2,k3
   real :: tinyD, tinyP


  ! for 3d only --------------
  real, allocatable,dimension(:,:,:,:,:) :: sig
  real, allocatable,dimension(:,:,:,:,:) :: lambda
  real, allocatable,dimension(:,:,:,:,:,:) :: leig
  real, allocatable,dimension(:,:,:,:,:,:) :: reig
  real, dimension(HY_VARINUMMAX) :: TransFluxXY,TransFluxYZ,TransFluxZX,&
                                    TransFluxYX,TransFluxZY,TransFluxXZ
  real :: dt2dxdy6,dt2dydz6,dt2dzdx6,hdtdx,hdtdy,hdtdz

  logical :: cons=.false.
  real    :: cs,ca,cf,as,af,uN
  real, dimension(MDIM) :: beta
  integer :: kGrav
  logical :: TransX_updateOnly,TransY_updateOnly,TransZ_updateOnly
  ! for 3d only --------------

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(i0,j0,k0,imax,jmax,kmax,k2,k3,scrch_Ctr,V0,normalFieldUpdate,&
  !$omp blockID,blkLimitsGC,dt,del,ogravX,ogravY,ogravZ,FlatCoeff,&
  !$omp hy_useGravity,hy_useGravPotUpdate,U,hy_gravConsv,hy_useGravHalfUpdate,&
  !$omp hy_order,hy_flattening,FlatTilde,blkLimits,hdt,hy_shockDetectOn,dataSize,&
  !$omp hgravX,hgravY,hgravZ,Bx,By,Bz,hy_forceHydroLimit,&
  !$omp sig,lambda,reig,leig,scrch_gpot,hy_killdivB,&
  !$omp gravXP,gravXN,gravYP,gravYN,gravZP,gravZN,&
  !$omp hy_use3dFullCTU,hy_transOrder,hy_geometry,hy_useHybridOrder,&
  !$omp xCenter,DivU,kGrav,soundSpeed) &
  !$omp private(i,j,k,Vxppp,Vxnnn,Vyp,Vyn,Vypp,Vynn,Vyppp,Vynnn,Vzp,&
  !$omp Vzn,Vzpp,Vznn,Vzppp,Vznnn,Wxp,Wxn,Wyp,Wyn,Wzp,Wzn,nOneHalfDens,&
  !$omp dp1,dp2,dv1,presL,presR,Sp,magPhi,magZ,sGeo_magp,sGeo_magz,&
  !$omp consVar,Vxpp,Vxnn,TransX_updateOnly,TransY_updateOnly,TransZ_updateOnly,&
  !$omp dt2dxdy6,dt2dydz6,dt2dzdx6,TransFluxXY,TransFluxYZ,TransFluxZX,&
  !$omp TransFluxYX,TransFluxZY,TransFluxXZ,Rinv,velPhi,velTht,&
  !$omp sGeo_dens,sGeo_velx,sGeo_velp,sGeo_pres,sGeo_trans,sGeo_eint,&
  !$omp geoFac,cs,eta,enth,dir,tinyD,tinyP)

  !$omp single

  kGrav=0
#ifdef GRAVITY
  kGrav=1
#endif

  k2=0
  k3=0

#if NDIM == 1
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = 3
  jmax =-1
  k0   = 3
  kmax =-1
#elif NDIM == 2
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = 3
  kmax =-1
  k2=1
#elif NDIM == 3
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = blkLimits(LOW, KAXIS)
  kmax = blkLimits(HIGH,KAXIS)
  k2=1
  k3=1
#endif

  hdt = 0.5*dt

!! DEV-dongwook: This call now is moved to the begining of hydro,
!! in order to provide an easier gc fill.
!! We noticed that divB blows up on AMR with BBT (e.g., vulcan3d)
!! because the tagging information is not available 
!! across the neighboring blocks at same refinement levels.
!! The fix was to call GC for RHCD_VAR after tagging, which is now
!! done at the very begining of hydro/MHD (see hy_uhd_unsplit).
!!
!!$  if (hy_shockDetectOn) then
!!$     !! Call shock detect algorithm to determine how to calculate 
!!$     !! transversal differencings in hy_dataReconstOnestep.
!!$     call hy_uhd_shockDetect(blockID)
!!$  endif
  !$omp end single


! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
  if (hy_order > 1) then
#if NFACE_VARS > 0
#if NDIM > 1
     !$omp single
     call Grid_getBlkPtr(blockID,Bx,FACEX)
     call Grid_getBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
     call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
     !$omp end single
#endif
#endif
  endif
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------
     !$omp single
     call Grid_getBlkPtr(blockID,U,CENTER)

     if (hy_geometry /= CARTESIAN) then
        ! Grab cell x-coords for this block  
        call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
     endif

     !! Allocate a temporary cell-centered array at each local block
     dataSize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
     allocate(V0(HY_VARINUMMAX,dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))

     !! For 3D transverse fluxes
     allocate(   sig(HY_VARINUMMAX,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(lambda(HY_WAVENUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(  leig(HY_WAVENUM,HY_VARINUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(  reig(HY_VARINUM,HY_WAVENUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate(  divU(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
     allocate( soundSpeed(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
     !$omp end single

     !$omp workshare
     V0(HY_DENS,:,:,:) = U(DENS_VAR,:,:,:)
     V0(HY_VELX:HY_VELZ,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)
     V0(HY_EINT,:,:,:) = U(EINT_VAR,:,:,:)
     !$omp end workshare nowait

#ifdef FLASH_UHD_3T
     !$omp workshare
     V0(HY_EELE,:,:,:) = U(EELE_VAR,:,:,:)
     V0(HY_EION,:,:,:) = U(EION_VAR,:,:,:)
     V0(HY_ERAD,:,:,:) = U(ERAD_VAR,:,:,:)
     !$omp end workshare nowait
#endif



! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
     if (.not. hy_forceHydroLimit) then
        !$omp workshare
        V0(HY_MAGX:HY_MAGZ,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)
        !$omp end workshare nowait
     else
        !$omp workshare
        V0(HY_MAGX:HY_MAGZ,:,:,:) = 0.
        !$omp end workshare nowait
     endif
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------
     !$omp workshare
     V0(HY_PRES,:,:,:) = U(PRES_VAR,:,:,:)
     V0(HY_GAMC,:,:,:) = U(GAMC_VAR,:,:,:)
     V0(HY_GAME,:,:,:) = U(GAME_VAR,:,:,:)
     !$omp end workshare

     !! pointers
     !$omp single
     call Grid_getBlkPtr(blockID,scrch_gpot,SCRATCH_CTR)
     ! We can store dens in unk to scratch_XP07 in hydro and scratch_XP10 in MHD
     if (hy_useGravity .and. hy_useGravPotUpdate) then
        scrch_gpot(VAR2_SCRATCH_CENTER_VAR,:,:,:)= U(DENS_VAR,:,:,:)
     endif
     !! Release pointers
     call Grid_releaseBlkPtr(blockID,scrch_gpot,SCRATCH_CTR)
     !$omp end single


!! DEV-Dongwook: I don't know how it ever worked with this zeroing out of
!! RHCD_VAR here, after calling shockDetect in the above???
!!$#ifdef RHCD_VAR
!!$     !$omp workshare
!!$     U(RHCD_VAR,:,:,:) = 0.
!!$     !$omp end workshare
!!$#endif

     ! call Grid_releaseBlkPtr(blockID,U,CENTER) ! move this - SMC


     !! -----------------------------------------------------------------------!
     !! Compute divergence of velocity fields and (magneto)sonic speed --------!
     !! -----------------------------------------------------------------------!
     if (hy_useHybridOrder) then
        !$omp single
        DivU = 0.
        soundSpeed = 0.

        do k=blkLimitsGC(LOW,KAXIS)+k3,blkLimitsGC(HIGH,KAXIS)-k3
           do j=blkLimitsGC(LOW,JAXIS)+k2,blkLimitsGC(HIGH,JAXIS)-k2
              do i=blkLimitsGC(LOW,IAXIS)+1,blkLimitsGC(HIGH,IAXIS)-1
                 ! (1) compute local (magneto)sonic speed. 
                 soundSpeed(i,j,k) = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)
#ifdef FLASH_USM_MHD
                 soundSpeed(i,j,k) = soundSpeed(i,j,k) &
                      + dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),U(MAGX_VAR:MAGZ_VAR,i,j,k))
#endif
                 soundSpeed(i,j,k) = soundSpeed(i,j,k)/U(DENS_VAR,i,j,k)
                 soundSpeed(i,j,k) = sqrt(soundSpeed(i,j,k))

                 ! (2) compute (undivided) divergence of velocity fields
                 DivU(i,j,k) = (U(VELX_VAR,i+1,j,k) - U(VELX_VAR,i-1,j,k))
                 if (NDIM > 1) then
                    DivU(i,j,k) = DivU(i,j,k)+ (U(VELY_VAR,i,j+1,k) - U(VELY_VAR,i,j-1,k))
                    if (NDIM == 3) then
                       DivU(i,j,k) = DivU(i,j,k)+ (U(VELZ_VAR,i,j,k+1) - U(VELZ_VAR,i,j,k-1))
                    endif
                 endif

              enddo
           enddo
        enddo
        DivU = 0.5*DivU
        !$omp end single
     endif

     !! -----------------------------------------------------------------------!
     !! Save old velocities ---------------------------------------------------!
     !! -----------------------------------------------------------------------!     
     !$omp workshare

#ifdef VOLX_VAR
              U(VOLX_VAR,:,:,:) = U(VELX_VAR,:,:,:)
#endif
#ifdef VOLY_VAR
              U(VOLY_VAR,:,:,:) = U(VELY_VAR,:,:,:)
#endif
#ifdef VOLZ_VAR
              U(VOLZ_VAR,:,:,:) = U(VELZ_VAR,:,:,:)
#endif

     !$omp end workshare

     !! -----------------------------------------------------------------------!
     !! Flattening begins here ------------------------------------------------!
     !! -----------------------------------------------------------------------!
     if (hy_flattening) then
        !$omp workshare
        FlatTilde = 0.0
        FlatCoeff = 0.0
        !$omp end workshare

        ! Flat tilde
#if NDIM == 3
        !$omp do schedule(static)
#endif
        do k=k0-2,kmax+2
#if NDIM == 2
           !$omp do schedule(static)
#endif
           do j=j0-2,jmax+2
#if NDIM == 1
              !$omp do schedule(static)
#endif
              do i=i0-2,imax+2
                 do dir=1,NDIM

                    select case (dir)
                    case (DIR_X)
                       dp1   = (V0(HY_PRES,i+1,j,k)-V0(HY_PRES,i-1,j,k))
                       dp2   = (V0(HY_PRES,i+2,j,k)-V0(HY_PRES,i-2,j,k))
                       dv1   =  V0(HY_VELX,i+1,j,k)-V0(HY_VELX,i-1,j,k)
                       presL = V0(HY_PRES,i-1,j,k)
                       presR = V0(HY_PRES,i+1,j,k)
#if NDIM > 1
                    case (DIR_Y)
                       dp1   = (V0(HY_PRES,i,j+1,k)-V0(HY_PRES,i,j-1,k))
                       dp2   = (V0(HY_PRES,i,j+2,k)-V0(HY_PRES,i,j-2,k))
                       dv1   =  V0(HY_VELY,i,j+1,k)-V0(HY_VELY,i,j-1,k)
                       presL = V0(HY_PRES,i,j-1,k)
                       presR = V0(HY_PRES,i,j+1,k)
#if NDIM > 2
                    case (DIR_Z)
                       dp1   = (V0(HY_PRES,i,j,k+1)-V0(HY_PRES,i,j,k-1))
                       dp2   = (V0(HY_PRES,i,j,k+2)-V0(HY_PRES,i,j,k-2))
                       dv1   =  V0(HY_VELZ,i,j,k+1)-V0(HY_VELZ,i,j,k-1)
                       presL = V0(HY_PRES,i,j,k-1)
                       presR = V0(HY_PRES,i,j,k+1)
#endif
#endif
                    end select

                    if (abs(dp2) > 1.e-15) then
                       Sp = dp1/dp2 - 0.75
                    else
                       Sp = 0.
                    endif

                    FlatTilde(dir,i,j,k) = max(0.0, min(1.0,10.0*Sp))
                    if ((abs(dp1)/min(presL,presR) < 1./3.) .or. dv1 > 0. ) then
                       FlatTilde(dir,i,j,k) = 0.
                    endif


                 enddo
              enddo
#if NDIM == 1
              !$omp end do
#endif
           enddo
#if NDIM == 2
           !$omp end do
#endif
        enddo
#if NDIM == 3
        !$omp end do
#endif
        !CD: Implict barrier needed because there are references to e.g k-1 and
        !k+1 elements of FlatTilde in the next nested loop.


        ! Flat coefficient
#if NDIM == 3
        !$omp do schedule(static)
#endif
        do k=k0-2,kmax+2
#if NDIM == 2
           !$omp do schedule(static)
#endif
           do j=j0-2,jmax+2
#if NDIM == 1
              !$omp do schedule(static)
#endif
              do i=i0-2,imax+2
                 do dir=1,NDIM

                    select case (dir)
                    case (DIR_X)
                       dp1   = (V0(HY_PRES,i+1,j,k)-V0(HY_PRES,i-1,j,k))

                       if ( dp1 < 0.0 ) then
                          FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i+1,j,k))
                       elseif (dp1 == 0.) then
                          FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                       else
                          FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i-1,j,k))
                       endif
#if NDIM > 1
                    case (DIR_Y)
                       dp1   = (V0(HY_PRES,i,j+1,k)-V0(HY_PRES,i,j-1,k))

                       if ( dp1 < 0.0 ) then
                          FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j+1,k))
                       elseif (dp1 == 0.) then
                          FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                       else
                          FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j-1,k))
                       endif
#if NDIM > 2
                    case (DIR_Z)
                       dp1   = (V0(HY_PRES,i,j,k+1)-V0(HY_PRES,i,j,k-1))

                       if ( dp1 < 0.0 ) then
                          FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j,k+1))
                       elseif (dp1 == 0.) then
                          FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                       else
                          FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j,k-1))
                       endif
#endif
#endif
                    end select

                 enddo
              enddo
#if NDIM == 1
              !$omp end do
#endif
           enddo
#if NDIM == 2
           !$omp end do
#endif
        enddo
#if NDIM == 3
        !$omp end do
#endif
     endif
     !! -----------------------------------------------------------------------!
     !! End of flattening -----------------------------------------------------!
     !! -----------------------------------------------------------------------!

     !! -----------------------------------------------------------------------!
     !! Arrays for gravity components -----------------------------------------!
     !! -----------------------------------------------------------------------!
#ifdef GRAVITY
     !$omp single
     allocate(gravXP(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
     allocate(gravXN(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
     if (NDIM > 1) then
        allocate(gravYP(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
        allocate(gravYN(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
        if (NDIM > 2) then
           allocate(gravZP(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
           allocate(gravZN(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)))
        endif
     endif
     !$omp end single
#endif


  !! Get block pointer for storages of Riemann states
#if defined(FLASH_UHD_NEED_SCRATCHVARS) || defined(FLASH_USM_MHD)
  !$omp single
  call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
  !$omp end single
#endif


! MHD only-------------------------------------------------------------------------------
!!$#ifdef FLASH_USM_MHD
!!$  endif !! End of if (normalFieldUpdate .and. (hy_order >= 2)) then
!!$#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

  !! Initialize arrays
  Vxpp = 0.
  Vxnn = 0.
  Vxppp= 0.
  Vxnnn= 0.

  Vyp  = 0.
  Vyn  = 0.
  Vypp = 0.
  Vynn = 0.
  Vyppp= 0.
  Vynnn= 0.

  Vzp  = 0.
  Vzn  = 0.
  Vzpp = 0.
  Vznn = 0.
  Vzppp= 0.
  Vznnn= 0.


  !! Compute Riemann states at each cell
#if NDIM < 3
  do k=k0-2,kmax+2
#if NDIM == 2
     !$omp do schedule(static)
#endif
     do j=j0-2,jmax+2
#if NDIM == 1
        !$omp do schedule(static)
#endif
        do i=i0-2,imax+2
#else
  ! Extra stencil is needed for 3D to correctly calculate transverse fluxes 
  !(i.e., cross derivatives in x,y, & z)
  !$omp do schedule(static)
  do k=k0-3,kmax+3
     do j=j0-3,jmax+3
        do i=i0-3,imax+3
#endif


! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
           if (.not. normalFieldUpdate) then
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

              TransX_updateOnly = .false.
              TransY_updateOnly = .false.
              TransZ_updateOnly = .false.

              if (i+2 .le. blkLimitsGC(HIGH,IAXIS)) then
                 Vxpp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i+2,j, k )
              else
                 TransX_updateOnly = .true.
              endif

              if (i-2 .ge. blkLimitsGC(LOW, IAXIS)) then
                 Vxnn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i-2,j, k )
              else
                 TransX_updateOnly = .true.
              endif

#if NGUARD > 4
              Vxppp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i+3,j, k )
              Vxnnn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i-3,j, k )
#endif
#if NDIM >= 2
              Vyp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j+1, k )
              Vyn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j-1, k )

              if (j+2 .le. blkLimitsGC(HIGH,JAXIS)) then
                 Vypp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j+2, k )
              else
                 TransY_updateOnly = .true.
              endif

              if (j-2 .ge. blkLimitsGC(LOW, JAXIS)) then
                 Vynn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j-2, k )
              else
                 TransY_updateOnly = .true.
              endif


#if NGUARD > 4
              Vyppp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j+3, k )
              Vynnn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j-3, k )
#endif

#if NDIM == 3
              Vzp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i, j, k+1)
              Vzn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i, j, k-1)

              if (k+2 .le. blkLimitsGC(HIGH,KAXIS)) then
                 Vzpp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i, j, k+2)
              else
                 TransZ_updateOnly = .true.
              endif

              if (k-2 .ge. blkLimitsGC(LOW, KAXIS)) then
                 Vznn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i, j, k-2)
              else
                 TransZ_updateOnly = .true.
              endif
#if NGUARD > 4
              Vzppp(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j, k+3 )
              Vznnn(HY_DENS:HY_END_VARS)=V0(HY_DENS:HY_END_VARS,i,j, k-3 )
#endif

#endif
#endif


              !! Left and right Riemann state reconstructions
              call hy_uhd_dataReconstOnestep&
                   (blockID,blkLimitsGC,i,j,k,dt,del,&
                    ogravX,ogravY,ogravZ,&
                    DivU,soundSpeed,&
                    V0(HY_DENS:HY_GRAV, i, j,k), &
                    V0(HY_DENS:HY_GRAV,i+1,j,k), &
                    V0(HY_DENS:HY_GRAV,i-1,j,k), &
                    Vyp,  Vyn,  &
                    Vzp,  Vzn,  & 
                    Vxpp, Vxnn, &
                    Vypp, Vynn, &
                    Vzpp, Vznn, &
                    Vxppp,Vxnnn,&
                    Vyppp,Vynnn,&
                    Vzppp,Vznnn,&
                    FlatCoeff,  &
                    TransX_updateOnly,&
                    TransY_updateOnly,&
                    TransZ_updateOnly,&
                    Wxp,Wxn,Wyp,Wyn,Wzp,Wzn,&
                    sig   (1,1,i,j,k), &
                    lambda(1,1,i,j,k),  &
                    leig  (1,1,1,i,j,k),&
                    reig  (1,1,1,i,j,k) )

              if (hy_geometry /= CARTESIAN) then
                 !! **************************************************************
                 !! Add geometric source terms in left and Right States          *
                 !! **************************************************************
                 Rinv = 1./xCenter(i)
                 select case (hy_geometry)
                 case (CYLINDRICAL)
                    velPhi = HY_VELZ
#ifdef FLASH_USM_MHD
                    magPhi = HY_MAGZ
                    magZ   = HY_MAGY
#endif
                    geoFac = Rinv
                 case (POLAR)
                    velPhi = HY_VELY
#ifdef FLASH_USM_MHD
                    magPhi = HY_MAGY
                    magZ = HY_MAGZ
#endif
                    geoFac = Rinv
                 case (SPHERICAL)
                    velPhi = HY_VELZ
                    velTht = HY_VELY
                    geoFac = 2.*Rinv
                 end select
                 
                 cs = sqrt(V0(HY_GAMC,i,j,k)*V0(HY_PRES,i,j,k)/V0(HY_DENS,i,j,k))
                 eta = (abs(V0(HY_VELX,i,j,k)) + cs) * dt/del(DIR_X)
                 eta = (1.-eta) / (cs*dt*abs(geoFac))
                 eta = min(1.,eta)
                 !! comment this line not to use the Colella hack
!!$                 geoFac = eta * geoFac
                 !! end of the Colella hack
                 enth = ((V0(HY_DENS,i,j,k)*V0(HY_EINT,i,j,k) + V0(HY_PRES,i,j,k))/V0(HY_DENS,i,j,k))/cs**2
                 !! right/left state source terms



                 !! right/left state source terms
                 if ( hy_geometry == CYLINDRICAL .OR. &
                      hy_geometry == POLAR .OR. &
                      hy_geometry == SPHERICAL ) then

                    sGeo_dens = -V0(HY_DENS,i,j,k) * V0(HY_VELX,i,j,k) * geoFac !src[DN]
                    sGeo_velx = (V0(velPhi,i,j,k)**2) * geoFac                  !src[VR]  
#ifdef FLASH_USM_MHD
                    sGeo_velx = sGeo_velx - (V0(magPhi,i,j,k)**2) * geoFac / V0(HY_DENS,i,j,k)
#endif                  
                    sGeo_velp = -V0(velPhi,i,j,k) * V0(HY_VELX,i,j,k) * geoFac !src[Vphi]
#ifdef FLASH_USM_MHD
                    sGeo_velp = sGeo_velp + V0(magPhi,i,j,k) * V0(HY_MAGX,i,j,k) * geoFac / V0(HY_DENS,i,j,k)
#endif                  
!!$                    sGeo_pres = -V0(HY_GAMC,i,j,k) * V0(HY_PRES,i,j,k) * V0(HY_VELX,i,j,k) * geoFac
                    sGeo_pres = sGeo_dens * cs**2                              !src[PR]
                    sGeo_eint = (sGeo_pres*enth)/V0(HY_DENS,i,j,k)  
#ifdef FLASH_USM_MHD
                    sGeo_magp = - V0(velPhi,i,j,k) * V0(HY_MAGX,i,j,k) * geoFac 
!                    sGeo_magp = 0.0 
!                    sGeo_magz = (V0(HY_MAGX,i,j,k) * V0(HY_VELY,i,j,k) &       !src[Bz]
!                               - V0(HY_VELX,i,j,k) * V0(magZ,i,j,k) ) * geoFac 
                    sGeo_magz =  - V0(HY_VELX,i,j,k) * V0(magZ,i,j,k)  * geoFac 

                    
                    

#endif                  
!!$                 else
!!$                    sGeo_dens = -2.*V0(HY_DENS,i,j,k) * V0(HY_VELX,i,j,k) * Rinv
!!$                    sGeo_velx =  2.*(V0(HY_PRES,i,j,k) + V0(velPhi,i,j,k)**2) * Rinv
!!$                    sGeo_velp = -2.*V0(velPhi,i,j,k) * V0(HY_VELX,i,j,k) * Rinv
!!$                    sGeo_pres = -2.*V0(HY_GAMC,i,j,k) * V0(HY_PRES,i,j,k) * V0(HY_VELX,i,j,k) * Rinv
                 endif

                 !! Add sources terms for n+1/2 Left state
                 Wxn(HY_DENS) = Wxn(HY_DENS) + hdt * sGeo_dens
                 Wxn(HY_VELX) = Wxn(HY_VELX) + hdt * sGeo_velx
                 Wxn(velPhi)  = Wxn(velPhi)  + hdt * sGeo_velp
                 Wxn(HY_PRES) = Wxn(HY_PRES) + hdt * sGeo_pres
                 Wxn(HY_EINT) = Wxn(HY_EINT) + hdt * sGeo_eint
#ifdef FLASH_USM_MHD
                 Wxn(magPhi)  = Wxn(magPhi)  + hdt * sGeo_magp
                 Wxn(magZ)    = Wxn(magZ)    + hdt * sGeo_magz
#endif                 
                 !! Add source terms for n+1/2 Right state
                 Wxp(HY_DENS) = Wxp(HY_DENS) + hdt * sGeo_dens
                 Wxp(HY_VELX) = Wxp(HY_VELX) + hdt * sGeo_velx
                 Wxp(velPhi)  = Wxp(velPhi)  + hdt * sGeo_velp
                 Wxp(HY_PRES) = Wxp(HY_PRES) + hdt * sGeo_pres
                 Wxp(HY_EINT) = Wxp(HY_EINT) + hdt * sGeo_eint
#ifdef FLASH_USM_MHD
                 Wxp(magPhi)  = Wxp(magPhi)  + hdt * sGeo_magp
                 Wxp(magZ)    = Wxp(magZ)    + hdt * sGeo_magz
#endif                 
                
                 if (xCenter(i) - 0.5*del(DIR_X) == 0.) then 
                    ! the velocity should be zero at r=0.
                    Wxn(HY_VELX) = 0.0
                    Wxn(velPhi)  = 0.0
#ifdef FLASH_USM_MHD
                    Wxn(HY_MAGX) = 0.0
                    Wxn(magPhi)  = 0.0
#endif                 
                 elseif (xCenter(i) + 0.5*del(DIR_X) == 0.) then
                    Wxp(HY_VELX) = 0.0
                    Wxp(velPhi)  = 0.0
#ifdef FLASH_USM_MHD
                    Wxp(HY_MAGX) = 0.0
                    Wxp(magPhi)  = 0.0
#endif                 
                 endif
                 !! Calculate R-momentum geometric source term for transverse fluxes
                 !! We will use the cell-centered states at t^n
                 sGeo_trans = (V0(velPhi,i,j,k)**2)*Rinv
#ifdef FLASH_USM_MHD
                 sGeo_trans = sGeo_trans - (V0(magPhi,i,j,k)**2)*Rinv / V0(HY_DENS,i,j,k)
#endif
                 if (hy_geometry == SPHERICAL) &
                      sGeo_trans = sGeo_trans + V0(velTht,i,j,k)**2*Rinv
#if NDIM > 1
                 Wyn(HY_VELX) = Wyn(HY_VELX) + hdt * sGeo_trans
                 Wyp(HY_VELX) = Wyp(HY_VELX) + hdt * sGeo_trans
#if NDIM == 3
                 Wzn(HY_VELX) = Wzn(HY_VELX) + hdt * sGeo_trans
                 Wzp(HY_VELX) = Wzp(HY_VELX) + hdt * sGeo_trans
#endif
#endif


                 !! Check positivity of density and pressure
                 tinyD=0. !hy_eswitch*V0(HY_DENS,i,j,k)
                 tinyP=0. !hy_eswitch*V0(HY_PRES,i,j,k)
                 if (Wxn(HY_DENS) < tinyD .or. Wxp(HY_DENS) < tinyD .or. &
                     Wxn(HY_PRES) < tinyP .or. Wxp(HY_PRES) < tinyP ) then
                    Wxn(HY_DENS:HY_PRES) = V0(HY_DENS:HY_PRES,i,j,k)
                    Wxp(HY_DENS:HY_PRES) = V0(HY_DENS:HY_PRES,i,j,k)
                 endif

#if NDIM >= 2
                 if (Wyn(HY_DENS) < tinyD .or. Wyp(HY_DENS) < tinyD .or. &
                     Wyn(HY_PRES) < tinyP .or. Wyp(HY_PRES) < tinyP ) then
                    Wyn(HY_DENS:HY_PRES) = V0(HY_DENS:HY_PRES,i,j,k)
                    Wyp(HY_DENS:HY_PRES) = V0(HY_DENS:HY_PRES,i,j,k)
                 endif
#if NDIM == 3
                 if (Wzn(HY_DENS) < tinyD .or. Wzp(HY_DENS) < tinyD .or. &
                     Wzn(HY_PRES) < tinyP .or. Wzp(HY_PRES) < tinyP ) then
                    Wzn(HY_DENS:HY_PRES) = V0(HY_DENS:HY_PRES,i,j,k)
                    Wzp(HY_DENS:HY_PRES) = V0(HY_DENS:HY_PRES,i,j,k)
                 endif
#endif
#endif

              endif

!!$!$                 if (hy_geometry == CYLINDRICAL .OR. hy_geometry == POLAR) then
!!$                    sGeo_dens = -V0(HY_DENS,i,j,k) * V0(HY_VELX,i,j,k) * geoFac
!!$                    sGeo_velx = 0. !(V0(HY_PRES,i,j,k) + V0(velPhi,i,j,k)**2) * geoFac
!!$                    sGeo_velp = 0. !-V0(velPhi,i,j,k) * V0(HY_VELX,i,j,k) * geoFac
!!$!$                    sGeo_pres = -V0(HY_GAMC,i,j,k) * V0(HY_PRES,i,j,k) * V0(HY_VELX,i,j,k) * geoFac
!!$                    sGeo_pres = sGeo_dens * cs**2
!!$                    sGeo_eint = (sGeo_pres*enth)/V0(HY_DENS,i,j,k)
!!$!$                 else
!!$!$                    sGeo_dens = -2.*V0(HY_DENS,i,j,k) * V0(HY_VELX,i,j,k) * Rinv
!!$!$                    sGeo_velx =  2.*(V0(HY_PRES,i,j,k) + V0(velPhi,i,j,k)**2) * Rinv
!!$!$                    sGeo_velp = -2.*V0(velPhi,i,j,k) * V0(HY_VELX,i,j,k) * Rinv
!!$!$                    sGeo_pres = -2.*V0(HY_GAMC,i,j,k) * V0(HY_PRES,i,j,k) * V0(HY_VELX,i,j,k) * Rinv
!!$!$                 endif
!!$
!!$                 !! Add sources terms for n+1/2 Left state
!!$                 Wxn(HY_DENS) = Wxn(HY_DENS) + hdt * sGeo_dens
!!$                 Wxn(HY_VELX) = Wxn(HY_VELX) + hdt * sGeo_velx
!!$                 Wxn(velPhi)  = Wxn(velPhi)  + hdt * sGeo_velp
!!$                 Wxn(HY_PRES) = Wxn(HY_PRES) + hdt * sGeo_pres
!!$                 Wxn(HY_EINT) = Wxn(HY_EINT) + hdt * sGeo_eint
!!$                 
!!$                 !! Add source terms for n+1/2 Right state
!!$                 Wxp(HY_DENS) = Wxp(HY_DENS) + hdt * sGeo_dens
!!$                 Wxp(HY_VELX) = Wxp(HY_VELX) + hdt * sGeo_velx
!!$                 Wxp(velPhi)  = Wxp(velPhi)  + hdt * sGeo_velp
!!$                 Wxp(HY_PRES) = Wxp(HY_PRES) + hdt * sGeo_pres
!!$                 Wxp(HY_EINT) = Wxp(HY_EINT) + hdt * sGeo_eint
!!$                
!!$                 if (xCenter(i) - 0.5*del(DIR_X) == 0.) then 
!!$                    ! the velocity should be zero at r=0.
!!$                    Wxn(HY_VELX) = 0.
!!$                 elseif (xCenter(i) + 0.5*del(DIR_X) == 0.) then
!!$                    Wxp(HY_VELX) = 0.
!!$                 endif
!!$                 !! Calculate R-momentum geometric source term for transverse fluxes
!!$                 !! We will use the cell-centered states at t^n
!!$                 sGeo_trans = (V0(velPhi,i,j,k)**2)*Rinv
!!$                 if (hy_geometry == SPHERICAL) &
!!$                      sGeo_trans = sGeo_trans + V0(velTht,i,j,k)**2*Rinv
!!$#if NDIM > 1
!!$                 Wyn(HY_VELX) = Wyn(HY_VELX) + hdt * sGeo_trans
!!$                 Wyp(HY_VELX) = Wyp(HY_VELX) + hdt * sGeo_trans
!!$#if NDIM == 3
!!$                 Wzn(HY_VELX) = Wzn(HY_VELX) + hdt * sGeo_trans
!!$                 Wzp(HY_VELX) = Wzp(HY_VELX) + hdt * sGeo_trans
!!$#endif
!!$#endif
!!$              endif

#ifdef GRAVITY
                 if (hy_useGravity .and. hy_useGravPotUpdate) then
                    ! we don't need the below at all now
                    ! here we can take averages of scratch Density and store them to unk

                    nOneHalfDens = (Wxp(HY_DENS) + Wxn(HY_DENS)) / 2.0
#if NDIM > 1
                    nOneHalfDens = (2.0 * nOneHalfDens + Wyp(HY_DENS) + Wyn(HY_DENS)) / 4.0
#if NDIM > 2
                    nOneHalfDens = (4.0 * nOneHalfDens + Wzp(HY_DENS) + Wzn(HY_DENS)) / 6.0
#endif
#endif
                    ! Now temporarily replace unk density with averaged n+1/2 density
                    U(DENS_VAR,i,j,k) = nOneHalfDens
                 endif
#endif



#ifdef GRAVITY
              if (hy_useGravity .and. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate) then
                 gravXP(i,j,k) = Wxp(HY_GRAV)
                 gravXN(i,j,k) = Wxn(HY_GRAV)
              endif
              if (hy_useGravity .and. .NOT. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate)then
                 Wxp(HY_VELX:HY_VELZ)=Wxp(HY_VELX:HY_VELZ)+hdt*(/Wxp(HY_GRAV),ogravY(i,j,k),ogravZ(i,j,k)/)
                 Wxn(HY_VELX:HY_VELZ)=Wxn(HY_VELX:HY_VELZ)+hdt*(/Wxn(HY_GRAV),ogravY(i,j,k),ogravZ(i,j,k)/)
              endif
#endif  
              !! Store Riemann states to scratch arrays
              scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)= Wxp(HY_DENS:HY_END_VARS-kGrav)
              scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)= Wxn(HY_DENS:HY_END_VARS-kGrav) 

#if NDIM >= 2
#ifdef GRAVITY
              if (hy_useGravity .and. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate) then
                 gravYP(i,j,k) = Wyp(HY_GRAV)
                 gravYN(i,j,k) = Wyn(HY_GRAV)   
              endif
              if (hy_useGravity .and. .NOT. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate)then
                 Wyp(HY_VELX:HY_VELZ)=Wyp(HY_VELX:HY_VELZ)+hdt*(/ogravX(i,j,k),Wyp(HY_GRAV),ogravZ(i,j,k)/)
                 Wyn(HY_VELX:HY_VELZ)=Wyn(HY_VELX:HY_VELZ)+hdt*(/ogravX(i,j,k),Wyn(HY_GRAV),ogravZ(i,j,k)/)
              endif
#endif
              scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=Wyp(HY_DENS:HY_END_VARS-kGrav)
              scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=Wyn(HY_DENS:HY_END_VARS-kGrav)
#if NDIM == 3
#ifdef GRAVITY
              if (hy_useGravity .and. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate) then
                 gravZP(i,j,k) = Wzp(HY_GRAV)
                 gravZN(i,j,k) = Wzn(HY_GRAV)
              endif
              if (hy_useGravity .and. .NOT. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate)then
                 Wzp(HY_VELX:HY_VELZ)=Wzp(HY_VELX:HY_VELZ)+hdt*(/ogravX(i,j,k),ogravY(i,j,k),Wzp(HY_GRAV)/)
                 Wzn(HY_VELX:HY_VELZ)=Wzn(HY_VELX:HY_VELZ)+hdt*(/ogravX(i,j,k),ogravY(i,j,k),Wzn(HY_GRAV)/)
              endif
#endif
              scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=Wzp(HY_DENS:HY_END_VARS-kGrav)
              scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=Wzn(HY_DENS:HY_END_VARS-kGrav)
#endif
#endif


#if NFACE_VARS > 0
#if NDIM > 1
              if (hy_order > 1 .and. hy_killdivB .and. (.not. hy_forceHydroLimit)) then
                 scrch_ctr(XP06_SCRATCH_CENTER_VAR,i,j,k)= Bx(MAG_FACE_VAR,i+1, j,   k  )
                 scrch_Ctr(XN06_SCRATCH_CENTER_VAR,i,j,k)= Bx(MAG_FACE_VAR,i,   j,   k  )

                 scrch_Ctr(YP07_SCRATCH_CENTER_VAR,i,j,k)= By(MAG_FACE_VAR,i,   j+1, k  )
                 scrch_Ctr(YN07_SCRATCH_CENTER_VAR,i,j,k)= By(MAG_FACE_VAR,i,   j,   k  )
#if NDIM == 3
                 scrch_Ctr(ZP08_SCRATCH_CENTER_VAR,i,j,k)= Bz(MAG_FACE_VAR,i,   j,   k+1)
                 scrch_Ctr(ZN08_SCRATCH_CENTER_VAR,i,j,k)= Bz(MAG_FACE_VAR,i,   j,   k  )
#endif
              endif
#endif
#endif

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
           else ! else of if (.not. normalFieldUpdate) then
                !! Updates of magnetic fields in normal direction
#if NFACE_VARS > 0
#if NDIM > 1
              if (hy_order > 1 .and. hy_killdivB .and. (.not. hy_forceHydroLimit)) then
                 scrch_Ctr(XP06_SCRATCH_CENTER_VAR,i,j,k)= Bx(MAGI_FACE_VAR,i+1, j,   k  )
                 scrch_Ctr(XN06_SCRATCH_CENTER_VAR,i,j,k)= Bx(MAGI_FACE_VAR,i,   j,   k  )
                 scrch_Ctr(YP07_SCRATCH_CENTER_VAR,i,j,k)= By(MAGI_FACE_VAR,i,   j+1, k  )
                 scrch_Ctr(YN07_SCRATCH_CENTER_VAR,i,j,k)= By(MAGI_FACE_VAR,i,   j,   k  )
#if NDIM == 3
                 scrch_Ctr(ZP08_SCRATCH_CENTER_VAR,i,j,k)= Bz(MAGI_FACE_VAR,i,   j,   k+1)
                 scrch_Ctr(ZN08_SCRATCH_CENTER_VAR,i,j,k)= Bz(MAGI_FACE_VAR,i,   j,   k  )
#endif
              endif
#endif
#endif
           endif ! end of if (.not. normalFieldUpdate) then
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------
        enddo
#if NDIM == 1
        !$omp end do
#endif
     enddo
#if NDIM == 2
     !$omp end do
#endif
  enddo
#if NDIM == 3
  !$omp end do
#endif


  !! Transverse correction terms for 3D
#if NDIM == 3
#ifdef FLASH_USM_MHD
  if (.not. normalFieldUpdate) then
#endif

     ! hy_use3dFullCTU=.true. will provide CFL <= 1; otherwise, CFL<0.5.
     if (hy_use3dFullCTU) then
        ! Set dt,dx,dy,dz factors
        dt2dxdy6=dt*dt/(6.*del(DIR_X)*del(DIR_Y))
        dt2dydz6=dt*dt/(6.*del(DIR_Y)*del(DIR_Z))
        dt2dzdx6=dt*dt/(6.*del(DIR_Z)*del(DIR_X))


        ! Now let's compute cross derivatives for CTU
#if NDIM == 3
        !$omp do schedule(static)
#endif
        do k=k0-2,kmax+2
#if NDIM == 2
           !$omp do schedule(static)
#endif
           do j=j0-2,jmax+2
#if NDIM == 1
              !$omp do schedule(static)
#endif
              do i=i0-2,imax+2

                 !! ============ x-direction ==================================================================
                 ! YZ cross derivatives for X states
                 call  upwindTransverseFlux&
                      (hy_transOrder,sig(:,DIR_Z,i,j-2:j+2,k),lambda(1,DIR_Y,i,j,k),leig(1,1,DIR_Y,i,j,k),&
                      reig(1,1,DIR_Y,i,j,k),TransFluxYZ(:))

                 ! ZY cross derivatives for X states
                 call upwindTransverseFlux&
                      (hy_transOrder,sig(:,DIR_Y,i,j,k-2:k+2),lambda(1,DIR_Z,i,j,k),leig(1,1,DIR_Z,i,j,k),&
                      reig(1,1,DIR_Z,i,j,k),TransFluxZY(:))

#ifdef FLASH_USM_MHD
                 TransFluxYZ(HY_MAGX) = 0.
                 TransFluxZY(HY_MAGX) = 0.
#endif
                 scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      +(TransFluxYZ(HY_DENS:HY_END_VARS-kGrav)+TransFluxZY(HY_DENS:HY_END_VARS-kGrav))*dt2dydz6

                 scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      +(TransFluxYZ(HY_DENS:HY_END_VARS-kGrav)+TransFluxZY(HY_DENS:HY_END_VARS-kGrav))*dt2dydz6

                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN X-DIRECTION
                 IF (scrch_Ctr(XP01_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(XN01_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(XP05_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(XN05_SCRATCH_CENTER_VAR,i,j,k) .le. 0 ) THEN

                    scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         V0(HY_DENS:HY_END_VARS-kGrav,i,j,k)

                    scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         V0(HY_DENS:HY_END_VARS-kGrav,i,j,k)
                 ENDIF

                 !! ============ y-direction ==================================================================
                 ! ZX cross derivatives for Y states
                 call  upwindTransverseFlux&
                      (hy_transOrder,sig(:,DIR_X,i,j,k-2:k+2),lambda(1,DIR_Z,i,j,k),leig(1,1,DIR_Z,i,j,k),&
                      reig(1,1,DIR_Z,i,j,k),TransFluxZX(:))

                 ! XZ cross derivatives for Y states
                 call  upwindTransverseFlux&
                      (hy_transOrder,sig(:,DIR_Z,i-2:i+2,j,k),lambda(1,DIR_X,i,j,k),leig(1,1,DIR_X,i,j,k),&
                      reig(1,1,DIR_X,i,j,k),TransFluxXZ(:))

#ifdef FLASH_USM_MHD
                 TransFluxZX(HY_MAGY) = 0.
                 TransFluxXZ(HY_MAGY) = 0.
#endif
                 scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxZX(HY_DENS:HY_END_VARS-kGrav)+TransFluxXZ(HY_DENS:HY_END_VARS-kGrav))*dt2dzdx6

                 scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxZX(HY_DENS:HY_END_VARS-kGrav)+TransFluxXZ(HY_DENS:HY_END_VARS-kGrav))*dt2dzdx6

                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN Y-DIRECTION
                 IF (scrch_Ctr(YP01_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(YN01_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(YP05_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(YN05_SCRATCH_CENTER_VAR,i,j,k) .le. 0 ) THEN

                    scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         V0(HY_DENS:HY_END_VARS-kGrav,i,j,k)

                    scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         V0(HY_DENS:HY_END_VARS-kGrav,i,j,k)
                 ENDIF


                 !! ============ z-direction ==================================================================
                 ! XY cross derivatives for Z states
                 call  upwindTransverseFlux&
                      (hy_transOrder,sig(:,DIR_Y,i-2:i+2,j,k),lambda(1,DIR_X,i,j,k),leig(1,1,DIR_X,i,j,k),&
                      reig(1,1,DIR_X,i,j,k),TransFluxXY(:))

                 ! YX cross derivatives for Z states
                 call  upwindTransverseFlux&
                      (hy_transOrder,sig(:,DIR_X,i,j-2:j+2,k),lambda(1,DIR_Y,i,j,k),leig(1,1,DIR_Y,i,j,k),&
                      reig(1,1,DIR_Y,i,j,k),TransFluxYX(:))

#ifdef FLASH_USM_MHD
                 TransFluxXY(HY_MAGZ) = 0.
                 TransFluxYX(HY_MAGZ) = 0.
#endif

                 scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxXY(HY_DENS:HY_END_VARS-kGrav)+TransFluxYX(HY_DENS:HY_END_VARS-kGrav))*dt2dxdy6

                 scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                      scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                      + (TransFluxXY(HY_DENS:HY_END_VARS-kGrav)+TransFluxYX(HY_DENS:HY_END_VARS-kGrav))*dt2dxdy6


                 !! CHECK FOR NEGATIVITY OF DENSITY AND PRESSURE IN Z-DIRECTION
                 IF (scrch_Ctr(ZP01_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(ZN01_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(ZP05_SCRATCH_CENTER_VAR,i,j,k) .le. 0. .or. &
                     scrch_Ctr(ZN05_SCRATCH_CENTER_VAR,i,j,k) .le. 0 ) THEN

                    scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         V0(HY_DENS:HY_END_VARS-kGrav,i,j,k)

                    scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         V0(HY_DENS:HY_END_VARS-kGrav,i,j,k)
                 ENDIF
                 

              enddo
#if NDIM == 1
        !$omp end do
#endif
           enddo
#if NDIM == 2
     !$omp end do
#endif
        enddo
#if NDIM == 3
  !$omp end do
#endif
        !CD: Implicit barrier to avoid a race condition on the release of U.

     end if !end of if (hy_use3dFullCTU) then

#ifdef FLASH_USM_MHD
  endif ! End of if (.not. normalFieldUpdate) then
#endif

#endif 
! end of #if NDIM == 3

  !$omp single
  call Grid_releaseBlkPtr(blockID,U,CENTER) ! moved from above - SMC
  !$omp end single nowait

  !! Add gravity source terms to "conservative" variables of Riemann states if gravity is used.
  !! A primitive variable method is done in the other routine, hy_uhd_addGravityUnsplit.F90.
  !! Note that if hy_useGravPotUpdate=.true. is used, then this conservative update is not
  !! available and adding gravitational source term should be used in primitive update.
#ifdef GRAVITY
  if (hy_useGravity .and. hy_gravConsv .and. hy_useGravHalfUpdate .and. .not. hy_useGravPotUpdate) then

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
     if (.not. normalFieldUpdate) then
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

#if NDIM == 3
        !$omp do schedule(static)
#endif
        do k=k0-2,kmax+2
#if NDIM == 2
           !$omp do schedule(static)
#endif
           do j=j0-2,jmax+2
#if NDIM == 1
              !$omp do schedule(static)
#endif
           do i=i0-2,imax+2

              !! Add interpolated gravity source terms to the *normal components* values at i+1/2 and i-1/2
              !! interfaces in a conservative way. Cell-centered, extrapolated gravity source terms are
              !! added to transverse components.
              !! Calling this conservative update may slow down the code due to conversions 
              !! between primitive and conservative variables.

              !! Xp --------------------------------------------
              call hy_uhd_prim2con(scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   consVar(HY_DENS:HY_VARINUM))

              ! (1) Need to update total energy before momenta update
              consVar(HY_ENER) = consVar(HY_ENER) + &
                   hdt*dot_product(consVar(HY_XMOM:HY_ZMOM),(/gravXP(i,j,k),hgravY(i,j,k),hgravZ(i,j,k)/))

              ! (2) Momenta update
              consVar(HY_XMOM)=consVar(HY_XMOM)+hdt*consVar(HY_DENS)*gravXP(i,j,k)
              consVar(HY_YMOM)=consVar(HY_YMOM)+hdt*consVar(HY_DENS)*hgravY(i,j,k)
              consVar(HY_ZMOM)=consVar(HY_ZMOM)+hdt*consVar(HY_DENS)*hgravZ(i,j,k)

              call hy_uhd_con2prim(consVar(HY_DENS:HY_VARINUM),&
                                   scrch_Ctr(XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k))

              !! Xn --------------------------------------------
              call hy_uhd_prim2con(scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   consVar(HY_DENS:HY_VARINUM))

              ! (1) Need to update total energy before momenta update
              consVar(HY_ENER) = consVar(HY_ENER) + &
                   hdt*dot_product(consVar(HY_XMOM:HY_ZMOM),(/gravXN(i,j,k),hgravY(i,j,k),hgravZ(i,j,k)/))

              ! (2) Momenta update
              consVar(HY_XMOM)=consVar(HY_XMOM)+hdt*consVar(HY_DENS)*gravXN(i,j,k)
              consVar(HY_YMOM)=consVar(HY_YMOM)+hdt*consVar(HY_DENS)*hgravY(i,j,k)
              consVar(HY_ZMOM)=consVar(HY_ZMOM)+hdt*consVar(HY_DENS)*hgravZ(i,j,k)

              call hy_uhd_con2prim(consVar(HY_DENS:HY_VARINUM),&
                                   scrch_Ctr(XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k))

#if NDIM >= 2

              !! Yp --------------------------------------------
              call hy_uhd_prim2con(scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   consVar(HY_DENS:HY_VARINUM))

              ! (1) Need to update total energy before momenta update
              consVar(HY_ENER) = consVar(HY_ENER) + &
                   hdt*dot_product(consVar(HY_XMOM:HY_ZMOM),(/hgravX(i,j,k),gravYP(i,j,k),hgravZ(i,j,k)/))

              ! (2) Momenta update
              consVar(HY_XMOM)=consVar(HY_XMOM)+hdt*consVar(HY_DENS)*hgravX(i,j,k)
              consVar(HY_YMOM)=consVar(HY_YMOM)+hdt*consVar(HY_DENS)*gravYP(i,j,k)
              consVar(HY_ZMOM)=consVar(HY_ZMOM)+hdt*consVar(HY_DENS)*hgravZ(i,j,k)

              call hy_uhd_con2prim(consVar(HY_DENS:HY_VARINUM),&
                                   scrch_Ctr(YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k))

              !! Yn --------------------------------------------
              call hy_uhd_prim2con(scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   consVar(HY_DENS:HY_VARINUM))

              ! (1) Need to update total energy before momenta update
              consVar(HY_ENER) = consVar(HY_ENER) + &
                   hdt*dot_product(consVar(HY_XMOM:HY_ZMOM),(/hgravX(i,j,k),gravYN(i,j,k),hgravZ(i,j,k)/))

              ! (2) Momenta update
              consVar(HY_XMOM)=consVar(HY_XMOM)+hdt*consVar(HY_DENS)*hgravX(i,j,k)
              consVar(HY_YMOM)=consVar(HY_YMOM)+hdt*consVar(HY_DENS)*gravYN(i,j,k)
              consVar(HY_ZMOM)=consVar(HY_ZMOM)+hdt*consVar(HY_DENS)*hgravZ(i,j,k)

              call hy_uhd_con2prim(consVar(HY_DENS:HY_VARINUM),&
                                   scrch_Ctr(YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k))

#if NDIM == 3

              !! Zp --------------------------------------------
              call hy_uhd_prim2con(scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   consVar(HY_DENS:HY_VARINUM))

              ! (1) Need to update total energy before momenta update
              consVar(HY_ENER) = consVar(HY_ENER) + &
                   hdt*dot_product(consVar(HY_XMOM:HY_ZMOM),(/hgravX(i,j,k),hgravY(i,j,k),gravZP(i,j,k)/))

              ! (2) Momenta update
              consVar(HY_XMOM)=consVar(HY_XMOM)+hdt*consVar(HY_DENS)*hgravX(i,j,k)
              consVar(HY_YMOM)=consVar(HY_YMOM)+hdt*consVar(HY_DENS)*hgravY(i,j,k)
              consVar(HY_ZMOM)=consVar(HY_ZMOM)+hdt*consVar(HY_DENS)*gravZP(i,j,k)

              call hy_uhd_con2prim(consVar(HY_DENS:HY_VARINUM),&
                                   scrch_Ctr(ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k))

              !! Zn --------------------------------------------
              call hy_uhd_prim2con(scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   consVar(HY_DENS:HY_VARINUM))

              ! (1) Need to update total energy before momenta update
              consVar(HY_ENER) = consVar(HY_ENER) + &
                   hdt*dot_product(consVar(HY_XMOM:HY_ZMOM),(/hgravX(i,j,k),hgravY(i,j,k),gravZN(i,j,k)/))

              ! (2) Momenta update
              consVar(HY_XMOM)=consVar(HY_XMOM)+hdt*consVar(HY_DENS)*hgravX(i,j,k)
              consVar(HY_YMOM)=consVar(HY_YMOM)+hdt*consVar(HY_DENS)*hgravY(i,j,k)
              consVar(HY_ZMOM)=consVar(HY_ZMOM)+hdt*consVar(HY_DENS)*gravZN(i,j,k)

              call hy_uhd_con2prim(consVar(HY_DENS:HY_VARINUM),&
                                   scrch_Ctr(ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k),&
                                   scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k))
#endif
#endif

           enddo
#if NDIM == 1
        !$omp end do
#endif
        enddo
#if NDIM == 2
     !$omp end do
#endif
     enddo
#if NDIM == 3
  !$omp end do
#endif

     !CD: Implicit barrier to avoid a race condition on the release of
     !scrch_Ctr, V0, Bx, By, Bz.

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
     endif !end of if (.not. normalFieldUpdate)) then
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

  endif !end of if (hy_useGravity .and. hy_useGravHalfUpdate .and. hy_gravConsv .and. (.not. hy_useGravPotUpdate)) then
#endif


  !! Release pointers
  !$omp single
!!$  call Grid_releaseBlkPtr(blockID,scrch_gpot,SCRATCH_CTR)
#if defined(FLASH_UHD_NEED_SCRATCHVARS) || defined(FLASH_USM_MHD)
  call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
#endif

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD

  if (hy_order > 1) then
#if NFACE_VARS > 0
#if NDIM > 1
     call Grid_releaseBlkPtr(blockID,Bx,FACEX)
     call Grid_releaseBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
     call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)
#endif
#endif
  else

#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

     !! Deallocate array
     deallocate(V0)
     deallocate(sig)
     deallocate(lambda)
     deallocate(leig)
     deallocate(reig)
     deallocate(divU)
     deallocate(soundSpeed)

#ifdef GRAVITY
     deallocate(gravXP)
     deallocate(gravXN)
     if (NDIM > 1) then
        deallocate(gravYP)
        deallocate(gravYN)
        if (NDIM > 2) then
           deallocate(gravZP)
           deallocate(gravZN)
        endif
     endif
#endif

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
  endif
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

  !$omp end single nowait
  !$omp end parallel


End Subroutine hy_uhd_getRiemannState


Subroutine upwindTransverseFlux(order,sig,lambda,leig,reig,TransFlux)

  implicit none

#include "UHD.h"
#include "Flash.h"

  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: order
  real,intent(IN),dimension(HY_WAVENUM) :: lambda
  real,intent(IN),dimension(HY_WAVENUM,HY_VARINUM) :: leig
  real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
  real,intent(IN),dimension(HY_VARINUMMAX,-2:2)   :: sig
  real,intent(OUT), dimension(HY_VARINUMMAX) :: TransFlux
  !!------------------------------------------------------------------------

  integer :: n,hyVarEnd
  real, dimension(HY_VARINUMMAX) :: vec

!!$#ifndef GRAVITY
!!$  HY_END_VARS = HY_EINT
!!$#else
!!$  HY_END_VARS = HY_GRAV
!!$#endif

#ifdef FLASH_USM_MHD
  hyVarEnd = HY_MAGZ
#else
  hyVarEnd = HY_PRES
#endif

  TransFlux = 0.

  do n=1,HY_WAVENUM

     if (lambda(n) < 0.) then

        select case (order)
        case(0)
           vec = 0.
        case(1)
           vec(HY_DENS:HY_END_VARS) = sig(HY_DENS:HY_END_VARS,1)-sig(HY_DENS:HY_END_VARS,0)
        case(3)
           !! This 3rd order breaks divB=0
!!$           vec(HY_DENS:HY_END_VARS) =(   -sig(HY_DENS:HY_END_VARS, 2)&  
!!$                                +6.*sig(HY_DENS:HY_END_VARS, 1)&  
!!$                                -3.*sig(HY_DENS:HY_END_VARS, 0)&  
!!$                                -2.*sig(HY_DENS:HY_END_VARS,-1))/6.
           vec(HY_DENS:HY_END_VARS) = sig(HY_DENS:HY_END_VARS,1)-sig(HY_DENS:HY_END_VARS,0)
        end select

     else
        select case (order)
        case(0)
           vec = 0.
        case(1)
           vec(HY_DENS:HY_END_VARS) = sig(HY_DENS:HY_END_VARS,0)-sig(HY_DENS:HY_END_VARS,-1)
        case(3)
           !! This 3rd order breaks divB=0
!!$           vec(HY_DENS:HY_END_VARS)= ( 2.*sig(HY_DENS:HY_END_VARS, 1)&
!!$                                +3.*sig(HY_DENS:HY_END_VARS, 0)&
!!$                                -6.*sig(HY_DENS:HY_END_VARS,-1)&
!!$                                   +sig(HY_DENS:HY_END_VARS,-2))/6.
           vec(HY_DENS:HY_END_VARS) = sig(HY_DENS:HY_END_VARS,0)-sig(HY_DENS:HY_END_VARS,-1)
        end select
     endif

     vec(HY_DENS:hyVarEnd)  = lambda(n)*reig(HY_DENS:hyVarEnd,n)&
          *dot_product(leig(n,HY_DENS:hyVarEnd),vec(HY_DENS:hyVarEnd))

     TransFlux(HY_DENS:hyVarEnd) = TransFlux(HY_DENS:hyVarEnd) + vec(HY_DENS:hyVarEnd)


  end do !end of do n=1,HY_WAVENUM

  !! Transverse flux for gamc, game
  TransFlux(HY_GAMC:HY_END_VARS) = lambda(HY_ENTROPY)*vec(HY_GAMC:HY_END_VARS)



End Subroutine upwindTransverseFlux
