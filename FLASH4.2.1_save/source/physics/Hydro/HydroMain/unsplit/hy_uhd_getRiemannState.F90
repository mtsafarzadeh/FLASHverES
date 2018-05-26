!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_getRiemannState
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
!!                          real(IN)    :: ogravX(:,:,:),
!!                          real(IN)    :: ogravY(:,:,:),
!!                          real(IN)    :: ogravZ(:,:,:),
!!                          real, pointer, dimension(:,:,:,:) :: scrch_Ctr,
!!                          real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig,
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
!!  ogravX       - gravity component in x-direction at n step
!!  ogravY       - gravity component in y-direction at n step
!!  ogravZ       - gravity component in z-direction at n step
!!  scrch_Ctr         - Pointer to the scrch array
!!  hy_SpcR,hy_SpcL,hy_SpcSig - Pointers for Species and mass scalar recon.
!!  normalFieldUpdate - a logical switch to choose normal magnetic fields updates only
!!                      (needed for MHD only)
!!   
!!***

!!REORDER(4):U, V0, scrch_Ctr, scrch_gpot, B[xyz]

Subroutine hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del,&
                                  ogravX,ogravY,ogravZ,&
                                  scrch_Ctr,hy_SpcR,hy_SpcL,hy_SpcSig,&
                                  normalFieldUpdate)
#include "Flash.h"
  use Hydro_data,           ONLY : hy_shockDetectOn,     &
                                   hy_order,             &
                                   hy_flattening,        &
                                   hy_useGravity,        &
                                   hy_useGravHalfUpdate, &
                                   hy_use3dFullCTU,      &
                                   hy_geometry,          &
                                   hy_useHybridOrder,    &
                                   hy_eswitch,           &
                                   hy_upwindTVD

#ifdef FLASH_USM_MHD
  use Hydro_data,           ONLY : hy_killDivB,   &
                                   hy_forceHydroLimit
#endif

  use hy_uhd_slopeLimiters, ONLY : mc
  use hy_uhd_interface,     ONLY : hy_uhd_dataReconstOnestep, &
                                   hy_uhd_shockDetect,        &
                                   hy_uhd_prim2con,           &
                                   hy_uhd_con2prim,           &
                                   hy_uhd_eigenParameters,    &
                                   hy_uhd_eigenValue,         &
                                   hy_uhd_eigenVector,        &
                                   hy_uhd_upwindTransverseFlux
  use Grid_interface,       ONLY : Grid_getBlkPtr,     &
                                   Grid_releaseBlkPtr, &
                                   Grid_getCellCoords

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
       ogravX,ogravY,ogravZ
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: ogravX,ogravY,ogravZ
#endif
  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
  real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
  logical, intent(IN), optional :: normalFieldUpdate
  !! ---------------------------------------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax, i, j, k
  integer,dimension(MDIM) :: dataSize
  real, pointer, dimension(:,:,:,:) :: U, scrch_gpot
  integer :: dir

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
#if NFACE_VARS > 0
#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: Bx,By,Bz
#endif
#endif
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

#ifdef FIXEDBLOCKSIZE
  real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: FlatCoeff,FlatTilde
  real, dimension(     GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: DivU
#else
  real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),&
                       blkLimitsGC(HIGH,JAXIS),&
                       blkLimitsGC(HIGH,KAXIS)) :: FlatCoeff, FlatTilde
  real, dimension(     blkLimitsGC(HIGH,IAXIS),&
                       blkLimitsGC(HIGH,JAXIS),&
                       blkLimitsGC(HIGH,KAXIS)) :: DivU
#endif
  real :: Sp, dv1, dp1, dp2, presL,presR,hdt

#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_IHI_GC) :: xCenter  
#else  
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter  
#endif

  integer :: k2,k3,kGrav,kHydro,order

  ! cylindrical geometry
  integer :: velPhi, velTht, magPhi, magZ
  integer :: HY_velPhi, HY_velTht, H_magPhi, H_magZ

  real :: Rinv, geoFac, eta, enth   
  real :: sGeo_dens, sGeo_velx, sGeo_velp, sGeo_pres
  real :: sGeo_trans, sGeo_eint, sGeo_magz,sGeo_magp

  ! for 3d only --------------
  real, dimension(HY_SPEC_END) :: TransFluxXY,TransFluxYZ,TransFluxZX,&
                                  TransFluxYX,TransFluxZY,TransFluxXZ
  real :: dt2dxdy6,dt2dydz6,dt2dzdx6,hdtdx,hdtdy,hdtdz
  logical :: cons=.false.
  real    :: cs,ca,cf,as,af,uN
  integer :: transOrder3D
  real, dimension(MDIM) :: beta
  logical :: TransX_updateOnly,TransY_updateOnly,TransZ_updateOnly
  ! for 3d only --------------

  integer :: iDim
  real, dimension(HY_VARINUMMAX,NDIM) :: Wp, Wn
  real, dimension(HY_END_VARS) :: Vc
  real, pointer, dimension(:)   :: SigmPtr,SigcPtr,SigpPtr

  real, allocatable,dimension(:,:,:,:,:),target :: sig
  real, allocatable,dimension(:,:,:,:,:) :: lambda
  real, allocatable,dimension(:,:,:,:,:,:) :: leftEig
  real, allocatable,dimension(:,:,:,:,:,:) :: rghtEig





  !! Set transverse flux interpolation order 
  !! (Do not change this!)
  transOrder3D = 1

  !! index for gravity
  kGrav = 0
#ifdef GRAVITY
  kGrav = 1
#endif

  !! index for pure Hydro
  kHydro = 1
#ifdef FLASH_USM_MHD
  kHydro = 0
#endif

  !! indices for various purposes
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

  ! half delta t
  hdt = 0.5*dt

  ! MHD only-------------------------------------------------------------------------------
#if defined(FLASH_USM_MHD) && (NFACE_VARS > 0) && (NDIM > 1)
  if (hy_order > 1) then
     call Grid_getBlkPtr(blockID,Bx,FACEX)
     call Grid_getBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
  endif
#endif /* endif of if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1 */
  ! MHD only-------------------------------------------------------------------------------


  !! Get block pointers for computing & storing Riemann states
  call Grid_getBlkPtr(blockID,U,CENTER)

#ifdef GPRO_VAR
  U(GPRO_VAR,:,:,:) = 0.
#endif


  if (hy_geometry /= CARTESIAN) then
     ! Grab cell x-coords for this block  
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
  endif

  dataSize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
  allocate( sig(HY_VARINUMMAX,NDIM,           dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
  allocate( lambda(HY_WAVENUM,NDIM,           dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
  allocate(leftEig(HY_VARINUM,HY_WAVENUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))
  allocate(rghtEig(HY_VARINUM,HY_WAVENUM,NDIM,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)))


  !! -----------------------------------------------------------------------!
  !! (1) Compute eigen structures once and for all -------------------------!
  !! (2) Compute divergence of velocity fields and (magneto)sonic speed ----!
  !! -----------------------------------------------------------------------!
  if (hy_order > 2) kHydro = 0
  if (hy_useHybridOrder) then
     DivU = 0.
     do k=k0-2-k3+kHydro*k3,kmax+2+k3-kHydro*k3
        do j=j0-2-k3+kHydro*k2,jmax+2+k3-kHydro*k2
           do i=i0-2-k3+kHydro,imax+2+k3-kHydro

              !! Compute undivided divergence of velocity fields and (magneto)sonic speed
              !! and store local (magneto)sonic speeds for hybrid order
              DivU(i,j,k) = U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k)
              if (NDIM > 1) then
                 DivU(i,j,k) = DivU(i,j,k) &
                      +U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k)
                 if (NDIM > 3) then
                    DivU(i,j,k) = DivU(i,j,k) &
                         +U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1)
                 endif
              endif
              DivU(i,j,k) = 0.5*DivU(i,j,k)

           enddo ! do i-loop
        enddo ! do j-loop
     enddo ! do k-loop
  endif !end of hybridOrder
  
  !! -----------------------------------------------------------------------!
  !! Save old velocities ---------------------------------------------------!
  !! -----------------------------------------------------------------------!     
#ifdef VOLX_VAR
  U(VOLX_VAR,:,:,:) = U(VELX_VAR,:,:,:)
#endif
#ifdef VOLY_VAR
  U(VOLY_VAR,:,:,:) = U(VELY_VAR,:,:,:)
#endif
#ifdef VOLZ_VAR
  U(VOLZ_VAR,:,:,:) = U(VELZ_VAR,:,:,:)
#endif


  if (hy_order > 2) kHydro=0
  !! -----------------------------------------------------------------------!
  !! (3) Flattening begins here --------------------------------------------!
  !! -----------------------------------------------------------------------!
  if (hy_flattening) then

     ! Initialize with zero
     FlatTilde = 0.
     FlatCoeff = 0.

     ! Flat tilde
     do k=k0-2+kHydro*k3,kmax+2-kHydro*k3
        do j=j0-2+kHydro*k2,jmax+2-kHydro*k2
           do i=i0-2+kHydro,imax+2-kHydro
              do dir=1,NDIM

                 select case (dir)
                 case (DIR_X)
                    dp1   = (U(PRES_VAR,i+1,j,k)-U(PRES_VAR,i-1,j,k))
                    dp2   = (U(PRES_VAR,i+2,j,k)-U(PRES_VAR,i-2,j,k))
                    dv1   =  U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k)
                    presL =  U(PRES_VAR,i-1,j,k)
                    presR =  U(PRES_VAR,i+1,j,k)
#if NDIM > 1
                 case (DIR_Y)
                    dp1   = (U(PRES_VAR,i,j+1,k)-U(PRES_VAR,i,j-1,k))
                    dp2   = (U(PRES_VAR,i,j+2,k)-U(PRES_VAR,i,j-2,k))
                    dv1   =  U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k)
                    presL =  U(PRES_VAR,i,j-1,k)
                    presR =  U(PRES_VAR,i,j+1,k)
#if NDIM > 2
                 case (DIR_Z)
                    dp1   = (U(PRES_VAR,i,j,k+1)-U(PRES_VAR,i,j,k-1))
                    dp2   = (U(PRES_VAR,i,j,k+2)-U(PRES_VAR,i,j,k-2))
                    dv1   =  U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1)
                    presL =  U(PRES_VAR,i,j,k-1)
                    presR =  U(PRES_VAR,i,j,k+1)
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
        enddo
     enddo

     ! Flat coefficient
     do k=k0-2+kHydro*k3,kmax+2-kHydro*k3
        do j=j0-2+kHydro*k2,jmax+2-kHydro*k2
           do i=i0-2+kHydro,imax+2-kHydro
              do dir=1,NDIM

                 select case (dir)
                 case (DIR_X)
                    dp1   = (U(PRES_VAR,i+1,j,k)-U(PRES_VAR,i-1,j,k))

                    if ( dp1 < 0.0 ) then
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i+1,j,k))
                    elseif (dp1 == 0.) then
                       FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                    else
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i-1,j,k))
                    endif
#if NDIM > 1
                 case (DIR_Y)
                    dp1   = (U(PRES_VAR,i,j+1,k)-U(PRES_VAR,i,j-1,k))

                    if ( dp1 < 0.0 ) then
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j+1,k))
                    elseif (dp1 == 0.) then
                       FlatCoeff(dir,i,j,k) = FlatTilde(dir,i,j,k)
                    else
                       FlatCoeff(dir,i,j,k) = max(FlatTilde(dir,i,j,k),FlatTilde(dir,i,j-1,k))
                    endif
#if NDIM > 2
                 case (DIR_Z)
                    dp1   = (U(PRES_VAR,i,j,k+1)-U(PRES_VAR,i,j,k-1))

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
        enddo
     enddo
  endif

  !! -----------------------------------------------------------------------!
  !! (4) Start calculating Riemann states ----------------------------------!
  !! -----------------------------------------------------------------------!
  !! Compute Riemann states at each cell
  do k=k0-2-k3+kHydro*k3,kmax+2+k3-kHydro*k3
     do j=j0-2-k3+kHydro*k2,jmax+2+k3-kHydro*k2
        do i=i0-2-k3+kHydro,imax+2+k3-kHydro
           ! Extra stencil is needed for 3D to correctly calculate transverse fluxes 
           !(i.e., cross derivatives in x,y, & z)

! MHD only-------------------------------------------------------------------------------
#ifdef FLASH_USM_MHD
           if (.not. normalFieldUpdate) then
#endif
! end of ifdef FLASH_USM_MHD
! MHD only-------------------------------------------------------------------------------

              !! save the cell center values for later use
              Vc(HY_DENS:HY_END_VARS-kGrav) = &
                      (/U(DENS_VAR,i,j,k)&
                       ,U(VELX_VAR:VELZ_VAR,i,j,k)&
                       ,U(PRES_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) 
                       ,U(MAGX_VAR:MAGZ_VAR,i,j,k)&
#endif
                       ,U(GAMC_VAR,i,j,k) &
                       ,U(GAME_VAR,i,j,k) &
                       ,U(EINT_VAR,i,j,k) &
#ifdef FLASH_UHD_3T
                       ,U(EELE_VAR,i,j,k) &
                       ,U(EION_VAR,i,j,k) &
                       ,U(ERAD_VAR,i,j,k) &
#endif
                       /)

              order = hy_order
#ifdef BDRY_VAR
              !! Reduce order in fluid cell near solid boundary if defined:
              !! Reduce order of spatial reconstruction depending on the distance to the solid boundary
              if (order > 2) then
                 if (maxval(U(BDRY_VAR,i-2:i+2,j-2*k2:j+2*k2,k-2*k3:k+2*k3)) < 0.) then !everyone is fluid
                    order = 3
                 else
                    order = 2
                 endif
              endif
              if ((U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i-1, j,   k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i+1, j,   k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j-k2,k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j+k2,k   ) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j,   k-k3) < 0.0) .or. &
                  (U(BDRY_VAR,i,j,k)*U(BDRY_VAR,i,   j,   k+k3) < 0.0)) then
                 order = 1
              endif
#endif

              if (order == 1) then
                 !! NOTE: The first-order scheme (or Donor-cell method) does not include
                 !! transverse fluxes, and therefore the CFL stability limit should satisfy the
                 !! upper limit, 1/NDIM. 
                 do iDim = 1,NDIM
                    Wn(HY_DENS:HY_END_VARS-kGrav,iDim) = Vc(HY_DENS:HY_END_VARS-kGrav)
                    Wp(HY_DENS:HY_END_VARS-kGrav,iDim) = Vc(HY_DENS:HY_END_VARS-kGrav)
                 enddo

              else ! for high-order schemes

                 !! Flag for tranverse update 
                 TransX_updateOnly = .false.
                 TransY_updateOnly = .false.
                 TransZ_updateOnly = .false.

                 if (i+2 > blkLimitsGC(HIGH,IAXIS) .or. &
                     i-2 < blkLimitsGC(LOW, IAXIS)) then
                    TransX_updateOnly = .true.
                 endif
#if NDIM > 1
                 if (j+2 > blkLimitsGC(HIGH,JAXIS) .or. &
                     j-2 < blkLimitsGC(LOW, JAXIS)) then
                    TransY_updateOnly = .true.
                 endif

#if NDIM > 2
                 if (k+2 > blkLimitsGC(HIGH,KAXIS) .or. & 
                     k-2 < blkLimitsGC(LOW, KAXIS)) then
                    TransZ_updateOnly = .true.
                 endif
#endif
#endif

                 !! Left and right Riemann state reconstructions
                 call hy_uhd_dataReconstOnestep&
                      (blockID,blkLimitsGC,    &
                       order,i,j,k,dt,del,     &
                       ogravX,ogravY,ogravZ,   &
                       DivU,FlatCoeff,         &
                       TransX_updateOnly,      &
                       TransY_updateOnly,      &
                       TransZ_updateOnly,      &
                       Wp,Wn,                  &
                       sig    (  1,1,i,j,k),   &
                       lambda (  1,1,i,j,k),   &
                       leftEig(1,1,1,i,j,k),   &
                       rghtEig(1,1,1,i,j,k),   &
                       hy_SpcR,hy_SpcL,hy_SpcSig)

              endif ! end of high-order reconstruction schemes

              if (hy_geometry /= CARTESIAN) then
                 !! **************************************************************
                 !! Add geometric source terms in left and Right States          *
                 !! **************************************************************
                 !Initialize geometric source terms
                 sGeo_dens = 0.
                 sGeo_velx = 0. 
                 sGeo_velp = 0.
                 sGeo_pres = 0.
                 sGeo_eint = 0.
                 sGeo_magz = 0.
                 sGeo_magp = 0.
                 sGeo_trans= 0.

                 Rinv = 1./xCenter(i)
                 select case (hy_geometry)
                 case (CYLINDRICAL)
                    velPhi    = VELZ_VAR
                    HY_velPhi = HY_VELZ
#if defined(FLASH_USM_MHD) 
                    magPhi    = MAGZ_VAR
                    magZ      = MAGY_VAR
                    H_magPhi  = HY_MAGZ
                    H_magZ    = HY_MAGY
#endif
                    geoFac    = Rinv
                 case (POLAR)
                    velPhi    = VELY_VAR
                    HY_velPhi = HY_VELY
#if defined(FLASH_USM_MHD) 
                    magPhi    = MAGY_VAR
                    magZ      = MAGZ_VAR
                    H_magPhi  = HY_MAGY
                    H_magZ    = HY_MAGZ
#endif
                    geoFac    = Rinv
                 case (SPHERICAL)
                    velPhi    = VELZ_VAR
                    velTht    = VELY_VAR
                    HY_velPhi = HY_VELZ
                    HY_velTht = HY_VELY
                    geoFac    = 2.*Rinv
                 end select

                 cs  = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
                 eta = (abs(U(VELX_VAR,i,j,k)) + cs) * dt/del(DIR_X)
                 eta = (1.-eta) / (cs*dt*abs(geoFac))
                 eta = min(1.,eta)
                 !! comment this line not to use the axis hack
                 !!   geoFac = eta * geoFac
                 !! end of the axis hack
                 enth = (U(EINT_VAR,i,j,k) + U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))/cs**2

                 !! right/left state source terms
                 sGeo_dens = -U(DENS_VAR,i,j,k) * U(VELX_VAR,i,j,k) * geoFac !src[DN]
                 sGeo_velx = (U(velPhi,i,j,k)**2) * geoFac                   !src[VR]  
#if defined(FLASH_USM_MHD) 
                 sGeo_velx = sGeo_velx - (U(magPhi,i,j,k)**2) * geoFac / U(DENS_VAR,i,j,k)
#endif                  
                 sGeo_velp = -U(velPhi,i,j,k) * U(VELX_VAR,i,j,k) * geoFac !src[Vphi]
#if defined(FLASH_USM_MHD) 
                 sGeo_velp = sGeo_velp + U(magPhi,i,j,k) * U(MAGX_VAR,i,j,k) * geoFac / U(DENS_VAR,i,j,k)
#endif                  
                 sGeo_pres = sGeo_dens * cs**2                              !src[PR]

                 sGeo_eint = (sGeo_pres*enth)/U(DENS_VAR,i,j,k)  
#if defined(FLASH_USM_MHD) 
                 sGeo_magp = - U(velPhi,  i,j,k) * U(MAGX_VAR,i,j,k) * geoFac 
                 sGeo_magz = - U(VELX_VAR,i,j,k) * U(magZ,    i,j,k) * geoFac 
#endif


                 !! Add sources terms for n+1/2 Left state
                 Wn(HY_DENS,  DIR_X) = Wn(HY_DENS,  DIR_X) + hdt * sGeo_dens
                 Wn(HY_VELX,  DIR_X) = Wn(HY_VELX,  DIR_X) + hdt * sGeo_velx
                 Wn(HY_velPhi,DIR_X) = Wn(HY_velPhi,DIR_X) + hdt * sGeo_velp
                 Wn(HY_PRES,  DIR_X) = Wn(HY_PRES,  DIR_X) + hdt * sGeo_pres
                 Wn(HY_EINT,  DIR_X) = Wn(HY_EINT,  DIR_X) + hdt * sGeo_eint
#if defined(FLASH_USM_MHD) 
                 Wn(H_magPhi, DIR_X) = Wn(H_magPhi, DIR_X) + hdt * sGeo_magp
                 Wn(H_magZ,   DIR_X) = Wn(H_magZ,   DIR_X) + hdt * sGeo_magz
#endif                 
                 !! Add source terms for n+1/2 Right state
                 Wp(HY_DENS,  DIR_X) = Wp(HY_DENS,  DIR_X) + hdt * sGeo_dens
                 Wp(HY_VELX,  DIR_X) = Wp(HY_VELX,  DIR_X) + hdt * sGeo_velx
                 Wp(HY_velPhi,DIR_X) = Wp(HY_velPhi,DIR_X) + hdt * sGeo_velp
                 Wp(HY_PRES,  DIR_X) = Wp(HY_PRES,  DIR_X) + hdt * sGeo_pres
                 Wp(HY_EINT,  DIR_X) = Wp(HY_EINT,  DIR_X) + hdt * sGeo_eint
#if defined(FLASH_USM_MHD) 
                 Wp(H_magPhi, DIR_X) = Wp(H_magPhi, DIR_X) + hdt * sGeo_magp
                 Wp(H_magZ,   DIR_X) = Wp(H_magZ,   DIR_X) + hdt * sGeo_magz
#endif                 


                 if (xCenter(i) - 0.5*del(DIR_X) == 0.) then 
                    ! the velocity should be zero at r=0.
                    Wn(HY_VELX,  DIR_X) = 0.0
                    Wn(HY_velPhi,DIR_X) = 0.0
#if defined(FLASH_USM_MHD) 
                    Wn(HY_MAGX,  DIR_X) = 0.0
                    Wn(H_magPhi, DIR_X) = 0.0
#endif                 
                 elseif (xCenter(i) + 0.5*del(DIR_X) == 0.) then
                    Wp(HY_VELX,  DIR_X) = 0.0
                    Wp(HY_velPhi,DIR_X) = 0.0
#if defined(FLASH_USM_MHD) 
                    Wp(HY_MAGX,  DIR_X) = 0.0
                    Wp(H_magPhi, DIR_X) = 0.0
#endif                 
                 endif
                 !! Calculate R-momentum geometric source term for transverse fluxes
                 !! We will use the cell-centered states at t^n
                 sGeo_trans = (U(velPhi,i,j,k)**2) * geoFac
#if defined(FLASH_USM_MHD) 
                 sGeo_trans = sGeo_trans - (U(magPhi,i,j,k)**2) * geoFac/U(DENS_VAR,i,j,k)
#endif
                 if (hy_geometry == SPHERICAL) &
                      sGeo_trans = sGeo_trans + U(velTht,i,j,k)**2 * geoFac
#if NDIM > 1
                 Wn(HY_VELX,DIR_Y) = Wn(HY_VELX,DIR_Y) + hdt * sGeo_trans
                 Wp(HY_VELX,DIR_Y) = Wp(HY_VELX,DIR_Y) + hdt * sGeo_trans
#if NDIM == 3
                 Wn(HY_VELX,DIR_Z) = Wn(HY_VELX,DIR_Z) + hdt * sGeo_trans
                 Wp(HY_VELX,DIR_Z) = Wp(HY_VELX,DIR_Z) + hdt * sGeo_trans
#endif
#endif
                 !! Check positivity of density and pressure
                 if (Wn(HY_DENS,DIR_X) < 0. .or. Wp(HY_DENS,DIR_X) < 0. .or. &
                     Wn(HY_PRES,DIR_X) < 0. .or. Wp(HY_PRES,DIR_X) < 0.) then
                    Wn(HY_DENS:HY_END_VARS-kGrav,DIR_X) = Vc(HY_DENS:HY_END_VARS-kGrav)
                    Wp(HY_DENS:HY_END_VARS-kGrav,DIR_X) = Vc(HY_DENS:HY_END_VARS-kGrav)
                 end if
#if NDIM > 1
                 if (Wn(HY_DENS,DIR_Y) < 0. .or. Wp(HY_DENS,DIR_Y) < 0. .or. &
                     Wn(HY_PRES,DIR_Y) < 0. .or. Wp(HY_PRES,DIR_Y) < 0.) then
                    Wn(HY_DENS:HY_END_VARS-kGrav,DIR_Y) = Vc(HY_DENS:HY_END_VARS-kGrav)
                    Wp(HY_DENS:HY_END_VARS-kGrav,DIR_Y) = Vc(HY_DENS:HY_END_VARS-kGrav)
                 end if
#if NDIM ==3
                 if (Wn(HY_DENS,DIR_Z) < 0. .or. Wp(HY_DENS,DIR_Z) < 0. .or. &
                     Wn(HY_PRES,DIR_Z) < 0. .or. Wp(HY_PRES,DIR_Z) < 0.) then
                    Wn(HY_DENS:HY_END_VARS-kGrav,DIR_Z) = Vc(HY_DENS:HY_END_VARS-kGrav)
                    Wp(HY_DENS:HY_END_VARS-kGrav,DIR_Z) = Vc(HY_DENS:HY_END_VARS-kGrav)
                 end if
#endif
#endif
              endif !end if of if (hy_geometry .ne. CARTESIAN)


#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 Wp(HY_VELX:HY_VELZ,DIR_X)=Wp(HY_VELX:HY_VELZ,DIR_X)&
                      +hdt*(/Wp(HY_GRAV,DIR_X),ogravY(i,j,k),ogravZ(i,j,k)/)

                 Wn(HY_VELX:HY_VELZ,DIR_X)=Wn(HY_VELX:HY_VELZ,DIR_X)&
                      +hdt*(/Wn(HY_GRAV,DIR_X),ogravY(i,j,k),ogravZ(i,j,k)/)
              endif
#endif  
              !! Store Riemann states to scratch arrays
              scrch_Ctr(XP01_SCRATCH_CENTER_VAR:XP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wp(HY_DENS:HY_END_VARS-kGrav,DIR_X)

              scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wn(HY_DENS:HY_END_VARS-kGrav,DIR_X)

#if NDIM >= 2
#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 Wp(HY_VELX:HY_VELZ,DIR_Y)=Wp(HY_VELX:HY_VELZ,DIR_Y)&
                      +hdt*(/ogravX(i,j,k),Wp(HY_GRAV,DIR_Y),ogravZ(i,j,k)/)

                 Wn(HY_VELX:HY_VELZ,DIR_Y)=Wn(HY_VELX:HY_VELZ,DIR_Y)&
                      +hdt*(/ogravX(i,j,k),Wn(HY_GRAV,DIR_Y),ogravZ(i,j,k)/)
              endif
#endif
              !! Store Riemann states to scratch arrays
              scrch_Ctr(YP01_SCRATCH_CENTER_VAR:YP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wp(HY_DENS:HY_END_VARS-kGrav,DIR_Y)

              scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wn(HY_DENS:HY_END_VARS-kGrav,DIR_Y)
#if NDIM == 3
#ifdef GRAVITY
              if (hy_useGravity .and. hy_useGravHalfUpdate)then
                 Wp(HY_VELX:HY_VELZ,DIR_Z)=Wp(HY_VELX:HY_VELZ,DIR_Z)&
                      +hdt*(/ogravX(i,j,k),ogravY(i,j,k),Wp(HY_GRAV,DIR_Z)/)

                 Wn(HY_VELX:HY_VELZ,DIR_Z)=Wn(HY_VELX:HY_VELZ,DIR_Z)&
                      +hdt*(/ogravX(i,j,k),ogravY(i,j,k),Wn(HY_GRAV,DIR_Z)/)
              endif
#endif
              !! Store Riemann states to scratch arrays
              scrch_Ctr(ZP01_SCRATCH_CENTER_VAR:ZP01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wp(HY_DENS:HY_END_VARS-kGrav,DIR_Z)

              scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)&
                   = Wn(HY_DENS:HY_END_VARS-kGrav,DIR_Z)
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
                 scrch_Ctr(XP06_SCRATCH_CENTER_VAR,i,j,k) = Bx(MAGI_FACE_VAR,i+1, j,   k  )
                 scrch_Ctr(XN06_SCRATCH_CENTER_VAR,i,j,k) = Bx(MAGI_FACE_VAR,i,   j,   k  )
                 scrch_Ctr(YP07_SCRATCH_CENTER_VAR,i,j,k) = By(MAGI_FACE_VAR,i,   j+1, k  )
                 scrch_Ctr(YN07_SCRATCH_CENTER_VAR,i,j,k) = By(MAGI_FACE_VAR,i,   j,   k  )
#if NDIM == 3
                 scrch_Ctr(ZP08_SCRATCH_CENTER_VAR,i,j,k) = Bz(MAGI_FACE_VAR,i,   j,   k+1)
                 scrch_Ctr(ZN08_SCRATCH_CENTER_VAR,i,j,k) = Bz(MAGI_FACE_VAR,i,   j,   k  )
#endif
              endif
#endif
#endif
           endif ! end of if (.not. normalFieldUpdate) then
#endif /* end of ifdef FLASH_USM_MHD */
! MHD only-------------------------------------------------------------------------------
        enddo
     enddo
  enddo


  !!---------------------------------------------------------------!
  !! (5) Transverse correction terms for 3D -----------------------!
  !!---------------------------------------------------------------!
#if NDIM == 3
!#if defined(FLASH_USM_MHD) 
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
        do k=k0-2,kmax+2
           do j=j0-2,jmax+2
              do i=i0-2,imax+2

                 !! save the cell center values for later use
                 Vc(HY_DENS:HY_END_VARS-kGrav) = &
                      (/U(DENS_VAR,i,j,k)&
                       ,U(VELX_VAR:VELZ_VAR,i,j,k)&
                       ,U(PRES_VAR,i,j,k)&
#if defined(FLASH_USM_MHD) 
                       ,U(MAGX_VAR:MAGZ_VAR,i,j,k)&
#endif
                       ,U(GAMC_VAR:GAME_VAR,i,j,k) &
                       ,U(EINT_VAR,i,j,k) &
#ifdef FLASH_UHD_3T
                       ,U(EELE_VAR,i,j,k) &
                       ,U(EION_VAR,i,j,k) &
                       ,U(ERAD_VAR,i,j,k) &
#endif
                       /)

                 !! ============ X-direction ==================================================================
                 If (i .ge. i0-1 .and. i .le. imax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((j .ge. j0   .and. j .le. jmax  ) .and. (k .ge. k0 .and. k .le. kmax)) then
#endif
                 ! YZ cross dervatives for X states
                 SigmPtr => sig(:,DIR_Z,i,j-1,k)
                 SigcPtr => sig(:,DIR_Z,i,j  ,k)
                 SigpPtr => sig(:,DIR_Z,i,j+1,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_END_VARS,TransFluxYZ(1))

                 ! ZY cross derivatives for X states
                 SigmPtr => sig(:,DIR_Y,i,j,k-1)
                 SigcPtr => sig(:,DIR_Y,i,j,k  )
                 SigpPtr => sig(:,DIR_Y,i,j,k+1)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_END_VARS,TransFluxZY(1))

#if defined(FLASH_USM_MHD) 
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
                         Vc(HY_DENS:HY_END_VARS-kGrav)

                    scrch_Ctr(XN01_SCRATCH_CENTER_VAR:XN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         Vc(HY_DENS:HY_END_VARS-kGrav)
                 ENDIF


#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
                 ! YZ cross dervatives for X states
                 SigmPtr => hy_SpcSig(:,i,j-1,k,DIR_Z)
                 SigcPtr => hy_SpcSig(:,i,j  ,k,DIR_Z)
                 SigpPtr => hy_SpcSig(:,i,j+1,k,DIR_Z)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_NSPEC,TransFluxYZ(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 ! ZY cross derivatives for X states
                 SigmPtr => hy_SpcSig(:,i,j,k-1,DIR_Y)
                 SigcPtr => hy_SpcSig(:,i,j,k  ,DIR_Y)
                 SigpPtr => hy_SpcSig(:,i,j,k+1,DIR_Y)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_NSPEC,TransFluxZY(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 hy_SpcR(1:HY_NSPEC,i,j,k,DIR_X) = hy_SpcR(1:HY_NSPEC,i,j,k,DIR_X)&
                      +(TransFluxYZ(HY_SPEC_BEG:HY_SPEC_END)+TransFluxZY(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6

                 hy_SpcL(1:HY_NSPEC,i,j,k,DIR_X) = hy_SpcL(1:HY_NSPEC,i,j,k,DIR_X)&
                      +(TransFluxYZ(HY_SPEC_BEG:HY_SPEC_END)+TransFluxZY(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6
#endif
#endif
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif

                 !! ============ y-direction ==================================================================
                 If (j .ge. j0-1 .and. j .le. jmax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((i .ge. i0   .and. i .le. imax  ) .and. (k .ge. k0 .and. k .le. kmax)) then
#endif
                 ! ZX cross derivatives for Y states
                 SigmPtr => sig(:,DIR_X,i,j,k-1)
                 SigcPtr => sig(:,DIR_X,i,j,k  )
                 SigpPtr => sig(:,DIR_X,i,j,k+1)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_END_VARS,TransFluxZX(1))

                 ! XZ cross derivatives for Y states
                 SigmPtr => sig(:,DIR_Z,i-1,j,k)
                 SigcPtr => sig(:,DIR_Z,i  ,j,k)
                 SigpPtr => sig(:,DIR_Z,i+1,j,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_END_VARS,TransFluxXZ(1))

#if defined(FLASH_USM_MHD) 
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
                         Vc(HY_DENS:HY_END_VARS-kGrav)

                    scrch_Ctr(YN01_SCRATCH_CENTER_VAR:YN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         Vc(HY_DENS:HY_END_VARS-kGrav)
                 ENDIF


#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
                 !ZX cross dervatives for X states
                 SigmPtr => hy_SpcSig(:,i,j,k-1,DIR_X)
                 SigcPtr => hy_SpcSig(:,i,j,k  ,DIR_X)
                 SigpPtr => hy_SpcSig(:,i,j,k+1,DIR_X)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Z,i,j,k),&
                       leftEig(1,1,DIR_Z,i,j,k),&
                       rghtEig(1,1,DIR_Z,i,j,k),&
                       HY_NSPEC,TransFluxZX(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 ! XZ cross derivatives for Y states
                 SigmPtr => hy_SpcSig(:,i-1,j,k,DIR_Z)
                 SigcPtr => hy_SpcSig(:,i  ,j,k,DIR_Z)
                 SigpPtr => hy_SpcSig(:,i+1,j,k,DIR_Z)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_NSPEC,TransFluxXZ(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Y) = hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Y)&
                      +(TransFluxZX(HY_SPEC_BEG:HY_SPEC_END)+TransFluxXZ(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6

                 hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Y) = hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Y)&
                      +(TransFluxZX(HY_SPEC_BEG:HY_SPEC_END)+TransFluxXZ(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6
#endif
#endif
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif


                 !! ============ z-direction ==================================================================
                 If (k .ge. k0-1 .and. k .le. kmax+1) then
#ifdef FLASH_UHD_HYDRO
                 If ((i .ge. i0   .and. i .le. imax  ) .and. (j .ge. j0 .and. j .le. jmax)) then
#endif
                 ! XY cross derivatives for Z states
                 SigmPtr => sig(:,DIR_Y,i-1,j,k)
                 SigcPtr => sig(:,DIR_Y,i  ,j,k)
                 SigpPtr => sig(:,DIR_Y,i+1,j,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_END_VARS,TransFluxXY(1))

                 ! YX cross derivatives for Z states
                 SigmPtr => sig(:,DIR_X,i,j-1,k)
                 SigcPtr => sig(:,DIR_X,i,j  ,k)
                 SigpPtr => sig(:,DIR_X,i,j+1,k)

                 call  hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_END_VARS,TransFluxYX(1))

#if defined(FLASH_USM_MHD) 
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
                         Vc(HY_DENS:HY_END_VARS-kGrav)

                    scrch_Ctr(ZN01_SCRATCH_CENTER_VAR:ZN01_SCRATCH_CENTER_VAR+HY_SCRATCH_NUM-1,i,j,k)=&
                         Vc(HY_DENS:HY_END_VARS-kGrav)
                 ENDIF



#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
                 ! XY cross dervatives for X states
                 SigmPtr => hy_SpcSig(:,i-1,j,k,DIR_Y)
                 SigcPtr => hy_SpcSig(:,i,  j,k,DIR_Y)
                 SigpPtr => hy_SpcSig(:,i+1,j,k,DIR_Y)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_X,i,j,k),&
                       leftEig(1,1,DIR_X,i,j,k),&
                       rghtEig(1,1,DIR_X,i,j,k),&
                       HY_NSPEC,TransFluxXY(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 ! YX cross derivatives for Y states
                 SigmPtr => hy_SpcSig(:,i,j-1,k,DIR_X)
                 SigcPtr => hy_SpcSig(:,i,j,  k,DIR_X)
                 SigpPtr => hy_SpcSig(:,i,j+1,k,DIR_X)

                 call hy_uhd_upwindTransverseFlux&
                      (dir,transOrder3D,&
                       SigmPtr,SigcPtr,SigpPtr,&
                       lambda   (1,DIR_Y,i,j,k),&
                       leftEig(1,1,DIR_Y,i,j,k),&
                       rghtEig(1,1,DIR_Y,i,j,k),&
                       HY_NSPEC,TransFluxYX(HY_SPEC_BEG),&
                       speciesScalar=.true.)

                 hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Z) = hy_SpcR(1:HY_NSPEC,i,j,k,DIR_Z)&
                      +(TransFluxXY(HY_SPEC_BEG:HY_SPEC_END)+TransFluxYX(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6

                 hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Z) = hy_SpcL(1:HY_NSPEC,i,j,k,DIR_Z)&
                      +(TransFluxXY(HY_SPEC_BEG:HY_SPEC_END)+TransFluxYX(HY_SPEC_BEG:HY_SPEC_END))*dt2dydz6
#endif
#endif
#ifdef FLASH_UHD_HYDRO
                 Endif
#endif
                 Endif

              enddo ! i-loop
           enddo ! j-loop
        enddo ! k-loop
     end if !end of if (hy_use3dFullCTU) then

!#if defined(FLASH_USM_MHD) 
#ifdef FLASH_USM_MHD
  endif ! End of if (.not. normalFieldUpdate) then
#endif

#endif /* end of #if NDIM == 3 */


  !! Release pointers
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  ! MHD only-------------------------------------------------------------------------------
#if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1
  if (hy_order > 1) then
     call Grid_releaseBlkPtr(blockID,Bx,FACEX)
     call Grid_releaseBlkPtr(blockID,By,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
  endif ! if (hy_order > 1) then
#endif /* endif of if defined(FLASH_USM_MHD) && NFACE_VARS > 0 && NDIM > 1 */
  ! MHD only-------------------------------------------------------------------------------


  !! Deallocate arrays
  deallocate(sig)
  deallocate(lambda)
  deallocate(leftEig)
  deallocate(rghtEig)

End Subroutine hy_uhd_getRiemannState
