!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_dataReconstOneStep
!!
!! NAME
!!
!!  hy_uhd_dataReconstOneStep
!!
!! SYNOPSIS
!!
!!  hy_uhd_dataReconstOneStep(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(:,:),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(IN)    :: dt,
!!                            real(IN)    :: del(MDIM),
!!                            real(IN)    :: ogravX(:,:,:),
!!                            real(IN)    :: ogravY(:,:,:),
!!                            real(IN)    :: ogravZ(:,:,:),
!!                            real(IN)    :: divV(:,:,:),
!!                            real(IN)    :: soundSpeed(:,:,:),
!!                            real(IN)    :: V0   (HY_VARINUMMAX),
!!                            real(IN)    :: Vxp  (HY_VARINUMMAX),
!!                            real(IN)    :: Vxn  (HY_VARINUMMAX),
!!                            real(IN)    :: Vyp  (HY_VARINUMMAX),
!!                            real(IN)    :: Vyn  (HY_VARINUMMAX),
!!                            real(IN)    :: Vzp  (HY_VARINUMMAX),
!!                            real(IN)    :: Vzn  (HY_VARINUMMAX),
!!                            real(IN)    :: Vxpp (HY_VARINUMMAX),
!!                            real(IN)    :: Vxnn (HY_VARINUMMAX),
!!                            real(IN)    :: Vypp (HY_VARINUMMAX),
!!                            real(IN)    :: Vynn (HY_VARINUMMAX),
!!                            real(IN)    :: Vzpp (HY_VARINUMMAX),
!!                            real(IN)    :: Vznn (HY_VARINUMMAX),
!!                            real(IN)    :: Vxnnn(HY_VARINUMMAX),
!!                            real(IN)    :: Vynnn(HY_VARINUMMAX),
!!                            real(IN)    :: Vznnn(HY_VARINUMMAX),
!!                            real(IN)    :: FlatCoeff(:,:,:,:),
!!                            logical(IN) :: TransX_updateOnly,
!!                            logical(IN) :: TransY_updateOnly,
!!                            logical(IN) :: TransZ_updateOnly,
!!                            real(OUT)   :: Wxp(HY_VARINUMMAX),
!!                            real(OUT)   :: Wxn(HY_VARINUMMAX),
!!                            real(OUT)   :: Wyp(HY_VARINUMMAX),
!!                            real(OUT)   :: Wyn(HY_VARINUMMAX),
!!                            real(OUT)   :: Wzp(HY_VARINUMMAX),
!!                            real(OUT)   :: Wzn(HY_VARINUMMAX),
!!                            real(OUT)   :: sig(HY_VARINUMMAX,NDIM),
!!                            real(OUT)   :: lambda(HY_WAVENUM,NDIM),
!!                            real(OUT)   :: leig(HY_WAVENUM,HY_VARINUM,NDIM),
!!                            real(OUT)   :: reig(HY_VARINUM,HY_WAVENUM,NDIM))
!!
!! ARGUMENTS
!!
!!  blockID     - local block ID
!!  blkLimitsGC - block limits including guardcells
!!  ix,iy,iz    - local indices
!!  dt          - timestep
!!  del         - deltas in each {x,y,z} direction
!!  ogravX,ogravY,ogravZ - gravity components at n step in each direction
!!  divV        - divergence of velocity fields
!!  soundSpeed  - local (magneto)sonic speed
!!  V0          - array containing primitive variables + gamc + game at cell center
!!  Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Vxpp,Vxnn,Vypp,Vynn,Vzpp,Vznn,Vxnnn,Vynnn,Vznnn - data in neighboring cells
!!  FlatCoeff   - flattening parameters
!!  TransX_updateOnly - a switch for a selective transverse update in x direction
!!  TransY_updateOnly - a switch for a selective transverse update in y direction
!!  TransZ_updateOnly - a switch for a selective transverse update in z direction
!!  Wxp,Wxn,Wyp,Wyn,Wzp,Wzn - left(m) and right(p) states
!!  sig    - a transverse flux term
!!  lambda - eigenvalue
!!  leig   - left eigenvector
!!  reig   - right eigenvector
!!
!!
!! DESCRIPTION
!!
!!  This onestep data reconstruction routine evolves the cell centered 
!!  values by dt/2 time step at each cell interface using characteristic 
!!  tracing method based on 2nd order MUSCL-Hancock, 3rd order PPM,
!!  or 5th order WENO schemes.
!!  For MHD, the cell interface values are reconstructed in multidimensional
!!  fashion with inclusions of MHD multidimensional terms that are
!!  proportional to divB for the unsplit MHD implementation.
!!
!! REFERENCES
!!
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Stone, Gardiner, Teuben, Hawley, Simon, "Athena: A new code for astrophysical MHD"
!!    arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
!!  * Colella and Woodward, 54, 174 (1984), JCP
!!  * Colella, 87, 171-200 (1990), JCP
!!
!!***



!!REORDER(4): B[xyz]


Subroutine hy_uhd_dataReconstOneStep(blockID,blkLimitsGC,ix,iy,iz, &
                                     dt,del,ogravX,ogravY,ogravZ,DivU,soundSpeed,V0,&
                                     Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn,  &
                                     Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn, &
                                     Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn,&
                                     FlatCoeff, &
                                     TransX_updateOnly,&
                                     TransY_updateOnly,&
                                     TransZ_updateOnly,&
                                     Wxp, Wxn, Wyp, Wyn, Wzp, Wzn, &
                                     sig,lambda,leig,reig )


  use Hydro_data,        ONLY : hy_eswitch, hy_killdivb, &
                                hy_forceHydroLimit, hy_order,&
                                hy_cfl, hy_cfl_original,&
                                hy_useHybridOrder,&
                                hy_shockDetectOn,&
                                hy_hybridOrderKappa,&
                                hy_geometry
  use Grid_interface,    ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                                Grid_getCellCoords

  use hy_uhd_interface,  ONLY : hy_uhd_checkRHjumpCond

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"


  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: blockID
  integer,intent(IN),dimension(LOW:HIGH,MDIM):: blkLimitsGC
  integer,intent(IN) :: ix,iy,iz
  real,   intent(IN) :: dt
  real,   intent(IN), dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: ogravX,ogravY,ogravZ
  real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: FlatCoeff
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: DivU,soundSpeed
  real, dimension(GRID_IHI_GC) :: xCenter
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: ogravX,ogravY,ogravZ
  real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: FlatCoeff
  real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: DivU,soundSpeed
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter
#endif
  logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly
  real, intent(INOUT),  dimension(HY_VARINUMMAX) :: V0, Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn, &
                                                  Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn,&
                                                  Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn
  real, intent(OUT), dimension(HY_VARINUMMAX)    :: Wxp, Wxn, Wyp, Wyn, Wzp, Wzn
  real, intent(OUT), dimension(HY_VARINUMMAX,NDIM) :: sig
  real, intent(OUT), dimension(HY_WAVENUM,NDIM) :: lambda
  real, intent(OUT), dimension(HY_WAVENUM,HY_VARINUM,NDIM) :: leig
  real, intent(OUT), dimension(HY_VARINUM,HY_WAVENUM,NDIM) :: reig
  !!------------------------------------------------------------------------

  real :: dx,dy,dz,hdtx,hdty,hdtz,hdt,uB,idx,idy,idz
  real :: dnBx,dnBy,dnBz
  real, dimension(HY_VARINUM,NDIM)  :: aBn
  integer :: order
  real, pointer, dimension(:,:,:,:) :: U
  integer :: k2, k3, k4
  logical :: SWxp,SWxn, SWyp,SWyn,SWzp,SWzn
  real :: epsilon != 1.e-8
  real :: xRight, xLeft
  real :: tinyD,tinyP

#if NDIM > 1
  real, dimension(:,:,:,:), pointer :: Bx, By, Bz
#if NDIM < 3
  !!NAG complains if this isn't here, since Bz is uninitialized.
  nullify(Bz)
#endif
!NDIM < 3 
#endif
!NDIM > 1


!!$#ifndef GRAVITY
!!$  HY_END_VARS = HY_EINT
!!$#else
!!$  HY_END_VARS = HY_GRAV
!!$#endif
!!$
!!$#ifdef FLASH_UHD_3T  
!!$  HY_END_VARS = HY_ERAD
!!$#endif


!!$  SWxp = .false.
!!$  SWyp = .false.
!!$  SWzp = .false.
!!$  SWxn = .false.
!!$  SWyn = .false.
!!$  SWzn = .false.


  ! dimensionality tuning parameters
  k2 = 0
  k3 = 0
  k4 = hy_order - 1

  hdt=.5*dt

  dx  =del(DIR_X)
  hdtx=0.5*dt/dx
  idx = 1./dx
#if NDIM >= 2
  dy  =del(DIR_Y)
  hdty=0.5*dt/dy
  idy = 1./dy
  k2 = 1
#if NDIM == 3
  dz  =del(DIR_Z)
  hdtz=0.5*dt/dz
  idz = 1./dz
  k3 = 1
#endif
#endif


  if (hy_geometry /= CARTESIAN) then
     !get coord info will use this for Areas and Volumes
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
  endif

  !! Array bound check for hybrid order in checking DivU
  if (ix-k4    < blkLimitsGC( LOW,IAXIS) .or. &
      ix+k4    > blkLimitsGC(HIGH,IAXIS) .or. &
      iy-k4*k2 < blkLimitsGC( LOW,JAXIS) .or. &
      iy+k4*k2 > blkLimitsGC(HIGH,JAXIS) .or. &
      iz-k4*k3 < blkLimitsGC( LOW,KAXIS) .or. &
      iz+k4*k3 > blkLimitsGC(HIGH,KAXIS)) then
     k4 = k4-1
  endif

  Wxp = 0.
  Wxn = 0.
  Wyp = 0.
  Wyn = 0.
  Wzp = 0.
  Wzn = 0.

#if NFACE_VARS > 0
#if NDIM >=2
  call Grid_getBlkPtr(blockID,Bx,FACEX)
  call Grid_getBlkPtr(blockID,By,FACEY)
#if NDIM == 3
  call Grid_getBlkPtr(blockID,Bz,FACEZ)
#endif
#endif
#endif

  !dot product of U & B fields
  uB = dot_product(V0(HY_VELX:HY_VELZ),V0(HY_MAGX:HY_MAGZ))

  !Multidimensional MHD terms
  if (.not. hy_forceHydroLimit) then
     aBn(:,DIR_X)= (/0.            , &
           -V0(HY_MAGX)/V0(HY_DENS), &
           -V0(HY_MAGY)/V0(HY_DENS), &
           -V0(HY_MAGZ)/V0(HY_DENS), &
           (V0(HY_GAME)-1.)*uB     , &     
            0.                     , &
           -V0(HY_VELY)            , &
           -V0(HY_VELZ)             /)
#if NDIM >= 2
     aBn(:,DIR_Y)=aBn(:,DIR_X)
     aBn(HY_MAGX,DIR_Y)=-V0(HY_VELX)
     aBn(HY_MAGY,DIR_Y)=0.

#if NDIM == 3
     aBn(:,DIR_Z)=aBn(:,DIR_X)
     aBn(HY_MAGX,DIR_Z)=-V0(HY_VELX)
     aBn(HY_MAGZ,DIR_Z)=0.
#endif
#endif
  else
     aBn = 0.
  endif


  !! Define the derivative of the normal B field
  dnBx = 0.
  dnBy = 0.
  dnBz = 0.

#if NFACE_VARS > 0
#if NDIM > 1
  if (hy_killdivb) then
     dnBx = Bx(MAG_FACE_VAR,ix+1,iy,iz)-Bx(MAG_FACE_VAR,ix,iy,iz)
     if (hy_geometry == CYLINDRICAL) then
       xRight = xCenter(ix) + 0.5*dx
       xLeft  = xCenter(ix) - 0.5*dx
       dnBx = (xRight*Bx(MAG_FACE_VAR,ix+1,iy,iz)-xLeft*Bx(MAG_FACE_VAR,ix,iy,iz))/xCenter(ix)
     endif     
     if (NDIM >= 2) then
        dnBy = By(MAG_FACE_VAR,ix,iy+1,iz)-By(MAG_FACE_VAR,ix,iy,iz)
        if (NDIM == 3) then
           dnBz = Bz(MAG_FACE_VAR,ix,iy,iz+1)-Bz(MAG_FACE_VAR,ix,iy,iz)
        endif
     endif
  else
     dnBx = .5*(Vxp(HY_MAGX)-Vxn(HY_MAGX))
     if (NDIM >= 2) then
        dnBy = .5*(Vyp(HY_MAGY)-Vyn(HY_MAGY))
        if (NDIM == 3) then
           dnBz = .5*(Vzp(HY_MAGZ)-Vzn(HY_MAGZ))
        endif
     endif
  endif
#endif
#endif

  !! *************************************************************
  !! (I) Perform data reconstruction in each normal direction    *
  !! *************************************************************
  !! Initialize data at the most distant neighboring location i+3 and i-3, which
  !! are only available for NGUARD > 4


#if NGUARD > 4
  Vxppp(HY_GRAV) = ogravX(ix+3,iy,iz)
  Vxnnn(HY_GRAV) = ogravX(ix-3,iy,iz)
#endif

  V0  (HY_GRAV) = ogravX(ix,  iy,iz)
  Vxp (HY_GRAV) = ogravX(ix+1,iy,iz)
  Vxn (HY_GRAV) = ogravX(ix-1,iy,iz)
  if (.not. TransX_updateOnly) then
     Vxpp(HY_GRAV) = ogravX(ix+2,iy,iz)
     Vxnn(HY_GRAV) = ogravX(ix-2,iy,iz)
  endif


  !! Reduce order in fluid cell near solid boundary if defined
  order = hy_order

#ifdef BDRY_VAR
  call Grid_getBlkPtr(blockID,U,CENTER)
  ! Reduce order of spatial reconstruction depending on the distance to the solid boundary
  if (order > 2) then
     if (maxval(U(BDRY_VAR,ix-2:ix+2,iy-2*k2:iy+2*k2,iz-2*k3:iz+2*k3)) < 0.) then !everyone is fluid
        order = 3
     else
        order = 2
     endif
  endif
  if ((U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix-1, iy,   iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix+1, iy,   iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy-k2,iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy+k2,iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy,   iz-k3) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy,   iz+k3) < 0.0)) then
     order = 1
  endif
#endif


  call DataReconstructNormalDir(DIR_X,order,dt,dx,V0,Vxp,Vxn,Vxpp,Vxnn,Vxppp,Vxnnn,&
                                dnBx,aBn(:,DIR_X),&
                                FlatCoeff(DIR_X,ix,iy,iz), &
                                TransX_updateOnly,&
                                lambda(1:HY_WAVENUM,DIR_X),&
                                Wxp,Wxn,&
                                sig (1:HY_VARINUMMAX,DIR_X),   &
                                leig(1:HY_WAVENUM,1:HY_VARINUM,DIR_X),&
                                reig(1:HY_VARINUM,1:HY_WAVENUM,DIR_X))
#if NDIM >= 2
  if (NDIM>=2) then
#if NGUARD > 4
     Vyppp(HY_GRAV) = ogravY(ix,iy+3,iz)
     Vynnn(HY_GRAV) = ogravY(ix,iy-3,iz)
#endif
     V0  (HY_GRAV) = ogravY(ix,iy,  iz)
     Vyp (HY_GRAV) = ogravY(ix,iy+1,iz)
     Vyn (HY_GRAV) = ogravY(ix,iy-1,iz)
     if (.not. TransY_updateOnly) then
        Vypp(HY_GRAV) = ogravY(ix,iy+2,iz)
        Vynn(HY_GRAV) = ogravY(ix,iy-2,iz)
     endif


  !! Reduce order in fluid cell near solid boundary if defined
  order = hy_order

#ifdef BDRY_VAR
  ! Reduce order of spatial reconstruction depending on the distance to the solid boundary
  if (order > 2) then
     if (maxval(U(BDRY_VAR,ix-2:ix+2,iy-2*k2:iy+2*k2,iz-2*k3:iz+2*k3)) < 0.) then !everyone is fluid
        order = 3
     else
        order = 2
     endif
  endif
  if ((U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix-1, iy,   iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix+1, iy,   iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy-k2,iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy+k2,iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy,   iz-k3) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy,   iz+k3) < 0.0)) then
     order = 1

     !! Reduce CFL for stability
     !$omp critical (Update_cfl)
     if (hy_cfl > 0.5) hy_cfl = 0.4
     !$omp end critical (Update_cfl)
  endif
#endif

     call DataReconstructNormalDir(DIR_Y,order,dt,dy,V0,Vyp,Vyn,Vypp,Vynn,Vyppp,Vynnn,&
                                   dnBy,aBn(:,DIR_Y),&
                                   FlatCoeff(DIR_Y,ix,iy,iz), &
                                   TransY_updateOnly,&
                                   lambda(1:HY_WAVENUM,DIR_Y),&
                                   Wyp,Wyn,&
                                   sig (1:HY_VARINUMMAX,DIR_Y),   &
                                   leig(1:HY_WAVENUM,1:HY_VARINUM,DIR_Y),&
                                   reig(1:HY_VARINUM,1:HY_WAVENUM,DIR_Y))
#if NDIM == 3
     if (NDIM == 3) then
#if NGUARD > 4
        Vzppp(HY_GRAV) = ogravZ(ix,iy,iz+3)
        Vznnn(HY_GRAV) = ogravZ(ix,iy,iz-3)
#endif
        V0  (HY_GRAV) = ogravZ(ix,iy,iz  )
        Vzp (HY_GRAV) = ogravZ(ix,iy,iz+1)
        Vzn (HY_GRAV) = ogravZ(ix,iy,iz-1)
        if (.not. TransZ_updateOnly) then
           Vzpp(HY_GRAV) = ogravZ(ix,iy,iz+2)
           Vznn(HY_GRAV) = ogravZ(ix,iy,iz-2)
        end if
        !!NAG complains about the array bound for FlatCoeff.
        !!since this is not executed except in 3d this should be ok
        !!for now. --PMR


  !! Reduce order in fluid cell near solid boundary if defined
  order = hy_order

#ifdef BDRY_VAR
  ! Reduce order of spatial reconstruction depending on the distance to the solid boundary
  if (order > 2) then
     if (maxval(U(BDRY_VAR,ix-2:ix+2,iy-2*k2:iy+2*k2,iz-2*k3:iz+2*k3)) < 0.) then !everyone is fluid
        order = 3
     else
        order = 2
     endif
  endif
  if ((U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix-1, iy,   iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix+1, iy,   iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy-k2,iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy+k2,iz   ) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy,   iz-k3) < 0.0) .or. &
      (U(BDRY_VAR,ix,iy,iz)*U(BDRY_VAR,ix,   iy,   iz+k3) < 0.0)) then
     order = 1

     !! Reduce CFL for stability
     !$omp critical (Update_cfl)
     if (hy_cfl > 0.33) hy_cfl = 0.3
     !$omp end critical (Update_cfl)
  endif

  call Grid_releaseBlkPtr(blockID,U,CENTER)

#endif

        call DataReconstructNormalDir(DIR_Z,order,dt,dz,V0,Vzp,Vzn,Vzpp,Vznn,Vzppp,Vznnn,&
                                      dnBz,aBn(:,DIR_Z),&
                                      FlatCoeff(DIR_Z,ix,iy,iz), &
                                      TransZ_updateOnly,&
                                      lambda(1:HY_WAVENUM,DIR_Z),&
                                      Wzp,Wzn,&
                                      sig (1:HY_VARINUMMAX,DIR_Z),   &
                                      leig(1:HY_WAVENUM,1:HY_VARINUM,DIR_Z),&
                                      reig(1:HY_VARINUM,1:HY_WAVENUM,DIR_Z))
    endif
#endif
!NDIM == 3
  endif
#endif
!NDIM >= 2


  !! *************************************************************
  !! (Ia) Avoid any possible negative states:                    *
  !!   - Use first-order Godunov scheme if they occur            *
  !! *************************************************************
  if (hy_useHybridOrder) then
     ! kappa = 0. is default; kappa = 0.4 is what Balsara uses.
     ! kappa = -0.1 for mh gives almost same bw shock profile without hybrid order
     ! kappa = -0.1 for ppm more overshoots than without it; -0.05 is bit more overshooting but almost identical
     !! epsilon = hy_hybridOrderKappa*soundSpeed(ix,iy,iz)
     epsilon = hy_hybridOrderKappa*soundSpeed(ix,iy,iz)

     if (minval(DivU(ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3)) < epsilon ) then
        call hy_uhd_checkRHjumpCond(DIR_X,idx,idy,idz,V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wxp,Wxn,SWxp,SWxn)
        if (NDIM >= 2) then
           call hy_uhd_checkRHjumpCond(DIR_Y,idx,idy,idz,V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wyp,Wyn,SWyp,SWyn)
           if (NDIM == 3) then
              call hy_uhd_checkRHjumpCond(DIR_Z,idx,idy,idz,V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wzp,Wzn,SWzp,SWzn)
           endif
        endif
     endif
  endif



#if NDIM >= 2
  !! *************************************************************
  !! (II) Perform to make contributions from transverse fluxes   *
  !! *************************************************************
  !! This is needed only for 2D & 3D

  !! For 2D
  Wxn(HY_DENS:HY_MAGZ) = Wxn(HY_DENS:HY_MAGZ) - hdty*sig(HY_DENS:HY_MAGZ,DIR_Y)     &
                                              - hdty*aBn(HY_DENS:HY_MAGZ,DIR_Y)*dnBy
  Wxp(HY_DENS:HY_MAGZ) = Wxp(HY_DENS:HY_MAGZ) - hdty*sig(HY_DENS:HY_MAGZ,DIR_Y)     &
                                              - hdty*aBn(HY_DENS:HY_MAGZ,DIR_Y)*dnBy
  ! gamc, game, eint, and gravity
  Wxn(HY_GAMC:HY_END_VARS)   = Wxn(HY_GAMC:HY_END_VARS)   - hdty*sig(HY_GAMC:HY_END_VARS,DIR_Y)
  Wxp(HY_GAMC:HY_END_VARS)   = Wxp(HY_GAMC:HY_END_VARS)   - hdty*sig(HY_GAMC:HY_END_VARS,DIR_Y)


  Wyn(HY_DENS:HY_MAGZ) = Wyn(HY_DENS:HY_MAGZ) - hdtx*sig(HY_DENS:HY_MAGZ,DIR_X)     &
                                              - hdtx*aBn(HY_DENS:HY_MAGZ,DIR_X)*dnBx
  Wyp(HY_DENS:HY_MAGZ) = Wyp(HY_DENS:HY_MAGZ) - hdtx*sig(HY_DENS:HY_MAGZ,DIR_X)     &
                                              - hdtx*aBn(HY_DENS:HY_MAGZ,DIR_X)*dnBx
  ! gamc, game, eint, and gravity
  Wyn(HY_GAMC:HY_END_VARS)   = Wyn(HY_GAMC:HY_END_VARS)   - hdtx*sig(HY_GAMC:HY_END_VARS,DIR_X)
  Wyp(HY_GAMC:HY_END_VARS)   = Wyp(HY_GAMC:HY_END_VARS)   - hdtx*sig(HY_GAMC:HY_END_VARS,DIR_X)

#if NDIM == 3
  Wxn(HY_DENS:HY_MAGZ) = Wxn(HY_DENS:HY_MAGZ) - hdtz*sig(HY_DENS:HY_MAGZ,DIR_Z)     &
                                              - hdtz*aBn(HY_DENS:HY_MAGZ,DIR_Z)*dnBz
  Wxp(HY_DENS:HY_MAGZ) = Wxp(HY_DENS:HY_MAGZ) - hdtz*sig(HY_DENS:HY_MAGZ,DIR_Z)     &
                                              - hdtz*aBn(HY_DENS:HY_MAGZ,DIR_Z)*dnBz
  ! gamc, game, and gravity
  Wxn(HY_GAMC:HY_END_VARS)   = Wxn(HY_GAMC:HY_END_VARS)   - hdtz*sig(HY_GAMC:HY_END_VARS,DIR_Z)
  Wxp(HY_GAMC:HY_END_VARS)   = Wxp(HY_GAMC:HY_END_VARS)   - hdtz*sig(HY_GAMC:HY_END_VARS,DIR_Z)


  Wyn(HY_DENS:HY_MAGZ) = Wyn(HY_DENS:HY_MAGZ) - hdtz*sig(HY_DENS:HY_MAGZ,DIR_Z)     &
                                              - hdtz*aBn(HY_DENS:HY_MAGZ,DIR_Z)*dnBz
  Wyp(HY_DENS:HY_MAGZ) = Wyp(HY_DENS:HY_MAGZ) - hdtz*sig(HY_DENS:HY_MAGZ,DIR_Z)     &
                                              - hdtz*aBn(HY_DENS:HY_MAGZ,DIR_Z)*dnBz
  ! gamc, game, and gravity
  Wyn(HY_GAMC:HY_END_VARS)   = Wyn(HY_GAMC:HY_END_VARS)   - hdtz*sig(HY_GAMC:HY_END_VARS,DIR_Z)
  Wyp(HY_GAMC:HY_END_VARS)   = Wyp(HY_GAMC:HY_END_VARS)   - hdtz*sig(HY_GAMC:HY_END_VARS,DIR_Z)


  Wzn(HY_DENS:HY_MAGZ) = Wzn(HY_DENS:HY_MAGZ) - hdtx*sig(HY_DENS:HY_MAGZ,DIR_X)     &
                                              - hdtx*aBn(HY_DENS:HY_MAGZ,DIR_X)*dnBx &
                                              - hdty*sig(HY_DENS:HY_MAGZ,DIR_Y)     &
                                              - hdty*aBn(HY_DENS:HY_MAGZ,DIR_Y)*dnBy

  Wzp(HY_DENS:HY_MAGZ) = Wzp(HY_DENS:HY_MAGZ) - hdtx*sig(HY_DENS:HY_MAGZ,DIR_X)     &
                                              - hdtx*aBn(HY_DENS:HY_MAGZ,DIR_X)*dnBx &
                                              - hdty*sig(HY_DENS:HY_MAGZ,DIR_Y)     &
                                              - hdty*aBn(HY_DENS:HY_MAGZ,DIR_Y)*dnBy

  ! gamc, game, eint, and gravity
  Wzn(HY_GAMC:HY_END_VARS)   = Wzn(HY_GAMC:HY_END_VARS)   - hdtx*sig(HY_GAMC:HY_END_VARS,DIR_X) &
                                                          - hdty*sig(HY_GAMC:HY_END_VARS,DIR_Y)

  Wzp(HY_GAMC:HY_END_VARS)   = Wzp(HY_GAMC:HY_END_VARS)   - hdtx*sig(HY_GAMC:HY_END_VARS,DIR_X) &
                                                          - hdty*sig(HY_GAMC:HY_END_VARS,DIR_Y)


  ! corrections for 3D MHD source terms
  sig(HY_DENS:HY_MAGZ,DIR_X) = sig(HY_DENS:HY_MAGZ,DIR_X) + aBn(HY_DENS:HY_MAGZ,DIR_X)*dnBx
  sig(HY_DENS:HY_MAGZ,DIR_Y) = sig(HY_DENS:HY_MAGZ,DIR_Y) + aBn(HY_DENS:HY_MAGZ,DIR_Y)*dnBy
  sig(HY_DENS:HY_MAGZ,DIR_Z) = sig(HY_DENS:HY_MAGZ,DIR_Z) + aBn(HY_DENS:HY_MAGZ,DIR_Z)*dnBz
  
#endif
#endif



  !! ***************************************************************
  !! (III) Use the continuous normal fields at each cell interface *
  !! ***************************************************************

!!$#if NFACE_VARS > 0
!!$  ! 2D and 3D MHD runs
!!$  if (hy_killdivb) then
!!$     Wxn(HY_MAGX)   = Bx(MAG_FACE_VAR,ix,   iy,  iz  )
!!$     Wxp(HY_MAGX)   = Bx(MAG_FACE_VAR,ix+1, iy,  iz  )
!!$
!!$     Wyn(HY_MAGY)   = By(MAG_FACE_VAR,ix,  iy,   iz  )
!!$     Wyp(HY_MAGY)   = By(MAG_FACE_VAR,ix,  iy+1, iz  )
!!$
!!$#if NDIM == 3
!!$     Wzn(HY_MAGZ)   = Bz(MAG_FACE_VAR,ix,  iy,   iz  )
!!$     Wzp(HY_MAGZ)   = Bz(MAG_FACE_VAR,ix,  iy,   iz+1)
!!$#endif
!!$
!!$  endif
!!$#else
!!$  ! 1D MHD or pure hydro mode
!!$  Wxn(HY_MAGX) = V0(HY_MAGX)
!!$  Wxp(HY_MAGX) = V0(HY_MAGX)
!!$  Wyn(HY_MAGY) = V0(HY_MAGY)
!!$  Wyp(HY_MAGY) = V0(HY_MAGY)
!!$  Wzn(HY_MAGZ) = V0(HY_MAGZ)
!!$  Wzp(HY_MAGZ) = V0(HY_MAGZ)
!!$#endif

#if NFACE_VARS > 0
  ! 2D and 3D MHD runs
  if (hy_killdivb) then
     Wxn(HY_MAGX)   = Bx(MAG_FACE_VAR,ix,   iy,  iz  )
     Wxp(HY_MAGX)   = Bx(MAG_FACE_VAR,ix+1, iy,  iz  )

     Wyn(HY_MAGY)   = By(MAG_FACE_VAR,ix,  iy,   iz  )
     Wyp(HY_MAGY)   = By(MAG_FACE_VAR,ix,  iy+1, iz  )

#if NDIM == 3
     Wzn(HY_MAGZ)   = Bz(MAG_FACE_VAR,ix,  iy,   iz  )
     Wzp(HY_MAGZ)   = Bz(MAG_FACE_VAR,ix,  iy,   iz+1)
#endif

  else
     Wxn(HY_MAGX)   = V0(HY_MAGX)
     Wxp(HY_MAGX)   = V0(HY_MAGX)

     Wyn(HY_MAGY)   = V0(HY_MAGY)
     Wyp(HY_MAGY)   = V0(HY_MAGY)

#if NDIM == 3
     Wzn(HY_MAGZ)   = V0(HY_MAGZ)
     Wzp(HY_MAGZ)   = V0(HY_MAGZ)
#endif

  endif
#else
  ! 1D MHD or pure hydro mode
  Wxn(HY_MAGX) = V0(HY_MAGX)
  Wxp(HY_MAGX) = V0(HY_MAGX)
  Wyn(HY_MAGY) = V0(HY_MAGY)
  Wyp(HY_MAGY) = V0(HY_MAGY)
  Wzn(HY_MAGZ) = V0(HY_MAGZ)
  Wzp(HY_MAGZ) = V0(HY_MAGZ)
#endif


  !! *************************************************************
  !! (IV) Avoid any possible negative states:                    *
  !!    - Use first-order Godunov scheme if they occur           *
  !! *************************************************************
  if (NDIM > 1) then
     if (hy_useHybridOrder) then
        !Note :epsilon = hy_hybridOrderKappa*soundSpeed(ix,iy,iz)
        if (minval(DivU(ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3)) < epsilon ) then
!!$#ifdef RHCD_VAR
!!$                    U(RHCD_VAR,ix,iy,iz) = 1.
!!$#endif
           call hy_uhd_checkRHjumpCond(DIR_X,idx,idy,idz,V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wxp,Wxn,SWxp,SWxn)
           if (NDIM >= 2) then
              call hy_uhd_checkRHjumpCond(DIR_Y,idx,idy,idz,V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wyp,Wyn,SWyp,SWyn)
              if (NDIM == 3) then
                 call hy_uhd_checkRHjumpCond(DIR_Z,idx,idy,idz,V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wzp,Wzn,SWzp,SWzn)
              endif
           endif

#ifdef RHCD_VAR
           call Grid_getBlkPtr(blockID,U,CENTER)
#endif

           ! set a lower cfl
           if (NDIM == 2) then
              if (SWxp .or. SWxn .or. SWyp .or. SWyn) then
#ifdef RHCD_VAR
                    U(RHCD_VAR,ix,iy,iz) = 1.
#endif
                 !$omp critical (Update_cfl)
                 if (hy_cfl > 0.5) then
                    hy_cfl = 0.4
                 endif
                 !$omp end critical (Update_cfl)
              endif
           elseif (NDIM == 3) then
              if (SWxp .or. SWxn .or. SWyp .or. SWyn .or. SWzp .or. SWzn) then

                 !$omp critical (Update_cfl)
                 if (hy_cfl > 0.33) then
                    hy_cfl = 0.3
#ifdef RHCD_VAR
                    U(RHCD_VAR,ix,iy,iz) = 1.
#endif
                 endif
                 !$omp end critical (Update_cfl)
              endif
           endif ! end if of (NDIM == 2) then

        endif !end if of (minval(DivU(ix-k4:ix+k4,iy-k4*k2:iy+k4*k2,iz-k4*k3:iz+k4*k3)) < 0. ) 
     endif ! end if of if (hy_useHybridOrder) then

  endif ! end if of if (NDIM > 1) then


  !! *************************************************************************
  !! (V) The last check for negative pressure and density in reconstruction  *
  !! *************************************************************************
  tinyD=0. !hy_eswitch*V0(HY_DENS)
  tinyP=0. !hy_eswitch*V0(HY_PRES)

  if (Wxn(HY_DENS) < tinyD .or. Wxp(HY_DENS) < tinyD .or. &
      Wxn(HY_PRES) < tinyP .or. Wxp(HY_PRES) < tinyP ) then
     Wxn(HY_DENS:HY_MAGZ) = V0(HY_DENS:HY_MAGZ)
     Wxp(HY_DENS:HY_MAGZ) = V0(HY_DENS:HY_MAGZ)
   endif
#if NDIM >= 2
  if (Wyn(HY_DENS) < tinyD .or. Wyp(HY_DENS) < tinyD .or. &
      Wyn(HY_PRES) < tinyP .or. Wyp(HY_PRES) < tinyP ) then
     Wyn(HY_DENS:HY_MAGZ) = V0(HY_DENS:HY_MAGZ)
     Wyp(HY_DENS:HY_MAGZ) = V0(HY_DENS:HY_MAGZ)
  endif
#if NDIM == 3
  if (Wzn(HY_DENS) < tinyD .or. Wzp(HY_DENS) < tinyD .or. &
      Wzn(HY_PRES) < tinyP .or. Wzp(HY_PRES) < tinyP ) then
     Wzn(HY_DENS:HY_MAGZ) = V0(HY_DENS:HY_MAGZ)
     Wzp(HY_DENS:HY_MAGZ) = V0(HY_DENS:HY_MAGZ)
  endif
#endif
#endif



#ifdef RHCD_VAR
  call Grid_releaseBlkPtr(blockID,U,CENTER)
#endif

  !! Release pointers for facevars
#if NFACE_VARS > 0
#if NDIM >= 2
  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
  call Grid_releaseBlkPtr(blockID,By,FACEY)
#if NDIM == 3
  call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
#endif
#endif
#endif

End Subroutine hy_uhd_dataReconstOneStep




Subroutine DataReconstructNormalDir&
     (dir,order,dt,delta,Vc,Vp,Vm,Vpp,Vmm,Vppp,Vmmm,        &
      dnBn,aBn,FlatCoeff,TransUpdateOnly,lambda,Wp,Wm,sig,leig,reig)

  use Hydro_data,        ONLY : hy_charLimiting,   &
                                hy_eswitch,        &
                                hy_RiemannSolver,  &
                                hy_tiny,           &
                                hy_entropy,        &
                                hy_flattening,     &
                                hy_transOrder,     &
                                hy_ContactSteepening,&
                                hy_upwindTVD,      &
                                hy_3Torder

  use hy_uhd_interface,     ONLY : hy_uhd_TVDslope,        &
                                   hy_uhd_TVDslopeUpwind,  &
                                   hy_uhd_eigenParameters, &
                                   hy_uhd_eigenValue,      &
                                   hy_uhd_eigenVector

  use hy_uhd_slopeLimiters, ONLY : checkMedian, mc, minmod, vanLeer

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: dir,order
  real,intent(IN) :: dt,delta
  real,intent(IN), dimension(HY_VARINUMMAX):: Vc,Vp,Vm,Vpp,Vmm,Vppp,Vmmm
  real,intent(IN) :: dnBn
  real,intent(IN), dimension(HY_VARINUM)  :: aBn
  real,intent(IN)  :: FlatCoeff
  logical, intent(IN) :: TransUpdateOnly
  real,intent(OUT),dimension(HY_WAVENUM)  :: lambda
  real,intent(OUT),dimension(HY_VARINUMMAX):: Wp,Wm
  real,intent(OUT),dimension(HY_VARINUMMAX):: sig
  real,intent(OUT),dimension(HY_WAVENUM,HY_VARINUM) :: leig
  real,intent(OUT),dimension(HY_VARINUM,HY_WAVENUM) :: reig
  !!------------------------------------------------------------------------
  integer :: n
  real    :: dtn,hdtn,factor
  real    :: constA,constB,constC,constD,lambdaMax,lambdaMin
  real    :: temp1, temp2, temp3, eta_steep, del2rhoR, del2rhoL
  real, dimension(HY_VARINUM)  :: vec,sigL,sigR
  real, dimension(HY_VARINUMMAX) :: vecL,vecR,delW,W6,delbar,delbarP,delbarN
  real, PARAMETER :: eta1=20.E0, eta2=0.05E0,epsln=0.01E0,K0=0.1E0

  ! EIG SYSTEM
  logical :: cons=.false.
  real    :: cs,ca,cf,as,af,uN
  real, dimension(MDIM) :: beta
  real, dimension(HY_WAVENUM)   :: lambdaP,lambdaN,lambdaPP,lambdaNN
  real, dimension(HY_WAVENUM,HY_VARINUM) :: leigP,leigN
  real, dimension(HY_VARINUM,HY_WAVENUM) :: reigP,reigN

  ! WENO parameters
  integer :: NVAR_DO
  real, dimension(HY_VARINUMMAX) ::  Vc0,Vp0,Vpp0,Vm0,Vmm0,Vppp0,Vmmm0,& ! primitive variables
                                     Wc0,Wp0,Wpp0,Wm0,Wmm0,Wppp0,Wmmm0,& ! characteristic variables
                                     vecL_temp, vecR_temp, &
                                     vecL_temp1,vecR_temp1,&
                                     vecL_temp2,vecR_temp2
  real, dimension(3,3) :: coeff1m,coeff1p
  real, dimension(3):: coeff2m,coeff2p
  real, dimension(3):: W5p,W5m,betaWeno,Alpha5p,Alpha5m,omega,omegaBar,V5p,V5m
  real, dimension(6):: Intw
  real :: sumAlpha
  real :: errorCheck,errorCheck1
  real :: Flattening
  real,dimension(HY_VARINUMMAX) :: TransFlux
!!$  integer :: HY_END_VARS
  integer :: nVar

!!$#ifndef GRAVITY
!!$  HY_END_VARS = HY_EINT
!!$#else
!!$  HY_END_VARS = HY_GRAV
!!$#endif
!!$
!!$#ifdef FLASH_UHD_3T  
!!$  HY_END_VARS = HY_ERAD
!!$#endif

  dtn=dt/delta
  hdtn=0.5*dtn

  vecL = 0.
  vecR = 0.
  vec  = 0.
  sigL = 0.
  sigR = 0.
  sig  = 0.


  !! Initialize in x-direction
  !! Cell center value
  Vc0(HY_DENS:HY_END_VARS) = Vc(HY_DENS:HY_END_VARS)

  if (dir==DIR_X) then
     Vc0(HY_MAGX)   = 0.
  elseif (dir==DIR_Y) then
     Vc0(HY_MAGY)   = 0.
  elseif (dir==DIR_Z) then
     Vc0(HY_MAGZ)   = 0.
  endif

  !! Prepare for 3rd order PPM or 5th order WENO
  if (order >= 3) then
     Vp0 (HY_DENS:HY_END_VARS) = Vp (HY_DENS:HY_END_VARS)
     Vpp0(HY_DENS:HY_END_VARS) = Vpp(HY_DENS:HY_END_VARS)
     Vm0 (HY_DENS:HY_END_VARS) = Vm (HY_DENS:HY_END_VARS)
     Vmm0(HY_DENS:HY_END_VARS) = Vmm(HY_DENS:HY_END_VARS)
     if (order == 5) then
        Vmmm0(HY_DENS:HY_END_VARS) = Vmmm(HY_DENS:HY_END_VARS)
        Vppp0(HY_DENS:HY_END_VARS) = Vppp(HY_DENS:HY_END_VARS)
     endif

     if (dir==DIR_X) then
        Vp0 (HY_MAGX) = 0.
        Vpp0(HY_MAGX) = 0.
        Vm0 (HY_MAGX) = 0.
        Vmm0(HY_MAGX) = 0.
        if (order == 5) then
           Vmmm0(HY_MAGX) = 0.
           Vppp0(HY_MAGX) = 0.
        endif
     elseif (dir==DIR_Y) then
        Vp0 (HY_MAGY) = 0.
        Vpp0(HY_MAGY) = 0.
        Vm0 (HY_MAGY) = 0.
        Vmm0(HY_MAGY) = 0.
        if (order == 5) then
           Vmmm0(HY_MAGY) = 0.
           Vmmm0(HY_MAGY) = 0.
        endif
     elseif (dir==DIR_Z) then
        Vp0 (HY_MAGZ) = 0.
        Vpp0(HY_MAGZ) = 0.
        Vm0 (HY_MAGZ) = 0.
        Vmm0(HY_MAGZ) = 0.
        if (order == 5) then
           Vmmm0(HY_MAGZ) = 0.
           Vppp0(HY_MAGZ) = 0.
        endif
     endif

  endif

  !! -------------------------------------------------------------------------------------!
  !! Calculate eigensystem  --------------------------------------------------------------!
  !! -------------------------------------------------------------------------------------!
  call hy_uhd_eigenParameters(Vc(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
  call hy_uhd_eigenValue(lambda,uN,cf,C_alfn=ca,C_slow=cs)
  call hy_uhd_eigenVector(leig,reig,Vc(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)


  IF (.not. TransUpdateOnly) THEN
     if (order > 2) then
        call hy_uhd_eigenParameters(Vp(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
        call hy_uhd_eigenValue(lambdaP,uN,cf,C_alfn=ca,C_slow=cs)
        call hy_uhd_eigenVector(leigP,reigP,Vp(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)

        call hy_uhd_eigenParameters(Vm(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
        call hy_uhd_eigenValue(lambdaN,uN,cf,C_alfn=ca,C_slow=cs)
        call hy_uhd_eigenVector(leigN,reigN,Vm(HY_DENS:HY_GAME),dir,cons,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)

        if (hy_entropy) then
           ! Entropy fix for low density region - this treatment is not needed for the data reconstruction
           ! in general, but users can turn it on when needed.
           call hy_uhd_entropyFix(lambda,lambdaN,lambdaP)
        endif
     endif
  ENDIF !end of IF (.not. TransUpdateOnly) THEN

  IF ((.not. TransUpdateOnly) .and. (order > 1)) THEN
     !! -------------------------------------------------------------------------------------!
     !! (1) Apply TVD slope limiter for normal gradients ------------------------------------!
     !! -------------------------------------------------------------------------------------!
     if (.not. hy_upwindTVD) then
        !! Original
        call hy_uhd_TVDslope&
             (dir,Vmm,Vm,Vc,Vp,Vpp,lambdaN,lambda,lambdaP,leig,delbar)

        if (order >= 3) then
           call hy_uhd_TVDslope&
                (dir,Vm,Vc,Vp,Vpp,Vppp,lambda,lambdaP,lambdaPP,leigP,delbarP)

           call hy_uhd_TVDslope&
                (dir,Vmmm,Vmm,Vm,Vc,Vp,lambdaNN,lambdaN,lambda,leigN,delbarN)
        endif
     else

        !! Upwinded TVD for PPM ---
        if ((lambdaN(HY_FASTRGHT) > 0. .and. lambdaP(HY_FASTLEFT) < 0.) .or. &
            (lambdaN(HY_FASTRGHT) > 0. .and. lambda (HY_FASTLEFT) < 0.) .or. &
            (lambda (HY_FASTRGHT) > 0. .and. lambdaP(HY_FASTLEFT) < 0.)) then
           call hy_uhd_TVDslopeUpwind&
                (dir,Vmm,Vm,Vc,Vp,Vpp,lambdaN,lambda,lambdaP,leig,delbar)
        else
           call hy_uhd_TVDslope&
                (dir,Vmm,Vm,Vc,Vp,Vpp,lambdaN,lambda,lambdaP,leig,delbar)
        endif

        if (order >= 3) then

           call hy_uhd_eigenParameters(Vpp(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
           call hy_uhd_eigenValue(lambdaPP,uN,cf,C_alfn=ca,C_slow=cs)

           call hy_uhd_eigenParameters(Vmm(HY_DENS:HY_GAME),dir,cons,uN,cf,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta)
           call hy_uhd_eigenValue(lambdaNN,uN,cf,C_alfn=ca,C_slow=cs)

           if (hy_entropy) then
              call hy_uhd_entropyFix(lambdaP,lambda,lambdaPP)
              call hy_uhd_entropyFix(lambdaN,lambdaNN,lambda)
           endif

           !! Upwinded TVD for PPM ---
           if ((lambda (HY_FASTRGHT) > 0. .and. lambdaPP(HY_FASTLEFT) < 0.) .or. &
               (lambda (HY_FASTRGHT) > 0. .and. lambdaP (HY_FASTLEFT) < 0.) .or. &
               (lambdaP(HY_FASTRGHT) > 0. .and. lambdaPP(HY_FASTLEFT) < 0.)) then
              call hy_uhd_TVDslopeUpwind&
                   (dir,Vm,Vc,Vp,Vpp,Vppp,lambda,lambdaP,lambdaPP,leigP,delbarP)
           else
              call hy_uhd_TVDslope&
                   (dir,Vm,Vc,Vp,Vpp,Vppp,lambda,lambdaP,lambdaPP,leigP,delbarP)
           endif


           !! Upwinded TVD for PPM ---
           if ((lambdaNN(HY_FASTRGHT) > 0. .and. lambda (HY_FASTLEFT) < 0.) .or. &
               (lambdaNN(HY_FASTRGHT) > 0. .and. lambdaN(HY_FASTLEFT) < 0.) .or. &
               (lambdaN (HY_FASTRGHT) > 0. .and. lambda (HY_FASTLEFT) < 0.)) then
              call hy_uhd_TVDslopeUpwind&
                   (dir,Vmmm,Vmm,Vm,Vc,Vp,lambdaNN,lambdaN,lambda,leigN,delbarN)
           else
              call hy_uhd_TVDslope&
                   (dir,Vmmm,Vmm,Vm,Vc,Vp,lambdaNN,lambdaN,lambda,leigN,delbarN)
           endif
        endif
     endif !end of if (.not. hy_upwindTVD) then


     !! -------------------------------------------------------------------------------------!
     !! (2) Begin high-order polynomial interpolation in normal direction -------------------!
     !! -------------------------------------------------------------------------------------!
     !! Calculate interface values at i+1/2 and i-1/2 using high-order polynomial

     !! First initialize flattening coefficients
     if (hy_flattening) then
        Flattening = FlatCoeff
     else
        Flattening = 0.
     endif


     if (order == 2) then
        ! We advect GAMC, GAME & 3T variables here.
        Wp(HY_GAMC:HY_END_VARS)=Vc(HY_GAMC:HY_END_VARS)&
             +0.5*delbar(HY_GAMC:HY_END_VARS)*(1.-Flattening)

        Wm(HY_GAMC:HY_END_VARS)=Vc(HY_GAMC:HY_END_VARS)&
             -0.5*delbar(HY_GAMC:HY_END_VARS)*(1.-Flattening)

        ! Ensure that the interpolated values lie between the cell-centered values
        do nVar=HY_GAMC,HY_END_VARS
           Wm(nVar) = max(min(Vc(nVar),Wm(nVar)),min(max(Vc(nVar),Wm(nVar)),Wm(nVar)))
           Wp(nVar) = max(min(Vc(nVar),Wp(nVar)),min(max(Vc(nVar),Wp(nVar)),Wp(nVar)))
        enddo

        Wp(HY_GAMC:HY_END_VARS)=Wp(HY_GAMC:HY_END_VARS)&
             -max(lambda(HY_ENTROPY),0.)*hdtn*delbar(HY_GAMC:HY_END_VARS)*(1.-Flattening)

        Wm(HY_GAMC:HY_END_VARS)=Wm(HY_GAMC:HY_END_VARS)&
             -min(lambda(HY_ENTROPY),0.)*hdtn*delbar(HY_GAMC:HY_END_VARS)*(1.-Flattening)

#ifdef GRAVITY
        ! No time advancement for gravity
        Wp(HY_GRAV)=Vc(HY_GRAV)+0.5*delbar(HY_GRAV)*(1.-Flattening)
        Wm(HY_GRAV)=Vc(HY_GRAV)-0.5*delbar(HY_GRAV)*(1.-Flattening)
#endif

!!$        ! Ensure that the interpolated values lie between the cell-centered values
!!$        ! before performing characteristic tracing.
!!$        ! At this point, we only applied slope limiters to 
!!$        ! gamc, game, gravity, and 3T variables. Other varialbes are done later
!!$        ! using characteristic tracing.
!!$        do n=HY_GAMC,HY_END_VARS
!!$           Wm(n) = max(min(Vc(n),Wm(n)),min(max(Vc(n),Wm(n)),Wm(n)))
!!$           Wp(n) = max(min(Vc(n),Wp(n)),min(max(Vc(n),Wp(n)),Wp(n)))
!!$        enddo


     elseif (order >= 3) then
        ! delbar contains characteristic variables that are limited.
        ! Project delbar to primitive variables now.
        if (hy_charLimiting) then
           delbar(HY_DENS:HY_MAGZ)  = reig(HY_DENS:HY_MAGZ,HY_FASTLEFT)*delbar(HY_FASTLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ALFNLEFT)*delbar(HY_ALFNLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_SLOWLEFT)*delbar(HY_SLOWLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ENTROPY )*delbar(HY_ENTROPY )+&
                                      reig(HY_DENS:HY_MAGZ,HY_SLOWRGHT)*delbar(HY_SLOWRGHT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ALFNRGHT)*delbar(HY_ALFNRGHT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_FASTRGHT)*delbar(HY_FASTRGHT)

           delbarP(HY_DENS:HY_MAGZ) = reigP(HY_DENS:HY_MAGZ,HY_FASTLEFT)*delbarP(HY_FASTLEFT)+&
                                      reigP(HY_DENS:HY_MAGZ,HY_ALFNLEFT)*delbarP(HY_ALFNLEFT)+&
                                      reigP(HY_DENS:HY_MAGZ,HY_SLOWLEFT)*delbarP(HY_SLOWLEFT)+&
                                      reigP(HY_DENS:HY_MAGZ,HY_ENTROPY )*delbarP(HY_ENTROPY )+&
                                      reigP(HY_DENS:HY_MAGZ,HY_SLOWRGHT)*delbarP(HY_SLOWRGHT)+&
                                      reigP(HY_DENS:HY_MAGZ,HY_ALFNRGHT)*delbarP(HY_ALFNRGHT)+&
                                      reigP(HY_DENS:HY_MAGZ,HY_FASTRGHT)*delbarP(HY_FASTRGHT)

           delbarN(HY_DENS:HY_MAGZ) = reigN(HY_DENS:HY_MAGZ,HY_FASTLEFT)*delbarN(HY_FASTLEFT)+&
                                      reigN(HY_DENS:HY_MAGZ,HY_ALFNLEFT)*delbarN(HY_ALFNLEFT)+&
                                      reigN(HY_DENS:HY_MAGZ,HY_SLOWLEFT)*delbarN(HY_SLOWLEFT)+&
                                      reigN(HY_DENS:HY_MAGZ,HY_ENTROPY )*delbarN(HY_ENTROPY )+&
                                      reigN(HY_DENS:HY_MAGZ,HY_SLOWRGHT)*delbarN(HY_SLOWRGHT)+&
                                      reigN(HY_DENS:HY_MAGZ,HY_ALFNRGHT)*delbarN(HY_ALFNRGHT)+&
                                      reigN(HY_DENS:HY_MAGZ,HY_FASTRGHT)*delbarN(HY_FASTRGHT)
        endif


        if (order == 3) then
           !! -------------------------------------------------------------------------------!
           !! (2)-a: Begining of polynomial interpolation for PPM ---------------------------!
           !! -------------------------------------------------------------------------------!
           !! Parabolic interpolation at the left and right cell interfaces
           !! Colella-Woodward Eqn 1.9, Sekora-Colella Eqn 7, Stone et al Eqn 46
           vecL(HY_DENS:HY_END_VARS) = 0.5*(Vc0(HY_DENS:HY_END_VARS)+Vm0(HY_DENS:HY_END_VARS)) &
                - (delbar(HY_DENS:HY_END_VARS)-delbarN(HY_DENS:HY_END_VARS))/6.

           vecR(HY_DENS:HY_END_VARS) = 0.5*(Vc0(HY_DENS:HY_END_VARS)+Vp0(HY_DENS:HY_END_VARS)) &
                - (delbarP(HY_DENS:HY_END_VARS)-delbar(HY_DENS:HY_END_VARS))/6.

           !! End of polynomial interpolation for PPM


        elseif (order == 5) then
           !! -------------------------------------------------------------------------------!
           !! (2)-b: Begining of polynomial interpolation for WENO --------------------------!
           !! -------------------------------------------------------------------------------!
           if (hy_charLimiting) then
              ! Initialize characteristic variables W*
              Wc0  = 0.
              Wp0  = 0.
              Wm0  = 0.
              Wpp0 = 0.
              Wmm0 = 0.
              Wppp0= 0.
              Wmmm0= 0.
              NVAR_DO = HY_WAVENUM
           else
              Wc0  = Vc0
              Wp0  = Vp0
              Wm0  = Vm0
              Wpp0 = Vpp0
              Wmm0 = Vmm0
              Wppp0= Vppp0
              Wmmm0= Vmmm0
              NVAR_DO = HY_VARINUM3
           endif


           do n=1,NVAR_DO

              if (hy_charLimiting) then
                 ! Note: It is important to use eigenvectors of the cell Vc at which the reconstruction of
                 !       the two left and right states, Wm and Wp, are considered.
                 Wc0(n)   = dot_product(leig(n,HY_DENS:HY_PRES), Vc  (HY_DENS:HY_PRES))

                 Wp0(n)   = dot_product(leig(n,HY_DENS:HY_PRES), Vp  (HY_DENS:HY_PRES))
                 Wpp0(n)  = dot_product(leig(n,HY_DENS:HY_PRES), Vpp (HY_DENS:HY_PRES))
                 Wppp0(n) = dot_product(leig(n,HY_DENS:HY_PRES), Vppp(HY_DENS:HY_PRES))

                 Wm0(n)   = dot_product(leig(n,HY_DENS:HY_PRES), Vm  (HY_DENS:HY_PRES))
                 Wmm0(n)  = dot_product(leig(n,HY_DENS:HY_PRES), Vmm (HY_DENS:HY_PRES))
                 Wmmm0(n) = dot_product(leig(n,HY_DENS:HY_PRES), Vmmm(HY_DENS:HY_PRES))
              endif

              coeff1p(1,1:3) = (/ 2., -7., 11./) !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
              coeff1p(2,1:3) = (/-1.,  5.,  2./) !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
              coeff1p(3,1:3) = (/ 2.,  5., -1./) !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
              coeff1p        = coeff1p/6.
              coeff2p(1:3)   = (/0.1, 0.6, 0.3/) !=(gamma1,gamma2,gamma3)

              coeff1m(1,1:3) = (/-1.,  5.,  2./) !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
              coeff1m(2,1:3) = (/ 2.,  5., -1./) !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
              coeff1m(3,1:3) = (/ 11.,-7.,  2./) !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
              coeff1m        = coeff1m/6.
              coeff2m(1:3)   = (/0.3, 0.6, 0.1/) !=(gamma1,gamma2,gamma3)

              ! Interpolation stencil for weno
              Intw(1:5) = (/Wmm0(n),Wm0(n),Wc0(n),Wp0(n),Wpp0(n)/)


              !! Calculate interface values at i+1/2
              W5p(1) = dot_product(coeff1p(1,1:3),Intw(1:3))
              W5p(2) = dot_product(coeff1p(2,1:3),Intw(2:4))
              W5p(3) = dot_product(coeff1p(3,1:3),Intw(3:5))

              !! Calculate interface values at i-1/2
              W5m(1) = dot_product(coeff1m(1,1:3),Intw(1:3))
              W5m(2) = dot_product(coeff1m(2,1:3),Intw(2:4))
              W5m(3) = dot_product(coeff1m(3,1:3),Intw(3:5))

              !! Calculate smoothness indicators at i+1/2
              betaWeno(1) = 13./12.*(Intw(1)-2.*Intw(2)+Intw(3))**2 + 0.25*(   Intw(1)-4.*Intw(2)+3.*Intw(3))**2
              betaWeno(2) = 13./12.*(Intw(2)-2.*Intw(3)+Intw(4))**2 + 0.25*(   Intw(2)              -Intw(4))**2
              betaWeno(3) = 13./12.*(Intw(3)-2.*Intw(4)+Intw(5))**2 + 0.25*(3.*Intw(3)-4.*Intw(4)   +Intw(5))**2

              !! Calculate weights at i+1/2
              Alpha5p(1) = coeff2p(1)/(hy_tiny+betaWeno(1))!**2
              Alpha5p(2) = coeff2p(2)/(hy_tiny+betaWeno(2))!**2
              Alpha5p(3) = coeff2p(3)/(hy_tiny+betaWeno(3))!**2

              !! Normalize weights at i+1/2
              sumAlpha = Alpha5p(1)+Alpha5p(2)+Alpha5p(3)
              omega(1) = Alpha5p(1)/sumAlpha
              omega(2) = Alpha5p(2)/sumAlpha
              omega(3) = Alpha5p(3)/sumAlpha

!!$           ! optimize weights
!!$           Alpha5p(1:3) = (Alpha5p(1:3)*(omega(1:3)+coeff2(1:3))**2&
!!$                            -3.*coeff2(1:3)*Alpha5p(1:3)+Alpha5p(1:3)**2)/&
!!$                            (coeff2(1:3)**2+Alpha5p(1:3)*(1.-2.*coeff2(1:3)))
!!$
!!$           ! normalize again
!!$           sumAlpha=Alpha5p(1)+Alpha5p(2)+Alpha5p(3)
!!$           omega(1)=Alpha5p(1)/sumAlpha
!!$           omega(2)=Alpha5p(2)/sumAlpha
!!$           omega(3)=Alpha5p(3)/sumAlpha

              !! Compute interface value at i+1/2
              vecR(n) = dot_product(omega(1:3), W5p(1:3))

              !! Calculate weights at i-1/2
              Alpha5m(1) = coeff2m(1)/(hy_tiny+betaWeno(1))!**2
              Alpha5m(2) = coeff2m(2)/(hy_tiny+betaWeno(2))!**2
              Alpha5m(3) = coeff2m(3)/(hy_tiny+betaWeno(3))!**2

              !! Normalize weights at i-1/2
              sumAlpha = Alpha5m(1)+Alpha5m(2)+Alpha5m(3)
              omega(1) = Alpha5m(1)/sumAlpha
              omega(2) = Alpha5m(2)/sumAlpha
              omega(3) = Alpha5m(3)/sumAlpha

!!$           ! optimize weights
!!$           Alpha5m(1:3) = (Alpha5m(1:3)*(omega(1:3)+coeff2(1:3))**2&
!!$                            -3.*coeff2(1:3)*Alpha5m(1:3)+Alpha5m(1:3)**2)/&
!!$                            (coeff2(1:3)**2+Alpha5m(1:3)*(1.-2.*coeff2(1:3)))
!!$
!!$           ! normalize again
!!$           sumAlpha=Alpha5m(1)+Alpha5m(2)+Alpha5m(3)
!!$           omega(1)=Alpha5p(1)/sumAlpha
!!$           omega(2)=Alpha5p(2)/sumAlpha
!!$           omega(3)=Alpha5p(3)/sumAlpha

              !! Compute interface value at i+1/2
              vecL(n) = dot_product(omega(1:3), W5m(1:3))
           enddo


           if (hy_charLimiting) then
              vecR(HY_DENS:HY_MAGZ) = reig(HY_DENS:HY_MAGZ,HY_FASTLEFT)*vecR(HY_FASTLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ALFNLEFT)*vecR(HY_ALFNLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_SLOWLEFT)*vecR(HY_SLOWLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ENTROPY )*vecR(HY_ENTROPY )+&
                                      reig(HY_DENS:HY_MAGZ,HY_SLOWRGHT)*vecR(HY_SLOWRGHT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ALFNRGHT)*vecR(HY_ALFNRGHT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_FASTRGHT)*vecR(HY_FASTRGHT)

              vecL(HY_DENS:HY_MAGZ) = reig(HY_DENS:HY_MAGZ,HY_FASTLEFT)*vecL(HY_FASTLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ALFNLEFT)*vecL(HY_ALFNLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_SLOWLEFT)*vecL(HY_SLOWLEFT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ENTROPY )*vecL(HY_ENTROPY )+&
                                      reig(HY_DENS:HY_MAGZ,HY_SLOWRGHT)*vecL(HY_SLOWRGHT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_ALFNRGHT)*vecL(HY_ALFNRGHT)+&
                                      reig(HY_DENS:HY_MAGZ,HY_FASTRGHT)*vecL(HY_FASTRGHT)

              ! gamc, game, eint and grav
              vecL(HY_GAMC:HY_END_VARS) = 0.5*(Vc0(HY_GAMC:HY_END_VARS)+Vm0(HY_GAMC:HY_END_VARS)) &
                   - (delbar(HY_GAMC:HY_END_VARS)-delbarN(HY_GAMC:HY_END_VARS))/6.

              vecR(HY_GAMC:HY_END_VARS) = 0.5*(Vc0(HY_GAMC:HY_END_VARS)+Vp0(HY_GAMC:HY_END_VARS)) &
                   - (delbarP(HY_GAMC:HY_END_VARS)-delbar(HY_GAMC:HY_END_VARS))/6.

           endif
           !! End of polynomial interpolation for WENO

        endif !! End of interpolation for PPM and WENO

        !! -------------------------------------------------------------------------------!
        !! (2)-c: Contact steepening for PPM or WENO -------------------------------------!
        !! -------------------------------------------------------------------------------!
        temp1 = Vp0(HY_DENS) - Vm0(HY_DENS)
        if (hy_ContactSteepening .and. abs(temp1) > hy_tiny) then
           ! Eqn 1.17 : Second derivatives
           del2rhoR = (Vpp0(HY_DENS)-2.*Vp0(HY_DENS)+ Vc0(HY_DENS))/(6.*delta*delta)
           del2rhoL = ( Vc0(HY_DENS)-2.*Vm0(HY_DENS)+Vmm0(HY_DENS))/(6.*delta*delta)

           ! Third derivative
           eta_steep = (del2rhoL-del2rhoR)*delta**2/temp1
           if (del2rhoR*del2rhoL >= 0.) then
              eta_steep = 0.
           endif
           if (epsln*min(Vp0(HY_DENS),Vm0(HY_DENS))-abs(Vp0(HY_DENS)-Vm0(HY_DENS)) >= 0.) then
              eta_steep = 0.
           endif

           ! Eqn 1.16
           eta_steep = max(0., min(1., eta1*(eta_steep - eta2)))

           ! Eqn 3.2
           temp2 = abs(Vp0(HY_PRES)-Vm0(HY_PRES))/min(Vp0(HY_PRES),Vm0(HY_PRES))
           temp3 = abs(Vp0(HY_DENS)-Vm0(HY_DENS))/min(Vp0(HY_DENS),Vm0(HY_DENS))

           if (Vc(HY_GAME)*K0*temp3-temp2 < 0.0) then
              eta_steep = 0.
           endif

           ! Eqn 1.15
           vecL(HY_DENS) = vecL(HY_DENS)*(1.-eta_steep) + (Vm0(HY_DENS)+0.5*delbarN(HY_DENS))*eta_steep
           vecR(HY_DENS) = vecR(HY_DENS)*(1.-eta_steep) + (Vp0(HY_DENS)-0.5*delbarP(HY_DENS))*eta_steep
        endif
        !! End of Contact steepening for PPM or WENO


        !! -------------------------------------------------------------------------------!
        !! (2)-d: Flattening for PPM or WENO ---------------------------------------------!
        !! -------------------------------------------------------------------------------!
        vecL(:) = Flattening*Vc(:) + (1.0-Flattening)*vecL(:)
        vecR(:) = Flattening*Vc(:) + (1.0-Flattening)*vecR(:)

        !! -------------------------------------------------------------------------------!
        !! (2)-e: Monotonicity constraint for PPM or WENO --------------------------------!
        !! -------------------------------------------------------------------------------!
!!$     errorCheck = 0.
!!$     do n=1,HY_END_VARS
!!$        vecL_temp(n) = checkMedian(vecL(n),Vc0(n),Vm0(n))
!!$        vecR_temp(n) = checkMedian(vecR(n),Vc0(n),Vp0(n))
!!$
!!$        vecL_temp1(n) = checkMedian(vecL_temp(n),Vc0(n),3.*Vc0(n)-2.*vecR_temp(n))
!!$        vecR_temp1(n) = checkMedian(vecR_temp(n),Vc0(n),3.*Vc0(n)-2.*vecL_temp(n))
!!$
!!$        errorCheck = errorCheck+(vecL(n) - vecL_temp1(n))*(vecR(n)-vecR_temp1(n))**2
!!$     enddo
!!$
!!$     if (errorCheck >= 1.E-12) then
!!$        ! if two states get modified then we check two possibilities:
!!$        ! (1) existence of local extremum (i.e., vecL = vecR = Vc)
!!$        errorCheck1 = 0.
!!$        do n=1,HY_END_VARS
!!$           ! check if vecL = vecR = Vc
!!$           errorCheck1 = errorCheck1 + abs(vecL(n) + vecR(n) - 2.*Vc0(n))
!!$        enddo
!!$        if (errorCheck1 < 1.E-12) then
!!$           ! we know there is a local extremum in this case and bound WENO values
!!$           vecL_temp2(n) = checkMedian(vecL_temp1(n),Vc0(n),3.*Vc0(n)-2.*vecR_temp(n))
!!$           vecR_temp2(n) = checkMedian(vecR_temp1(n),Vc0(n),3.*Vc0(n)-2.*vecL_temp(n))
!!$
!!$           vecL(n) = checkMedian(vecL_temp(n),vecL_temp1(n),vecL_temp2(n))
!!$           vecR(n) = checkMedian(vecR_temp(n),vecR_temp1(n),vecR_temp2(n))
!!$        endif
!!$     endif


        ! Ensure that the interpolated values lie between the cell-centered values
        ! Limit according to Colella-Woodward Eqn 1.10
        do n=1,HY_END_VARS
           vecL(n) = max(min(Vc0(n),Vm0(n)),min(max(Vc0(n),Vm0(n)),vecL(n)))
           vecR(n) = max(min(Vc0(n),Vp0(n)),min(max(Vc0(n),Vp0(n)),vecR(n)))
           if (order == 3) then !Note that WENO5 doesn't need this check
              if ( (vecR(n) - Vc0(n))*(Vc0(n)-vecL(n)) <= 0.) then
                 vecL(n) = Vc0(n)
                 vecR(n) = Vc0(n)
              endif
           endif
           if ( 6.*(vecR(n)-vecL(n))*(Vc0(n)-0.5*(vecL(n)+vecR(n))) > (vecR(n) - vecL(n))**2  ) then
              vecL(n) = 3.*Vc0(n) - 2.*vecR(n)
           endif
           if ( 6.*(vecR(n)-vecL(n))*(Vc0(n)-0.5*(vecL(n)+vecR(n))) < -(vecR(n) - vecL(n))**2  ) then
              vecR(n) = 3.*Vc0(n) - 2.*vecL(n)
           endif
        enddo
        !! End of Contact steepeing, Flattening, and Monotonicity constraint for PPM and WENO


        !! -------------------------------------------------------------------------------!
        !! (2)-f: Take initial guesses for the left and right states ---------------------!
        !! -------------------------------------------------------------------------------!
        !! PPM coefficients for parabolic interpolations
        delW(HY_DENS:HY_END_VARS) = vecR(HY_DENS:HY_END_VARS)-vecL(HY_DENS:HY_END_VARS)
        W6(HY_DENS:HY_END_VARS)   = 6.*(Vc0(HY_DENS:HY_END_VARS)-0.5*(vecR(HY_DENS:HY_END_VARS)+vecL(HY_DENS:HY_END_VARS)))

        !! Right states ===================================================================
        lambdaMax =max(lambda(HY_FASTRGHT),0.)
        Wp(HY_DENS:HY_MAGZ) = vecR(HY_DENS:HY_MAGZ) - lambdaMax*hdtn &
             *(delW(HY_DENS:HY_MAGZ) - (1.0 - lambdaMax*hdtn*4./3.)*W6(HY_DENS:HY_MAGZ) )

        !! gamc, game, eint, and 3T vars (eele, eion, erad)
        lambdaMax = max(lambda(HY_ENTROPY),0.)
        Wp(HY_GAMC:HY_END_VARS) = vecR(HY_GAMC:HY_END_VARS) - lambdaMax*hdtn &
             *(delW(HY_GAMC:HY_END_VARS) - (1.0 - lambdaMax*hdtn*4./3.)*W6(HY_GAMC:HY_END_VARS) )

        !! Left states  ===================================================================
        lambdaMin = -min(lambda(HY_FASTLEFT),0.)
        Wm(HY_DENS:HY_MAGZ) = vecL(HY_DENS:HY_MAGZ) + lambdaMin*hdtn &
             *(delW(HY_DENS:HY_MAGZ) + (1.0 - lambdaMin*hdtn*4./3.)*W6(HY_DENS:HY_MAGZ) )

        !! gamc, game, eint, and 3T vars (eele, eion, erad)
        lambdaMin = -min(lambda(HY_ENTROPY),0.)
        Wm(HY_GAMC:HY_END_VARS) = vecL(HY_GAMC:HY_END_VARS) + lambdaMin*hdtn &
             *(delW(HY_GAMC:HY_END_VARS) + (1.0 - lambdaMin*hdtn*4./3.)*W6(HY_GAMC:HY_END_VARS) )

#ifdef GRAVITY
        !! gravity component includes only spatial interpolation
        Wp(HY_GRAV) = vecR(HY_GRAV)
        Wm(HY_GRAV) = vecL(HY_GRAV)
#endif
        !! End of initial guesses
     endif
     !! End of high-order polynomial interpolations (MUSCL-Hancock, PPM or WENO) for interface values

  ENDIF  !end of IF ((.not. TransUpdateOnly) .and. (order > 1)) THEN

  !! -------------------------------------------------------------------------------!
  !! (3) Advance the above interpolated interface values by 1/2 time step using     !
  !!     characteristic tracing method for MUSCL-Hancock, PPM or WENO               !
  !! -------------------------------------------------------------------------------!


  do n=1,HY_WAVENUM

     IF ((.not. TransUpdateOnly) .and. (order > 1)) THEN

        if (hy_charLimiting) then
           !! Apply slope limiter on characteristic variables

           if (hy_RiemannSolver == ROE) then
              if (lambda(n) < 0.) then
                 if (order == 2) then
                    !! Apply monotone slope limiting for normal flux
                    vecL(HY_DENS:HY_MAGZ) = .5*(-1.-dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)*delbar(n)*(1.-Flattening)
                    sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)
                 elseif (order >= 3) then
                    ! PPM step 10
                    constC = hdtn*(lambda(HY_FASTLEFT) - lambda(n))
                    constD = 1./3. *(dtn)**2 * (lambda(HY_FASTLEFT)**2 - lambda(n)**2)
                    vecL(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                         (constC*(delW(HY_DENS:HY_MAGZ)+W6(HY_DENS:HY_MAGZ))+constD*W6(HY_DENS:HY_MAGZ)))&
                         *reig(HY_DENS:HY_MAGZ,n)
                    sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)
                 endif

              elseif (lambda(n) > 0.) then
                 if (order == 2) then
                    !! Apply monotone slope limiting for normal flux term
                    vecR(HY_DENS:HY_MAGZ) = .5*(1.-dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)*delbar(n)*(1.-Flattening)
                    sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
                 elseif (order >= 3) then
                    constA = hdtn*(lambda(HY_FASTRGHT) - lambda(n))
                    constB = 1./3. *(dtn)**2 * (lambda(HY_FASTRGHT)**2 - lambda(n)**2 )
                    vecR(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                         (constA*(delW(HY_DENS:HY_MAGZ)-W6(HY_DENS:HY_MAGZ))+constB*W6(HY_DENS:HY_MAGZ)))&
                         *reig(HY_DENS:HY_MAGZ,n)
                    sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
                 endif
              endif
           else
              !! Left states for HLL* type solvers
              !! For more detail, see "Athena: A new code for astrophysical MHD"
              !! by Stone, Gardiner, Teuben, Hawley, Simon, arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
              !! Apply monotone slope limiting for normal flux
              if (order == 2) then
                 vecL(HY_DENS:HY_MAGZ) = .5*(-1.-dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)*delbar(n)*(1.-Flattening)
                 sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)

                 !! Right states for HLL* type solvers
                 !! Apply monotone slope limiting for normal flux term
                 vecR(HY_DENS:HY_MAGZ) = .5*(1.-dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)*delbar(n)*(1.-Flattening)
                 sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
              elseif (order >= 3) then
                 ! PPM step 10
                 constC = hdtn*(lambda(HY_FASTLEFT) - lambda(n))
                 constD = 1./3. *(dtn)**2 * (lambda(HY_FASTLEFT)**2 - lambda(n)**2)
                 vecL(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                      (constC*(delW(HY_DENS:HY_MAGZ)+W6(HY_DENS:HY_MAGZ))+constD*W6(HY_DENS:HY_MAGZ)))&
                      *reig(HY_DENS:HY_MAGZ,n)
                 sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)

                 constA = hdtn*(lambda(HY_FASTRGHT) - lambda(n))
                 constB = 1./3. *(dtn)**2 * (lambda(HY_FASTRGHT)**2 - lambda(n)**2 )
                 vecR(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                      (constA*(delW(HY_DENS:HY_MAGZ)-W6(HY_DENS:HY_MAGZ))+constB*W6(HY_DENS:HY_MAGZ)))&
                      *reig(HY_DENS:HY_MAGZ,n)
                 sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
              endif
           endif


        else
           !! Apply slope limiter on primitive variables

           if (hy_RiemannSolver == ROE) then
              if (lambda(n) < 0.) then
                 if (order == 2) then
                    !! Apply monotone slope limiting for normal flux
                    vecL(HY_DENS:HY_MAGZ) = .5*(-1.- dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)&
                         *dot_product(leig(n,HY_DENS:HY_MAGZ), delbar(HY_DENS:HY_MAGZ)*(1.-Flattening))
                    sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)
                 elseif (order >= 3) then
                    ! PPM
                    ! some coefficients
                    constC = hdtn*(lambda(HY_FASTLEFT) - lambda(n))
                    constD = 1./3. *(dtn)**2 * (lambda(HY_FASTLEFT)**2 - lambda(n)**2)
                    vecL(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                         (constC*(delW(HY_DENS:HY_MAGZ)+W6(HY_DENS:HY_MAGZ))+constD*W6(HY_DENS:HY_MAGZ)))&
                         *reig(HY_DENS:HY_MAGZ,n)
                    sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)
                 endif

              elseif (lambda(n) > 0.) then
                 if (order == 2) then
                    !! Apply monotone slope limiting for normal flux
                    vecR(HY_DENS:HY_MAGZ) = .5*(1.-dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)&
                         *dot_product(leig(n,HY_DENS:HY_MAGZ), delbar(HY_DENS:HY_MAGZ)*(1.-Flattening))
                    sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
                 elseif (order >= 3) then
                    ! PPM
                    ! some coefficients
                    constA = hdtn*(lambda(HY_FASTRGHT) - lambda(n))
                    constB = 1./3. *(dtn)**2 * (lambda(HY_FASTRGHT)**2 - lambda(n)**2 )
                    vecR(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                         (constA*(delW(HY_DENS:HY_MAGZ)-W6(HY_DENS:HY_MAGZ))+constB*W6(HY_DENS:HY_MAGZ)))&
                         *reig(HY_DENS:HY_MAGZ,n)
                    sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
                 endif

              endif
           else
              !! Left states for HLL* type solvers
              !! Apply monotone slope limiting for normal flux
              if (order == 2) then
                 vecL(HY_DENS:HY_MAGZ) = .5*(-1.- dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)&
                      *dot_product(leig(n,HY_DENS:HY_MAGZ), delbar(HY_DENS:HY_MAGZ)*(1.-Flattening))
                 sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)

                 !! Right states for HLL* type solvers
                 !! Apply monotone slope limiting for normal flux
                 vecR(HY_DENS:HY_MAGZ) = .5*(1.-dtn*lambda(n))*reig(HY_DENS:HY_MAGZ,n)&
                      *dot_product(leig(n,HY_DENS:HY_MAGZ), delbar(HY_DENS:HY_MAGZ)*(1.-Flattening))
                 sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)

              elseif (order >= 3) then
                 ! PPM step 10
                 constC = hdtn*(lambda(HY_FASTLEFT) - lambda(n))
                 constD = 1./3. *(dtn)**2 * (lambda(HY_FASTLEFT)**2 - lambda(n)**2)
                 vecL(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                      (constC*(delW(HY_DENS:HY_MAGZ)+W6(HY_DENS:HY_MAGZ))+constD*W6(HY_DENS:HY_MAGZ)))&
                      *reig(HY_DENS:HY_MAGZ,n)
                 sigL(HY_DENS:HY_MAGZ) = sigL(HY_DENS:HY_MAGZ) + vecL(HY_DENS:HY_MAGZ)

                 constA = hdtn*(lambda(HY_FASTRGHT) - lambda(n))
                 constB = 1./3. *(dtn)**2 * (lambda(HY_FASTRGHT)**2 - lambda(n)**2 )
                 vecR(HY_DENS:HY_MAGZ) = dot_product(leig(n,HY_DENS:HY_MAGZ),&
                      (constA*(delW(HY_DENS:HY_MAGZ)-W6(HY_DENS:HY_MAGZ))+constB*W6(HY_DENS:HY_MAGZ)))&
                      *reig(HY_DENS:HY_MAGZ,n)
                 sigR(HY_DENS:HY_MAGZ) = sigR(HY_DENS:HY_MAGZ) + vecR(HY_DENS:HY_MAGZ)
              endif

           endif

        endif ! End of if (hy_charLimiting)
     ENDIF ! end of IF ((.not. TransUpdateOnly) .and. (order > 1)) THEN

     !! Transverse flux -------------------------------------------
     !! Apply upwinding differencing for transverse flux
     TransFlux = 0.
     if (lambda(n) < 0.) then
        select case (hy_transOrder)
        case(0)
           TransFlux = 0.
        case(1)
           !! First order upwind for transverse flux derivative works fine but diffusive
           TransFlux(HY_DENS:HY_END_VARS)=Vp(HY_DENS:HY_END_VARS)-Vc(HY_DENS:HY_END_VARS)
        case(3)
           TransFlux(HY_DENS:HY_END_VARS)= ( -Vpp(HY_DENS:HY_END_VARS)&  
                                      +6.*Vp(HY_DENS:HY_END_VARS)&   
                                      -3.*Vc(HY_DENS:HY_END_VARS)&   
                                      -2.*Vm(HY_DENS:HY_END_VARS))/6.
        case(4)
           do nVar = HY_DENS,HY_END_VARS
              TransFlux(nVar)=minmod((-Vpp(nVar)+6.*Vp(nVar)-3.*Vc(nVar)-2.*Vm(nVar))/6.,&
                                       Vp(nVar)-Vc(nVar))
           enddo
        end select
    else
       select case (hy_transOrder)
       case(0)
          TransFlux = 0.
       case(1)
          TransFlux(HY_DENS:HY_END_VARS)=Vc(HY_DENS:HY_END_VARS)-Vm(HY_DENS:HY_END_VARS)
       case(3)
          TransFlux(HY_DENS:HY_END_VARS)= (2.*Vp(HY_DENS:HY_END_VARS)&   
                                    +3.*Vc(HY_DENS:HY_END_VARS)&   
                                    -6.*Vm(HY_DENS:HY_END_VARS)&   
                                      +Vmm(HY_DENS:HY_END_VARS))/6.
       case(4)
           do nVar = HY_DENS,HY_END_VARS
              TransFlux(nVar)=minmod(( 2.*Vp(nVar)+3.*Vc(nVar)-6.*Vm(nVar)+Vmm(nVar))/6.,&
                                       Vc(nVar)-Vm(nVar))
           enddo
        end select
     endif

     ! usual transeverse terms
     vec(HY_DENS:HY_MAGZ)  = lambda(n)*reig(HY_DENS:HY_MAGZ,n)&
          *dot_product(leig(n,HY_DENS:HY_MAGZ),TransFlux(HY_DENS:HY_MAGZ))


     !! Transverse fluxes for density, velocities, pressure, and magnetic fields.
     sig(HY_DENS:HY_MAGZ) = sig(HY_DENS:HY_MAGZ) + vec(HY_DENS:HY_MAGZ)

  end do ! End of do n=1,HY_WAVENUM



  !! Transverse flux for gamc, game
  !! NOTE: Not including transverse flux terms for 3T variables & eint (also for gravity??) works better
  !!       for preserving positivity in HEDP simulations.
  !! Before we had: sig(HY_GAMC:HY_END_VARS) = lambda(HY_ENTROPY)*TransFlux(HY_GAMC:HY_END_VARS)
  sig(HY_GAMC:HY_GAME) = lambda(HY_ENTROPY)*TransFlux(HY_GAMC:HY_GAME)


#ifdef FLASH_UHD_3T
  IF (.not. TransUpdateOnly) THEN
     IF (hy_3Torder .ne. order) THEN

        if (hy_3Torder == 1) then
           Wp(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)
           Wm(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)

        elseif (hy_3Torder == 2) then
           Wp(HY_EINT:HY_ERAD)=Vc(HY_EINT:HY_ERAD)&
                +0.5*delbar(HY_EINT:HY_ERAD)*(1.-Flattening)

           Wm(HY_EINT:HY_ERAD)=Vc(HY_EINT:HY_ERAD)&
                -0.5*delbar(HY_EINT:HY_ERAD)*(1.-Flattening)

           do nVar=HY_EINT,HY_ERAD
              Wm(nVar) = max(min(Vc(nVar),Wm(nVar)),min(max(Vc(nVar),Wm(nVar)),Wm(nVar)))
              Wp(nVar) = max(min(Vc(nVar),Wp(nVar)),min(max(Vc(nVar),Wp(nVar)),Wp(nVar)))
           enddo

           Wp(HY_EINT:HY_ERAD)=Wp(HY_EINT:HY_ERAD)&
                -max(lambda(HY_ENTROPY),0.)*hdtn*delbar(HY_EINT:HY_ERAD)*(1.-Flattening)

           Wm(HY_EINT:HY_ERAD)=Wm(HY_EINT:HY_ERAD)&
                -min(lambda(HY_ENTROPY),0.)*hdtn*delbar(HY_EINT:HY_ERAD)*(1.-Flattening)
        endif
     ENDIF
  ENDIF
#endif


  IF (.not. TransUpdateOnly) THEN
     !! Riemann states in normal direction
     if (order == 1) then
        Wm(HY_DENS:HY_END_VARS) = Vc(HY_DENS:HY_END_VARS)!-hdtn*aBn(HY_DENS:HY_MAGZ)*dnBn
        Wp(HY_DENS:HY_END_VARS) = Vc(HY_DENS:HY_END_VARS)!-hdtn*aBn(HY_DENS:HY_MAGZ)*dnBn
     elseif (order == 2) then
        Wm(HY_DENS:HY_MAGZ) = Vc0(HY_DENS:HY_MAGZ)+sigL(HY_DENS:HY_MAGZ)-hdtn*aBn(HY_DENS:HY_MAGZ)*dnBn
        Wp(HY_DENS:HY_MAGZ) = Vc0(HY_DENS:HY_MAGZ)+sigR(HY_DENS:HY_MAGZ)-hdtn*aBn(HY_DENS:HY_MAGZ)*dnBn
     elseif (order >= 3) then
        Wm(HY_DENS:HY_MAGZ) = Wm(HY_DENS:HY_MAGZ)+sigL(HY_DENS:HY_MAGZ)-hdtn*aBn(HY_DENS:HY_MAGZ)*dnBn
        Wp(HY_DENS:HY_MAGZ) = Wp(HY_DENS:HY_MAGZ)+sigR(HY_DENS:HY_MAGZ)-hdtn*aBn(HY_DENS:HY_MAGZ)*dnBn
     endif

     !! End of advancing the interpolated interface values by 1/2 time step for PPM or WENO
  ENDIF ! end of IF (.not. TransUpdateOnly) THEN

end Subroutine DataReconstructNormalDir
