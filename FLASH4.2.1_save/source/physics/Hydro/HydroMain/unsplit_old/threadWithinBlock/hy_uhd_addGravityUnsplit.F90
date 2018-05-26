!!****if* source/physics/Hydro/HydroMain/unsplit_old/threadWithinBlock/hy_uhd_addGravityUnsplit
!!
!!
!! NAME
!!
!!  hy_uhd_addGravityUnsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_addGravityUnsplit( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimits(LOW:HIGH,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real(IN)    :: dt,
!!                            real(IN)    :: gravX,
!!                            real(IN)    :: gravY,
!!                            real(IN)    :: gravZ )
!!
!!
!! DESCRIPTION
!!
!!  Adds the second part of the gravitational force to the momenta and energy.
!!  The first half is added in unsplitUpdate.  This centers the gravitational
!!  force appropriately on the current timestep rather than extrapolating.
!!
!! ARGUMENTS
!!
!!  blockID   - local block ID
!!  blkLimits - block limits 
!!  dataSize  - array size of gravX,Y,Z in non-fixed block size mode in UG
!!  dt        - timestep
!!  gravX     - gravity source term in x-direction
!!  gravY     - gravity source term in y-direction
!!  gravZ     - gravity source term in z-direction
!!
!!
!!***

!!REORDER(4):scrch_Ctr

Subroutine hy_uhd_addGravityUnsplit&
     (blockID,blkLimits,dataSize,dt,gravX,gravY,gravZ,halfTimeGravUpdate)

  use Hydro_data,      ONLY : hy_useGravity,        &
                              hy_useGravHalfUpdate, &
                              hy_threadWithinBlock

  use Grid_interface,  ONLY : Grid_getBlkPtr,       &
                              Grid_releaseBlkPtr

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  integer, dimension(MDIM), intent(IN) :: dataSize
  real,    intent(IN) :: dt

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
       gravX,gravY,gravZ
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(IN) :: &
       gravX,gravY,gravZ
#endif
  logical, intent(IN) :: halfTimeGravUpdate
  !! -----------------------------------------------------
  real :: hdt
  integer :: i,j,k
  real, pointer, dimension(:,:,:,:) :: U
  real, dimension(3) :: momentaOld, momentaNew



  integer, dimension(2) :: gravPos
  integer :: i0,imax,j0,jmax,k0,kmax
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: grav
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: grav
#endif

  real, pointer, dimension(:,:,:,:) :: scrch_Ctr
  real, dimension(MDIM)        :: Gxp,Gxn,Gyp,Gyn,Gzp,Gzn
  real,allocatable,dimension(:,:,:) :: gravXP,gravXN,gravYP,gravYN,gravZP,gravZN

  hdt = 0.5 * dt

  if (halfTimeGravUpdate) then

!! Execute this routine only when Gravity potential is used and 
!! therefore needing scratchvar setup!
#ifdef FLASH_UHD_NEED_SCRATCHVARS

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
#elif NDIM == 3
  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = blkLimits(LOW, KAXIS)
  kmax = blkLimits(HIGH,KAXIS)
#endif


  !! Get block pointer for storages of Riemann states
  call Grid_getBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)

  do k=k0-2,kmax+2
     do j=j0-2,jmax+2
        do i=i0-2,imax+2
           call gravReconOneZone(blockID,blkLimits,dataSize,i,j,k,&
                                 gravX,gravY,gravZ, & 
                                 Gxp,Gxn,Gyp,Gyn,Gzp,Gzn)

           scrch_Ctr(XP02_SCRATCH_CENTER_VAR:XP04_SCRATCH_CENTER_VAR,i,j,k)=&
                scrch_Ctr(XP02_SCRATCH_CENTER_VAR:XP04_SCRATCH_CENTER_VAR,i,j,k)+hdt*Gxp(1:3)

           scrch_Ctr(XN02_SCRATCH_CENTER_VAR:XN04_SCRATCH_CENTER_VAR,i,j,k)=&
                scrch_Ctr(XN02_SCRATCH_CENTER_VAR:XN04_SCRATCH_CENTER_VAR,i,j,k)+hdt*Gxn(1:3)
#if NDIM >= 2
           scrch_Ctr(YP02_SCRATCH_CENTER_VAR:YP04_SCRATCH_CENTER_VAR,i,j,k)=&
                scrch_Ctr(YP02_SCRATCH_CENTER_VAR:YP04_SCRATCH_CENTER_VAR,i,j,k)+hdt*Gyp(1:3)

           scrch_Ctr(YN02_SCRATCH_CENTER_VAR:YN04_SCRATCH_CENTER_VAR,i,j,k)=&
                scrch_Ctr(YN02_SCRATCH_CENTER_VAR:YN04_SCRATCH_CENTER_VAR,i,j,k)+hdt*Gyn(1:3)
#if NDIM == 3
           scrch_Ctr(ZP02_SCRATCH_CENTER_VAR:ZP04_SCRATCH_CENTER_VAR,i,j,k)=&
                scrch_Ctr(ZP02_SCRATCH_CENTER_VAR:ZP04_SCRATCH_CENTER_VAR,i,j,k)+hdt*Gzp(1:3)

           scrch_Ctr(ZN02_SCRATCH_CENTER_VAR:ZN04_SCRATCH_CENTER_VAR,i,j,k)=&
                scrch_Ctr(ZN02_SCRATCH_CENTER_VAR:ZN04_SCRATCH_CENTER_VAR,i,j,k)+hdt*Gzn(1:3)
#endif
#endif
        enddo
     enddo
  enddo
  
  !! Release block pointer for storages of Riemann states
  call Grid_releaseBlkPtr(blockID,scrch_Ctr,SCRATCH_CTR)



#endif
! endif of #ifdef FLASH_UHD_NEED_SCRATCHVARS




  else

  !! Get block pointer for storages of Riemann states
  call Grid_getBlkPtr(blockID,U,CENTER)

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared(U,blkLimits,hdt,gravX,gravY,gravZ) &
  !$omp private(i,j,k,momentaOld,momentaNew)

  !$omp do schedule(static) 
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           momentaOld(1:3) = U(DENS_VAR,i,j,k)*U(VELX_VAR:VELZ_VAR,i,j,k)

           momentaNew(1:3) = momentaOld(1:3)&
                + hdt*U(DENS_VAR,i,j,k)*(/gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k)/)

           U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k) &
                + hdt*dot_product(momentaNew(1:3),(/gravX(i,j,k),gravY(i,j,k),gravZ(i,j,k)/))

           U(VELX_VAR:VELZ_VAR,i,j,k) = momentaNew(1:3)/U(DENS_VAR,i,j,k)
        enddo
     enddo
  enddo
  !$omp end do
  !$omp end parallel

  !! Release block pointer for storages of Riemann states
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  endif

end Subroutine hy_uhd_addGravityUnsplit


#define GRAV_X 1
#define GRAV_Y 2
#define GRAV_Z 3

subroutine gravReconOneZone(blockID,blkLimits,dataSize,i,j,k, &
                           gravX,gravY,gravZ, & 
                           Gxp,Gxn,Gyp,Gyn,Gzp,Gzn)

  use Grid_interface, ONLY: Grid_getDeltas
  use hy_uhd_slopeLimiters, ONLY : minmod,mc,vanLeer

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID, i, j, k
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  integer, dimension(MDIM), intent(IN) :: dataSize

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: &
       gravX,gravY,gravZ
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(IN) :: &
       gravX,gravY,gravZ
#endif

  real, dimension(MDIM), intent(OUT) :: Gxp,Gxn,Gyp,Gyn,Gzp,Gzn

  !! -----------------------------------------------------

  real, dimension(MDIM) :: delta
  real :: grad

  !! Note: It works better with using slope limiters although gpol is usually smooth.
  !!       Also, vanLeer seems most suitable among the choices of slope limiters, and
  !!       mc is not too good in the dustcollapse problem in 2D cylindrical.
  logical :: oldMethod = .false.

  Gxp = 0.
  Gxn = 0.
  Gyp = 0.
  Gyn = 0.
  Gzp = 0.
  Gzn = 0.

  call Grid_getDeltas(blockID,delta)

  !! Note: grav*(:,:,:) are gravity at n+1/2 state, where *=X,Y, & Z.
  !(1)-x
  ! x-component grav at (i+1/2,j,k) & (i-1/2,j,k)
  if (oldMethod) then
     grad = (-gravX(i+2,j,k)+8.*gravX(i+1,j,k)-8.*gravX(i-1,j,k)+gravX(i-2,j,k))/12.
  else
     grad = vanLeer(gravX(i+1,j,k)-gravX(i,j,k),gravX(i,j,k)-gravX(i-1,j,k))
  endif
  grad = 0.5*grad

  Gxp(GRAV_X) = gravX(i,j,k) + grad
  Gxn(GRAV_X) = gravX(i,j,k) - grad

#if NDIM >= 2
  !(2)-x
  ! y-component grav at (i+1/2,j,k) & (i-1/2,j,k)
  if (oldMethod) then
     grad = (-gravY(i+2,j,k)+8.*gravY(i+1,j,k)-8.*gravY(i-1,j,k)+gravY(i-2,j,k))/12.
  else
     grad = vanLeer(gravY(i+1,j,k)-gravY(i,j,k),gravY(i,j,k)-gravY(i-1,j,k))
  endif
  grad = 0.5*grad

  Gxp(GRAV_Y) = gravY(i,j,k) + grad
  Gxn(GRAV_Y) = gravY(i,j,k) - grad


  !(1)-y
  ! x-component grav at (i,j+1/2,k) & (i,j-1/2,k)
  if (oldMethod) then
     grad = (-gravX(i,j+2,k)+8.*gravX(i,j+1,k)-8.*gravX(i,j-1,k)+gravX(i,j-2,k))/12.
  else
     grad = vanLeer(gravX(i,j+1,k)-gravX(i,j,k),gravX(i,j,k)-gravX(i,j-1,k))
  endif
  grad = 0.5*grad

  Gyp(GRAV_X) = gravX(i,j,k) + grad
  Gyn(GRAV_X) = gravX(i,j,k) - grad

  !(2)-y
  ! y-component grav at (i,j+1/2,k) & (i,j-1/2,k)
  if (oldMethod) then
     grad = (-gravY(i,j+2,k)+8.*gravY(i,j+1,k)-8.*gravY(i,j-1,k)+gravY(i,j-2,k))/12.
  else
     grad = vanLeer(gravY(i,j+1,k)-gravY(i,j,k),gravY(i,j,k)-gravY(i,j-1,k))
  endif
  grad = 0.5*grad

  Gyp(GRAV_Y) = gravY(i,j,k) + grad
  Gyn(GRAV_Y) = gravY(i,j,k) - grad

#if NDIM == 3

  !(3)-x
  ! z-component grav at (i+1/2,j,k) & (i-1/2,j,k)
  if (oldMethod) then
     grad = (-gravZ(i+2,j,k)+8.*gravZ(i+1,j,k)-8.*gravZ(i-1,j,k)+ gravZ(i-2,j,k))/12.
  else
     grad = vanLeer(gravZ(i+1,j,k)-gravZ(i,j,k),gravZ(i,j,k)-gravZ(i-1,j,k))
  endif
  grad = 0.5*grad

  Gxp(GRAV_Z) = gravZ(i,j,k) + grad
  Gxn(GRAV_Z) = gravZ(i,j,k) - grad


  !(3)-y
  ! z-component grav at (i,j+1/2,k) & (i,j-1/2,k)
  if (oldMethod) then
     grad = (-gravZ(i,j+2,k)+8.*gravZ(i,j+1,k)-8.*gravZ(i,j-1,k)+gravZ(i,j-2,k))/12.
  else
     grad = vanLeer(gravZ(i,j+1,k)-gravZ(i,j,k),gravZ(i,j,k)-gravZ(i,j-1,k))
  endif
  grad = 0.5*grad

  Gyp(GRAV_Z) = gravZ(i,j,k) + grad
  Gyn(GRAV_Z) = gravZ(i,j,k) - grad

  !(1)-z
  ! x-component grav at (i,j,k+1/2) & (i,j,k-1/2)
  if (oldMethod) then
     grad = (-gravX(i,j,k+2)+8.*gravX(i,j,k+1)-8.*gravX(i,j,k-1)+ gravX(i,j,k-2))/12.
  else
     grad = vanLeer(gravX(i,j,k+1)-gravX(i,j,k),gravX(i,j,k)-gravX(i,j,k-1))
  endif
  grad = 0.5*grad

  Gzp(GRAV_X) = gravX(i,j,k) + grad
  Gzn(GRAV_X) = gravX(i,j,k) - grad

  !(2)-z
  ! y-component grav at (i,j,k+1/2) & (i,j,k-1/2)
  if (oldMethod) then
     grad = (-gravY(i,j,k+2)+8.*gravY(i,j,k+1)-8.*gravY(i,j,k-1)+gravY(i,j,k-2))/12.
  else
     grad = vanLeer(gravY(i,j,k+1)-gravY(i,j,k),gravY(i,j,k)-gravY(i,j,k-1))
  endif
  grad = 0.5*grad

  Gzp(GRAV_Y) = gravY(i,j,k) + grad
  Gzn(GRAV_Y) = gravY(i,j,k) - grad

  !(3)-z
  ! z-component grav at (i,j,k+1/2) & (i,j,k-1/2)
  if (oldMethod) then
     grad = (-gravZ(i,j,k+2)+8.*gravZ(i,j,k+1)-8.*gravZ(i,j,k-1)+gravZ(i,j,k-2))/12.
  else
     grad = vanLeer(gravZ(i,j,k+1)-gravZ(i,j,k),gravZ(i,j,k)-gravZ(i,j,k-1))
  endif
  grad = 0.5*grad

  Gzp(GRAV_Z) = gravZ(i,j,k) + grad
  Gzn(GRAV_Z) = gravZ(i,j,k) - grad

#endif
#endif

  return

end subroutine gravReconOneZone
