!!****if* source/physics/Hydro/HydroMain/unsplit_old/threadWithinBlock/hy_uhd_putGravityUnsplit
!!
!! NAME
!!
!!  hy_uhd_PutGravityUnsplit
!!
!! SYNOPSIS
!!
!!  hy_uhd_putGravityUnsplit( integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(2,MDIM),
!!                            integer(IN) :: dataSize(MDIM),
!!                            real   (IN) :: dt,
!!                            real   (IN) :: dtOld,
!!                            real(OUT)   :: gravX(3,:,:,:),
!!                            real(OUT)   :: gravY(3,:,:,:),
!!                            real(OUT)   :: gravZ(3,:,:,:))
!!
!! ARGUMENTS
!!
!!  blockID     - a current block ID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  dataSize    - dimensions for gravX, gravY and gravZ arrays
!!  dt          - timestep
!!  dtOld       - old timestep (needed for temporal extrapolations of gravity)
!!  gravX       - gravity components in x-direcition at time steps n, n+1/2 and n+1
!!  gravY       - gravity components in y-direcition at time steps n, n+1/2 and n+1
!!  gravZ       - gravity components in z-direcition at time steps n, n+1/2 and n+1
!!
!! DESCRIPTION
!!
!!  This routine puts gravity components to arrays gravX, gravY and gravZ:
!!  gravX(1,:,:,:) includes gravity components at time step n,
!!  gravX(2,:,:,:) includes extrapolated gravity components at time step n+1/2,
!!  gravX(3,:,:,:) includes extrapolated gravity components at time step n+1.
!!
!!*** 

Subroutine hy_uhd_putGravityUnsplit&
     (blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)

  use Gravity_interface, ONLY : Gravity_accelOneRow
  use Hydro_data, ONLY : hy_threadWithinBlock

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"
#include "Flash_omp.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
  integer, dimension(MDIM), intent(IN) :: dataSize
  real,    intent(IN) :: dt, dtOld

#ifdef FIXEDBLOCKSIZE
  real, dimension(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(OUT) :: &
       gravX,gravY,gravZ
#else
  real, dimension(3,dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(OUT) :: &
       gravX,gravY,gravZ
#endif
  !! -----------------------------------------------------

  integer, dimension(2) :: gravPos
  integer :: ix, iy, iz
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: &
       grav,gravOld
#else
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)) :: &
       grav,gravOld
#endif
  real    :: dtFactor


  ! time factor to extrapolate gravity components to n+1/2 & n+1 steps from using
  ! gravity components at n-1 and n steps.
  dtFactor = dt/dtOld

  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp shared( gravX,gravY,gravZ,blockID,dataSize,blkLimitsGC,dtFactor) &
  !$omp private(grav,gravOld,ix,iy,iz,gravPos)

  ! initialize arrays
  grav = 0.
  gravOld = 0.

#if NDIM == 3
  !$omp do schedule(static)
#endif  
  do iz=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
#if NDIM == 2
     !$omp do schedule(static)
#endif    
     do iy=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        gravPos(1)=iy
        gravPos(2)=iz

#if defined(GPOT_VAR) && defined(GPOL_VAR) && defined(FLASH_GRAVITY_TIMEDEP)

        ! Gravity implementation defines FLASH_GRAVITY_TIMEDEP -> time-dependent gravity field,
        ! interpolate the acceleration linearly in time (pointwise) - KW
        call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),   grav(:,iy,iz),GPOT_VAR)

        ! gravity at time step n
        gravX(1,:,iy,iz) = grav(:,iy,iz)

        ! gravity at time step n+1/2
        gravX(2,:,iy,iz) = grav(:,iy,iz) + 0.5*dtFactor * (grav(:,iy,iz) - gravOld(:,iy,iz))

        ! gravity at time step n+1
        gravX(3,:,iy,iz) = grav(:,iy,iz) +     dtFactor * (grav(:,iy,iz) - gravOld(:,iy,iz))

#else
        ! FLASH_GRAVITY_TIMEDEP not defined -> assume time-independent gravity field.
        ! Also if GPOT_VAR or GPOL_VAR defined -> use current accel without time
        ! interpolation, i.e., handle like time-independent gravity field - KW
        call Gravity_accelOneRow(gravPos,DIR_X,blockID,dataSize(IAXIS),   gravX(1,:,iy,iz))
        gravX(2,:,iy,iz)=gravX(1,:,iy,iz)
        gravX(3,:,iy,iz)=gravX(1,:,iy,iz)
#endif



     enddo
#if NDIM == 2
     !$omp end do
#endif    
  enddo
#if NDIM == 3
  !$omp end do
#endif    

  if (NDIM >= 2) then
     ! initialize arrays
     grav = 0.
     gravOld = 0.

#if NDIM == 3     
     !$omp do schedule(static)
#endif    
     do iz=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
#if NDIM == 2
        !$omp do schedule(static)
#endif         
        do ix=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           gravPos(1)=ix
           gravPos(2)=iz

#if defined(GPOT_VAR) && defined(GPOL_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
           call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),   grav(ix,:,iz),GPOT_VAR)

           ! gravity at time step n
           gravY(1,ix,:,iz) = grav(ix,:,iz)

           ! gravity at time step n+1/2 (extrapolation)
           gravY(2,ix,:,iz) = grav(ix,:,iz) + 0.5*dtFactor * (grav(ix,:,iz) - gravOld(ix,:,iz))

           ! gravity at time step n+1 (extrapolation)
           gravY(3,ix,:,iz) = grav(ix,:,iz) +     dtFactor * (grav(ix,:,iz) - gravOld(ix,:,iz))
#else
           call Gravity_accelOneRow(gravPos,DIR_Y,blockID,dataSize(JAXIS),gravY(1,ix,:,iz))
           gravY(2,ix,:,iz) = gravY(1,ix,:,iz)
           gravY(3,ix,:,iz) = gravY(1,ix,:,iz)
#endif
        enddo
#if NDIM == 2
       !$omp end do
#endif                
     enddo
#if NDIM == 3
     !$omp end do
#endif                

     if (NDIM == 3) then
        ! initialize arrays
        grav = 0.
        gravOld = 0.

        !$omp do schedule(static)
        do iy=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do ix=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              gravPos(1)=ix
              gravPos(2)=iy

#if defined(GPOT_VAR) && defined(GPOL_VAR) && defined(FLASH_GRAVITY_TIMEDEP)
              call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),   grav(ix,iy,:),GPOT_VAR)

              ! gravity at time step n
              gravZ(1,ix,iy,:) = grav(ix,iy,:)

              ! gravity at time step n+1/2 (extrapolation)
              gravZ(2,ix,iy,:) = grav(ix,iy,:) + 0.5*dtFactor * (grav(ix,iy,:) - gravOld(ix,iy,:))

              ! gravity at time step n+1 (extrapolation)
              gravZ(3,ix,iy,:) = grav(ix,iy,:) +     dtFactor * (grav(ix,iy,:) - gravOld(ix,iy,:))
#else
              call Gravity_accelOneRow(gravPos,DIR_Z,blockID,dataSize(KAXIS),gravZ(1,ix,iy,:))
              gravZ(2,ix,iy,:) = gravZ(1,ix,iy,:)
              gravZ(3,ix,iy,:) = gravZ(1,ix,iy,:)
#endif
           enddo
        enddo
        !$omp end do
     endif
  endif
  !$omp end parallel
End Subroutine hy_uhd_putGravityUnsplit
