!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_shockDetect
!!
!! NAME
!!
!!  hy_uhd_shockDetect
!!
!! SYNOPSIS
!!
!!  hy_uhd_shockDetect( integer (IN) :: blockID )
!!
!! DESCRIPTION
!!
!!  This routine detects strongly compressive motions in simulation
!!  by calculating undivided pressure gradients and divergence of
!!  velocity fields. Two parameters beta and delta have been set 
!!  to detect strong shocks. If such shocks exist then the unsplit
!!  scheme applies its robust flux differencings using limited slopes
!!  in data reconstruction step (see hy_uhd_dataReconstruct.F90).
!!  Different shock strengths can also be detected by lowering/increasing
!!  beta and delta values.
!!
!! ARGUMENTS
!!
!!  blockID  - local block ID
!!
!! REFERENCE 
!!
!!  Balsara and Spicer, JCP, 149:270--292, 1999.
!!
!!***

!!REORDER(4): U, scrch_Ptr

!#define DEBUG_HY_SHOCK_DETECT

Subroutine hy_uhd_shockDetect( blockID )


  use Grid_interface,    ONLY : Grid_getBlkIndexLimits, &
                                Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getDeltas
  use Hydro_data,        ONLY : hy_cfl, hy_cfl_original,&
                                hy_meshMe,              &
                                hy_RiemannSolver,       &
                                hy_hybridRiemannOnly,   &
                                hy_geometry
  use Driver_data,       ONLY : dr_nStep
  use Logfile_interface, ONLY : Logfile_open,Logfile_close
  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockID
  !! -----------------------------------------------------

  integer :: i,j,k
  logical :: SW1, SW2
  integer :: istat,k2,k3
  integer, dimension(MDIM) :: dataSize
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  real :: divv,gradPx,gradPy,gradPz
  real :: minP,minC,beta,delta
  real, dimension(:,:,:), allocatable :: Vc
  real, dimension(:,:,:,:), pointer   :: U


  ! Two parameters that can be adjusted to detect shocks
  ! with different strengths:
  ! (a) The lower the values the weaker shocks detected 
  !     (lower the values to detect more shock regions)
  ! (b) The larger the values the stronger shocks detected
  !     (increase the values to detect less shock regions)
  beta = 0.5 !0.5 !10. ! gradP
  delta= 0.1  !0.1 !2. ! divV
!!$  beta  = 0.1 !0.1
!!$  delta = 0.01

  k2=0
  k3=0
  ! Set dimensional indices
  if (NDIM > 1) k2=1
  if (NDIM > 2) k3=1


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,U,CENTER)

#ifdef SHOK_VAR
     U(SHOK_VAR,:,:,:)=0.
#else
  if (hy_RiemannSolver == HYBR) then
     call Driver_abortFlash&
          ("[hy_uhd_shockDetect]: SHOK_VAR has not been defined for shock detection")
  endif
#endif


  !! Allocate a temporary cell-centered array for sound speed
  dataSize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
  allocate(Vc(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)), stat=istat)

  !! Compute sound speed
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           Vc(i,j,k) = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
        enddo
     enddo
  enddo

  ! Always revert back to the original cfl
  hy_cfl = hy_cfl_original

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

           minP = minval(U(PRES_VAR,i-1:i+1,j-k2:j+k2,k-k3:k+k3))
           minC = minval(        Vc(i-1:i+1,j-k2:j+k2,k-k3:k+k3))

           !! We do not need to include non-Cartesian geom factors here.
           !! Undivided divV
           divv =        U(VELX_VAR,i+1,j,  k  ) - U(VELX_VAR,i-1,j,  k  )
#if NDIM > 1
           divv = divv + U(VELY_VAR,i,  j+1,k  ) - U(VELY_VAR,i,  j-1,k  )
#if NDIM == 3
           divv = divv + U(VELZ_VAR,i,  j,  k+1) - U(VELZ_VAR,i,  j,  k-1)
#endif
#endif
           divv = 0.5*divv  

           !! Undivided grad pres
           gradPx = 0.5*(U(PRES_VAR,i+1,j,  k  ) - U(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.
#if NDIM > 1
           gradPy = 0.5*(U(PRES_VAR,i,  j+1,k  ) - U(PRES_VAR,i,  j-1,k  ))
#if NDIM == 3
           gradPz = 0.5*(U(PRES_VAR,i,  j,  k+1) - U(PRES_VAR,i,  j,  k-1))
#endif
#endif
           if ( abs(gradPx)+abs(gradPy)+abs(gradPz) .ge. beta*minP ) then
              SW1 = .true.
           endif

           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif


           if (SW1 .and. SW2) then
              !local hybrid method which applies (a diffusive) HLL solver when SHOK_VAR = 1.
#ifdef SHOK_VAR
              U(SHOK_VAR,i,j,k) = 1.
#endif

              if (.not. hy_hybridRiemannOnly) then !a case with shock detection is needed
                 !$omp critical (Update_cfl)
                 !! DL: turn these cfl reduction on now
                 if (NDIM == 1) then
                    if (hy_cfl > 0.60) hy_cfl = 0.60
                 elseif (NDIM == 2) then
                    if (hy_cfl > 0.45) hy_cfl = 0.45
                 elseif (NDIM == 3) then
                    if (hy_cfl > 0.25) hy_cfl = 0.25
                 endif
                 !$omp end critical (Update_cfl)

              endif ! endif (.not. hy_hybridRiemannOnly) then

           endif !endif (SW1 .and. SW2) then

        enddo !enddo i-loop
     enddo !enddo j-loop
  enddo !enddo k-loop

  ! Release block pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  ! Deallocate sound speed array
  deallocate(Vc)

End Subroutine hy_uhd_shockDetect
