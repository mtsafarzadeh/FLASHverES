!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_shockDetect
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
!!  to detect very strong shocks. If such shocks exit then the unsplit
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
                                hy_RiemannSolver,       &
                                hy_hybridRiemannOnly
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
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real :: divv,gradPx,gradPy,gradPz
  real, pointer, dimension(:,:,:,:) :: U
  logical :: SW1, SW2
  real :: minP,minC,beta,delta
  real, dimension(MDIM) :: del
  integer,dimension(MDIM) :: dataSize
  real, allocatable, dimension(:,:,:) :: Vc
  integer :: istat

  ! Two parameters that can be adjusted to detect shocks
  ! with different strengths:
  ! (a) The lower the values the weaker shocks detected and
  ! (b) The larger the values the stronger shocks detected.
  beta = 0.5 !10.
  delta= 0.1 !2.

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,U,CENTER)
  call Grid_getDeltas(blockID,del)
  if (hy_RiemannSolver == HYBR) then
#ifdef SHOK_VAR
     U(SHOK_VAR,:,:,:)=0.
#else
     call Driver_abortFlash&
          ("[hy_uhd_shockDetect]: SHOK_VAR has not been defined for shock detection")
#endif
  endif


  !! Allocate a temporary cell-centered array for sound speed
  dataSize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
  allocate(Vc(dataSize(IAXIS), dataSize(JAXIS), dataSize(KAXIS)), stat=istat)

  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           Vc(i,j,k) = sqrt(U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k))
        enddo
     enddo
  enddo


  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

#if NDIM == 1
           minP = minval(U(PRES_VAR,i-1:i+1,j,k))
           minC = minval(Vc(i-1:i+1,j,k))

           divv   = 0.5*(U(VELX_VAR,i+1,j,  k  ) - U(VELX_VAR,i-1,j,  k  ))

           gradPx = 0.5*(U(PRES_VAR,i+1,j,  k  ) - U(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.

           if ( abs(gradPx) .ge. beta*minP ) then
              SW1 = .true.
           endif
           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif

           if (SW1 .and. SW2) then
              if (hy_RiemannSolver == HYBR) then
                 !local hybrid method which applies (a diffusive) HLL solver when SHOK_VAR = 1.
#ifdef SHOK_VAR
                 U(SHOK_VAR,i,j,k) = 1.
#endif
              endif

              if (.not. hy_hybridRiemannOnly) then !a case with shock detection is needed
                 if (hy_cfl > 0.6) hy_cfl = 0.6
              endif
           endif



#elif NDIM == 2
           minP = minval(U(PRES_VAR,i-1:i+1,j-1:j+1,k))
           minC = minval(Vc(i-1:i+1,j-1:j+1,k))

           divv   = 0.5*(U(VELX_VAR,i+1,j,  k  ) - U(VELX_VAR,i-1,j,  k  )&
                        +U(VELY_VAR,i,  j+1,k  ) - U(VELY_VAR,i,  j-1,k  ))

           gradPx = 0.5*(U(PRES_VAR,i+1,j,  k  ) - U(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.5*(U(PRES_VAR,i,  j+1,k  ) - U(PRES_VAR,i,  j-1,k  ))
           gradPz = 0.

           if ( abs(gradPx)+abs(gradPy) .ge. beta*minP ) then
              SW1 = .true.
           endif
           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif

           if (SW1 .and. SW2) then
              if (hy_RiemannSolver == HYBR) then
                 !local hybrid method which applies (a diffusive) HLL solver when SHOK_VAR = 1.
#ifdef SHOK_VAR
                 U(SHOK_VAR,i,j,k) = 1.
#endif
              endif

              if (.not. hy_hybridRiemannOnly) then
                 if (hy_cfl > 0.45) hy_cfl = 0.45
              endif
           endif


#elif NDIM == 3
           minP = minval(U(PRES_VAR,i-1:i+1,j-1:j+1,k-1:k+1))
           minC = minval(Vc(i-1:i+1,j-1:j+1,k-1:k+1))

           divv   = 0.5*(U(VELX_VAR,i+1,j,  k  ) - U(VELX_VAR,i-1,j,  k  )&
                        +U(VELY_VAR,i,  j+1,k  ) - U(VELY_VAR,i,  j-1,k  )&
                        +U(VELZ_VAR,i,  j,  k+1) - U(VELZ_VAR,i,  j,  k-1))

           gradPx = 0.5*(U(PRES_VAR,i+1,j,  k  ) - U(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.5*(U(PRES_VAR,i,  j+1,k  ) - U(PRES_VAR,i,  j-1,k  ))
           gradPz = 0.5*(U(PRES_VAR,i,  j,  k+1) - U(PRES_VAR,i,  j,  k-1))

           if ( abs(gradPx)+abs(gradPy)+abs(gradPz) .ge. beta*minP ) then
              SW1 = .true.
           endif
           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif

           if (SW1 .and. SW2) then
              if (hy_RiemannSolver == HYBR) then
                 !local hybrid method which applies (a diffusive) HLL solver when SHOK_VAR = 1.
#ifdef SHOK_VAR
                 U(SHOK_VAR,i,j,k) = 1.
#endif
              endif

              if (.not. hy_hybridRiemannOnly) then
                 if (hy_cfl > 0.25) hy_cfl = 0.25
              endif
           endif

#endif

        enddo
     enddo
  enddo

  ! Release block pointer
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  ! Deallocate sound speed array
  deallocate(Vc)

End Subroutine hy_uhd_shockDetect
