!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_fluxLimiter
!!
!!  NAME 
!!
!!  Diffuse_fluxLimiter
!!
!!  SYNOPSIS
!!
!! 
!!  call Diffuse_fluxLimiter(integer, intent(in) :: idcoeff, 
!!                           integer, intent(in) :: ifunc,
!!                           integer, intent(in) :: ifl,
!!                           integer, intent(in) :: mode,
!!                           integer, intent(IN) :: blkcnt,
!!                           integer, intent(IN) :: blklst(blkcnt))
!!
!!  DESCRIPTION 
!!      This routine modifes the diffusion coefficient (idcoeff) 
!!      and applies limiting.
!!
!! ARGUMENTS
!!
!!   idcoeff : index into solution vector, indicating a variable that holds
!!             the coefficient to which the limiter is applied.
!!   ifunc   : index into solution vector giving the quantity whose flux
!!             is to be limited.
!!   ifl     : index into solution vector giving flux imiter variable.
!!   mode    : Flux limiter mode.
!!   blkcnt  : The number of blocks in the list
!!   blklst  : The list of blocks on which the solution must be updated.
!!
!! SIDE EFFECTS
!!
!! Modifies the diffusion coefficient variable, indicted by the
!! argument idcoef
!!
!!***
subroutine Diffuse_fluxLimiter(idcoef, ifunc, ifl, mode, blkcnt, blklst)
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getDeltas, Grid_fillGuardCells
  use Diffuse_data, ONLY: diff_meshMe, diff_geometry
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
  
  ! Arguments:
  integer, intent(in) :: idcoef
  integer, intent(in) :: ifunc
  integer, intent(in) :: ifl
  integer, intent(in) :: mode
  integer, intent(IN) :: blkcnt
  integer, intent(IN) :: blklst(blkcnt)

  integer :: lb, i, j, k
  real, pointer :: blkPtr(:,:,:,:)
  integer :: blklim(2,MDIM), blklimgc(2,MDIM)
  real :: dcoef
  real :: fl, fp1, fm1
  real :: delta(MDIM) ! The cell width
  real :: maggrad ! The magnitude of the gradient
  real, allocatable :: xcent(:)
  real, allocatable :: ycent(:)
  real, allocatable :: zcent(:)
  real :: gradi, gradj, gradk
  
    
  if(mode == FL_NONE) return

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb), blklim, blklimgc)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     call Grid_getDeltas(blklst(lb), delta)
     
     allocate(xcent(blklimgc(HIGH, IAXIS)))
     call Grid_getCellCoords(IAXIS, blklst(lb), CENTER, .true., &
          xcent, blklimgc(HIGH, IAXIS))

     allocate(ycent(blklimgc(HIGH, JAXIS)))
     call Grid_getCellCoords(JAXIS, blklst(lb), CENTER, .true., &
          ycent, blklimgc(HIGH, JAXIS))

     allocate(zcent(blklimgc(HIGH, KAXIS)))          
     call Grid_getCellCoords(KAXIS, blklst(lb), CENTER, .true., &
          zcent, blklimgc(HIGH, KAXIS))

     do k = blklim(LOW,KAXIS)-K3D, blklim(HIGH,KAXIS)+K3D
        do j = blklim(LOW,JAXIS)-K2D, blklim(HIGH,JAXIS)+K2D
           do i = blklim(LOW,IAXIS)-1, blklim(HIGH,IAXIS)+1

              ! Compute the magnitude of the gradient:
              fm1 = blkPtr(ifunc, i-1, j, k)
              fp1 = blkPtr(ifunc, i+1, j, k)
              gradi = ((fp1-fm1)/(2.0*delta(IAXIS)))**2
              maggrad = gradi

#if NDIM >= 2
              fm1 = blkPtr(ifunc, i, j-1, k)
              fp1 = blkPtr(ifunc, i, j+1, k)
              if(diff_geometry == POLAR) then
                 gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2 / xcent(i)
              else
                 gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2
              end if
              maggrad = maggrad + gradj
#endif

#if NDIM == 3
              fm1 = blkPtr(ifunc, i, j, k-1)
              fp1 = blkPtr(ifunc, i, j, k+1)
              gradk = ((fp1-fm1)/(2.0*delta(KAXIS)))**2
              maggrad = maggrad + gradk
#endif

              maggrad = sqrt(maggrad)

              ! Now compute the new diffusion coefficient:
              dcoef = blkPtr(idcoef, i, j, k)
              fl    = blkPtr(ifl, i, j, k)

              select case(mode)
              case(FL_HARMONIC)
                 dcoef = 1.0 / (1.0/dcoef + maggrad/(fl + 1.0D-100))

              case(FL_MINMAX)
                 dcoef = min(dcoef, fl/(maggrad + 1.0e-10))

              case(FL_LARSEN)
                 dcoef = 1.0 / sqrt((1.0/dcoef)**2 + (maggrad/(fl + 1.0D-100))**2)

              case DEFAULT
                 call Driver_abortFlash("Invalid Flux limiter type")
              end select

              blkPtr(idcoef, i, j, k) = dcoef
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     
     deallocate(xcent)
     deallocate(ycent)
     deallocate(zcent)
     
  end do

end subroutine Diffuse_fluxLimiter
