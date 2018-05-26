!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsCheck1DRec
!!
!! NAME
!!
!!  ed_beamsCheck1DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck1DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 1D rectangular
!!  grids (cartesian + spherical). All checks which depend on domain grid details should go
!!  in here. Currently it contains the following:
!!
!!         1) Check, if all beam target points are completely within the domain.
!!         2) Check, if all beam lens points are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck1DRec ()

  use Driver_interface,         ONLY : Driver_abortFlash

  use EnergyDeposition_data,    ONLY : ed_beams,         &
                                       ed_numberOfBeams, &
                                       ed_xminDomain,    &
                                       ed_xmaxDomain
  
  implicit none

#include "Flash.h"
#include "EnergyDeposition.h"
#include "constants.h"

  logical :: inDomain
  logical :: outOfDomain

  integer :: beam

  real    :: lensX
  real    :: targetX
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     lensX   = ed_beams (beam) % lensX
     targetX = ed_beams (beam) % targetX
!
!
!     ...Check, if beam target location is completely within domain.
!
!
     outOfDomain = (targetX < ed_xminDomain) .or. (targetX > ed_xmaxDomain)

     if (outOfDomain) then
         call Driver_abortFlash ("ed_beamsCheck1DRec: Beam target outside of domain!")
     end if
!
!
!     ...Repeat the same procedure with the lens location.
!
!
     inDomain = (lensX >= ed_xminDomain) .and. (lensX <= ed_xmaxDomain)

     if (inDomain) then
         call Driver_abortFlash ("ed_beamsCheck1DRec: Beam lens inside the domain!")
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsCheck1DRec
