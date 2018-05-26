!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamsCheck3DRec
!!
!! NAME
!!
!!  ed_beamsCheck3DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck3DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 3D rectangular
!!  grids (cartesian). All checks which depend on domain grid details should go in here.
!!  Currently it contains the following:
!!
!!         1) Check, if all beam elliptical target areas are completely within the domain.
!!         2) Check, if all beam elliptical lens areas are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck3DRec ()

  use Driver_interface,         ONLY : Driver_abortFlash

  use EnergyDeposition_data,    ONLY : ed_beams,         &
                                       ed_numberOfBeams, &
                                       ed_xminDomain,    &
                                       ed_xmaxDomain,    &
                                       ed_yminDomain,    &
                                       ed_ymaxDomain,    &
                                       ed_zminDomain,    &
                                       ed_zmaxDomain
  
  implicit none

#include "Flash.h"
#include "EnergyDeposition.h"
#include "constants.h"

  logical :: inDomain
  logical :: outOfDomain

  integer :: beam

  real    :: ex, ey, ez
  real    :: lensX, lensY, lensZ
  real    :: s1, s2
  real    :: targetX, targetY, targetZ
  real    :: u1x, u1y, u1z
  real    :: u2x, u2y, u2z
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, ed_numberOfBeams

     u1x     = ed_beams (beam) % semiAxisUnitMajorX
     u1y     = ed_beams (beam) % semiAxisUnitMajorY
     u1z     = ed_beams (beam) % semiAxisUnitMajorZ
     u2x     = ed_beams (beam) % semiAxisUnitMinorX
     u2y     = ed_beams (beam) % semiAxisUnitMinorY
     u2z     = ed_beams (beam) % semiAxisUnitMinorZ
     s1      = ed_beams (beam) % targetSemiAxisMajor
     s2      = ed_beams (beam) % targetSemiAxisMinor
     targetX = ed_beams (beam) % targetX
     targetY = ed_beams (beam) % targetY
     targetZ = ed_beams (beam) % targetZ
!
!
!     ...The implicit form of the elliptical target boundary curve is:
!
!              e = s1 * u1 * cos (t)  +  s2 * u2 * sin (t)
!
!        where u1,u2 stand for the two semiaxis unit vectors in the local target coordinate
!        system, s1,s2 are the two semiaxes lengths, 'e' is a vector on the elliptical boundary
!        and 't' is the implicit parameter ranging from 0 to 2pi. Differentiating this
!        equation with respect to 't' we obtain:
!
!                     t = arctan ( [s2 * u2] / [s1 * u1] )
!
!        which, when simplifying, leads to the simple minimax equation:
!
!                  e (min,max) = +/- sqrt ([s1 *u1]^2 + [s2 *u2]^2)
!
!        or, in component form:
!
!                 ex (min,max) = +/- sqrt ([s1 *u1x]^2 + [s2 *u2x]^2)
!                 ey (min,max) = +/- sqrt ([s1 *u1y]^2 + [s2 *u2y]^2)
!                 ez (min,max) = +/- sqrt ([s1 *u1z]^2 + [s2 *u2z]^2)
!
!        Since the domain boundaries are in terms of the global coordinate system, we need
!        to convert the elliptical curve points from the local to the global coordinate
!        system via:
!
!                                E = e + T
!
!        where T is the position of the target center. The results of E will then
!        be tested against the domain boundaries.
!
!
     ex = sqrt ( (s1 * u1x) ** 2 + (s2 * u2x) ** 2 )
     ey = sqrt ( (s1 * u1y) ** 2 + (s2 * u2y) ** 2 )
     ez = sqrt ( (s1 * u1z) ** 2 + (s2 * u2z) ** 2 )
!
!
!     ...Check, if beam target area is incident completely within domain.
!
!
     outOfDomain =     (targetX - ex < ed_xminDomain) &
                  .or. (targetY - ey < ed_yminDomain) &
                  .or. (targetZ - ez < ed_zminDomain) &
                  .or. (targetX + ex > ed_xmaxDomain) &
                  .or. (targetY + ey > ed_ymaxDomain) &
                  .or. (targetZ + ez > ed_zmaxDomain)

     if (outOfDomain) then
         call Driver_abortFlash ("ed_beamsCheck3DRec: Beam target (partially) outside of domain!")
     end if
!
!
!     ...Repeat the same procedure with the lens.
!
!
     s1    = ed_beams (beam) % lensSemiAxisMajor
     s2    = ed_beams (beam) % lensSemiAxisMinor
     lensX = ed_beams (beam) % lensX
     lensY = ed_beams (beam) % lensY
     lensZ = ed_beams (beam) % lensZ

     ex = sqrt ( (s1 * u1x) ** 2 + (s2 * u2x) ** 2 )
     ey = sqrt ( (s1 * u1y) ** 2 + (s2 * u2y) ** 2 )
     ez = sqrt ( (s1 * u1z) ** 2 + (s2 * u2z) ** 2 )

     inDomain =      (lensX - ex >= ed_xminDomain) &
               .and. (lensY - ey >= ed_yminDomain) &
               .and. (lensZ - ez >= ed_zminDomain) &
               .and. (lensX + ex <= ed_xmaxDomain) &
               .and. (lensY + ey <= ed_ymaxDomain) &
               .and. (lensZ + ez <= ed_zmaxDomain)

     if (inDomain) then
         call Driver_abortFlash ("ed_beamsCheck3DRec: Beam lens (partially) inside the domain!")
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beamsCheck3DRec
