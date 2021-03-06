!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beam3DGridWeightSquare
!!
!! NAME
!!
!!  ed_beam3DGridWeightSquare
!!
!! SYNOPSIS
!!
!!  call ed_beam3DGridWeightSquare (real,              intent (in)  :: semiAxisMajor,
!!                                  real,              intent (in)  :: semiAxisMinor,
!!                                  character (len=*), intent (in)  :: crossSectionFunctionType,
!!                                  real,              intent (in)  :: gaussianExponent,
!!                                  real,              intent (in)  :: gaussianRadiusMajor,
!!                                  real,              intent (in)  :: gaussianRadiusMinor,
!!                                  real,              intent (in)  :: gaussianCenterMajor,
!!                                  real,              intent (in)  :: gaussianCenterMinor,
!!                                  integer,           intent (in)  :: nTicsSemiAxisMajor,
!!                                  integer,           intent (in)  :: nTicsSemiAxisMinor,
!!                                  real,              intent (in)  :: delta,
!!                                  integer,           intent (in)  :: totalGridPoints,
!!                                  real,              intent (out) :: gridWeight)
!!
!! DESCRIPTION
!!
!!  Calculates the weight for the 3D beam square grid, which is defined as the sum of the individual
!!  square grid point weights. The square grid must have been set up for the grid points to be ready
!!  for retrieval. If the number of grid points processed mismatches the total number of grid points
!!  of the beam grid, the program is aborted with a message.
!!
!! ARGUMENTS
!!
!!  semiAxisMajor            : the elliptical major semiaxis of the 3D beam
!!  semiAxisMinor            : the elliptical minor semiaxis of the 3D beam 
!!  crossSectionFunctionType : the cross section function type defining the grid point weigths
!!  gaussianExponent         : the gaussian exponent
!!  gaussianRadiusMajor      : the gaussian radius along the major semiaxis
!!  gaussianRadiusMinor      : the gaussian radius along the minor semiaxis
!!  gaussianCenterMajor      : the gaussian center location along the major semiaxis
!!  gaussianCenterMinor      : the gaussian center location along the minor semiaxis
!!  nTicsSemiAxisMajor       : # of grid positions along the major semiaxis
!!  nTicsSemiAxisMinor       : # of grid positions along the minor semiaxis
!!  delta                    : the separation distance between two consecutive tics
!!  totalGridPoints          : the total number of grid points of the beam grid
!!  gridWeight               : the value of the grid weight
!!
!! NOTES
!!
!!***

subroutine ed_beam3DGridWeightSquare (semiAxisMajor,                    &
                                      semiAxisMinor,                    &
                                      crossSectionFunctionType,         &
                                      gaussianExponent,                 &
                                      gaussianRadiusMajor,              &
                                      gaussianRadiusMinor,              &
                                      gaussianCenterMajor,              &
                                      gaussianCenterMinor,              &
                                      nTicsSemiAxisMajor,               &
                                      nTicsSemiAxisMinor,               &
                                      delta,                            &
                                      totalGridPoints,                  &
                                                             gridWeight )

  use Driver_interface,             ONLY : Driver_abortFlash
  use ed_interface,                 ONLY : ed_beam3DGridPointsSquare
  use ed_beamCrossSectionFunctions, ONLY : ed_beamCrossSectionWeight

  implicit none

#include "EnergyDeposition.h"

  real,              intent (in)  :: semiAxisMajor
  real,              intent (in)  :: semiAxisMinor
  character (len=*), intent (in)  :: crossSectionFunctionType
  real,              intent (in)  :: gaussianExponent
  real,              intent (in)  :: gaussianRadiusMajor
  real,              intent (in)  :: gaussianRadiusMinor
  real,              intent (in)  :: gaussianCenterMajor
  real,              intent (in)  :: gaussianCenterMinor
  integer,           intent (in)  :: nTicsSemiAxisMajor
  integer,           intent (in)  :: nTicsSemiAxisMinor
  real,              intent (in)  :: delta
  integer,           intent (in)  :: totalGridPoints
  real,              intent (out) :: gridWeight

  logical :: moreGridPoints
  logical :: startGrid

  integer :: maxGridPoints
  integer :: n
  integer :: nGridPoints
  integer :: usedGridPoints

  real :: weightCrossSection
  real :: x,y

  real, allocatable :: xGrid (:)
  real, allocatable :: yGrid (:)
!
!
!     ...Calculate the grid weight.
!
!
  gridWeight     =  0.0

  startGrid      = .true.
  moreGridPoints = .true.
  usedGridPoints =  0
  maxGridPoints  =  BEAM_GRID_ARRAYSIZE

  allocate (xGrid (1:maxGridPoints))
  allocate (yGrid (1:maxGridPoints))

  do while (moreGridPoints)

     call ed_beam3DGridPointsSquare (semiAxisMajor,                 &
                                     semiAxisMinor,                 &
                                     nTicsSemiAxisMajor,            &
                                     nTicsSemiAxisMinor,            &
                                     delta,                         &
                                     startGrid,                     &
                                     maxGridPoints,                 &
                                                    moreGridPoints, &
                                                    nGridPoints,    &
                                                    xGrid,          &
                                                    yGrid           )

     do n = 1, nGridPoints

        x = xGrid (n)
        y = yGrid (n)

        weightCrossSection = ed_beamCrossSectionWeight (crossSectionFunctionType,       &
                                                        x        = x,                   &
                                                        y        = y,                   &
                                                        Cx       = gaussianCenterMajor, &
                                                        Cy       = gaussianCenterMinor, &
                                                        Rx       = gaussianRadiusMajor, &
                                                        Ry       = gaussianRadiusMinor, &
                                                        Exponent = gaussianExponent     )

        gridWeight = gridWeight + weightCrossSection

     end do

     usedGridPoints = usedGridPoints + nGridPoints

  end do

  if (usedGridPoints /= totalGridPoints) then
      call Driver_abortFlash ("ed_beam3DGridWeightSquare: # of used and total grid points mismatch !")
  end if

  deallocate (xGrid)
  deallocate (yGrid)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_beam3DGridWeightSquare
