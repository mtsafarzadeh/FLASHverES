!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_traceBlockRays2DCyl3D
!!
!! NAME
!!
!!  ed_traceBlockRays2DCyl3D
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays2DCyl3D (real    (in)    :: timeStep,
!!                                 integer (in)    :: rayFirst
!!                                 integer (in)    :: rayLast,
!!                                 integer (in)    :: iminBlock,
!!                                 integer (in)    :: imaxBlock,
!!                                 integer (in)    :: jminBlock,
!!                                 integer (in)    :: jmaxBlock,
!!                                 real    (in)    :: xminBlock,
!!                                 real    (in)    :: xmaxBlock,
!!                                 real    (in)    :: zminBlock,
!!                                 real    (in)    :: zmaxBlock,
!!                                 real    (in)    :: deltaInvX,
!!                                 real    (in)    :: deltaInvZ,
!!                                 logical (in)    :: blockReflectMinX,
!!                                 logical (in)    :: blockReflectMaxX,
!!                                 logical (in)    :: blockReflectMinZ,
!!                                 logical (in)    :: blockReflectMaxZ,
!!                                 real    (inout) :: wedgeEnergyDepot (:,:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active 3D cartesian rays through
!!  one block for 2D cylindrical geometries. On exit, each ray has either:
!!
!!            i)  reached a different (yet unknown) 2D cylindrical block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!!  The 3-dimensional shape of the 2D cylindrical geometry is approximated as a 3D wedge,
!!  whose width extends in the y-direction. The y-direction constitutes the linear angular
!!  approximation. The smaller the width of the wedge, the more accurate the cylindrical
!!  representation becomes. Rays hitting on one of the wedge's y-directional boundaries
!!  stay in the same 2D cylindrical block and automatically jump to the other y-directional
!!  boundary, thus mimicking travel inside a cylindrical shell. The y-directional boundaries
!!  are defined by two lines 'y = m * x' and 'y = - m * x' with opposite slope, whose value
!!  depends on the wedges angle of aperture. The following picture shows the wedge with
!!  origin at the 2D cylindrical domain origin:
!!
!!
!!     y axis
!!       |                                                     *  <--- y = + m * x wedge boundary
!!       |                                            *        |
!!       |                                   *        |        |
!!       |                          *        |        |        |
!!       |                 *        |        |        |        |
!!       |        *        |        |        |        |         
!!       O--------|--------|--------|--------|--------|------- X ---------------------> R or x
!!       |        *        |        |        |        |         
!!       |                 *        |        |        |        |  X = 2D cylindrical
!!       |                          *        |        |        |      domain boundary
!!       |                                   *        |        |
!!       |                                            *        |
!!       |                                                     *  <--- y = - m * x wedge boundary
!!
!!
!!
!!  The z-direction is the same for either the 2D cylindrical or the 3D cartesian picture.
!!
!! ARGUMENTS
!!
!!  timeStep         : current timestep value
!!  rayFirst         : first ray index to be considered
!!  rayLast          : last ray index to be considered
!!  iminBlock        : minimum wedge i-index limit defining the interior block
!!  imaxBlock        : maximum wedge i-index limit defining the interior block
!!  jminBlock        : minimum wedge j-index limit defining the interior block
!!  jmaxBlock        : maximum wedge j-index limit defining the interior block
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  zminBlock        : minimum z-coordinate limit of the block
!!  zmaxBlock        : maximum z-coordinate limit of the block
!!  deltaInvX        : inverse of the wedge's x-dimension
!!  deltaInvZ        : inverse of the wedge's z-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  blockReflectMinZ : is the block boundary on the minimum z-side reflective ?
!!  blockReflectMaxZ : is the block boundary on the maximum z-side reflective ?
!!  wedgeEnergyDepot : array collecting the ray energy deposition for each wedge
!!
!! NOTES
!!
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation. 
!!  Inside each wedge, the paths of the rays are evaluated stepwise, using the
!!  bicubic expansions of the number electron density and electron temperature
!!  grid.
!!
!!***

subroutine ed_traceBlockRays2DCyl3D (timeStep,                          &
                                     rayFirst,  rayLast,                &
                                     iminBlock, imaxBlock,              &
                                     jminBlock, jmaxBlock,              &
                                     xminBlock, xmaxBlock,              &
                                     zminBlock, zmaxBlock,              &
                                     deltaInvX,                         &
                                     deltaInvZ,                         &
                                     blockReflectMinX,                  &
                                     blockReflectMaxX,                  &
                                     blockReflectMinZ,                  &
                                     blockReflectMaxZ,                  &
                                                       wedgeEnergyDepot ) 

  use EnergyDeposition_data,  ONLY : ed_Boltzmann,                   &
                                     ed_cellCubicNele,               &
                                     ed_cellCubicTele,               &
                                     ed_cellDensity,                 &
                                     ed_cellEdges,                   &
                                     ed_cellRayDifferentialStep,     &
                                     ed_cellVolume,                  &
                                     ed_cellWallThickness,           &
                                     ed_cellZbar,                    &
                                     ed_depoVarIsPerMass,            &
                                     ed_electronMass,                &
                                     ed_electronCharge,              &
                                     ed_energyOutTimestep,           &
                                     ed_infiniteTime,                &
                                     ed_laser3Din2DwedgeCosine,      &
                                     ed_laser3Din2DwedgeSine,        &
                                     ed_laser3Din2DwedgeSlope,       &
                                     ed_laserIOMaxNumberOfPositions, &
                                     ed_laserIOMaxNumberOfRays,      &
                                     ed_laserIONumberOfPositions,    &
                                     ed_laserIONumberOfRaysWritten,  &
                                     ed_laserIORayPositions,         &
                                     ed_laserIORayPower,             &
                                     ed_laserIORayFrequency,         &
                                     ed_laserIORayTags,              &
                                     ed_laserIOWrite,                &
                                     ed_rays,                        &
                                     ed_raysMovedIntoDomain,         &
                                     ed_rayZeroPower,                &
                                     ed_speedOfLightSquared,         &
                                     ed_unitRoundoff,                &
                                     ed_xminDomain,                  &
                                     ed_xmaxDomain,                  &
                                     ed_yminDomain,                  &
                                     ed_ymaxDomain

  use Driver_interface,       ONLY : Driver_abortFlash

  use ed_interface,           ONLY : ed_CoulombFactor,             &
                                     ed_inverseBremsstrahlungRate

  implicit none

#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  real,    intent (in)    :: timeStep
  integer, intent (in)    :: rayFirst,  rayLast   
  integer, intent (in)    :: iminBlock, imaxBlock
  integer, intent (in)    :: jminBlock, jmaxBlock
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: zminBlock, zmaxBlock
  real,    intent (in)    :: deltaInvX
  real,    intent (in)    :: deltaInvZ
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinZ
  logical, intent (in)    :: blockReflectMaxZ
  real,    intent (inout) :: wedgeEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock)

  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: blockWedgeMinX, blockWedgeMaxX
  logical :: blockWedgeMinZ, blockWedgeMaxZ
  logical :: crossX, crossZ
  logical :: impossibleRay
  logical :: inDomain, inBlock
  logical :: newWedge, newWedgeIJ
  logical :: onBlockFaceBoundary, onBlockWedgeBoundary
  logical :: rayOutOfBlock
  logical :: reflectX, reflectZ
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: velYeq0, velYgt0, velYlt0
  logical :: velZeq0, velZgt0, velZlt0
  logical :: wedgeFaceMinX, wedgeFaceMaxX
  logical :: wedgeFaceMinY, wedgeFaceMaxY
  logical :: wedgeFaceMinZ, wedgeFaceMaxZ
  logical :: writeRay

  integer :: face
  integer :: i,j
  integer :: ip,jp
  integer :: nRayWritePositions
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: a,b,c,d,q
  real    :: accX,accZ
  real    :: c2div1nc, c2div2nc, c2div4nc
  real    :: distToFaceMinX, distToFaceMaxX
  real    :: distToFaceMinY, distToFaceMaxY
  real    :: distToFaceMinZ, distToFaceMaxZ
  real    :: gradNeleX, gradNeleZ
  real    :: gradTeleX, gradTeleZ
  real    :: integral
  real    :: lnLambda
  real    :: minFaceDistance
  real    :: Nele, Tele
  real    :: newX, newY, newZ
  real    :: nu
  real    :: nudgeX, nudgeZ
  real    :: powerLossFactor
  real    :: R,S,U,W
  real    :: rayCritDens, rayCritDensInv
  real    :: rayPower
  real    :: rayX, rayY, rayZ
  real    :: time2Face
  real    :: stepTime, stepTimeHalf
  real    :: stepVelocity
  real    :: velX, velY, velZ
  real    :: wedgeCosine, wedgeSine
  real    :: wedgeDensity
  real    :: wedgeEnergy
  real    :: wedgeMass, wedgeMassInv
  real    :: wedgePower
  real    :: wedgeSlope
  real    :: wedgeVolume, wedgeVolumeInv
  real    :: wedgeWallThicknessHalf
  real    :: wedgeZbar
  real    :: x,z
  real    :: xminWedge, yminWedge, zminWedge
  real    :: xmaxWedge, ymaxWedge, zmaxWedge

  real, parameter :: GaussianRoot1 = 1.577350269189626         ! is 1 + 1/sqrt(3) for integration limits [0,1]
  real, parameter :: GaussianRoot2 = 4.226497308103742e-01     ! is 1 - 1/sqrt(3) for integration limits [0,1]

  real    :: aCoeff (1:6)
  real    :: bCoeff (1:6)
  real    :: cCoeff (1:6)
  real    :: xz     (1:16)
!
!
!     ...Define some variables.
!
!
  wedgeWallThicknessHalf = 0.5 * ed_cellWallThickness

  wedgeCosine = ed_laser3Din2DwedgeCosine
  wedgeSine   = ed_laser3Din2DwedgeSine
  wedgeSlope  = ed_laser3Din2DwedgeSlope
!
!
!     ...Outer (threaded) loop over all rays associated with the current 2D cylindrical block.
!
!
!$omp do schedule (dynamic)

  do ray = rayFirst , rayLast
        
     rayTag         = int (ed_rays (RAY_TAGS,ray))
     rayX           =      ed_rays (RAY_POSX,ray)
     rayY           =      ed_rays (RAY_POSY,ray)
     rayZ           =      ed_rays (RAY_POSZ,ray)
     velX           =      ed_rays (RAY_VELX,ray)
     velY           =      ed_rays (RAY_VELY,ray)
     velZ           =      ed_rays (RAY_VELZ,ray)
     rayPower       =      ed_rays (RAY_POWR,ray)
     rayCritDens    =      ed_rays (RAY_DENC,ray)
     rayCritDensInv = 1.0 / rayCritDens

     c2div1nc = ed_speedOfLightSquared * rayCritDensInv
     c2div4nc = 0.25 * c2div1nc
     c2div2nc = c2div4nc + c2div4nc
!
!
!     ...Decide, if we should write this ray out. If the case, start the writeout procedure.
!        If threading is done, the ray writing index must be protected from incrementation
!        by another thread and is saved in a local thread variable.
!
!
     writeRay = ed_laserIOWrite

     if (writeRay) then
         writeRay = mod (rayTag, ed_laserIORayFrequency) == 0
         if (writeRay) then
             !$omp critical (WriteRayIndex)
                   ed_laserIONumberOfRaysWritten = ed_laserIONumberOfRaysWritten + 1
                   rayWriteIndex = ed_laserIONumberOfRaysWritten
             !$omp end critical (WriteRayIndex)
         end if
     end if

     if(writeRay) then
        nRayWritePositions = 0
        ed_laserIORayTags (rayWriteIndex) = rayTag
     end if
!
!
!     ...Find the indices (i,j) of the initial wedge through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the four possible 2D cylindrical wedge faces the ray currently is. The
!        current position of the ray is such that it is not exactly on the block
!        boundary but nudged into the block. The distance to the closest wedge face
!        is always less than the cell wall thickness.
!
!
     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )
     j = jminBlock + int ( (rayZ - zminBlock) * deltaInvZ )

     blockWedgeMinX = (i == iminBlock)
     blockWedgeMaxX = (i == imaxBlock)
     blockWedgeMinZ = (j == jminBlock)
     blockWedgeMaxZ = (j == jmaxBlock)

     onBlockWedgeBoundary = (     blockWedgeMinX &
                             .or. blockWedgeMaxX &
                             .or. blockWedgeMinZ &
                             .or. blockWedgeMaxZ )

     if (.not.onBlockWedgeBoundary) then
          call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray not in a block boundary wedge')
     end if

     xminWedge = ed_cellEdges (i  ,1)
     xmaxWedge = ed_cellEdges (i+1,1)
     zminWedge = ed_cellEdges (j  ,2)
     zmaxWedge = ed_cellEdges (j+1,2)

     distToFaceMinX = abs (xminWedge - rayX)
     distToFaceMaxX = abs (xmaxWedge - rayX)
     distToFaceMinZ = abs (zminWedge - rayZ)
     distToFaceMaxZ = abs (zmaxWedge - rayZ)

     minFaceDistance = min (distToFaceMinX, distToFaceMaxX, &
                            distToFaceMinZ, distToFaceMaxZ)

     if (minFaceDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray too far away from block boundary wedge face')
     end if

     wedgeFaceMinX = (distToFaceMinX <= ed_cellWallThickness)
     wedgeFaceMaxX = (distToFaceMaxX <= ed_cellWallThickness)
     wedgeFaceMinZ = (distToFaceMinZ <= ed_cellWallThickness)
     wedgeFaceMaxZ = (distToFaceMaxZ <= ed_cellWallThickness)

     impossibleRay =     (wedgeFaceMinX .and. wedgeFaceMaxX) &
                    .or. (wedgeFaceMinZ .and. wedgeFaceMaxZ)

     if (impossibleRay) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray cannot be on two opposite wedge faces')
     end if

     blockFaceMinX = (blockWedgeMinX .and. wedgeFaceMinX)
     blockFaceMaxX = (blockWedgeMaxX .and. wedgeFaceMaxX)
     blockFaceMinZ = (blockWedgeMinZ .and. wedgeFaceMinZ)
     blockFaceMaxZ = (blockWedgeMaxZ .and. wedgeFaceMaxZ)

     onBlockFaceBoundary =     blockFaceMinX &
                          .or. blockFaceMaxX &
                          .or. blockFaceMinZ &
                          .or. blockFaceMaxZ

     if (.not.onBlockFaceBoundary) then
          call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray not on block face boundary')
     end if

     velXgt0 = (velX  > 0.0)
     velXeq0 = (velX == 0.0)
     velXlt0 = (velX  < 0.0)
     velYgt0 = (velY  > 0.0)
     velYeq0 = (velY == 0.0)
     velYlt0 = (velY  < 0.0)
     velZgt0 = (velZ  > 0.0)
     velZeq0 = (velZ == 0.0)
     velZlt0 = (velZ  < 0.0)

     stationaryRay = (velXeq0 .and. velYeq0 .and. velZeq0)

     if (stationaryRay) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: stationary ray at a block face boundary')
     end if

     rayOutOfBlock =     (blockFaceMinX .and. velXlt0) &
                    .or. (blockFaceMaxX .and. velXgt0) &
                    .or. (blockFaceMinZ .and. velZlt0) &
                    .or. (blockFaceMaxZ .and. velZgt0)

     if (rayOutOfBlock .and. ed_raysMovedIntoDomain) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray moves out of block')
     end if
!
!
!     ...Calculate the electron number density and the electron temperature as well as their
!        gradients (d/dx and d/dy) at the current ray position in the (i,j) wedge. These
!        are calculated using the corresponding bicubic expansion coefficients. If either of
!        these values become negative, the calculation must be stopped.
!
!        Algorithm implemented to evaluate the function values and their derivatives
!        ---------------------------------------------------------------------------
!
!        The bicubic expansion coefficients a (i,j) for each wedge contain all the info that
!        is needed. For a point (X,Z) inside the wedge with rescaled [0,1] x,z coordinates,
!        we have:
!
!
!                              3   3            i j
!                   F (X,Z) = sum sum  a (i,j) x z
!                             i=0 j=0
!
!                              3   3                i-1 j
!                     dF/dX = sum sum  i * a (i,j) x   z  / (wedge x-dimension)
!                             i=1 j=0
!
!                              3   3                i j-1
!                     dF/dZ = sum sum  j * a (i,j) x z    / (wedge z-dimension)
!                             i=0 j=1
!
!
!        Since the coordinate monomial parts are all common to these expressions, we first
!        construct the monomial array 'w', which will contain all needed monomials in the
!        same i,j order as the order for the a (i,j) is defined:
!
!                          i j
!                         x z   --->  placed into  ---> xz (i + 4j + 1)
!
!
!        The function values are then given by the scalar product between the monomial
!        and the bicubic expansion coefficients vectors. For the derivative evaluations,
!        only a subset of the monimials and the a (i,j)'s are needed. Below is a picture
!        to show which sections of the two arrays are needed for each individual type
!        of derivative.
!
!
!                       a (i,j)
!                                   --
!                                  |  |
!                          --     /|  |        <--- shown is a section of 4 consecutive
!                         |XX|  1/ |--|             boxes each containing an equal amount
!                         |XX|  /  |  |             of elements. The slope lines indicate
!                         |--| /  /|  |             which boxes are multiplied together
!                         |  |/ 2/ |--|             (element by element). The numerical
!                         |  |  /  |  |             value on the slope lines indicate the
!                         |--| /  /|  |             multiplying factor. Boxes marked XX
!                         |  |/ 3/ |--|             are not needed.
!                         |  |  /  |XX|
!                         |--| /   |XX|
!                         |  |/     --
!                         |  |
!                          --
!                                 xz (i,j)
!
!
!        dF/dX  -->  The picture shows a particular j-index section with i-index range
!                    from 0 to 3. Each box contains only 1 element.
!
!        dF/dZ  -->  The picture shows the entire vectors. Each box corresponds to
!                    a particular j-index and contains 4 elements with i-index
!                    range from 0 to 3.
!
!
!        Accumulate first the monomial array. This array will be used for both electron number
!        density and electron temperature evaluations.
!
!
     x = (rayX - xminWedge) * deltaInvX            ! rescaled [0,1] ray x coordinate
     z = (rayZ - zminWedge) * deltaInvZ            ! rescaled [0,1] ray z coordinate

     xz ( 1)    = 1.0                              !
     xz ( 2)    = x                                !
     xz ( 3)    = x * xz (2)                       !               i j
     xz ( 4)    = x * xz (3)                       ! the monomial x z  array
     xz ( 5: 8) = z * xz (1:4)                     !
     xz ( 9:12) = z * xz (5:8)                     !
     xz (13:16) = z * xz (9:12)                    !
!
!
!     ...Electron number density + gradients.
!
!
     Nele = sum (xz (1:16) * ed_cellCubicNele (1:16,i,j,1))

     a = sum ( xz (1:13:4) * ed_cellCubicNele (2:14:4,i,j,1))
     b = sum ( xz (2:14:4) * ed_cellCubicNele (3:15:4,i,j,1))
     c = sum ( xz (3:15:4) * ed_cellCubicNele (4:16:4,i,j,1))

     gradNeleX = (a + b + b + c + c + c) * deltaInvX

     a = sum (  xz ( 1: 4) * ed_cellCubicNele ( 5: 8,i,j,1))
     b = sum (  xz ( 5: 8) * ed_cellCubicNele ( 9:12,i,j,1))
     c = sum (  xz ( 9:12) * ed_cellCubicNele (13:16,i,j,1))

     gradNeleZ = (a + b + b + c + c + c) * deltaInvZ
!
!
!     ...Electron temperature + gradients.
!
!
     Tele = sum (xz (1:16) * ed_cellCubicTele (1:16,i,j,1))

     a = sum ( xz (1:13:4) * ed_cellCubicTele (2:14:4,i,j,1))
     b = sum ( xz (2:14:4) * ed_cellCubicTele (3:15:4,i,j,1))
     c = sum ( xz (3:15:4) * ed_cellCubicTele (4:16:4,i,j,1))

     gradTeleX = (a + b + b + c + c + c) * deltaInvX

     a = sum (  xz ( 1: 4) * ed_cellCubicTele ( 5: 8,i,j,1))
     b = sum (  xz ( 5: 8) * ed_cellCubicTele ( 9:12,i,j,1))
     c = sum (  xz ( 9:12) * ed_cellCubicTele (13:16,i,j,1))

     gradTeleZ = (a + b + b + c + c + c) * deltaInvZ
!
!
!     ...Check, if obtained values of electron number density and electron temeperature
!        make sense.
!
!
     if (Nele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Nele <= 0 for a wedge (initial)')
     end if

     if (Tele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Tele <= 0 for a wedge (initial)')
     end if
!
!
!     ...Get extra needed info about the initial wedge (i,j).
!
!
     wedgeZbar      = ed_cellZbar    (i,j,1)
     wedgeDensity   = ed_cellDensity (i,j,1)
     wedgeVolume    = ed_cellVolume  (i,j,1)
     wedgeVolumeInv = 1.0 / wedgeVolume
     wedgeMass      = wedgeDensity * wedgeVolume
     wedgeMassInv   = 1.0 / wedgeMass
!
!
!     ...We are ready to follow the ray's path through all the wedges of the current
!        block. The current wedge indices (i,j) and the previous wedge indices (ip,jp)
!        will be updated as we move through the block. In case a laser IO is performed
!        on the rays, store the initial ray IO data.
!
!
     if (writeRay) then
         nRayWritePositions = nRayWritePositions + 1
         if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
            ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
            ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayZ
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
         end if
     end if
!
!
!-------------------- Loop following ray through wedges in block --------------------------------------------
!
!
     do                                ! indefinite loop through the block wedges
                                       ! will be broken (exit) by the various conditions
                                       ! of the ray (no power, out of domain, etc)
!
!
!     ...Determine the correct stepping time, either for the ray to stay within the current
!        wedge (i,j) or to one of the 6 wedge walls. After this section, we have the correct
!        stepping time and the ray's new position. The ray is considered to stay in the same
!        wedge, if it is not found within the wedge wall. If it is in the wedge wall, the ray
!        is formally between two wedges, i.e. crossing to the neighboring wedge is considered
!        a possibility. Note, that crossing the y-directional boundaries of the current wedge
!        is also considered to be a criterion for the ray to enter a new wedge.
!
!        The stepping time is determined from the current ray velocity and the (user's) given
!        ray stepping differential. This ray stepping size is a fraction of the cell's dimension
!        and controlls accuracy of the ray's path. The current implementation has the stepping
!        velocity set equal to the maximum absolute ray velocity component. This results in
!        different ray distances travelled, depending on how the ray velocity is oriented in
!        space. Maximum ray distance travelled (i.e. equal to the ray stepping size) will occur,
!        if the ray velocity vector is parallel to one of the coordinate axis. Minimum ray distance
!        travelled will be 1/sqrt(2) of the ray stepping size (when all ray velocity components
!        are equal). The major cost of this approach is the one division by the stepping
!        velocity. If, on the other hand, we want to enforce always the ray distance travelled
!        to be equal to the ray stepping size, we would have an additional square root cost.
!
!
        stepVelocity = max (abs (velX), abs (velZ))
        stepTime     = ed_cellRayDifferentialStep / stepVelocity
        stepTimeHalf = 0.5 * stepTime

        accX = - c2div2nc * gradNeleX               ! acceleration in x-direction
        accZ = - c2div2nc * gradNeleZ               ! acceleration in z-direction

        newX = rayX + (velX + accX * stepTimeHalf) * stepTime
        newY = rayY + (velY                      ) * stepTime
        newZ = rayZ + (velZ + accZ * stepTimeHalf) * stepTime

        newWedge =      (newX <           xminWedge + wedgeWallThicknessHalf) &
                   .or. (newX >           xmaxWedge - wedgeWallThicknessHalf) &
                   .or. (newZ <           zminWedge + wedgeWallThicknessHalf) &
                   .or. (newZ >           zmaxWedge - wedgeWallThicknessHalf) &
                   .or. (newY < - wedgeSlope * newX + wedgeWallThicknessHalf) &    ! test current yminWedge
                   .or. (newY > + wedgeSlope * newX - wedgeWallThicknessHalf)      ! test current ymaxWedge

        if (newWedge) then

            aCoeff (1:2) =   0.5 * accX
            aCoeff (3:3) =   0.5 * accX * wedgeSlope
            aCoeff (4:4) = - 0.5 * accX * wedgeSlope
            aCoeff (5:6) =   0.5 * accZ

            bCoeff (1:2) = velX
            bCoeff (3:3) = velY + wedgeSlope * velX
            bCoeff (4:4) = velY - wedgeSlope * velX
            bCoeff (5:6) = velZ

            cCoeff (1) = rayX - xminWedge
            cCoeff (2) = rayX - xmaxWedge
            cCoeff (3) = rayY + wedgeSlope * rayX     ! this is current - yminWedge
            cCoeff (4) = rayY - wedgeSlope * rayX     ! this is current - ymaxWedge
            cCoeff (5) = rayZ - zminWedge
            cCoeff (6) = rayZ - zmaxWedge

            stepTime = ed_infiniteTime

            do face = 1,6

               a = aCoeff (face)
               b = bCoeff (face)
               c = cCoeff (face)

               if (a == 0.0) then
                   if (b /= 0.0) then
                       time2Face = - c / b
                       if (time2Face > 0.0 .and. time2Face < stepTime) stepTime = time2Face
                   end if
               else
                   if (c == 0.0) then
                       time2Face = - b / a
                       if (time2Face > 0.0 .and. time2Face < stepTime) stepTime = time2Face
                   else
                       d = b * b - 4.0 * a * c
                       if ((d > 0.0)) then
                            q = -0.5 * (b + sign (1.0,b) * sqrt (d))
                            time2Face = q / a
                            if (time2Face > 0.0 .and. time2Face < stepTime) stepTime = time2Face
                            time2Face = c / q
                            if (time2Face > 0.0 .and. time2Face < stepTime) stepTime = time2Face
                       else if (d == 0.0) then
                            time2Face = -0.5 * b / a
                            if (time2Face > 0.0 .and. time2Face < stepTime) stepTime = time2Face
                       end if
                   end if
               end if

            end do

            if (stepTime == ed_infiniteTime) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: infinite wedge crossing time')
            end if

            stepTimeHalf = 0.5 * stepTime

            rayX = rayX + (velX + accX * stepTimeHalf) * stepTime
            rayY = rayY + (velY                      ) * stepTime
            rayZ = rayZ + (velZ + accZ * stepTimeHalf) * stepTime

        else

            rayX = newX
            rayY = newY
            rayZ = newZ

        end if
!
!
!     ...We determined the stepping time of the ray to be such that either the ray stays within
!        the current wedge or hits a particular wedge face plane. Calculate the power deposition as
!        the ray travels during this time. This is done by evaluating an integral using the initial
!        velocities and gradients. The method of integral evaluation is Gaussian Quadrature with weight
!        function equal to 1. The associated orthogonal polynomials are the Legendre Polynomials.
!        If the remaining ray power is considered to have reached a 'zero' value, mark the ray as
!        nonexistent and exit the indefinite block loop.
!
!
        lnLambda = ed_CoulombFactor (wedgeZbar,         &
                                     ed_electronCharge, &
                                     ed_Boltzmann,      &
                                     Tele,              &
                                     Nele               )

        nu = ed_inverseBremsstrahlungRate (wedgeZbar,         &
                                           ed_electronCharge, &
                                           ed_electronMass,   &
                                           ed_Boltzmann,      &
                                           Tele,              &
                                           Nele,              &
                                           rayCritDens,       &
                                           lnLambda           )

        U =   (     velX * gradNeleX +      velZ * gradNeleZ) / Nele
        W =   (     velX * gradTeleX +      velZ * gradTeleZ) / Tele
        R = - (gradNeleX * gradNeleX + gradNeleZ * gradNeleZ) * c2div4nc / Nele
        S = - (gradNeleX * gradTeleX + gradNeleZ * gradTeleZ) * c2div4nc / Tele

        a = GaussianRoot1 * stepTimeHalf
        b = a * a
        c = 1.0 + U * a + R * b
        d = sqrt (1.0 + W * a + S * b)

        integral = c * c / (d * d * d)

        a = GaussianRoot2 * stepTimeHalf
        b = a * a
        c = 1.0 + U * a + R * b
        d = sqrt (1.0 + W * a + S * b)

        integral = integral + c * c / (d * d * d)
        integral = integral * nu * stepTimeHalf

        powerLossFactor = exp (-integral)
        wedgePower      = rayPower * (1.0 - powerLossFactor)
        wedgeEnergy     = wedgePower * timeStep

        if (ed_depoVarIsPerMass) then
            wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy * wedgeMassInv
        else
            wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy * wedgeVolumeInv
        end if

        rayPower = rayPower * powerLossFactor

        if (rayPower <= ed_rayZeroPower) then
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            exit
        end if
!
!
!     ...Calculate the new ray velocities. If the ray is stationary (no movement), mark the
!        ray as nonexistent and exit the block loop. In case a laser IO is performed on the
!        rays, store the current ray IO data.
!
!
        velX = velX + accX * stepTime
        velZ = velZ + accZ * stepTime

        velXgt0 = (velX  > 0.0)
        velXeq0 = (velX == 0.0)
        velXlt0 = (velX  < 0.0)
        velYgt0 = (velY  > 0.0)
        velYeq0 = (velY == 0.0)
        velYlt0 = (velY  < 0.0)
        velZgt0 = (velZ  > 0.0)
        velZeq0 = (velZ == 0.0)
        velZlt0 = (velZ  < 0.0)

        stationaryRay  = velXeq0 .and. velYeq0 .and. velZeq0

        if (stationaryRay) then
            write (*,*) ' stationary ray detected! Removing it from list ... '
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            exit
        end if

        if (writeRay) then
            nRayWritePositions = nRayWritePositions + 1
            if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
               ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
               ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayZ
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new wedge, we have to determine: 1) which new
!        wedge (i,j) it is and 2) the appropriate nudging values on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original wedge. After handling the logistics inside the following 'if'
!        statement, the new wedge indices i,j are either the old ones or new ones. The keyword
!        'newWedgeIJ' will indicate, if the indices i,j will change.
!
!
        newWedgeIJ = .false.

        if (newWedge) then

            ymaxWedge = + wedgeSlope * rayX
            yminWedge = - ymaxWedge

            distToFaceMinX = abs (xminWedge - rayX)
            distToFaceMaxX = abs (xmaxWedge - rayX)
            distToFaceMinY = abs (yminWedge - rayY)
            distToFaceMaxY = abs (ymaxWedge - rayY)
            distToFaceMinZ = abs (zminWedge - rayZ)
            distToFaceMaxZ = abs (zmaxWedge - rayZ)

            minFaceDistance = min (distToFaceMinX, distToFaceMaxX, &
                                   distToFaceMinY, distToFaceMaxY, &
                                   distToFaceMinZ, distToFaceMaxZ)

            if (minFaceDistance > ed_cellWallThickness) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray too far away from wedge face')
            end if

            wedgeFaceMinX = (distToFaceMinX <= wedgeWallThicknessHalf)
            wedgeFaceMaxX = (distToFaceMaxX <= wedgeWallThicknessHalf)
            wedgeFaceMinY = (distToFaceMinY <= wedgeWallThicknessHalf) .and. (rayY < 0.0)
            wedgeFaceMaxY = (distToFaceMaxY <= wedgeWallThicknessHalf) .and. (rayY > 0.0)
            wedgeFaceMinZ = (distToFaceMinZ <= wedgeWallThicknessHalf)
            wedgeFaceMaxZ = (distToFaceMaxZ <= wedgeWallThicknessHalf)

            crossX = .false.
            crossZ = .false.

            nudgeX = 0.0
            nudgeZ = 0.0

            ip = i
            jp = j

            if (wedgeFaceMinX) then

                rayX   = xminWedge
                nudgeX = + wedgeWallThicknessHalf

                if (velXlt0) then
                    i = i - 1
                    crossX = .true.
                end if

            else if (wedgeFaceMaxX) then

                rayX   = xmaxWedge
                nudgeX = - wedgeWallThicknessHalf

                if (velXgt0) then
                    i = i + 1
                    crossX = .true.
                end if

            end if

            if (wedgeFaceMinY) then

                if (wedgeSlope * velX < - velY) then
                    rayY = ymaxWedge
                    velX = velX * wedgeCosine - velY * wedgeSine
                    velY = velX * wedgeSine   + velY * wedgeCosine
                else
                    rayY = yminWedge
                end if

            else if (wedgeFaceMaxY) then

                if (wedgeSlope * velX < velY) then
                    rayY = yminWedge
                    velX = velY * wedgeSine   + velX * wedgeCosine
                    velY = velY * wedgeCosine - velX * wedgeSine
                else
                    rayY = ymaxWedge
                end if

            end if

            if (wedgeFaceMinZ) then

                rayZ   = zminWedge
                nudgeZ = + wedgeWallThicknessHalf

                if (velZlt0) then
                    j = j - 1
                    crossZ = .true.
                end if

            else if (wedgeFaceMaxZ) then

                rayZ   = zmaxWedge
                nudgeZ = - wedgeWallThicknessHalf

                if (velZgt0) then
                    j = j + 1
                    crossZ = .true.
                end if

            end if

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)
            blockFaceMinZ = (rayZ == zminBlock)
            blockFaceMaxZ = (rayZ == zmaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
            reflectZ =     (blockFaceMinZ .and. blockReflectMinZ .and. velZlt0) &
                      .or. (blockFaceMaxZ .and. blockReflectMaxZ .and. velZgt0)

            if (reflectX) then
                i = ip
                velX = - velX
                crossX = .false.
            end if

            if (reflectZ) then
                j = jp
                velZ = - velZ
                crossZ = .false.
            end if

            if (crossX) then
                nudgeX = (i - ip) * wedgeWallThicknessHalf
            end if

            if (crossZ) then
                nudgeZ = (j - jp) * wedgeWallThicknessHalf
            end if

            rayX = rayX + nudgeX
            rayZ = rayZ + nudgeZ

            newWedgeIJ = crossX .or. crossZ

        end if
!
!
!     ...We are now sure about the target wedge. Check, if the target wedge (i,j) is still within the block.
!        If it is, we check if this is a wedge with new i,j indices, in which case we update the i,j info.
!        Next we have to calculate the new electron number density and temperature as well as their gradients
!        at the point where the ray is located in the target wedge. If the target wedge is not within the
!        block, check if the ray coordinates are still within the defined domain. If not, store its latest
!        data and mark it as nonexistent. If the ray is still within the domain boundaries, exit the current
!        block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock)

        if (inBlock) then

            if (newWedgeIJ) then

                xminWedge    = ed_cellEdges (i  ,1)
                xmaxWedge    = ed_cellEdges (i+1,1)
                zminWedge    = ed_cellEdges (j  ,2)
                zmaxWedge    = ed_cellEdges (j+1,2)

                wedgeZbar      = ed_cellZbar    (i,j,1)
                wedgeDensity   = ed_cellDensity (i,j,1)
                wedgeVolume    = ed_cellVolume  (i,j,1)
                wedgeVolumeInv = 1.0 / wedgeVolume
                wedgeMass      = wedgeDensity * wedgeVolume
                wedgeMassInv   = 1.0 / wedgeMass

            end if

            x = (rayX - xminWedge) * deltaInvX
            z = (rayZ - zminWedge) * deltaInvZ

            xz ( 1)    = 1.0
            xz ( 2)    = x
            xz ( 3)    = x * xz (2)
            xz ( 4)    = x * xz (3)
            xz ( 5: 8) = z * xz (1:4)
            xz ( 9:12) = z * xz (5:8)
            xz (13:16) = z * xz (9:12)

            Nele = sum (xz (1:16) * ed_cellCubicNele (1:16,i,j,1))

            a = sum ( xz (1:13:4) * ed_cellCubicNele (2:14:4,i,j,1))
            b = sum ( xz (2:14:4) * ed_cellCubicNele (3:15:4,i,j,1))
            c = sum ( xz (3:15:4) * ed_cellCubicNele (4:16:4,i,j,1))

            gradNeleX = (a + b + b + c + c + c) * deltaInvX

            a = sum (  xz ( 1: 4) * ed_cellCubicNele ( 5: 8,i,j,1))
            b = sum (  xz ( 5: 8) * ed_cellCubicNele ( 9:12,i,j,1))
            c = sum (  xz ( 9:12) * ed_cellCubicNele (13:16,i,j,1))

            gradNeleZ = (a + b + b + c + c + c) * deltaInvZ

            Tele = sum (xz (1:16) * ed_cellCubicTele (1:16,i,j,1))

            a = sum ( xz (1:13:4) * ed_cellCubicTele (2:14:4,i,j,1))
            b = sum ( xz (2:14:4) * ed_cellCubicTele (3:15:4,i,j,1))
            c = sum ( xz (3:15:4) * ed_cellCubicTele (4:16:4,i,j,1))

            gradTeleX = (a + b + b + c + c + c) * deltaInvX

            a = sum (  xz ( 1: 4) * ed_cellCubicTele ( 5: 8,i,j,1))
            b = sum (  xz ( 5: 8) * ed_cellCubicTele ( 9:12,i,j,1))
            c = sum (  xz ( 9:12) * ed_cellCubicTele (13:16,i,j,1))

            gradTeleZ = (a + b + b + c + c + c) * deltaInvZ

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Nele <= 0 for a wedge (target)')
            end if

            if (Tele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Tele <= 0 for a wedge (target)')
            end if

        else

            ed_rays (RAY_POSX,ray) = rayX
            ed_rays (RAY_POSY,ray) = rayY
            ed_rays (RAY_POSZ,ray) = rayZ
            ed_rays (RAY_VELX,ray) = velX
            ed_rays (RAY_VELY,ray) = velY
            ed_rays (RAY_VELZ,ray) = velZ
            ed_rays (RAY_POWR,ray) = rayPower

            inDomain =      (rayX > ed_xminDomain) &
                      .and. (rayX < ed_xmaxDomain) &
                      .and. (rayZ > ed_yminDomain) &
                      .and. (rayZ < ed_ymaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX             ! undo the nudging, since
                 ed_rays (RAY_POSZ,ray) = rayZ - nudgeZ             ! it is not needed anymore
                 ed_rays (RAY_BLCK,ray) = real (RAY_OUTDOMAIN)
                 ed_energyOutTimeStep   = ed_energyOutTimeStep + rayPower * timeStep
            end if

            exit

        end if
!
!
!-------------------- End loop following ray through wedges in block --------------------------------------------
!
!
     end do
!
!
!     ...Check to see if we ran out of laser IO buffer space
!
!
     if (writeRay .and. (nRayWritePositions > ed_laserIOMaxNumberOfPositions) ) then
         print *, "[ed_traceBlockRays2DCyl3D] Ray ", ray, &
                  " ran out of IO buffer space (", nRayWritePositions, ")"
     end if
!
!
!     ...Consider next ray.
!
!
  end do
!$omp end do nowait
!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceBlockRays2DCyl3D
