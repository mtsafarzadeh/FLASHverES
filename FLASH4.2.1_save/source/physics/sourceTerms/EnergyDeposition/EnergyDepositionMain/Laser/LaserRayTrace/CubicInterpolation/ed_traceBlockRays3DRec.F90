!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_traceBlockRays3DRec
!!
!! NAME
!!
!!  ed_traceBlockRays3DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays3DRec (real    (in)    :: timeStep,
!!                               integer (in)    :: rayFirst
!!                               integer (in)    :: rayLast,
!!                               integer (in)    :: iminBlock,
!!                               integer (in)    :: imaxBlock,
!!                               integer (in)    :: jminBlock,
!!                               integer (in)    :: jmaxBlock,
!!                               integer (in)    :: kminBlock,
!!                               integer (in)    :: kmaxBlock,
!!                               real    (in)    :: xminBlock,
!!                               real    (in)    :: xmaxBlock,
!!                               real    (in)    :: yminBlock,
!!                               real    (in)    :: ymaxBlock,
!!                               real    (in)    :: zminBlock,
!!                               real    (in)    :: zmaxBlock,
!!                               real    (in)    :: deltaInvX,
!!                               real    (in)    :: deltaInvY,
!!                               real    (in)    :: deltaInvZ,
!!                               logical (in)    :: blockReflectMinX,
!!                               logical (in)    :: blockReflectMaxX,
!!                               logical (in)    :: blockReflectMinY,
!!                               logical (in)    :: blockReflectMaxY,
!!                               logical (in)    :: blockReflectMinZ,
!!                               logical (in)    :: blockReflectMaxZ,
!!                               real    (inout) :: cellEnergyDepot (:,:,:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active rays through one block
!!  for those geometries consisting formally of 3D rectangular grids (cartesian).
!!  On exit, each ray has either:
!!
!!            i)  reached a different (yet unknown) block
!!           ii)  has been absorbed -> is nonexistent
!!          iii)  has reached the domain boundary and exited -> is nonexistent.
!!
!! ARGUMENTS
!!
!!  timeStep         : current timestep value
!!  rayFirst         : first ray index to be considered
!!  rayLast          : last ray index to be considered
!!  iminBlock        : minimum cell i-index limit defining the interior block
!!  imaxBlock        : maximum cell i-index limit defining the interior block
!!  jminBlock        : minimum cell j-index limit defining the interior block
!!  jmaxBlock        : maximum cell j-index limit defining the interior block
!!  kminBlock        : minimum cell k-index limit defining the interior block
!!  kmaxBlock        : maximum cell k-index limit defining the interior block
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  yminBlock        : minimum y-coordinate limit of the block
!!  ymaxBlock        : maximum y-coordinate limit of the block
!!  zminBlock        : minimum z-coordinate limit of the block
!!  zmaxBlock        : maximum z-coordinate limit of the block
!!  deltaInvX        : inverse of the cell's x-dimension
!!  deltaInvY        : inverse of the cell's y-dimension
!!  deltaInvZ        : inverse of the cell's z-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  blockReflectMinY : is the block boundary on the minimum y-side reflective ?
!!  blockReflectMaxY : is the block boundary on the maximum y-side reflective ?
!!  blockReflectMinZ : is the block boundary on the minimum z-side reflective ?
!!  blockReflectMaxZ : is the block boundary on the maximum z-side reflective ?
!!  cellEnergyDepot  : array collecting the ray energy deposition for each cell
!!
!! NOTES
!!        
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation.
!!  Inside each cell, the paths of the rays are evaluated stepwise, using the
!!  tricubic expansions of the number electron density and electron temperature
!!  grid.
!!
!!***

subroutine ed_traceBlockRays3DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   jminBlock, jmaxBlock,              &
                                   kminBlock, kmaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   yminBlock, ymaxBlock,              &
                                   zminBlock, zmaxBlock,              &
                                   deltaInvX,                         &
                                   deltaInvY,                         &
                                   deltaInvZ,                         &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                   blockReflectMinY,                  &
                                   blockReflectMaxY,                  &
                                   blockReflectMinZ,                  &
                                   blockReflectMaxZ,                  &
                                                      cellEnergyDepot ) 

  use EnergyDeposition_data,  ONLY : ed_Boltzmann,                   &
                                     ed_cellCubicNele,               &
                                     ed_cellCubicTele,               &
                                     ed_cellDensity,                 &
                                     ed_cellEdges,                   &
                                     ed_cellRayDifferentialStep,     &
                                     ed_cellWallThickness,           &
                                     ed_cellZbar,                    &
                                     ed_depoVarIsPerMass,            &
                                     ed_electronMass,                &
                                     ed_electronCharge,              &
                                     ed_energyOutTimestep,           &
                                     ed_infiniteTime,                &
                                     ed_laserIOMaxNumberOfPositions, &
                                     ed_laserIOMaxNumberOfRays,      &
                                     ed_laserIONumberOfPositions,    &
                                     ed_laserIONumberOfRaysWritten,  &
                                     ed_laserIORayFrequency,         &
                                     ed_laserIORayPositions,         &
                                     ed_laserIORayPower,             &
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
                                     ed_ymaxDomain,                  &
                                     ed_zminDomain,                  &
                                     ed_zmaxDomain

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
  integer, intent (in)    :: kminBlock, kmaxBlock
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: yminBlock, ymaxBlock
  real,    intent (in)    :: zminBlock, zmaxBlock
  real,    intent (in)    :: deltaInvX
  real,    intent (in)    :: deltaInvY
  real,    intent (in)    :: deltaInvZ
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  logical, intent (in)    :: blockReflectMinY
  logical, intent (in)    :: blockReflectMaxY
  logical, intent (in)    :: blockReflectMinZ
  logical, intent (in)    :: blockReflectMaxZ
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock,jminBlock:jmaxBlock,kminBlock:kmaxBlock)

  logical :: blockCellMinX, blockCellMaxX
  logical :: blockCellMinY, blockCellMaxY
  logical :: blockCellMinZ, blockCellMaxZ
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: blockFaceMinY, blockFaceMaxY
  logical :: blockFaceMinZ, blockFaceMaxZ
  logical :: cellFaceMinX,  cellFaceMaxX
  logical :: cellFaceMinY,  cellFaceMaxY
  logical :: cellFaceMinZ,  cellFaceMaxZ
  logical :: crossX, crossY, crossZ
  logical :: impossibleRay
  logical :: inDomain, inBlock
  logical :: newCell
  logical :: onBlockFaceBoundary, onBlockCellBoundary
  logical :: rayOutOfBlock
  logical :: reflectX, reflectY, reflectZ
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: velYeq0, velYgt0, velYlt0
  logical :: velZeq0, velZgt0, velZlt0
  logical :: writeRay

  integer :: face
  integer :: i,j,k
  integer :: ip,jp,kp
  integer :: nRayWritePositions
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: a,b,c,d,q
  real    :: accX, accY, accZ
  real    :: c2div1nc, c2div2nc, c2div4nc
  real    :: cellDensity
  real    :: cellEnergy
  real    :: cellMass, cellMassInv
  real    :: cellPower
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: cellZbar
  real    :: distToFaceMinX, distToFaceMinY, distToFaceMinZ
  real    :: distToFaceMaxX, distToFaceMaxY, distToFaceMaxZ
  real    :: gradNeleX, gradNeleY, gradNeleZ
  real    :: gradTeleX, gradTeleY, gradTeleZ
  real    :: integral
  real    :: lnLambda
  real    :: minFaceDistance
  real    :: Nele, Tele
  real    :: newX, newY, newZ
  real    :: nu
  real    :: nudgeX, nudgeY, nudgeZ
  real    :: powerLossFactor
  real    :: R,S,U,W
  real    :: rayCritDens, rayCritDensInv
  real    :: rayPower
  real    :: rayX, rayY, rayZ
  real    :: stepTime, stepTimeHalf
  real    :: stepVelocity
  real    :: time2Face
  real    :: velX, velY, velZ
  real    :: x, y, z
  real    :: xminCell, yminCell, zminCell
  real    :: xmaxCell, ymaxCell, zmaxCell

  real, parameter :: GaussianRoot1 = 1.577350269189626         ! is 1 + 1/sqrt(3) for integration limits [0,1]
  real, parameter :: GaussianRoot2 = 4.226497308103742e-01     ! is 1 - 1/sqrt(3) for integration limits [0,1]

  real    :: aCoeff (1:6)
  real    :: bCoeff (1:6)
  real    :: cCoeff (1:6)
  real    :: xyz    (1:64)
!
!
!     ...Define some variables. Fix the cell volume, which is independent of the cell location
!        inside the block.
!
!
  cellVolume            = 1.0 / (deltaInvX * deltaInvY * deltaInvZ)
  cellVolumeInv         = 1.0 / cellVolume
  cellWallThicknessHalf = 0.5 * ed_cellWallThickness
!
!
!     ...Outer (threaded) loop over all rays associated with the current block.
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
                   writeRay = ed_laserIONumberOfRaysWritten < ed_laserIOMaxNumberOfRays
                   if (writeRay) then
                       ed_laserIONumberOfRaysWritten = ed_laserIONumberOfRaysWritten + 1
                       rayWriteIndex = ed_laserIONumberOfRaysWritten
                   end if
             !$omp end critical (WriteRayIndex)
         end if
     end if

     if(writeRay) then
        nRayWritePositions = 0
        ed_laserIORayTags (rayWriteIndex) = rayTag
     end if
!
!
!     ...Find the indices (i,j,k) of the initial cell through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the four possible faces the ray currently is. The current position of the
!        ray is such that it is not exactly on the block boundary but nudged into
!        the block. The distance to the closest cell face is always less than
!        the cell wall thickness.
!
!
     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )
     j = jminBlock + int ( (rayY - yminBlock) * deltaInvY )
     k = kminBlock + int ( (rayZ - zminBlock) * deltaInvZ )

     blockCellMinX = (i == iminBlock)
     blockCellMaxX = (i == imaxBlock)
     blockCellMinY = (j == jminBlock)
     blockCellMaxY = (j == jmaxBlock)
     blockCellMinZ = (k == kminBlock)
     blockCellMaxZ = (k == kmaxBlock)

     onBlockCellBoundary = (     blockCellMinX &
                            .or. blockCellMaxX &
                            .or. blockCellMinY &
                            .or. blockCellMaxY &
                            .or. blockCellMinZ &
                            .or. blockCellMaxZ )

     if (.not.onBlockCellBoundary) then
          call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray not in a block boundary cell')
     end if

     xminCell = ed_cellEdges (i  ,1)
     xmaxCell = ed_cellEdges (i+1,1)
     yminCell = ed_cellEdges (j  ,2)
     ymaxCell = ed_cellEdges (j+1,2)
     zminCell = ed_cellEdges (k  ,3)
     zmaxCell = ed_cellEdges (k+1,3)

     distToFaceMinX = abs (xminCell - rayX)
     distToFaceMaxX = abs (xmaxCell - rayX)
     distToFaceMinY = abs (yminCell - rayY)
     distToFaceMaxY = abs (ymaxCell - rayY)
     distToFaceMinZ = abs (zminCell - rayZ)
     distToFaceMaxZ = abs (zmaxCell - rayZ)

     minFaceDistance = min (distToFaceMinX, distToFaceMaxX, &
                            distToFaceMinY, distToFaceMaxY, &
                            distToFaceMinZ, distToFaceMaxZ)

     if (minFaceDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray too far away from block boundary cell face')
     end if

     cellFaceMinX = (distToFaceMinX <= ed_cellWallThickness)
     cellFaceMaxX = (distToFaceMaxX <= ed_cellWallThickness)
     cellFaceMinY = (distToFaceMinY <= ed_cellWallThickness)
     cellFaceMaxY = (distToFaceMaxY <= ed_cellWallThickness)
     cellFaceMinZ = (distToFaceMinZ <= ed_cellWallThickness)
     cellFaceMaxZ = (distToFaceMaxZ <= ed_cellWallThickness)

     impossibleRay =     (cellFaceMinX .and. cellFaceMaxX) &
                    .or. (cellFaceMinY .and. cellFaceMaxY) &
                    .or. (cellFaceMinZ .and. cellFaceMaxZ)

     if (impossibleRay) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray cannot be on two opposite cell faces')
     end if

     blockFaceMinX = (blockCellMinX .and. cellFaceMinX)
     blockFaceMaxX = (blockCellMaxX .and. cellFaceMaxX)
     blockFaceMinY = (blockCellMinY .and. cellFaceMinY)
     blockFaceMaxY = (blockCellMaxY .and. cellFaceMaxY)
     blockFaceMinZ = (blockCellMinZ .and. cellFaceMinZ)
     blockFaceMaxZ = (blockCellMaxZ .and. cellFaceMaxZ)

     onBlockFaceBoundary =     blockFaceMinX &
                          .or. blockFaceMaxX &
                          .or. blockFaceMinY &
                          .or. blockFaceMaxY &
                          .or. blockFaceMinZ &
                          .or. blockFaceMaxZ

     if (.not.onBlockFaceBoundary) then
          call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray not on block face boundary')
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
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: stationary ray at a block face boundary')
     end if

     rayOutOfBlock =     (blockFaceMinX .and. velXlt0) &
                    .or. (blockFaceMaxX .and. velXgt0) &
                    .or. (blockFaceMinY .and. velYlt0) &
                    .or. (blockFaceMaxY .and. velYgt0) &
                    .or. (blockFaceMinZ .and. velZlt0) &
                    .or. (blockFaceMaxZ .and. velZgt0)

     if (rayOutOfBlock .and. ed_raysMovedIntoDomain) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray moves out of block')
     end if
!
!
!     ...Calculate the electron number density and the electron temperature as well as their
!        gradients (d/dx, d/dy and d/dz) at the current ray position in the (i,j,k) cell. These
!        are calculated using the corresponding tricubic expansion coefficients. If either of
!        these values become negative, the calculation must be stopped.
!
!        Algorithm implemented to evaluate the function values and their derivatives
!        ---------------------------------------------------------------------------
!
!        The tricubic expansion coefficients a (i,j,k) for each cell contain all the info that
!        is needed. For a point (X,Y,Z) inside the cell with rescaled [0,1] x,y,z coordinates,
!        we have:
!
!
!                         3   3   3              i j k
!            F (X,Y,Z) = sum sum sum  a (i,j,k) x y z
!                        i=0 j=0 k=0
!
!                         3   3   3                  i-1 j k
!                dF/dX = sum sum sum  i * a (i,j,k) x   y z  / (cell x-dimension)
!                        i=1 j=0 k=0
!
!                         3   3   3                  i j-1 k
!                dF/dY = sum sum sum  j * a (i,j,k) x y   z  / (cell y-dimension)
!                        i=0 j=1 k=0
!
!                         3   3   3                  i j k-1
!                dF/dZ = sum sum sum  k * a (i,j,k) x y z    / (cell z-dimension)
!                        i=0 j=0 k=1
!
!
!        Since the coordinate monomial parts are all common to these expressions, we first
!        construct the monomial array 'w', which will contain all needed monomials in the
!        same i,j,k order as the order for the a (i,j,k) is defined:
!
!                       i j k
!                      x y z   --->  placed into  ---> xyz (i + 4j + 16k + 1)
!
!
!        The function values are then given by the scalar product between the monomial
!        and the tricubic expansion coefficients vectors. For the derivative evaluations,
!        only a subset of the monimials and the a (i,j,k)'s are needed. Below is a picture
!        to show which sections of the two arrays are needed for each individual type
!        of derivative.
!
!
!                      a (i,j,k)
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
!                                 xyz (i,j,k)
!
!
!        dF/dX  -->  The picture shows a particular jk-index pair section with
!                    i-index range from 0 to 3. Each box contains only 1 element.
!
!        dF/dY  -->  The picture shows a particular k-index section with i- and
!                    j-index ranges from 0 to 3. Each box corresponds to a particular
!                    j-index and contains 4 element with i-index range from 0 to 3.
!
!        dF/dZ  -->  The picture shows the entire vectors. Each box corresponds to
!                    a particular k-index and contains 16 elements with i- and j-index
!                    ranges from 0 to 3.
!
!
!        Accumulate first the monomial array. This array will be used for both electron number
!        density and electron temperature evaluations.
!
!
     x = (rayX - xminCell) * deltaInvX            ! rescaled [0,1] ray x coordinate
     y = (rayY - yminCell) * deltaInvY            ! rescaled [0,1] ray y coordinate
     z = (rayZ - zminCell) * deltaInvZ            ! rescaled [0,1] ray z coordinate

     xyz ( 1)    = 1.0                            !
     xyz ( 2)    = x                              !
     xyz ( 3)    = x * xyz ( 2)                   !
     xyz ( 4)    = x * xyz ( 3)                   !               i j k
     xyz ( 5: 8) = y * xyz ( 1:4)                 ! the monomial x y z  array
     xyz ( 9:12) = y * xyz ( 5:8)                 !
     xyz (13:16) = y * xyz ( 9:12)                !
     xyz (17:32) = z * xyz ( 1:16)                !
     xyz (33:48) = z * xyz (17:32)                !
     xyz (49:64) = z * xyz (33:48)                !
!
!
!     ...Electron number density + gradients.
!
!
     Nele = sum (xyz (1:64) * ed_cellCubicNele (1:64,i,j,k))

     a = sum ( xyz (1:61:4) * ed_cellCubicNele (2:62:4,i,j,k))
     b = sum ( xyz (2:62:4) * ed_cellCubicNele (3:63:4,i,j,k))
     c = sum ( xyz (3:63:4) * ed_cellCubicNele (4:64:4,i,j,k))

     gradNeleX = (a + b + b + c + c + c) * deltaInvX

     a = sum (  xyz ( 1: 4) * ed_cellCubicNele ( 5: 8,i,j,k)  &
              + xyz (17:20) * ed_cellCubicNele (21:24,i,j,k)  &
              + xyz (33:36) * ed_cellCubicNele (37:40,i,j,k)  &
              + xyz (49:52) * ed_cellCubicNele (53:56,i,j,k))

     b = sum (  xyz ( 5: 8) * ed_cellCubicNele ( 9:12,i,j,k)  &
              + xyz (21:24) * ed_cellCubicNele (25:28,i,j,k)  &
              + xyz (37:40) * ed_cellCubicNele (41:44,i,j,k)  &
              + xyz (53:56) * ed_cellCubicNele (57:60,i,j,k))

     c = sum (  xyz ( 9:12) * ed_cellCubicNele (13:16,i,j,k)  &
              + xyz (25:28) * ed_cellCubicNele (29:32,i,j,k)  &
              + xyz (41:44) * ed_cellCubicNele (45:48,i,j,k)  &
              + xyz (57:60) * ed_cellCubicNele (61:64,i,j,k))

     gradNeleY = (a + b + b + c + c + c) * deltaInvY

     a = sum (  xyz ( 1:16) * ed_cellCubicNele (17:32,i,j,k))
     b = sum (  xyz (17:32) * ed_cellCubicNele (33:48,i,j,k))
     c = sum (  xyz (33:48) * ed_cellCubicNele (49:64,i,j,k))

     gradNeleZ = (a + b + b + c + c + c) * deltaInvZ
!
!
!     ...Electron temperature + gradients.
!
!
     Tele = sum (xyz (1:64) * ed_cellCubicTele (1:64,i,j,k))

     a = sum ( xyz (1:61:4) * ed_cellCubicTele (2:62:4,i,j,k))
     b = sum ( xyz (2:62:4) * ed_cellCubicTele (3:63:4,i,j,k))
     c = sum ( xyz (3:63:4) * ed_cellCubicTele (4:64:4,i,j,k))

     gradTeleX = (a + b + b + c + c + c) * deltaInvX

     a = sum (  xyz ( 1: 4) * ed_cellCubicTele ( 5: 8,i,j,k)  &
              + xyz (17:20) * ed_cellCubicTele (21:24,i,j,k)  &
              + xyz (33:36) * ed_cellCubicTele (37:40,i,j,k)  &
              + xyz (49:52) * ed_cellCubicTele (53:56,i,j,k))

     b = sum (  xyz ( 5: 8) * ed_cellCubicTele ( 9:12,i,j,k)  &
              + xyz (21:24) * ed_cellCubicTele (25:28,i,j,k)  &
              + xyz (37:40) * ed_cellCubicTele (41:44,i,j,k)  &
              + xyz (53:56) * ed_cellCubicTele (57:60,i,j,k))

     c = sum (  xyz ( 9:12) * ed_cellCubicTele (13:16,i,j,k)  &
              + xyz (25:28) * ed_cellCubicTele (29:32,i,j,k)  &
              + xyz (41:44) * ed_cellCubicTele (45:48,i,j,k)  &
              + xyz (57:60) * ed_cellCubicTele (61:64,i,j,k))

     gradTeleY = (a + b + b + c + c + c) * deltaInvY

     a = sum (  xyz ( 1:16) * ed_cellCubicTele (17:32,i,j,k))
     b = sum (  xyz (17:32) * ed_cellCubicTele (33:48,i,j,k))
     c = sum (  xyz (33:48) * ed_cellCubicTele (49:64,i,j,k))

     gradTeleZ = (a + b + b + c + c + c) * deltaInvZ
!
!
!     ...Check, if obtained values of electron number density and electron temeperature
!        make sense.
!
!
     if (Nele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Nele <= 0 for a cell (initial)')
     end if

     if (Tele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Tele <= 0 for a cell (initial)')
     end if
!
!
!     ...Get extra needed info about the initial cell (i,j,k).
!
!
     cellZbar    = ed_cellZbar    (i,j,k)
     cellDensity = ed_cellDensity (i,j,k)
     cellMass    = cellDensity * cellVolume
     cellMassInv = 1.0 / cellMass
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. At this point, the ray is in the initial cell (i,j,k) and is ready to move
!        through the block. In case a laser IO is performed on the rays, store the initial
!        ray IO data.
!
!
     if (writeRay) then
         nRayWritePositions = nRayWritePositions + 1
         if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
            ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
            ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayY
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = rayZ
         end if
     end if
!
!
!-------------------- Loop following ray through cells in block --------------------------------------------
!
!
     do                                ! indefinite loop through the block cells
                                       ! will be broken (exit) by the various conditions
                                       ! of the ray (no power, out of domain, etc)
!
!
!     ...Determine the correct stepping time, either for the ray to stay within the current
!        cell (i,j,k) or to one of the 6 cell walls. After this section, we have the correct
!        stepping time and the ray's new position. The ray is considered to stay in the same
!        cell, if it is not found within the cell wall. If it is in the cell wall, the ray
!        is formally between two cells, i.e. crossing to the neighboring cell is considered
!        a possibility.
!
!        The stepping time is determined from the current ray velocity and the (user's) given
!        ray stepping differential. This ray stepping size is a fraction of the cell's dimension
!        and controlls accuracy of the ray's path. The current implementation has the stepping
!        velocity set equal to the maximum absolute ray velocity component. This results in
!        different ray distances travelled, depending on how the ray velocity is oriented in
!        space. Maximum ray distance travelled (i.e. equal to the ray stepping size) will occur,
!        if the ray velocity vector is parallel to one of the coordinate axis. Minimum ray distance
!        travelled will be 1/sqrt(3) of the ray stepping size (when all ray velocity components
!        are equal). The major cost of this approach is the one division by the stepping
!        velocity. If, on the other hand, we want to enforce always the ray distance travelled
!        to be equal to the ray stepping size, we would have an additional square root cost.
!
!
        stepVelocity = max (abs (velX), abs (velY), abs (velZ))
        stepTime     = ed_cellRayDifferentialStep / stepVelocity
        stepTimeHalf = 0.5 * stepTime

        accX = - c2div2nc * gradNeleX               ! acceleration in x-direction
        accY = - c2div2nc * gradNeleY               ! acceleration in y-direction
        accZ = - c2div2nc * gradNeleZ               ! acceleration in z-direction

        newX = rayX + (velX + accX * stepTimeHalf) * stepTime
        newY = rayY + (velY + accY * stepTimeHalf) * stepTime
        newZ = rayZ + (velZ + accZ * stepTimeHalf) * stepTime



!    write (*,*) ' rayX     = ',rayX
!    write (*,*) ' rayY     = ',rayY
!    write (*,*) ' rayZ     = ',rayZ
!    write (*,*) ' velX     = ',velX
!    write (*,*) ' velY     = ',velY
!    write (*,*) ' velZ     = ',velZ
!    write (*,*) ' accX     = ',accX
!    write (*,*) ' accY     = ',accY
!    write (*,*) ' accZ     = ',accZ
!    write (*,*) ' newX     = ',newX
!    write (*,*) ' newY     = ',newY
!    write (*,*) ' newZ     = ',newZ
!    write (*,*) ' xminCell = ',xminCell
!    write (*,*) ' xmaxCell = ',xmaxCell
!    write (*,*) ' yminCell = ',yminCell
!    write (*,*) ' ymaxCell = ',ymaxCell
!    write (*,*) ' zminCell = ',zminCell
!    write (*,*) ' zmaxCell = ',zmaxCell
!    write (*,*) ' --------------- '



        newCell =      (newX < xminCell) &
                  .or. (newX > xmaxCell) &
                  .or. (newY < yminCell) &
                  .or. (newY > ymaxCell) &
                  .or. (newZ < zminCell) &
                  .or. (newZ > zmaxCell)

        if (newCell) then

            aCoeff (1:2) = 0.5 * accX
            aCoeff (3:4) = 0.5 * accY
            aCoeff (5:6) = 0.5 * accZ
            bCoeff (1:2) = velX
            bCoeff (3:4) = velY
            bCoeff (5:6) = velZ

            cCoeff (1) = rayX - xminCell
            cCoeff (2) = rayX - xmaxCell
            cCoeff (3) = rayY - yminCell
            cCoeff (4) = rayY - ymaxCell
            cCoeff (5) = rayZ - zminCell
            cCoeff (6) = rayZ - zmaxCell

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
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: infinite crossing time for a cell')
            end if

            stepTimeHalf = 0.5 * stepTime

            rayX = rayX + (velX + accX * stepTimeHalf) * stepTime
            rayY = rayY + (velY + accY * stepTimeHalf) * stepTime
            rayZ = rayZ + (velZ + accZ * stepTimeHalf) * stepTime

        else

            rayX = newX
            rayY = newY
            rayZ = newZ

        end if



!    write (*,*) ' rayX     = ',rayX
!    write (*,*) ' rayY     = ',rayY
!    write (*,*) ' rayZ     = ',rayZ
!
!
!     ...We determined the stepping time of the ray to be such that either the ray stays within
!        the current cell or hits a particular cell face plane. Calculate the power deposition as
!        the ray travels during this time. This is done by evaluating an integral using the initial
!        velocities and gradients. The method of integral evaluation is Gaussian Quadrature with weight
!        function equal to 1. The associated orthogonal polynomials are the Legendre Polynomials.
!        If the remaining ray power is considered to have reached a 'zero' value, mark the ray as
!        nonexistent and exit the indefinite block loop.
!
!
        lnLambda = ed_CoulombFactor (cellZbar,          &
                                     ed_electronCharge, &
                                     ed_Boltzmann,      &
                                     Tele,              &
                                     Nele               )

        nu = ed_inverseBremsstrahlungRate (cellZbar,          &
                                           ed_electronCharge, &
                                           ed_electronMass,   &
                                           ed_Boltzmann,      &
                                           Tele,              &
                                           Nele,              &
                                           rayCritDens,       &
                                           lnLambda           )

        U =   (     velX * gradNeleX +      velY * gradNeleY +      velZ * gradNeleZ) / Nele
        W =   (     velX * gradTeleX +      velY * gradTeleY +      velZ * gradTeleZ) / Tele
        R = - (gradNeleX * gradNeleX + gradNeleY * gradNeleY + gradNeleZ * gradNeleZ) * c2div4nc / Nele
        S = - (gradNeleX * gradTeleX + gradNeleY * gradTeleY + gradNeleZ * gradTeleZ) * c2div4nc / Tele

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
        cellPower       = rayPower * (1.0 - powerLossFactor)
        cellEnergy      = cellPower * timeStep

        if (ed_depoVarIsPerMass) then
            cellEnergyDepot (i,j,k) = cellEnergyDepot (i,j,k) + cellEnergy * cellMassInv
        else
            cellEnergyDepot (i,j,k) = cellEnergyDepot (i,j,k) + cellEnergy * cellVolumeInv
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
        velY = velY + accY * stepTime
        velZ = velZ + accZ * stepTime



!    write (*,*) ' velX     = ',velX
!    write (*,*) ' velY     = ',velY
!    write (*,*) ' velZ     = ',velZ



        velXeq0 = (velX == 0.0)
        velYeq0 = (velY == 0.0)
        velZeq0 = (velZ == 0.0)
 
        stationaryRay = velXeq0 .and. velYeq0 .and. velZeq0

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
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayY
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = rayZ
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new cell, we have to determine: 1) which new
!        cell (i,j,k) it is and 2) the appropriate nudging values on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original cell. After handling the logistics inside the following 'if'
!        statement, the new cell indices i,j,k are either the old ones or new ones.
!
!
        if (newCell) then

            distToFaceMinX = abs (xminCell - rayX)
            distToFaceMaxX = abs (xmaxCell - rayX)
            distToFaceMinY = abs (yminCell - rayY)
            distToFaceMaxY = abs (ymaxCell - rayY)
            distToFaceMinZ = abs (zminCell - rayZ)
            distToFaceMaxZ = abs (zmaxCell - rayZ)

            minFaceDistance = min (distToFaceMinX, distToFaceMaxX, &
                                   distToFaceMinY, distToFaceMaxY, &
                                   distToFaceMinZ, distToFaceMaxZ)

            if (minFaceDistance > ed_cellWallThickness) then
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: ray to far away from cell face')
            end if

            cellFaceMinX = (distToFaceMinX <= ed_cellWallThickness)
            cellFaceMaxX = (distToFaceMaxX <= ed_cellWallThickness)
            cellFaceMinY = (distToFaceMinY <= ed_cellWallThickness)
            cellFaceMaxY = (distToFaceMaxY <= ed_cellWallThickness)
            cellFaceMinZ = (distToFaceMinZ <= ed_cellWallThickness)
            cellFaceMaxZ = (distToFaceMaxZ <= ed_cellWallThickness)

            velXgt0 = (velX  > 0.0)
            velXlt0 = (velX  < 0.0)
            velYgt0 = (velY  > 0.0)
            velYlt0 = (velY  < 0.0)
            velZgt0 = (velZ  > 0.0)
            velZlt0 = (velZ  < 0.0)

            crossX = .false.
            crossY = .false.
            crossZ = .false.

            nudgeX = 0.0
            nudgeY = 0.0
            nudgeZ = 0.0

            ip = i
            jp = j
            kp = k

            if (cellFaceMinX) then

                rayX   = xminCell
                nudgeX = + cellWallThicknessHalf

                if (velXlt0) then
                    i = i - 1
                    crossX = .true.
                end if

            else if (cellFaceMaxX) then

                rayX   = xmaxCell
                nudgeX = - cellWallThicknessHalf

                if (velXgt0) then
                    i = i + 1
                    crossX = .true.
                end if

            end if

            if (cellFaceMinY) then

                rayY   = yminCell
                nudgeY = + cellWallThicknessHalf

                if (velYlt0) then
                    j = j - 1
                    crossY = .true.
                end if

            else if (cellFaceMaxY) then

                rayY   = ymaxCell
                nudgeY = - cellWallThicknessHalf

                if (velYgt0) then
                    j = j + 1
                    crossY = .true.
                end if

            end if

            if (cellFaceMinZ) then

                rayZ   = zminCell
                nudgeZ = + cellWallThicknessHalf

                if (velZlt0) then
                    k = k - 1
                    crossZ = .true.
                end if

            else if (cellFaceMaxZ) then

                rayZ   = zmaxCell
                nudgeZ = - cellWallThicknessHalf

                if (velZgt0) then
                    k = k + 1
                    crossZ = .true.
                end if

            end if

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)
            blockFaceMinY = (rayY == yminBlock)
            blockFaceMaxY = (rayY == ymaxBlock)
            blockFaceMinZ = (rayZ == zminBlock)
            blockFaceMaxZ = (rayZ == zmaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
            reflectY =     (blockFaceMinY .and. blockReflectMinY .and. velYlt0) &
                      .or. (blockFaceMaxY .and. blockReflectMaxY .and. velYgt0)
            reflectZ =     (blockFaceMinZ .and. blockReflectMinZ .and. velZlt0) &
                      .or. (blockFaceMaxZ .and. blockReflectMaxZ .and. velZgt0)

            if (reflectX) then
                i = ip
                velX = - velX
                crossX = .false.
            end if

            if (reflectY) then
                j = jp
                velY = - velY
                crossY = .false.
            end if

            if (reflectZ) then
                k = kp
                velZ = - velZ
                crossZ = .false.
            end if

            if (crossX) then
                nudgeX = (i - ip) * cellWallThicknessHalf
            end if

            if (crossY) then
                nudgeY = (j - jp) * cellWallThicknessHalf
            end if

            if (crossZ) then
                nudgeZ = (k - kp) * cellWallThicknessHalf
            end if

            rayX = rayX + nudgeX
            rayY = rayY + nudgeY
            rayZ = rayZ + nudgeZ

            newCell = crossX .or. crossY .or. crossZ


!    write (*,*) ' cellFaceMinX = ',cellFaceMinX
!    write (*,*) ' cellFaceMaxX = ',cellFaceMaxX
!    write (*,*) ' cellFaceMinY = ',cellFaceMinY
!    write (*,*) ' cellFaceMaxY = ',cellFaceMaxY
!    write (*,*) ' cellFaceMinZ = ',cellFaceMinZ
!    write (*,*) ' cellFaceMaxZ = ',cellFaceMaxZ
!    write (*,*) ' crossX       = ',crossX
!    write (*,*) ' crossY       = ',crossY
!    write (*,*) ' crossZ       = ',crossZ
!    write (*,*) ' rayX         = ',rayX
!    write (*,*) ' rayY         = ',rayY
!    write (*,*) ' rayZ         = ',rayZ


        end if

!    write (*,*) ' --------------- '
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i,j,k) is still within the block.
!        If it is, we check if this is a new cell, in which case we update the cell info. Next we have to
!        calculate the new electron number density and temperature as well as their gradients at the point
!        where the ray is located in the target cell. If the target cell is not within the block, check if
!        the ray coordinates are still within the defined domain. If not, store its latest data and mark
!        it as nonexistent. If the ray is still within the domain boundaries, exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock) &
                 .and. (k >= kminBlock) &
                 .and. (k <= kmaxBlock)

        if (inBlock) then

            if (newCell) then

                xminCell    = ed_cellEdges   (i  ,1)
                xmaxCell    = ed_cellEdges   (i+1,1)
                yminCell    = ed_cellEdges   (j  ,2)
                ymaxCell    = ed_cellEdges   (j+1,2)
                zminCell    = ed_cellEdges   (k  ,3)
                zmaxCell    = ed_cellEdges   (k+1,3)

                cellZbar    = ed_cellZbar    (i,j,k)
                cellDensity = ed_cellDensity (i,j,k)
                cellMass    = cellDensity * cellVolume
                cellMassInv = 1.0 / cellMass

            end if

            x = (rayX - xminCell) * deltaInvX
            y = (rayY - yminCell) * deltaInvY
            z = (rayZ - zminCell) * deltaInvZ

            xyz ( 1)    = 1.0
            xyz ( 2)    = x
            xyz ( 3)    = x * xyz ( 2)
            xyz ( 4)    = x * xyz ( 3)
            xyz ( 5: 8) = y * xyz ( 1:4)
            xyz ( 9:12) = y * xyz ( 5:8)
            xyz (13:16) = y * xyz ( 9:12)
            xyz (17:32) = z * xyz ( 1:16)
            xyz (33:48) = z * xyz (17:32)
            xyz (49:64) = z * xyz (33:48)

            Nele = sum (xyz (1:64) * ed_cellCubicNele (1:64,i,j,k))

            a = sum ( xyz (1:61:4) * ed_cellCubicNele (2:62:4,i,j,k))
            b = sum ( xyz (2:62:4) * ed_cellCubicNele (3:63:4,i,j,k))
            c = sum ( xyz (3:63:4) * ed_cellCubicNele (4:64:4,i,j,k))

            gradNeleX = (a + b + b + c + c + c) * deltaInvX

            a = sum (  xyz ( 1: 4) * ed_cellCubicNele ( 5: 8,i,j,k)  &
                     + xyz (17:20) * ed_cellCubicNele (21:24,i,j,k)  &
                     + xyz (33:36) * ed_cellCubicNele (37:40,i,j,k)  &
                     + xyz (49:52) * ed_cellCubicNele (53:56,i,j,k))

            b = sum (  xyz ( 5: 8) * ed_cellCubicNele ( 9:12,i,j,k)  &
                     + xyz (21:24) * ed_cellCubicNele (25:28,i,j,k)  &
                     + xyz (37:40) * ed_cellCubicNele (41:44,i,j,k)  &
                     + xyz (53:56) * ed_cellCubicNele (57:60,i,j,k))

            c = sum (  xyz ( 9:12) * ed_cellCubicNele (13:16,i,j,k)  &
                     + xyz (25:28) * ed_cellCubicNele (29:32,i,j,k)  &
                     + xyz (41:44) * ed_cellCubicNele (45:48,i,j,k)  &
                     + xyz (57:60) * ed_cellCubicNele (61:64,i,j,k))

            gradNeleY = (a + b + b + c + c + c) * deltaInvY

            a = sum (  xyz ( 1:16) * ed_cellCubicNele (17:32,i,j,k))
            b = sum (  xyz (17:32) * ed_cellCubicNele (33:48,i,j,k))
            c = sum (  xyz (33:48) * ed_cellCubicNele (49:64,i,j,k))

            gradNeleZ = (a + b + b + c + c + c) * deltaInvZ

            Tele = sum (xyz (1:64) * ed_cellCubicTele (1:64,i,j,k))

            a = sum ( xyz (1:61:4) * ed_cellCubicTele (2:62:4,i,j,k))
            b = sum ( xyz (2:62:4) * ed_cellCubicTele (3:63:4,i,j,k))
            c = sum ( xyz (3:63:4) * ed_cellCubicTele (4:64:4,i,j,k))

            gradTeleX = (a + b + b + c + c + c) * deltaInvX

            a = sum (  xyz ( 1: 4) * ed_cellCubicTele ( 5: 8,i,j,k)  &
                     + xyz (17:20) * ed_cellCubicTele (21:24,i,j,k)  &
                     + xyz (33:36) * ed_cellCubicTele (37:40,i,j,k)  &
                     + xyz (49:52) * ed_cellCubicTele (53:56,i,j,k))

            b = sum (  xyz ( 5: 8) * ed_cellCubicTele ( 9:12,i,j,k)  &
                     + xyz (21:24) * ed_cellCubicTele (25:28,i,j,k)  &
                     + xyz (37:40) * ed_cellCubicTele (41:44,i,j,k)  &
                     + xyz (53:56) * ed_cellCubicTele (57:60,i,j,k))

            c = sum (  xyz ( 9:12) * ed_cellCubicTele (13:16,i,j,k)  &
                     + xyz (25:28) * ed_cellCubicTele (29:32,i,j,k)  &
                     + xyz (41:44) * ed_cellCubicTele (45:48,i,j,k)  &
                     + xyz (57:60) * ed_cellCubicTele (61:64,i,j,k))

            gradTeleY = (a + b + b + c + c + c) * deltaInvY

            a = sum (  xyz ( 1:16) * ed_cellCubicTele (17:32,i,j,k))
            b = sum (  xyz (17:32) * ed_cellCubicTele (33:48,i,j,k))
            c = sum (  xyz (33:48) * ed_cellCubicTele (49:64,i,j,k))

            gradTeleZ = (a + b + b + c + c + c) * deltaInvZ

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Nele <= 0 for a cell (target)')
            end if

            if (Tele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Tele <= 0 for a cell (target)')
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
                      .and. (rayY > ed_yminDomain) &
                      .and. (rayY < ed_ymaxDomain) &
                      .and. (rayZ > ed_zminDomain) &
                      .and. (rayZ < ed_zmaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX             ! undo the nudging
                 ed_rays (RAY_POSY,ray) = rayY - nudgeY             ! since it is not
                 ed_rays (RAY_POSZ,ray) = rayZ - nudgeZ             ! needed anymore
                 ed_rays (RAY_BLCK,ray) = real (RAY_OUTDOMAIN)
                 ed_energyOutTimeStep   = ed_energyOutTimeStep + rayPower * timeStep
            end if

            exit

        end if
!
!
!-------------------- End loop following ray through cells in block --------------------------------------------
!
!
     end do
!
!
!     ...Check to see if we ran out of laser IO buffer space
!
!
     if(writeRay .and. (nRayWritePositions > ed_laserIOMaxNumberOfPositions) ) then
        print *, "[ed_traceBlockRays3DRec] Ray ", ray, &
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
end subroutine ed_traceBlockRays3DRec
