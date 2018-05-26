!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_traceBlockRays3DRec
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
                                     ed_cellCenters,                 &
                                     ed_cellDensity,                 &
                                     ed_cellEdges,                   &
                                     ed_cellGradNele,                &
                                     ed_cellGradTele,                &
                                     ed_cellNele,                    &
                                     ed_cellTele,                    &
                                     ed_cellZbar,                    &
                                     ed_cellWallThickness,           &
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

  use ed_commInterface, ONLY : ed_commProgressTransport, &
       ed_commHandleOffBlkRay, &
       ed_commIncrementDeadRays

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
  logical :: onBlockFaceBoundary, onBlockCellBoundary
  logical :: rayCrossesBoundaryX, rayCrossesBoundaryY, rayCrossesBoundaryZ
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
  real    :: cellMass
  real    :: cellPower
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: cellZbar
  real    :: centerCellX, centerCellY, centerCellZ
  real    :: centerNele, centerNeleOld, centerTele
  real    :: crossTime, crossTimeHalf
  real    :: distToFaceMinX, distToFaceMinY, distToFaceMinZ
  real    :: distToFaceMaxX, distToFaceMaxY, distToFaceMaxZ
  real    :: distX, distY, distZ
  real    :: distXOld, distYOld, distZOld
  real    :: gradNeleX, gradNeleY, gradNeleZ
  real    :: gradNeleXOld, gradNeleYOld, gradNeleZOld
  real    :: gradTeleX, gradTeleY, gradTeleZ
  real    :: integral
  real    :: lnLambda
  real    :: minFaceDistance
  real    :: Nele, NeleOld, Tele
  real    :: nu
  real    :: nudgeX, nudgeY, nudgeZ
  real    :: powerLossFactor
  real    :: R,S,U,W
  real    :: rayCritDens, rayCritDensInv
  real    :: rayPower
  real    :: rayX, rayY, rayZ
  real    :: time2Face
  real    :: velX, velY, velZ
  real    :: vnewSqr
  real    :: xminCell, yminCell, zminCell
  real    :: xmaxCell, ymaxCell, zmaxCell

  real, parameter :: GaussianRoot1 = 1.577350269189626         ! is 1 + 1/sqrt(3) for integration limits [0,1]
  real, parameter :: GaussianRoot2 = 4.226497308103742e-01     ! is 1 - 1/sqrt(3) for integration limits [0,1]

  real    :: aCoeff (1:6)
  real    :: bCoeff (1:6)
  real    :: cCoeff (1:6)

  integer :: numDeadRays
  numDeadRays = 0
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

     call ed_commProgressTransport()

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
!     ...Calculate the electron density and the electron temperature at the current ray
!        position in the (i,j,k) cell. These are calculated via interpolation from the
!        corresponding gradient values. If either of these values become negative, the
!        calculation must be stopped.
!
!
     centerCellX = ed_cellCenters (i,1)
     centerCellY = ed_cellCenters (j,2)
     centerCellZ = ed_cellCenters (k,3)
     centerNele  = ed_cellNele     (  i,j,k)
     centerTele  = ed_cellTele     (  i,j,k)
     gradNeleX   = ed_cellGradNele (1,i,j,k)
     gradNeleY   = ed_cellGradNele (2,i,j,k)
     gradNeleZ   = ed_cellGradNele (3,i,j,k)
     gradTeleX   = ed_cellGradTele (1,i,j,k)
     gradTeleY   = ed_cellGradTele (2,i,j,k)
     gradTeleZ   = ed_cellGradTele (3,i,j,k)

     distX   = rayX - centerCellX
     distY   = rayY - centerCellY
     distZ   = rayZ - centerCellZ
     Nele    = centerNele  +  gradNeleX * distX  +  gradNeleY * distY  +  gradNeleZ * distZ
     Tele    = centerTele  +  gradTeleX * distX  +  gradTeleY * distY  +  gradTeleZ * distZ

     if (Nele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Nele <= 0 for a cell (initial)')
     end if

     if (Tele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Tele <= 0 for a cell (initial)')
     end if
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. The current cell indices (i,j,k) and the previous cell indices (ip,jp,kp)
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

        xminCell    = ed_cellEdges   (i  ,1)
        xmaxCell    = ed_cellEdges   (i+1,1)
        yminCell    = ed_cellEdges   (j  ,2)
        ymaxCell    = ed_cellEdges   (j+1,2)
        zminCell    = ed_cellEdges   (k  ,3)
        zmaxCell    = ed_cellEdges   (k+1,3)

        cellZbar    = ed_cellZbar    (i,j,k)
        cellDensity = ed_cellDensity (i,j,k)
        cellMass    = cellDensity * cellVolume
!
!
!     ...The ray is being refracted through the cell. We have to follow its path through
!        the cell (i,j,k) and determine its exit cell. At that point the previous cell
!        becomes simply the current cell, that is cell (ip,jp,kp) = cell (i,j,k) and
!        the loop is closed.
!
!        Since the ray is moving classically, given an initial position 's0', an initial
!        velocity 'v0' and a constant acceleration 'a', at time 't' we have the new
!        ray position 's' and velocity 'v' given by the equations:
!
!                            s(t) = s0 + v0*t + (a/2)*t^2
!                            v(t) = v0 + a*t
!
!        Note, that s(t) is not the distance travelled by the ray. That would have to
!        be determined by a line integral over the individual s(t)'s. For our purposes
!        we need to find s(t) and v(t) at the position where the ray will exit cell (i,j,k).
!        We have to examine, which 't' leads to crossing of a cell face and pick the
!        smallest 't'. This can be done by looking at the quadratic position equation:
!
!                                   A*t^2 + B*t + C = s(t)
!
!        where in our case:
!
!                            A = a/2 = -c^2*grad(ne)/(4*nc)
!                            B = v0
!                            C = s0
!
!        and 'c' is the speed of light, 'ne' the electron density, 'nc' the critical
!        density, 'v0' the initial ray velocity and 's0' the initial ray position. The
!        three quadratic equation coefficients A,B,C stand for the three x,y,z components,
!        so in fact we have the 3 component vector equation:
!
!                       |Ax|           |Bx|         |Cx|     |sx(t)|
!                       |Ay| * t^2  +  |By| * t  +  |Cy|  =  |sy(t)|
!                       |Az|           |Bz|         |Cz|     |sz(t)|
!
!        The general equation of a plane in 3D space is (E,F,G,H = numerical coefficients):
!
!                                     Ex + Fy + Gz + H = 0
!
!        If the ray crosses that plane in time 't', then at position sx(t),sy(t) and sz(t)
!        the plane equation is fulfilled, that is:
!
!                             E * sx(t) + F * sy(t) + G * sz(t) + H = 0
!
!        To find the time it takes to cross that plane we now substitute the x,y,z
!        components of the quadratic time equation into the plane equation. We get
!
!           E[Ax*t^2 + Bx*t + Cx] + F[Ay*t^2 + By*t + Cy] + G[Az*t^2 + Bz*t + Cz] + H = 0
!
!        and collecting the 't' terms:
!
!        [E*Ax + F*Ay + G*Az] * t^2 + [E*Bx + F*By + G*Bz] * t + [E*Cx + F*Cy + G*Cz + H] = 0
!
!        Solving this equation gives two possible 't' values. The plane will be crossed
!        if 't' is >= 0. Note, that this does not mean that an actual quadrilateral (i.e.
!        a face) contained in the plane will be crossed. This has to be checked after one
!        determines the plane's crossing points.
!
!        In our present case the quadratic equation simplifies considerably because the
!        planes containing the cell's faces are all parallel to the coordinate system. The
!        quadratic equation represents actually a set of 6 different quadratic equations:
!        3 for each coordinate direction times 2 for each of the possible 2 cell faces in
!        each direction. We hence have to solve 6 quadratic equations, each giving a pair
!        of t's.
!
!        Significance of the type of t solutions:
!
!          1) t is negative/imaginary -> the ray is moving away from that particular
!                                        cell face plane and will not hit it.
!
!          2) t is positive           -> the ray will hit the cell face plane in t time.
!
!          3) t is zero               -> the ray is already on the cell face plane.
!
!        Among all 't' solutions we need the smallest positive 't' of all 12 possible
!        solutions. The order of the four cell faces is as follows:
!
!                  face 1 -> face in the plane with constant minimum x-axis value 
!                  face 2 -> face in the plane with constant maximum x-axis value 
!                  face 3 -> face in the plane with constant minimum y-axis value 
!                  face 4 -> face in the plane with constant maximum y-axis value 
!                  face 5 -> face in the plane with constant minimum z-axis value 
!                  face 6 -> face in the plane with constant maximum z-axis value 
!
!
        accX = - c2div2nc * gradNeleX               ! acceleration in x-direction
        accY = - c2div2nc * gradNeleY               ! acceleration in y-direction
        accZ = - c2div2nc * gradNeleZ               ! acceleration in z-direction

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

        crossTime = ed_infiniteTime

        do face = 1,6

           a = aCoeff (face)
           b = bCoeff (face)
           c = cCoeff (face)

           if (a == 0.0) then
               if (b /= 0.0) then
                   time2Face = - c / b
                   if (time2Face > 0.0 .and. time2Face < crossTime) crossTime = time2Face
               end if
           else
               if (c == 0.0) then
                   time2Face = - b / a
                   if (time2Face > 0.0 .and. time2Face < crossTime) crossTime = time2Face
               else
                   d = b * b - 4.0 * a * c
                   if ((d > 0.0)) then
                        q = -0.5 * (b + sign (1.0,b) * sqrt (d))
                        time2Face = q / a
                        if (time2Face > 0.0 .and. time2Face < crossTime) crossTime = time2Face
                        time2Face = c / q
                        if (time2Face > 0.0 .and. time2Face < crossTime) crossTime = time2Face
                   else if (d == 0.0) then
                        time2Face = -0.5 * b / a
                        if (time2Face > 0.0 .and. time2Face < crossTime) crossTime = time2Face
                   end if
               end if
           end if

        end do

        if (crossTime == ed_infiniteTime) then
            call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: infinite crossing time for a cell')
        end if

        crossTimeHalf = 0.5 * crossTime
!
!
!     ...We found a reasonable crossing time to a particular cell face plane. Calculate the
!        power deposition as the ray traverses the cell. This is done by evaluating an integral
!        using the initial (before crossing) velocities. The method of integral evaluation is
!        Gaussian Quadrature with weight function equal to 1. The associated orthogonal
!        polynomials are the Legendre Polynomials. If the remaining ray power is considered to
!        have reached a 'zero' value, mark the ray as nonexistent and exit the indefinite loop.
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

        a = GaussianRoot1 * crossTimeHalf
        b = a * a
        c = 1.0 + U * a + R * b
        d = sqrt (1.0 + W * a + S * b)

        integral = c * c / (d * d * d)

        a = GaussianRoot2 * crossTimeHalf
        b = a * a
        c = 1.0 + U * a + R * b
        d = sqrt (1.0 + W * a + S * b)

        integral = integral + c * c / (d * d * d)
        integral = integral * nu * crossTimeHalf

        powerLossFactor = exp (-integral)
        cellPower       = rayPower * (1.0 - powerLossFactor)
        cellEnergy      = cellPower * timeStep

        if (ed_depoVarIsPerMass) then
           cellEnergyDepot (i,j,k) = cellEnergyDepot (i,j,k) + cellEnergy / cellMass
        else
           cellEnergyDepot (i,j,k) = cellEnergyDepot (i,j,k) + cellEnergy * cellVolumeInv
        end if

        rayPower = rayPower * powerLossFactor

        if (rayPower <= ed_rayZeroPower) then
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            numDeadRays = numDeadRays + 1
            exit
        end if
!
!
!     ...Calculate the new ray position and velocities and check, if the ray crosses
!        a cell face. In case a laser IO is performed on the rays, store the current
!        ray IO data.
!
!
        rayX = rayX + (velX + accX * crossTimeHalf) * crossTime
        rayY = rayY + (velY + accY * crossTimeHalf) * crossTime
        rayZ = rayZ + (velZ + accZ * crossTimeHalf) * crossTime

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

        velX = velX + accX * crossTime
        velY = velY + accY * crossTime
        velZ = velZ + accZ * crossTime

        velXgt0 = (velX  > 0.0)
        velXeq0 = (velX == 0.0)
        velXlt0 = (velX  < 0.0)
        velYgt0 = (velY  > 0.0)
        velYeq0 = (velY == 0.0)
        velYlt0 = (velY  < 0.0)
        velZgt0 = (velZ  > 0.0)
        velZeq0 = (velZ == 0.0)
        velZlt0 = (velZ  < 0.0)

        stationaryRay = velXeq0 .and. velYeq0 .and. velZeq0

        if (stationaryRay) then
            write (*,*) ' stationary ray detected! Removing it from list ... '
            ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
            numDeadRays = numDeadRays + 1
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
!     ...Calculate the cell indices of the new cell where the ray moves into. Again, the ray
!        is (at least) on one of the six faces and moves into the next cell from the current
!        cell, which becomes cell (ip,jp,kp). The ray can be on:
!
!                       1 cell face  -> within the cell face
!                       2 cell faces -> on the edge between the cell faces
!                       3 cell faces -> in the corner defined by the three cell faces
!
!        The goal here is now to determine into which cell (i,j,k) the ray will go. Cell (i,j,k)
!        will in most of the cases be different from cell (ip,jp,kp), but occasionally both cells
!        can remain identical. If the ray's position is considered to be close to a face and it
!        is certain that the ray will cross that face, the ray is forced (temporarily) to be
!        exactly on that face. This is done in order to calculate the correct refraction
!        properties of the ray, if needed. The nudging value(s) and their directions to keep
!        the ray in the current cell are also set here.
!
!        The variables 'cross(X,Y,Z)' indicate, if the ray crosses the corresponding cell boundaries
!        and Snell's law has to be applied to the X,Y,Z components of the ray's velocity in case of
!        refraction. Note, that if at least one of the 'cross(X,Y,Z)' is true, then the ray crosses
!        a cell boundary from cell (ip,jp,kp) to cell (i,j,k).
!
!
        ip = i
        jp = j
        kp = k

        crossX = .false.
        crossY = .false.
        crossZ = .false.

        nudgeX = 0.0
        nudgeY = 0.0
        nudgeZ = 0.0

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
!
!
!     ...At this stage, the ray situation is as follows:
!
!                                     ---------- ----------
!                                   /|         /|         /|
!                                  / |        / |        / |
!                                  ---------- ----------   |          P = previous ray position
!                                 |  |       |  |       |  |
!                                 | P|    *  | N|    +  |  |          N = current ray position
!                                 |  /       |  /       |  /
!                                 | /        | /        | /         *,+ = cell centers
!                                 |/ip,jp,kp |/  i,j,k  |/
!                                  ---------- ----------
!
!        The new target cell (i,j,k) is not final yet. It depends on the outcome of Snell's law and
!        the block face boundary conditions. Check first, if a block face reflection needs to be honored.
!        After that, check if Snell's law needs to be applied to find the new velocity components.
!        For this we need to know the number of electrons in the previous old cell and the new cell
!        right at the current boundary. Again this is calculated via interpolation from the gradient
!        values. If the ray reflects, the new target cell is updated and the corresponding velocity
!        component is inverted. If the ray still crosses a cell boundary, determine the proper
!        nudging for the ray and calculate the new nudged ray positions.
!
!
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

        rayCrossesBoundaryX = (.not.reflectX .and. crossX)
        rayCrossesBoundaryY = (.not.reflectY .and. crossY)
        rayCrossesBoundaryZ = (.not.reflectZ .and. crossZ)

        if (rayCrossesBoundaryX .or. rayCrossesBoundaryY .or. rayCrossesBoundaryZ) then

            centerNeleOld = centerNele                          ! for cell (ip,jp,kp) at center *
            gradNeleXOld  = gradNeleX                           ! for cell (ip,jp,kp) at center *
            gradNeleYOld  = gradNeleY                           ! for cell (ip,jp,kp) at center *
            gradNeleZOld  = gradNeleZ                           ! for cell (ip,jp,kp) at center *
            distXOld      = rayX - centerCellX                  ! for cell (ip,jp,kp) at position N
            distYOld      = rayY - centerCellY                  ! for cell (ip,jp,kp) at position N
            distZOld      = rayZ - centerCellZ                  ! for cell (ip,jp,kp) at position N
            centerNele    = ed_cellNele     (  i,j,k)           ! for cell (i ,j ,k ) at center +
            gradNeleX     = ed_cellGradNele (1,i,j,k)           ! for cell (i ,j ,k ) at center +
            gradNeleY     = ed_cellGradNele (2,i,j,k)           ! for cell (i ,j ,k ) at center +
            gradNeleZ     = ed_cellGradNele (3,i,j,k)           ! for cell (i ,j ,k ) at center +
            distX         = rayX - ed_cellCenters (i,1)         ! for cell (i ,j ,k ) at position N
            distY         = rayY - ed_cellCenters (j,2)         ! for cell (i ,j ,k ) at position N
            distZ         = rayZ - ed_cellCenters (k,3)         ! for cell (i ,j ,k ) at position N

            NeleOld = centerNeleOld + gradNeleXOld * distXOld &
                                    + gradNeleYOld * distYOld &
                                    + gradNeleZOld * distZOld

            Nele    = centerNele    + gradNeleX    * distX    &
                                    + gradNeleY    * distY    &
                                    + gradNeleZ    * distZ

            if (NeleOld <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: NeleOld < 0 for a cell (Snell)')
            end if

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays3DRec] ERROR: Nele <= 0 for a cell (Snell)')
            end if

            if (rayCrossesBoundaryX) then
                vNewSqr  = velX * velX + (NeleOld - Nele) * c2div1nc
                reflectX = vNewSqr < 0.0
                if (.not.reflectX) then
                     velX = sign (1.,velX) * sqrt (vNewSqr)
                end if
            end if

            if (rayCrossesBoundaryY) then
                vNewSqr  = velY * velY + (NeleOld - Nele) * c2div1nc
                reflectY = vNewSqr < 0.0
                if (.not.reflectY) then
                     velY = sign (1.,velY) * sqrt (vNewSqr)
                end if
            end if

            if (rayCrossesBoundaryZ) then
                vNewSqr  = velZ * velZ + (NeleOld - Nele) * c2div1nc
                reflectZ = vNewSqr < 0.0
                if (.not.reflectZ) then
                     velZ = sign (1.,velZ) * sqrt (vNewSqr)
                end if
            end if

        end if

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
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i,j,k) is still within the block.
!        If it is, we have to calculate the new electron density and the new electron temperature where
!        the ray is located in the target cell. If the target cell is not within the block, check if the
!        ray coordinates are still within the defined domain. If not, store its latest data and mark it
!        as nonexistent. If the ray is still within the domain boundaries, exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock) &
                 .and. (k >= kminBlock) &
                 .and. (k <= kmaxBlock)

        if (inBlock) then

            centerCellX = ed_cellCenters (i,1)
            centerCellY = ed_cellCenters (j,2)
            centerCellZ = ed_cellCenters (k,3)
            centerNele  = ed_cellNele     (  i,j,k)
            centerTele  = ed_cellTele     (  i,j,k)
            gradNeleX   = ed_cellGradNele (1,i,j,k)
            gradNeleY   = ed_cellGradNele (2,i,j,k)
            gradNeleZ   = ed_cellGradNele (3,i,j,k)
            gradTeleX   = ed_cellGradTele (1,i,j,k)
            gradTeleY   = ed_cellGradTele (2,i,j,k)
            gradTeleZ   = ed_cellGradTele (3,i,j,k)

            distX   = rayX - centerCellX
            distY   = rayY - centerCellY
            distZ   = rayZ - centerCellZ
            Nele    = centerNele  +  gradNeleX * distX  +  gradNeleY * distY  +  gradNeleZ * distZ
            Tele    = centerTele  +  gradTeleX * distX  +  gradTeleY * distY  +  gradTeleZ * distZ

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
                 numDeadRays = numDeadRays + 1
            else
                 call ed_commHandleOffBlkRay(ray)
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

  if (numDeadRays > 0) then
     call ed_commIncrementDeadRays(numDeadRays)
  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceBlockRays3DRec
