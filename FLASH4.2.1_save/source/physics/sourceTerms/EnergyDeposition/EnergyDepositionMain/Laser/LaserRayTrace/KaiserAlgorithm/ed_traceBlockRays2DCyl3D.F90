!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_traceBlockRays2DCyl3D
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
                                     ed_cellCenters,                 &
                                     ed_cellDensity,                 &
                                     ed_cellEdges,                   &
                                     ed_cellGradNele,                &
                                     ed_cellGradTele,                &
                                     ed_cellNele,                    &
                                     ed_cellTele,                    &
                                     ed_cellVolume,                  &
                                     ed_cellZbar,                    &
                                     ed_cellWallThickness,           &
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
  logical :: onBlockFaceBoundary, onBlockWedgeBoundary
  logical :: rayCrossesBoundaryX, rayCrossesBoundaryZ
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
  real    :: centerWedgeX, centerWedgeZ
  real    :: centerNele, centerNeleOld, centerTele
  real    :: crossTime, crossTimeHalf
  real    :: distToFaceMinX, distToFaceMaxX
  real    :: distToFaceMinY, distToFaceMaxY
  real    :: distToFaceMinZ, distToFaceMaxZ
  real    :: distX, distZ
  real    :: distXOld, distZOld
  real    :: gradNeleX, gradNeleZ
  real    :: gradNeleXOld, gradNeleZOld
  real    :: gradTeleX, gradTeleZ
  real    :: integral
  real    :: lnLambda
  real    :: minFaceDistance
  real    :: Nele, NeleOld, Tele
  real    :: nu
  real    :: nudgeX, nudgeZ
  real    :: powerLossFactor
  real    :: R,S,U,W
  real    :: rayCritDens, rayCritDensInv
  real    :: rayPower
  real    :: rayX, rayY, rayZ
  real    :: time2Face
  real    :: velX, velY, velZ
  real    :: vnewSqr
  real    :: wedgeCosine
  real    :: wedgeDensity
  real    :: wedgeEnergy
  real    :: wedgeMass
  real    :: wedgePower
  real    :: wedgeSine
  real    :: wedgeSlope
  real    :: wedgeVolume
  real    :: wedgeWallThicknessHalf
  real    :: wedgeZbar
  real    :: xminWedge, yminWedge, zminWedge
  real    :: xmaxWedge, ymaxWedge, zmaxWedge

  real, parameter :: GaussianRoot1 = 1.577350269189626         ! is 1 + 1/sqrt(3) for integration limits [0,1]
  real, parameter :: GaussianRoot2 = 4.226497308103742e-01     ! is 1 - 1/sqrt(3) for integration limits [0,1]

  real    :: aCoeff (1:6)
  real    :: bCoeff (1:6)
  real    :: cCoeff (1:6)

  integer :: numDeadRays
  numDeadRays = 0
!
!
!     ...Define some variables.
!
!
  wedgeCosine            =       ed_laser3Din2DwedgeCosine
  wedgeSine              =       ed_laser3Din2DwedgeSine
  wedgeSlope             =       ed_laser3Din2DwedgeSlope
  wedgeWallThicknessHalf = 0.5 * ed_cellWallThickness
!
!
!     ...Outer (threaded) loop over all rays associated with the current 2D cylindrical block.
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

!     if (rayOutOfBlock .and. ed_raysMovedIntoDomain) then
!         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: ray moves out of block')
!     end if
!
!
!     ...Calculate the electron density and the electron temperature at the current ray
!        position in the (i,j) wedge. These are calculated via interpolation from the
!        corresponding 2D cylindrical gradient values. If either of these values become
!        negative, the calculation must be stopped.
!
!
     centerWedgeX = ed_cellCenters  (i,1)
     centerWedgeZ = ed_cellCenters  (j,2)
     centerNele   = ed_cellNele     (  i,j,1)
     centerTele   = ed_cellTele     (  i,j,1)
     gradNeleX    = ed_cellGradNele (1,i,j,1)
     gradNeleZ    = ed_cellGradNele (2,i,j,1)
     gradTeleX    = ed_cellGradTele (1,i,j,1)
     gradTeleZ    = ed_cellGradTele (2,i,j,1)

     distX = rayX - centerWedgeX
     distZ = rayZ - centerWedgeZ
     Nele  = centerNele  +  gradNeleX * distX  +  gradNeleZ * distZ
     Tele  = centerTele  +  gradTeleX * distX  +  gradTeleZ * distZ

     if (Nele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Nele <= 0 for a wedge (initial)')
     end if

     if (Tele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Tele <= 0 for a wedge (initial)')
     end if
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

        xminWedge    = ed_cellEdges   (i  ,1)
        xmaxWedge    = ed_cellEdges   (i+1,1)
        zminWedge    = ed_cellEdges   (j  ,2)
        zmaxWedge    = ed_cellEdges   (j+1,2)

        wedgeZbar    = ed_cellZbar    (i,j,1)
        wedgeVolume  = ed_cellVolume  (i,j,1)
        wedgeDensity = ed_cellDensity (i,j,1)
        wedgeMass    = wedgeDensity * wedgeVolume
!
!
!     ...The 3D ray is being refracted through the wedge (i,j). We have to follow its path
!        through the wedge and determine its exit wedge. At that point the previous wedge
!        becomes simply the current wedge, that is wedge (ip,jp) = wedge (i,j) and the loop
!        is closed.
!
!        Since the 3D ray is moving classically, given an initial position 's0', an initial
!        velocity 'v0' and a constant acceleration 'a', at time 't' we have the new
!        ray position 's' and velocity 'v' given by the equations:
!
!                            s(t) = s0 + v0*t + (a/2)*t^2
!                            v(t) = v0 + a*t
!
!        Note, that s(t) is not the distance travelled by the ray. That would have to
!        be determined by a line integral over the individual s(t)'s. For our purposes
!        we need to find s(t) and v(t) at the position where the ray will exit wedge (i,j).
!        We have to examine, which 't' leads to crossing of a wedge face and pick the
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
!        The 6 planes that delimit the 3D (i,j) wedge are:
!
!                    1) x = xminWedge                   E = 1 , H = - xminWedge
!                    2) x = xmaxWedge                   E = 1 , H = - xmaxWedge
!                    3) y = - wedgeSlope * x            E = + wedgeSlope , F = 1
!                    4) y = + wedgeSlope * x            E = - wedgeSlope , F = 1
!                    5) z = zminWedge                   G = 1 , H = - zminWedge
!                    6) z = zmaxWedge                   G = 1 , H = - zmaxWedge
!
!        The quadratic equation represents actually a set of 6 different quadratic equations:
!        3 for each coordinate direction times 2 for each of the possible 2 wedge faces in
!        each direction. We hence have to solve 6 quadratic equations, each giving a pair
!        of t's.
!
!        Significance of the type of t solutions:
!
!          1) t is negative/imaginary -> the ray is moving away from that particular
!                                        wedge face plane and will not hit it.
!
!          2) t is positive           -> the ray will hit the wedge face plane in t time.
!
!          3) t is zero               -> the ray is already on the wedge face plane.
!
!        Among all 't' solutions we need the smallest positive 't' of all 8 possible
!        solutions. The order of the six wedge faces is as follows:
!
!                  face 1 -> face in the plane with constant minimum x-axis value 
!                  face 2 -> face in the plane with constant maximum x-axis value 
!                  face 3 -> face in the plane y = - wedgeSlope * x
!                  face 4 -> face in the plane y = + wedgeSlope * x
!                  face 5 -> face in the plane with constant minimum z-axis value 
!                  face 6 -> face in the plane with constant maximum z-axis value 
!
!
        accX = - c2div2nc * gradNeleX           ! acceleration in x-direction (radial)
        accZ = - c2div2nc * gradNeleZ           ! acceleration in z-direction

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
        cCoeff (3) = rayY + wedgeSlope * rayX
        cCoeff (4) = rayY - wedgeSlope * rayX
        cCoeff (5) = rayZ - zminWedge
        cCoeff (6) = rayZ - zmaxWedge

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
            call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: infinite wedge crossing time')
        end if

        crossTimeHalf = 0.5 * crossTime
!
!
!     ...We found a reasonable crossing time to a particular wedge face plane. Calculate the
!        power deposition as the ray traverses the wedge. This is done by evaluating an integral
!        using the initial (before crossing) velocities. The method of integral evaluation is
!        Gaussian Quadrature with weight function equal to 1. The associated orthogonal
!        polynomials are the Legendre Polynomials. If the remaining ray power is considered to
!        have reached a 'zero' value, mark the ray as nonexistent and exit the indefinite loop.
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
        wedgePower      = rayPower * (1.0 - powerLossFactor)
        wedgeEnergy     = wedgePower * timeStep

        if (ed_depoVarIsPerMass) then
           wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy / wedgeMass
        else
           wedgeEnergyDepot (i,j) = wedgeEnergyDepot (i,j) + wedgeEnergy / wedgeVolume
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
!        a wedge face. In case a laser IO is performed on the rays, store the current
!        ray IO data.
!
!
        rayX = rayX + (velX + accX * crossTimeHalf) * crossTime
        rayY = rayY + (velY                       ) * crossTime
        rayZ = rayZ + (velZ + accZ * crossTimeHalf) * crossTime

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

        velX = velX + accX * crossTime
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

        stationaryRay  = velXeq0 .and. velYeq0 .and. velZeq0

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
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = rayZ
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
            end if
        end if
!
!
!     ...Calculate the wedge 2D cylindrical indices of the new wedge where the ray moves into.
!        If the ray crosses either the x-faces and/or the z-faces, then these indices will
!        change. If the ray crosses only the y-faces, then the indices stay the same.
!
!        The goal here is now to determine into which wedge (i,j) the ray will go. If the ray's
!        position is considered to be close to a 2D cylindrical face and it is certain that the
!        ray will cross that face, the ray is forced (temporarily) to be exactly on that face.
!        This is done in order to calculate the correct refraction properties of the ray, if needed.
!        The nudging value(s) for the 2D cylindrical faces and their directions to keep the ray
!        in the current wedge are also set here.
!
!        The variables 'cross(X,Z)' indicate, if the ray crosses the corresponding 2D cylindrical
!        wedge boundaries and Snell's law has to be applied to the X,Z components of the ray's
!        velocity in case of refraction. Note, that if at least one of the 'cross(X,Z)' is true, then
!        the ray crosses a 2D cylindrical wedge boundary from old wedge (ip,jp) to new wedge (i,j).
!
!        The crossing of a y-face is treated differently, because the y-face represents a jump
!        in cylindrical curvature. When a ray crosses any of the y-faces, then the ray enters
!        a new wedge in which the x,y coordinate system has been rotated when compared to the
!        original wedge. This rotation implies new rotated x,y velocity components for the ray
!        travelling now in the new wedge.
!
!
        ip = i
        jp = j

        crossX = .false.
        crossZ = .false.

        nudgeX = 0.0
        nudgeZ = 0.0

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
!
!
!     ...At this stage, the ray situation is as follows:
!
!
!                                  ---------- ----------
!                                 |          |          |
!                                 |          |          |          P = previous ray position
!                                 |     *    N     +    |          N = current ray position
!                                 P          |          |        *,+ = wedge centers
!                                 |   ip,jp  |    i,j   |
!                                  ---------- ----------
!
!        The new target wedge (i,j) is not final yet. It depends on the outcome of Snell's law and
!        the block face boundary conditions. Check first, if a block face reflection needs to be honored.
!        After that, check if Snell's law needs to be applied to find the new velocity components.
!        For this we need to know the number of electrons in the previous old wedge and the new wedge
!        right at the current boundary. Again this is calculated via interpolation from the gradient
!        values. If the ray reflects, the new target wedge is updated and the corresponding velocity
!        component is inverted. If the ray still crosses a 2D cylindrical wedge boundary, determine the
!        proper nudging for the ray and calculate the new nudged ray positions.
!
!
        blockFaceMinX = (rayX == xminBlock)
        blockFaceMaxX = (rayX == xmaxBlock)
        blockFaceMinZ = (rayZ == zminBlock)
        blockFaceMaxZ = (rayZ == zmaxBlock)

        reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                  .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)
        reflectZ =     (blockFaceMinZ .and. blockReflectMinZ .and. velZlt0) &
                  .or. (blockFaceMaxZ .and. blockReflectMaxZ .and. velZgt0)

        rayCrossesBoundaryX = (.not.reflectX .and. crossX)
        rayCrossesBoundaryZ = (.not.reflectZ .and. crossZ)

        if (rayCrossesBoundaryX .or. rayCrossesBoundaryZ) then

            centerNeleOld = centerNele                          ! for wedge (ip,jp) at center *
            gradNeleXOld  = gradNeleX                           ! for wedge (ip,jp) at center *
            gradNeleZOld  = gradNeleZ                           ! for wedge (ip,jp) at center *
            distXOld      = rayX - centerWedgeX                 ! for wedge (ip,jp) at position N
            distZOld      = rayZ - centerWedgeZ                 ! for wedge (ip,jp) at position N
            centerNele    = ed_cellNele       (  i,j,1)         ! for wedge (i ,j ) at center +
            gradNeleX     = ed_cellGradNele   (1,i,j,1)         ! for wedge (i ,j ) at center +
            gradNeleZ     = ed_cellGradNele   (2,i,j,1)         ! for wedge (i ,j ) at center +
            distX         = rayX - ed_cellCenters (i,1)         ! for wedge (i ,j ) at position N
            distZ         = rayZ - ed_cellCenters (j,2)         ! for wedge (i ,j ) at position N

            NeleOld = centerNeleOld + gradNeleXOld * distXOld + gradNeleZOld * distZOld
            Nele    = centerNele    + gradNeleX    * distX    + gradNeleZ    * distZ

            if (NeleOld <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: NeleOld < 0 for a wedge (Snell)')
            end if

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays2DCyl3D] ERROR: Nele <= 0 for a wedge (Snell)')
            end if

            if (rayCrossesBoundaryX) then
                vNewSqr  = velX * velX + (NeleOld - Nele) * c2div1nc
                reflectX = vNewSqr < 0.0
                if (.not.reflectX) then
                     velX = sign (1.,velX) * sqrt (vNewSqr)
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
!
!
!     ...We are now sure about the target wedge. Check, if the target wedge (i,j) is still within the block.
!        If it is, we have to calculate the new electron density and the new electron temperature where
!        the ray is located in the target wedge. If the target wedge is not within the block, check if the
!        ray coordinates are still within the defined domain. If not, store its latest data and mark it as
!        nonexistent. If the ray is still within the domain boundaries, exit the current block loop.
!
!
        inBlock =      (i >= iminBlock) &
                 .and. (i <= imaxBlock) &
                 .and. (j >= jminBlock) &
                 .and. (j <= jmaxBlock)

        if (inBlock) then

            centerWedgeX = ed_cellCenters  (i,1)
            centerWedgeZ = ed_cellCenters  (j,2)
            centerNele   = ed_cellNele     (  i,j,1)
            centerTele   = ed_cellTele     (  i,j,1)
            gradNeleX    = ed_cellGradNele (1,i,j,1)
            gradNeleZ    = ed_cellGradNele (2,i,j,1)
            gradTeleX    = ed_cellGradTele (1,i,j,1)
            gradTeleZ    = ed_cellGradTele (2,i,j,1)

            distX   = rayX - centerWedgeX
            distZ   = rayZ - centerWedgeZ
            Nele    = centerNele  +  gradNeleX * distX  +  gradNeleZ * distZ
            Tele    = centerTele  +  gradTeleX * distX  +  gradTeleZ * distZ

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
                 numDeadRays = numDeadRays + 1
            else
                 call ed_commHandleOffBlkRay(ray)
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

  if (numDeadRays > 0) then
     call ed_commIncrementDeadRays(numDeadRays)
  end if

!
!
!     ...Ready!
!
!
  return
end subroutine ed_traceBlockRays2DCyl3D
