!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_traceBlockRays1DRec
!!
!! NAME
!!
!!  ed_traceBlockRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_traceBlockRays1DRec (real    (in)    :: timeStep,
!!                               integer (in)    :: rayFirst
!!                               integer (in)    :: rayLast,
!!                               integer (in)    :: iminBlock,
!!                               integer (in)    :: imaxBlock,
!!                               real    (in)    :: xminBlock,
!!                               real    (in)    :: xmaxBlock,
!!                               real    (in)    :: deltaInvX,
!!                               logical (in)    :: blockReflectMinX,
!!                               logical (in)    :: blockReflectMaxX,
!!                               real    (inout) :: cellEnergyDepot (:))
!!
!! DESCRIPTION
!!
!!  Traces the movement of the current collection of active rays through one block for
!!  those geometries consisting formally of 1D rectangular grids (cartesian + spherical).
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
!!  xminBlock        : minimum x-coordinate limit of the block
!!  xmaxBlock        : maximum x-coordinate limit of the block
!!  deltaInvX        : inverse of the cell's x-dimension
!!  blockReflectMinX : is the block boundary on the minimum x-side reflective ?
!!  blockReflectMaxX : is the block boundary on the maximum x-side reflective ?
!!  cellEnergyDepot  : array collecting the ray energy deposition for each cell
!!
!! NOTES
!!
!!  The code allows for threading to be used on the outer ray trace loop.
!!  The paths of the rays are computed using the geometric optics approximation. 
!!  Inside each cell, the paths of the rays are evaluated stepwise, using the
!!  monocubic expansions of the number electron density and electron temperature
!!  grid.
!!
!!***

subroutine ed_traceBlockRays1DRec (timeStep,                          &
                                   rayFirst,  rayLast,                &
                                   iminBlock, imaxBlock,              &
                                   xminBlock, xmaxBlock,              &
                                   deltaInvX,                         &
                                   blockReflectMinX,                  &
                                   blockReflectMaxX,                  &
                                                      cellEnergyDepot ) 

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
                                     ed_xmaxDomain

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
  real,    intent (in)    :: xminBlock, xmaxBlock
  real,    intent (in)    :: deltaInvX
  logical, intent (in)    :: blockReflectMinX
  logical, intent (in)    :: blockReflectMaxX
  real,    intent (inout) :: cellEnergyDepot (iminBlock:imaxBlock)

  logical :: blockCellMinX, blockCellMaxX
  logical :: blockFaceMinX, blockFaceMaxX
  logical :: cellFaceMinX, cellFaceMaxX
  logical :: crossX
  logical :: impossibleRay
  logical :: inDomain, inBlock
  logical :: newCell
  logical :: onBlockFaceBoundary, onBlockCellBoundary
  logical :: rayOutOfBlock
  logical :: reflectX
  logical :: stationaryRay
  logical :: velXeq0, velXgt0, velXlt0
  logical :: writeRay

  integer :: face
  integer :: i
  integer :: ip
  integer :: nRayWritePositions
  integer :: ray
  integer :: rayTag
  integer :: rayWriteIndex

  real    :: a,b,c,d,q
  real    :: accX
  real    :: c2div1nc, c2div2nc, c2div4nc
  real    :: cellDensity
  real    :: cellEnergy
  real    :: cellMass, cellMassInv
  real    :: cellPower
  real    :: cellVolume, cellVolumeInv
  real    :: cellWallThicknessHalf
  real    :: cellZbar
  real    :: distToFaceMinX, distToFaceMaxX
  real    :: gradNeleX
  real    :: gradTeleX
  real    :: integral
  real    :: lnLambda
  real    :: minFaceDistance
  real    :: Nele, Tele
  real    :: newX
  real    :: nu
  real    :: nudgeX
  real    :: powerLossFactor
  real    :: R,S,U,W
  real    :: rayCritDens, rayCritDensInv
  real    :: rayPower
  real    :: rayX
  real    :: stepTime, stepTimeHalf
  real    :: stepVelocity
  real    :: time2Face
  real    :: velX
  real    :: x1, x2, x3
  real    :: xminCell, xmaxCell

  real, parameter :: GaussianRoot1 = 1.577350269189626         ! is 1 + 1/sqrt(3) for integration limits [0,1]
  real, parameter :: GaussianRoot2 = 4.226497308103742e-01     ! is 1 - 1/sqrt(3) for integration limits [0,1]

  real    :: aCoeff (1:2)
  real    :: bCoeff (1:2)
  real    :: cCoeff (1:2)
!
!
!     ...Define some variables.
!
!
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
     velX           =      ed_rays (RAY_VELX,ray)
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
!     ...Find the index (i) of the initial cell through which the ray will
!        enter the block. We know for sure that the ray enters the block, because
!        otherwise it would not be on the current block list. Check, on which of
!        the two possible faces the ray currently is. The current position of the
!        ray is such that it is not exactly on the block boundary but nudged into
!        the block. The distance to the closest cell face is always less than
!        the cell wall thickness.
!
!
     i = iminBlock + int ( (rayX - xminBlock) * deltaInvX )

     blockCellMinX = (i == iminBlock)
     blockCellMaxX = (i == imaxBlock)

     onBlockCellBoundary = blockCellMinX .or. blockCellMaxX

     if (.not.onBlockCellBoundary) then
          call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray not in a block boundary cell')
     end if

     xminCell = ed_cellEdges (i  ,1)
     xmaxCell = ed_cellEdges (i+1,1)

     distToFaceMinX = abs (xminCell - rayX)
     distToFaceMaxX = abs (xmaxCell - rayX)

     minFaceDistance = min (distToFaceMinX, distToFaceMaxX)

     if (minFaceDistance > ed_cellWallThickness) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray too far away from block boundary cell face')
     end if

     cellFaceMinX = (distToFaceMinX <= ed_cellWallThickness)
     cellFaceMaxX = (distToFaceMaxX <= ed_cellWallThickness)

     impossibleRay = cellFaceMinX .and. cellFaceMaxX

     if (impossibleRay) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray cannot be on two opposite cell faces')
     end if

     blockFaceMinX = (blockCellMinX .and. cellFaceMinX)
     blockFaceMaxX = (blockCellMaxX .and. cellFaceMaxX)

     onBlockFaceBoundary = blockFaceMinX .or. blockFaceMaxX

     if (.not.onBlockFaceBoundary) then
          call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray not on block face boundary')
     end if

     velXgt0 = (velX  > 0.0)
     velXeq0 = (velX == 0.0)
     velXlt0 = (velX  < 0.0)

     stationaryRay = velXeq0

     if (stationaryRay) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: stationary ray at a block face boundary')
     end if

     rayOutOfBlock =     (blockFaceMinX .and. velXlt0) &
                    .or. (blockFaceMaxX .and. velXgt0)

     if (rayOutOfBlock .and. ed_raysMovedIntoDomain) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray moves out of block')
     end if
!
!
!     ...Calculate the electron number density and the electron temperature as well as their
!        gradient (d/dx) at the current ray position in the (i) cell. These are calculated
!        using the corresponding monocubic expansion coefficients. If either of these values
!        become negative, the calculation must be stopped.
!
!        Algorithm implemented to evaluate the function values and their derivatives
!        ---------------------------------------------------------------------------
!
!        The monocubic expansion coefficients a (i) for each cell contain all the info that
!        is needed. For a point (X) inside the cell with rescaled [0,1] x-coordinate, we have:
!
!
!                              3          i
!                     F (X) = sum  a (i) x
!                             i=0
!
!                              3              i-1
!                     dF/dX = sum  i * a (i) x    / (cell x-dimension)
!                             i=1
!
!
     x1 = (rayX - xminCell) * deltaInvX           ! rescaled [0,1] ray x coordinate
     x2 = x1 * x1
     x3 = x1 * x2
!
!
!     ...Electron number density + gradient.
!
!
     Nele =        ed_cellCubicNele (1,i,1,1) &
            + x1 * ed_cellCubicNele (2,i,1,1) &
            + x2 * ed_cellCubicNele (3,i,1,1) &
            + x3 * ed_cellCubicNele (4,i,1,1)

     a =      ed_cellCubicNele (2,i,1,1)
     b = x1 * ed_cellCubicNele (3,i,1,1)
     c = x2 * ed_cellCubicNele (4,i,1,1)

     gradNeleX = (a + b + b + c + c + c) * deltaInvX
!
!
!     ...Electron temperature + gradient.
!
!
     Tele =        ed_cellCubicTele (1,i,1,1) &
            + x1 * ed_cellCubicTele (2,i,1,1) &
            + x2 * ed_cellCubicTele (3,i,1,1) &
            + x3 * ed_cellCubicTele (4,i,1,1)

     a =      ed_cellCubicTele (2,i,1,1)
     b = x1 * ed_cellCubicTele (3,i,1,1)
     c = x2 * ed_cellCubicTele (4,i,1,1)

     gradTeleX = (a + b + b + c + c + c) * deltaInvX
!
!
!     ...Check, if obtained values of electron number density and electron temeperature
!        make sense.
!
!
     if (Nele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: Nele <= 0 for a cell (initial)')
     end if

     if (Tele <= 0.0) then
         call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: Tele <= 0 for a cell (initial)')
     end if
!
!
!     ...Get extra needed info about the initial cell (i).
!
!
     cellZbar      = ed_cellZbar    (i,1,1)
     cellDensity   = ed_cellDensity (i,1,1)
     cellVolume    = ed_cellVolume  (i,1,1)
     cellVolumeInv = 1.0 / cellVolume
     cellMass      = cellDensity * cellVolume
     cellMassInv   = 1.0 / cellMass
!
!
!     ...We are ready to follow the ray's path through all the cells of the current
!        block. The current cell index (i) and the previous cell index (ip) will be
!        updated as we move through the block. In case a laser IO is performed on the
!        rays, store the initial ray IO data.
!
!
     if (writeRay) then
         nRayWritePositions = nRayWritePositions + 1
         if(nRayWritePositions <= ed_laserIOMaxNumberOfPositions) then
            ed_laserIONumberOfPositions (rayWriteIndex                           ) = nRayWritePositions
            ed_laserIORayPower          (rayWriteIndex, nRayWritePositions       ) = rayPower
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, IAXIS) = rayX
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = 0.0
            ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
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
!        cell (i) or to one of the 2 cell walls. After this section, we have the correct
!        stepping time and the ray's new position. The ray is considered to stay in the same
!        cell, if it is not found within the cell wall. If it is in the cell wall, the ray
!        is formally between two cells, i.e. crossing to the neighboring cell is considered
!        a possibility.
!
!        The stepping time is determined from the current ray velocity and the (user's) given
!        ray stepping differential. This ray stepping size is a fraction of the cell's dimension
!        and controlls accuracy of the ray's path.
!
!
        stepVelocity = abs (velX)
        stepTime     = ed_cellRayDifferentialStep / stepVelocity
        stepTimeHalf = 0.5 * stepTime

        accX = - c2div2nc * gradNeleX               ! acceleration in x-direction

        newX = rayX + (velX + accX * stepTimeHalf) * stepTime

        newCell =      (newX < xminCell + cellWallThicknessHalf) &
                  .or. (newX > xmaxCell - cellWallThicknessHalf)

        if (newCell) then

            aCoeff (1:2) = 0.5 * accX
            bCoeff (1:2) = velX

            cCoeff (1) = rayX - xminCell
            cCoeff (2) = rayX - xmaxCell

            stepTime = ed_infiniteTime

            do face = 1,2

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
                call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: infinite crossing time for a cell')
            end if

            stepTimeHalf = 0.5 * stepTime

            rayX = rayX + (velX + accX * stepTimeHalf) * stepTime
        else
            rayX = newX
        end if
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

        U =        velX * gradNeleX / Nele
        W =        velX * gradTeleX / Tele
        R = - gradNeleX * gradNeleX * c2div4nc / Nele
        S = - gradNeleX * gradTeleX * c2div4nc / Tele

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
            cellEnergyDepot (i) = cellEnergyDepot (i) + cellEnergy * cellMassInv
        else
            cellEnergyDepot (i) = cellEnergyDepot (i) + cellEnergy * cellVolumeInv
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

        velXeq0 = (velX == 0.0)
 
        stationaryRay = velXeq0

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
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, JAXIS) = 0.0
               ed_laserIORayPositions      (rayWriteIndex, nRayWritePositions, KAXIS) = 0.0
            end if
        end if
!
!
!     ...If, at the current stage, the ray enters a new cell, we have to determine: 1) which new
!        cell (i) it is and 2) the appropriate nudging value on the ray's position. Due to
!        possible reflective boundary conditions on the block faces, it can happen that the ray
!        stays in the original cell. After handling the logistics inside the following 'if'
!        statement, the new cell index i is either the old one or a new one.
!
!
        if (newCell) then

            distToFaceMinX = abs (xminCell - rayX)
            distToFaceMaxX = abs (xmaxCell - rayX)

            minFaceDistance = min (distToFaceMinX, distToFaceMaxX)

            if (minFaceDistance > ed_cellWallThickness) then
                call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: ray to far away from cell face')
            end if

            cellFaceMinX = (distToFaceMinX <= cellWallThicknessHalf)
            cellFaceMaxX = (distToFaceMaxX <= cellWallThicknessHalf)

            velXgt0 = (velX  > 0.0)
            velXlt0 = (velX  < 0.0)

            crossX = .false.
            nudgeX = 0.0

            ip = i

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

            blockFaceMinX = (rayX == xminBlock)
            blockFaceMaxX = (rayX == xmaxBlock)

            reflectX =     (blockFaceMinX .and. blockReflectMinX .and. velXlt0) &
                      .or. (blockFaceMaxX .and. blockReflectMaxX .and. velXgt0)

            if (reflectX) then
                i = ip
                velX = - velX
                crossX = .false.
            end if

            if (crossX) then
                nudgeX = (i - ip) * cellWallThicknessHalf
            end if

            rayX = rayX + nudgeX

            newCell = crossX

        end if
!
!
!     ...We are now sure about the target cell. Check, if the target cell (i) is still within the block.
!        If it is, we check if this is a new cell, in which case we update the cell info. Next we have to
!        calculate the new electron number density and temperature as well as their gradient at the point
!        where the ray is located in the target cell. If the target cell is not within the block, check if
!        the ray coordinates are still within the defined domain. If not, store its latest data and mark
!        it as nonexistent. If the ray is still within the domain boundaries, exit the current block loop.
!
!
        inBlock = (i >= iminBlock) .and. (i <= imaxBlock)

        if (inBlock) then

            if (newCell) then

                xminCell = ed_cellEdges   (i  ,1)
                xmaxCell = ed_cellEdges   (i+1,1)

                cellZbar      = ed_cellZbar    (i,1,1)
                cellDensity   = ed_cellDensity (i,1,1)
                cellVolume    = ed_cellVolume  (i,1,1)
                cellVolumeInv = 1.0 / cellVolume
                cellMass      = cellDensity * cellVolume
                cellMassInv   = 1.0 / cellMass

            end if

            x1 = (rayX - xminCell) * deltaInvX
            x2 = x1 * x1
            x3 = x1 * x2

            Nele =        ed_cellCubicNele (1,i,1,1) &
                   + x1 * ed_cellCubicNele (2,i,1,1) &
                   + x2 * ed_cellCubicNele (3,i,1,1) &
                   + x3 * ed_cellCubicNele (4,i,1,1)

            a =      ed_cellCubicNele (2,i,1,1)
            b = x1 * ed_cellCubicNele (3,i,1,1)
            c = x2 * ed_cellCubicNele (4,i,1,1)

            gradNeleX = (a + b + b + c + c + c) * deltaInvX

            Tele =        ed_cellCubicTele (1,i,1,1) &
                   + x1 * ed_cellCubicTele (2,i,1,1) &
                   + x2 * ed_cellCubicTele (3,i,1,1) &
                   + x3 * ed_cellCubicTele (4,i,1,1)

            a =      ed_cellCubicTele (2,i,1,1)
            b = x1 * ed_cellCubicTele (3,i,1,1)
            c = x2 * ed_cellCubicTele (4,i,1,1)

            gradTeleX = (a + b + b + c + c + c) * deltaInvX

            if (Nele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: Nele <= 0 for a cell (target)')
            end if

            if (Tele <= 0.0) then
                call Driver_abortFlash ('[ed_traceBlockRays1DRec] ERROR: Tele <= 0 for a cell (target)')
            end if

        else

            ed_rays (RAY_POSX,ray) = rayX
            ed_rays (RAY_VELX,ray) = velX
            ed_rays (RAY_POWR,ray) = rayPower

            inDomain = (rayX > ed_xminDomain) .and. (rayX < ed_xmaxDomain)

            if (.not.inDomain) then
                 ed_rays (RAY_POSX,ray) = rayX - nudgeX          ! undo the nudging (not needed anymore)
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
        print *, "[ed_traceBlockRays1DRec] Ray ", ray, &
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
end subroutine ed_traceBlockRays1DRec
