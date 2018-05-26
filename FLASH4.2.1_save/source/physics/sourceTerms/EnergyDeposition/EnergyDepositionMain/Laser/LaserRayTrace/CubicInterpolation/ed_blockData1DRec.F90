!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_blockData1DRec
!!
!! NAME
!!
!!  ed_blockData1DRec
!!
!! SYNOPSIS
!!
!!  call ed_blockData1DRec (integer (in) :: blockID,
!!                          integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: iminData,
!!                          integer (in) :: imaxData,
!!                          integer (in) :: iminDerv,
!!                          integer (in) :: imaxDerv,
!!                          real    (in) :: deltaInvI,
!!                          real    (in) :: blockData (:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally
!!  of 1D rectangular grids (cartesian + spherical). The block is specified through
!!  its number ID and the blockData array, which contains the needed data for the block.
!!  The following is computed and/or stored into specific arrays:
!!
!!     1) the cell Densities
!!     2) the cell Volumes
!!     3) the cell Zbar values
!!     4) the cell center Nele (electron number density) values
!!     5) the cell center Tele (electron temperature) values
!!     6) the cell vertex Nele values
!!     7) the cell vertex Tele values
!!     8) the cell vertex d/dx derivative Nele values
!!     9) the cell vertex d/dx derivative Tele values
!!    10) the cell monocubic expansion coefficients for the vertex Nele grid
!!    11) the cell monocubic expansion coefficients for the vertex Tele grid
!!
!!  The necessary arrays for storage must have been allocated before calling this routine.
!!  No checks are done on the passed cell index limits. This is done before calling this
!!  routine.
!!
!! ARGUMENTS
!!
!!  blockID   : the block ID number
!!  iminBlock : minimum cell i-index limit defining the interior block
!!  imaxBlock : maximum cell i-index limit defining the interior block
!!  iminData  : minimum cell i-index limit needed for evaluating center Nele and Tele values
!!  imaxData  : maximum cell i-index limit needed for evaluating center Nele and Tele values
!!  iminDerv  : minimum cell i-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  imaxDerv  : maximum cell i-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  deltaInvI : inverse of the cell's x-dimension
!!  blockData : two-dimensional array containing the block data
!!
!!***

subroutine ed_blockData1DRec (blockID,              &
                              iminBlock, imaxBlock, &
                              iminData , imaxData,  &
                              iminDerv , imaxDerv,  &
                              deltaInvI,            &
                              blockData             )

  use Driver_interface,                ONLY : Driver_abortFlash

  use Eos_interface,                   ONLY : Eos_getAbarZbar

  use Grid_interface,                  ONLY : Grid_getSingleCellVol

  use ut_cubicInterpolationInterface,  ONLY : ut_monoCubicCoefficients

  use ed_slopeLimiters,                ONLY : ed_mc

  use EnergyDeposition_data,           ONLY : ed_Avogadro,            &
                                              ed_cellDensity,         &
                                              ed_cellCubicNele,       &
                                              ed_cellCubicTele,       &
                                              ed_cellVolume,          &
                                              ed_cellZbar,            &
                                              ed_cellWallThickness

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: iminDerv , imaxDerv
  real,    intent (in) :: deltaInvI
  real,    intent (in) :: blockData (:,:)

  logical :: inBlock

  integer :: i
  integer :: numberOfCells

  integer :: cellIndex (1:3)

  real    :: abar, zbar
  real    :: cellDensity
  real    :: cellVolume
  real    :: centerNele, centerTele
  real    :: notusedDump
  real    :: vm, v0, vp

  real, parameter :: half = 0.5

  real, allocatable :: vertices    (:)
  real, allocatable :: derivatives (:)
!
!
!     ...Dump the arguments not used. This avoids compiler warnings about not using
!        these arguments.
!
!
  notusedDump = deltaInvI
!
!
!     ...Compute total number of interior cells in block.
!
!
  numberOfCells = imaxBlock - iminBlock + 1
!
!
!     ...Allocate the necessary intermediate arrays.
!
!
  allocate (vertices    (iminData : imaxData+1))
  allocate (derivatives (iminDerv : imaxDerv  ))
!
!
!     ...Compute all needed center cell quantities for the # of electrons.
!
!             In block cells: Density, Volume, Zbar and Nele.
!             In guard cells: Nele.
!
!        The center Nele values are not stored, but are immediately processed
!        to compute their contribution to the vertex Nele values. For each
!        cell (i) we compute its 2 vertex Nele values using linear interpolation
!        of the 2 surrounding center Nele values:
!
!               vertex Nele = 1/2 * sum of two surrounding center Nele
!
!        Note, that if all center Nele values are positive, so will be
!        all the vertex Nele values. The storage of the resulting vertex Nele
!        grid is such, that redundant values are avoided. This is achieved by
!        creating a 'vertices' array, which in its (i) position will contain
!        the lower vertex Nele value of the cell (i). Thus for each cell (i) we
!        have the complete set of 2 vertex Nele values stored as follows:
!
!                  (lower vertex of cell)      vertices (i, )
!                  (upper vertex of cell)      vertices (i+1)
!
!        The outermost values of the 'vertices' array will not be accurate, but
!        they are not needed for evaluating the derivatives below.
!
!        The index limiters for the center Nele values must be:
!
!                    iminData = iminBlock - 2
!                    imaxData = imaxBlock + 2
!
!
  vertices = 0.0

  do i = iminData,imaxData

     cellDensity = blockData (DENS_VAR,i)
     call Eos_getAbarZbar (blockData (:,i),    abar, zbar)

     centerNele = cellDensity * ed_Avogadro * zbar / abar

     vertices (i  ) = vertices (i  ) + centerNele
     vertices (i+1) = vertices (i+1) + centerNele

     inBlock = (i >= iminBlock) .and. (i <= imaxBlock)

     if (inBlock) then

         cellIndex (1) = i
         cellIndex (2) = 1
         cellIndex (3) = 1

         call Grid_getSingleCellVol (blockID, EXTERIOR, cellIndex, cellVolume)

         ed_cellDensity (i,1,1) = cellDensity
         ed_cellVolume  (i,1,1) = cellVolume
         ed_cellZbar    (i,1,1) = zbar

     end if
  enddo

  vertices = half * vertices
!
!
!     ...Compute all extra needed derivative values of the vertex Nele grid for all cells
!        of the block. The derivatives are computed using symmetric differences. For getting
!        the cubic expansion coefficients we need the derivatives rescaled to cell coordinates
!        between 0 and 1. Hence, the rescaled d/dx at grid point i is evaluated as
!        (central differencing, h = cell x-dimension = x(i+1)-x(i)):
!
!                  true d/dx at point i = [v(i+1) - v(i-1)] / 2h
!                                       = [v(i+1) - v(i-1)] / [x(i+1) - x(i-1)]
!              rescaled d/dx at point i = [v(i+1) - v(i-1)] / 2
!
!        The index limiters for the derivative values must be:
!
!                    iminDerv = iminBlock
!                    imaxDerv = imaxBlock + 1
!
!
  do i = iminDerv, imaxDerv

     v0 = vertices (i  )
     vm = vertices (i-1)
     vp = vertices (i+1)

     derivatives (i) = ed_mc (v0 - vm , vp - v0)   ! rescaled derivatives need 1/2 factor -> ed_mc has 1/2 factor
                                                   ! -> overall factor = 1
  enddo
!
!
!     ...Scatter the vertices/derivative data into the appropriate cell places for obtaining the
!        monocubic expansion coefficients. The order of the vertices for each cell will be such
!        that it corresponds to the order needed for evaluating the monocubic expansion
!        coefficients:
!
!                                1) lower vertex of cell
!                                2) upper vertex of cell
!
!
!                                  |------------| ---- i
!                                  1            2
!
!
  do i = iminBlock, imaxBlock

     ed_cellCubicNele (1,i,1,1) =    vertices (i  )   ! vertex values
     ed_cellCubicNele (2,i,1,1) =    vertices (i+1)   !
     ed_cellCubicNele (3,i,1,1) = derivatives (i  )   ! d/dx values
     ed_cellCubicNele (4,i,1,1) = derivatives (i+1)   !

  enddo
!
!
!     ...Calculate all Nele monocubic expansion coefficients in one shot.
!
!
  call  ut_monoCubicCoefficients (numberOfCells, ed_cellCubicNele)
!
!
!     ...Repeat the same procedure for the Tele values.
!
!
  vertices = 0.0

  do i = iminData,imaxData
     centerTele     = blockData (TELE_VAR,i)
     vertices (i  ) = vertices (i  ) + centerTele
     vertices (i+1) = vertices (i+1) + centerTele
  enddo

  vertices = half * vertices

  do i = iminDerv, imaxDerv
     v0 = vertices (i  )
     vm = vertices (i-1)
     vp = vertices (i+1)
     derivatives (i) = ed_mc (v0 - vm , vp - v0)
  enddo

  do i = iminBlock, imaxBlock
     ed_cellCubicTele (1,i,1,1) =    vertices (i  )
     ed_cellCubicTele (2,i,1,1) =    vertices (i+1)
     ed_cellCubicTele (3,i,1,1) = derivatives (i  )
     ed_cellCubicTele (4,i,1,1) = derivatives (i+1)
  enddo

  call  ut_monoCubicCoefficients (numberOfCells, ed_cellCubicTele)
!
!
!     ...Deallocate the intermediate arrays.
!
!
  deallocate (vertices)
  deallocate (derivatives)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData1DRec
