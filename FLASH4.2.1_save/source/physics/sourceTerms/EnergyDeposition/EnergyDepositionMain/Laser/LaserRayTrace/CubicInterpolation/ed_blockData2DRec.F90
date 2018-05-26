!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_blockData2DRec
!!
!! NAME
!!
!!  ed_blockData2DRec
!!
!! SYNOPSIS
!!
!!  call ed_blockData2DRec (integer (in) :: blockID,
!!                          integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: jminBlock,
!!                          integer (in) :: jmaxBlock,
!!                          integer (in) :: iminData,
!!                          integer (in) :: imaxData,
!!                          integer (in) :: jminData,
!!                          integer (in) :: jmaxData,
!!                          integer (in) :: iminDerv,
!!                          integer (in) :: imaxDerv,
!!                          integer (in) :: jminDerv,
!!                          integer (in) :: jmaxDerv,
!!                          real    (in) :: deltaI,
!!                          real    (in) :: deltaJ,
!!                          real    (in) :: deltaInvI,
!!                          real    (in) :: deltaInvJ,
!!                          real    (in) :: blockData (:,:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally
!!  of 2D rectangular grids (cartesian + cylindrical). The block is specified through
!!  its number ID and the blockData array, which contains the needed data for the block.
!!  The following is computed and/or stored into specific arrays:
!!
!!     1) the cell Densities
!!     2) the cell Zbar values
!!     3) the cell center Nele (electron number density) values
!!     4) the cell center Tele (electron temperature) values
!!     5) the cell vertex Nele values
!!     6) the cell vertex Tele values
!!     7) the cell vertex (mixed, up to 2nd order) derivative Nele values
!!     8) the cell vertex (mixed, up to 2nd order) derivative Tele values
!!     9) the cell bicubic expansion coefficients for the vertex Nele grid
!!    10) the cell bicubic expansion coefficients for the vertex Tele grid
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
!!  jminBlock : minimum cell j-index limit defining the interior block
!!  jmaxBlock : maximum cell j-index limit defining the interior block
!!  iminData  : minimum cell i-index limit needed for evaluating center Nele and Tele values
!!  imaxData  : maximum cell i-index limit needed for evaluating center Nele and Tele values
!!  jminData  : minimum cell j-index limit needed for evaluating center Nele and Tele values
!!  jmaxData  : maximum cell j-index limit needed for evaluating center Nele and Tele values
!!  iminDerv  : minimum cell i-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  imaxDerv  : maximum cell i-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  jminDerv  : minimum cell j-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  jmaxDerv  : maximum cell j-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  deltaI    : the cell's x-dimension
!!  deltaJ    : the cell's y-dimension
!!  deltaInvI : inverse of the cell's x-dimension
!!  deltaInvJ : inverse of the cell's y-dimension
!!  blockData : three-dimensional array containing the block data
!!
!!***

subroutine ed_blockData2DRec (blockID,              &
                              iminBlock, imaxBlock, &
                              jminBlock, jmaxBlock, &
                              iminData , imaxData,  &
                              jminData , jmaxData,  &
                              iminDerv , imaxDerv,  &
                              jminDerv , jmaxDerv,  &
                              deltaI   , deltaJ,    &
                              deltaInvI, deltaInvJ, &
                              blockData             )

  use Driver_interface,                ONLY : Driver_abortFlash

  use Eos_interface,                   ONLY : Eos_getAbarZbar

  use Grid_interface,                  ONLY : Grid_getSingleCellVol

  use ut_cubicInterpolationInterface,  ONLY : ut_biCubicCoefficients

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
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: jminData , jmaxData
  integer, intent (in) :: iminDerv , imaxDerv
  integer, intent (in) :: jminDerv , jmaxDerv
  real,    intent (in) :: deltaI   , deltaJ
  real,    intent (in) :: deltaInvI, deltaInvJ
  real,    intent (in) :: blockData (:,:,:)

  logical :: inBlock

  integer :: i, im, ip
  integer :: j, jm, jp
  integer :: numberOfCells

  integer :: cellIndex (1:3)

  real    :: abar, zbar
  real    :: cellDensity
  real    :: cellVolume
  real    :: centerNele, centerTele
  real    :: deltaInvIJ
  real    :: dx, dy, dxdy
  real    :: notusedDump
  real    :: vmm, vm0, vmp, &
             v0m, v00, v0p, &
             vpm, vp0, vpp

  real, parameter :: fourth = 0.25

  real, allocatable :: vertices    (:,:)
  real, allocatable :: derivatives (:,:,:)
!
!
!     ...Dump the arguments not used. This avoids compiler warnings about not using
!        these arguments.
!
!
  notusedDump = deltaI
  notusedDump = deltaJ
!
!
!     ...Compute necessary variables.
!
!
  deltaInvIJ  =  deltaInvI * deltaInvJ
!
!
!     ...Compute total number of interior cells in block.
!
!
  numberOfCells =   (imaxBlock - iminBlock + 1) &
                  * (jmaxBlock - jminBlock + 1)
!
!
!     ...Allocate the necessary intermediate arrays.
!
!
  allocate (vertices    (     iminData : imaxData+1 , jminData : jmaxData+1))
  allocate (derivatives (1:3, iminDerv : imaxDerv   , jminDerv : jmaxDerv  ))
!
!
!     ...Compute all needed center cell quantities for the # of electrons.
!
!             In block cells: Density, Volume, Zbar and Nele.
!             In guard cells: Nele.
!
!        The center Nele values are not stored, but are immediately processed
!        to compute their contribution to the vertex Nele values. For each
!        cell (i,j) we compute its 4 vertex Nele values using bilinear
!        interpolation of the 4 surrounding center Nele values:
!
!               vertex Nele = 1/4 * sum of all surrounding center Nele
!
!        Note, that if all center Nele values are positive, so will be
!        all the vertex Nele values. The storage of the resulting vertex Nele
!        grid is such, that redundant values are avoided. This is achieved by
!        creating a 'vertices' array, which in its (i,j) position will contain
!        the lower left corner vertex Nele value of the cel (i,j). Thus for
!        each cell (i,j) we have the complete set of 4 vertex Nele values
!        stored as follows:
!
!                0 0  (lower left vertex of cell)      vertices (i,  j, )
!                1 0  (vertex on i-axis of cell)       vertices (i+1,j, )
!                0 1  (vertex on j-axis of cell)       vertices (i,  j+1)
!                1 1  (vertex on ij-plane of cell)     vertices (i+1,j+1)
!
!        The outermost values of the 'vertices' array will not be accurate, but
!        they are not needed for evaluating the (mixed) derivatives below.
!
!        The index limiters for the center Nele values must be:
!
!                    iminData = iminBlock - 2
!                    imaxData = imaxBlock + 2
!                    jminData = jminBlock - 2
!                    jmaxData = jmaxBlock + 2
!
!
  vertices = 0.0

  do j = jminData,jmaxData
     jp = j + 1
     do i = iminData,imaxData
        ip = i + 1

        cellDensity = blockData (DENS_VAR,i,j)
        call Eos_getAbarZbar (blockData (:,i,j),    abar, zbar)

        centerNele = cellDensity * ed_Avogadro * zbar / abar

        vertices (i , j ) = vertices (i , j ) + centerNele
        vertices (ip, j ) = vertices (ip, j ) + centerNele
        vertices (i , jp) = vertices (i , jp) + centerNele
        vertices (ip, jp) = vertices (ip, jp) + centerNele

        inBlock =       (i >= iminBlock) &
                  .and. (i <= imaxBlock) &
                  .and. (j >= jminBlock) &
                  .and. (j <= jmaxBlock)

        if (inBlock) then

            cellIndex (1) = i
            cellIndex (2) = j
            cellIndex (3) = 1

            call Grid_getSingleCellVol (blockID, EXTERIOR, cellIndex, cellVolume)

            ed_cellDensity (i,j,1) = cellDensity
            ed_cellVolume  (i,j,1) = cellVolume
            ed_cellZbar    (i,j,1) = zbar
        end if

     enddo
  enddo

  vertices = fourth * vertices
!
!
!     ...Compute all extra needed derivative values of the vertex Nele grid for all cells
!        of the block. The derivatives are computed using symmetric differences. For getting
!        the cubic expansion coefficients we need the derivatives rescaled to cell coordinates
!        between 0 and 1. Hence, for example the rescaled d/dx at grid point i is evaluated
!        as (central differencing, h = cell x-dimension = x(i+1)-x(i)):
!
!                  true d/dx at point i = [v(i+1) - v(i-1)] / 2h
!                                       = [v(i+1) - v(i-1)] / [x(i+1) - x(i-1)]
!              rescaled d/dx at point i = [v(i+1) - v(i-1)] / 2
!
!        The index limiters for the derivative values must be:
!
!                    iminDerv = iminBlock
!                    imaxDerv = imaxBlock + 1
!                    jminDerv = jminBlock
!                    jmaxDerv = jmaxBlock + 1
!
!
  do j = jminDerv, jmaxDerv
     jm = j - 1
     jp = j + 1
     do i = iminDerv, imaxDerv
        im = i - 1
        ip = i + 1

        v00 = vertices (i , j )
        vm0 = vertices (im, j )
        vp0 = vertices (ip, j )
        v0m = vertices (i , jm)
        vmm = vertices (im, jm)
        vpm = vertices (ip, jm)
        v0p = vertices (i , jp)
        vmp = vertices (im, jp)
        vpp = vertices (ip, jp)

        dx = ed_mc (v00 - vm0 , vp0 - v00)                  ! rescaled derivatives need 1/2 factor ->
        dy = ed_mc (v00 - v0m , v0p - v00)                  ! ed_mc has 1/2 factor -> overall factor = 1

        dxdy = (  ed_mc (v0p - vmp , vpp - v0p) &           ! rescaled derivatives need 1/4
                - ed_mc (v0m - vmm , vpm - v0m) &           ! factor -> ed_mc has 1/2 factor
                + ed_mc (vp0 - vpm , vpp - vp0) &           ! and relevant vertices like vpp
                - ed_mc (vm0 - vmm , vmp - vm0) ) * fourth  ! occur 2x -> overall factor = 1/4

        derivatives (1,i,j) = dx
        derivatives (2,i,j) = dy
        derivatives (3,i,j) = dxdy

     enddo
  enddo
!
!
!     ...Scatter the vertices/derivative data into the appropriate cell places for obtaining the
!        bicubic expansion coefficients. The order of the vertices for each cell will be such
!        that it corresponds to the order needed for evaluating the bicubic expansion
!        coefficients:
!
!
!                     1)    0 0  (lower left vertex of cell)
!                     2)    1 0  (vertex on i-axis of cell)
!                     3)    0 1  (vertex on j-axis of cell)
!                     4)    1 1  (vertex on ij-plane of cell)
!
!
!                                   j
! 
!                                   |
!                                   |
!                                   |
!                                  3 ----------- 4
!                                   |           |
!                                   |           |
!                                   |           |
!                                   |           |
!                                   |           |
!                                    ----------- ---- i
!                                  1             2
!
!
  do j = jminBlock, jmaxBlock
     jp = j + 1
     do i = iminBlock, imaxBlock
        ip = i + 1

        ed_cellCubicNele ( 1,i,j,1) =    vertices (   i , j )   ! vertex values
        ed_cellCubicNele ( 2,i,j,1) =    vertices (   ip, j )   !
        ed_cellCubicNele ( 3,i,j,1) =    vertices (   i , jp)   !
        ed_cellCubicNele ( 4,i,j,1) =    vertices (   ip, jp)   !

        ed_cellCubicNele ( 5,i,j,1) = derivatives (1, i , j )   ! d/dx values
        ed_cellCubicNele ( 6,i,j,1) = derivatives (1, ip, j )   !
        ed_cellCubicNele ( 7,i,j,1) = derivatives (1, i , jp)   !
        ed_cellCubicNele ( 8,i,j,1) = derivatives (1, ip, jp)   !

        ed_cellCubicNele ( 9,i,j,1) = derivatives (2, i , j )   ! d/dy values
        ed_cellCubicNele (10,i,j,1) = derivatives (2, ip, j )   !
        ed_cellCubicNele (11,i,j,1) = derivatives (2, i , jp)   !
        ed_cellCubicNele (12,i,j,1) = derivatives (2, ip, jp)   !

        ed_cellCubicNele (13,i,j,1) = derivatives (3, i , j )   ! d2/dxdy values
        ed_cellCubicNele (14,i,j,1) = derivatives (3, ip, j )   !
        ed_cellCubicNele (15,i,j,1) = derivatives (3, i , jp)   !
        ed_cellCubicNele (16,i,j,1) = derivatives (3, ip, jp)   !

     enddo
  enddo
!
!
!     ...Calculate all Nele bicubic expansion coefficients in one shot.
!
!
  call  ut_biCubicCoefficients (numberOfCells, ed_cellCubicNele)
!
!
!     ...Repeat the same procedure for the Tele values.
!
!
  vertices = 0.0

  do j = jminData,jmaxData
     jp = j + 1
     do i = iminData,imaxData
        ip = i + 1

        centerTele = blockData (TELE_VAR,i,j)

        vertices (i , j ) = vertices (i , j ) + centerTele
        vertices (ip, j ) = vertices (ip, j ) + centerTele
        vertices (i , jp) = vertices (i , jp) + centerTele
        vertices (ip, jp) = vertices (ip, jp) + centerTele

     enddo
  enddo

  vertices = fourth * vertices

  do j = jminDerv, jmaxDerv
     jm = j - 1
     jp = j + 1
     do i = iminDerv, imaxDerv
        im = i - 1
        ip = i + 1

        v00 = vertices (i , j )
        vm0 = vertices (im, j )
        vp0 = vertices (ip, j )
        v0m = vertices (i , jm)
        vmm = vertices (im, jm)
        vpm = vertices (ip, jm)
        v0p = vertices (i , jp)
        vmp = vertices (im, jp)
        vpp = vertices (ip, jp)

        dx = ed_mc (v00 - vm0 , vp0 - v00)
        dy = ed_mc (v00 - v0m , v0p - v00)

        dxdy = (  ed_mc (v0p - vmp , vpp - v0p) &
                - ed_mc (v0m - vmm , vpm - v0m) &
                + ed_mc (vp0 - vpm , vpp - vp0) &
                - ed_mc (vm0 - vmm , vmp - vm0) ) * fourth

        derivatives (1,i,j) = dx
        derivatives (2,i,j) = dy
        derivatives (3,i,j) = dxdy

     enddo
  enddo

  do j = jminBlock, jmaxBlock
     jp = j + 1
     do i = iminBlock, imaxBlock
        ip = i + 1

        ed_cellCubicTele ( 1,i,j,1) =    vertices (   i , j )   ! vertex values
        ed_cellCubicTele ( 2,i,j,1) =    vertices (   ip, j )   !
        ed_cellCubicTele ( 3,i,j,1) =    vertices (   i , jp)   !
        ed_cellCubicTele ( 4,i,j,1) =    vertices (   ip, jp)   !

        ed_cellCubicTele ( 5,i,j,1) = derivatives (1, i , j )   ! d/dx values
        ed_cellCubicTele ( 6,i,j,1) = derivatives (1, ip, j )   !
        ed_cellCubicTele ( 7,i,j,1) = derivatives (1, i , jp)   !
        ed_cellCubicTele ( 8,i,j,1) = derivatives (1, ip, jp)   !

        ed_cellCubicTele ( 9,i,j,1) = derivatives (2, i , j )   ! d/dy values
        ed_cellCubicTele (10,i,j,1) = derivatives (2, ip, j )   !
        ed_cellCubicTele (11,i,j,1) = derivatives (2, i , jp)   !
        ed_cellCubicTele (12,i,j,1) = derivatives (2, ip, jp)   !

        ed_cellCubicTele (13,i,j,1) = derivatives (3, i , j )   ! d2/dxdy values
        ed_cellCubicTele (14,i,j,1) = derivatives (3, ip, j )   !
        ed_cellCubicTele (15,i,j,1) = derivatives (3, i , jp)   !
        ed_cellCubicTele (16,i,j,1) = derivatives (3, ip, jp)   !

     enddo
  enddo

  call  ut_biCubicCoefficients (numberOfCells, ed_cellCubicTele)
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
end subroutine ed_blockData2DRec
