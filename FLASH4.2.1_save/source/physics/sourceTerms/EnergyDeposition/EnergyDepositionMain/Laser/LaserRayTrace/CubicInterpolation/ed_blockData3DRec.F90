!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_blockData3DRec
!!
!! NAME
!!
!!  ed_blockData3DRec
!!
!! SYNOPSIS
!!
!!  call ed_blockData3DRec (integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: jminBlock,
!!                          integer (in) :: jmaxBlock,
!!                          integer (in) :: kminBlock,
!!                          integer (in) :: kmaxBlock,
!!                          integer (in) :: iminData,
!!                          integer (in) :: imaxData,
!!                          integer (in) :: jminData,
!!                          integer (in) :: jmaxData,
!!                          integer (in) :: kminData,
!!                          integer (in) :: kmaxData,
!!                          integer (in) :: iminDerv,
!!                          integer (in) :: imaxDerv,
!!                          integer (in) :: jminDerv,
!!                          integer (in) :: jmaxDerv,
!!                          integer (in) :: kminDerv,
!!                          integer (in) :: kmaxDerv,
!!                          real    (in) :: deltaI,
!!                          real    (in) :: deltaJ,
!!                          real    (in) :: deltaK,
!!                          real    (in) :: deltaInvI,
!!                          real    (in) :: deltaInvJ,
!!                          real    (in) :: deltaInvK,
!!                          real    (in) :: blockData (:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally of
!!  3D rectangular grids (cartesian). The block is specified through its number ID and the
!!  blockData array, which contains the needed data for the block. The following is computed
!!  (and some stored) into specific arrays:
!!
!!     1) the cell Densities                                                  (stored)
!!     2) the cell Zbar values                                                (stored)
!!     3) the cell center Nele (electron number density) values               (computed)
!!     4) the cell center Tele (electron temperature) values                  (computed)
!!     5) the cell vertex Nele values                                         (computed)
!!     6) the cell vertex Tele values                                         (computed)
!!     7) the cell vertex (mixed, up to 3rd order) derivative Nele values     (computed)
!!     8) the cell vertex (mixed, up to 3rd order) derivative Tele values     (computed)
!!     9) the cell tricubic expansion coefficients for the vertex Nele grid   (stored)
!!    10) the cell tricubic expansion coefficients for the vertex Tele grid   (stored)
!!
!!  The necessary arrays for storage must have been allocated before calling this routine.
!!  No checks are done on the passed cell index limits. This is done before calling this
!!  routine.
!!
!! ARGUMENTS
!!
!!  iminBlock : minimum cell i-index limit defining the interior block
!!  imaxBlock : maximum cell i-index limit defining the interior block
!!  jminBlock : minimum cell j-index limit defining the interior block
!!  jmaxBlock : maximum cell j-index limit defining the interior block
!!  kminBlock : minimum cell k-index limit defining the interior block
!!  kmaxBlock : maximum cell k-index limit defining the interior block
!!  iminData  : minimum cell i-index limit needed for evaluating center Nele and Tele values
!!  imaxData  : maximum cell i-index limit needed for evaluating center Nele and Tele values
!!  jminData  : minimum cell j-index limit needed for evaluating center Nele and Tele values
!!  jmaxData  : maximum cell j-index limit needed for evaluating center Nele and Tele values
!!  kminData  : minimum cell k-index limit needed for evaluating center Nele and Tele values
!!  kmaxData  : maximum cell k-index limit needed for evaluating center Nele and Tele values
!!  iminDerv  : minimum cell i-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  imaxDerv  : maximum cell i-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  jminDerv  : minimum cell j-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  jmaxDerv  : maximum cell j-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  kminDerv  : minimum cell k-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  kmaxDerv  : maximum cell k-index limit needed for evaluating derivatives of vertex Nele and Tele
!!  deltaI    : the cell's x-dimension
!!  deltaJ    : the cell's y-dimension
!!  deltaK    : the cell's z-dimension
!!  deltaInvI : inverse of the cell's x-dimension
!!  deltaInvJ : inverse of the cell's y-dimension
!!  deltaInvK : inverse of the cell's z-dimension
!!  blockData : four-dimensional array containing the block data
!!
!!***

subroutine ed_blockData3DRec (iminBlock, imaxBlock,            &
                              jminBlock, jmaxBlock,            &
                              kminBlock, kmaxBlock,            &
                              iminData , imaxData,             &
                              jminData , jmaxData,             &
                              kminData , kmaxData,             &
                              iminDerv , imaxDerv,             &
                              jminDerv , jmaxDerv,             &
                              kminDerv , kmaxDerv,             &
                              deltaI   , deltaJ   , deltaK,    &
                              deltaInvI, deltaInvJ, deltaInvK, &
                              blockData                        )

  use Driver_interface,                ONLY : Driver_abortFlash

  use Eos_interface,                   ONLY : Eos_getAbarZbar

  use ut_cubicInterpolationInterface,  ONLY : ut_triCubicCoefficients

  use ed_slopeLimiters,                ONLY : ed_mc

  use EnergyDeposition_data,           ONLY : ed_Avogadro,          &
                                              ed_cellCubicNele,     &
                                              ed_cellCubicTele,     &
                                              ed_cellDensity,       &
                                              ed_cellZbar

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: kminBlock, kmaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: jminData , jmaxData
  integer, intent (in) :: kminData , kmaxData
  integer, intent (in) :: iminDerv , imaxDerv
  integer, intent (in) :: jminDerv , jmaxDerv
  integer, intent (in) :: kminDerv , kmaxDerv
  real,    intent (in) :: deltaI   , deltaJ   , deltaK
  real,    intent (in) :: deltaInvI, deltaInvJ, deltaInvK
  real,    intent (in) :: blockData (:,:,:,:)

  logical :: inBlock

  integer :: i, im, ip
  integer :: j, jm, jp
  integer :: k, km, kp
  integer :: numberOfCells

  real    :: abar, zbar
  real    :: cellDensity
  real    :: centerNele, centerTele
  real    :: deltaInvIJ, deltaInvIK, deltaInvJK, deltaInvIJK
  real    :: dx, dy, dz, dxdy, dxdz, dydz, dxdydz
  real    :: notusedDump
  real    :: vmmm, vmm0, vmmp, vm0m, vm00, vm0p, vmpm, vmp0, vmpp, &
             v0mm, v0m0, v0mp, v00m, v000, v00p, v0pm, v0p0, v0pp, &
             vpmm, vpm0, vpmp, vp0m, vp00, vp0p, vppm, vpp0, vppp

  real, parameter :: fourth  = 0.25
  real, parameter :: eighth  = 0.125
  real, parameter :: twelfth = 1.0 / 12.0

  real, allocatable :: vertices    (:,:,:)
  real, allocatable :: derivatives (:,:,:,:)
!
!
!     ...Dump the arguments not used. This avoids compiler warnings about not using
!        these arguments.
!
!
  notusedDump = deltaI
  notusedDump = deltaJ
  notusedDump = deltaK
!
!
!     ...Compute necessary variables.
!
!
  deltaInvIJ  =  deltaInvI * deltaInvJ
  deltaInvIK  =  deltaInvI * deltaInvK
  deltaInvJK  =  deltaInvJ * deltaInvK
  deltaInvIJK = deltaInvIJ * deltaInvK
!
!
!     ...Compute total number of interior cells in block.
!
!
  numberOfCells =   (imaxBlock - iminBlock + 1) &
                  * (jmaxBlock - jminBlock + 1) &
                  * (kmaxBlock - kminBlock + 1)
!
!
!     ...Allocate the necessary intermediate arrays.
!
!
  allocate (vertices    (     iminData : imaxData+1 , jminData : jmaxData+1 , kminData : kmaxData+1))
  allocate (derivatives (1:7, iminDerv : imaxDerv   , jminDerv : jmaxDerv   , kminDerv : kmaxDerv  ))
!
!
!     ...Compute all needed center cell quantities for the # of electrons.
!
!             In block cells: Density, Zbar and Nele.
!             In guard cells: Nele.
!
!        The center Nele values are not stored, but are immediately processed
!        to compute their contribution to the vertex Nele values. For each
!        cell (i,j,k) we compute its 8 vertex Nele values using trilinear
!        interpolation of the 8 surrounding center Nele values:
!
!               vertex Nele = 1/8 * sum of all surrounding center Nele
!
!        Note, that if all center Nele values are positive, so will be
!        all the vertex Nele values. The storage of the resulting vertex Nele
!        grid is such, that redundant values are avoided. This is achieved by
!        creating a 'vertices' array, which in its (i,j,k) position will contain
!        the lower left corner vertex Nele value of the cel (i,j,k). Thus for
!        each cell (i,j,k) we have the complete set of 8 vertex Nele values
!        stored as follows:
!
!                0 0 0  (lower left vertex of cell)      vertices (i,  j,  k  )
!                1 0 0  (vertex on i-axis of cell)       vertices (i+1,j,  k  )
!                0 1 0  (vertex on j-axis of cell)       vertices (i,  j+1,k  )
!                0 0 1  (vertex on k-axis of cell)       vertices (i,  j,  k+1)
!                1 1 0  (vertex on ij-plane of cell)     vertices (i+1,j+1,k  )
!                1 0 1  (vertex on ik-plane of cell)     vertices (i+1,j,  k+1)
!                0 1 1  (vertex on jk-plane of cell)     vertices (i,  j+1,k+1)
!                1 1 1  (upper right vertex of cell)     vertices (i+1,j+1,k+1)
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
!                    kminData = kminBlock - 2
!                    kmaxData = kmaxBlock + 2
!
!
  vertices = 0.0

  do k = kminData,kmaxData
     kp = k + 1
     do j = jminData,jmaxData
        jp = j + 1
        do i = iminData,imaxData
           ip = i + 1

           cellDensity = blockData (DENS_VAR,i,j,k)
           call Eos_getAbarZbar (blockData (:,i,j,k),    abar, zbar)

           centerNele = cellDensity * ed_Avogadro * zbar / abar

           vertices (i , j , k ) = vertices (i , j , k ) + centerNele
           vertices (ip, j , k ) = vertices (ip, j , k ) + centerNele
           vertices (i , jp, k ) = vertices (i , jp, k ) + centerNele
           vertices (i , j , kp) = vertices (i , j , kp) + centerNele
           vertices (ip, jp, k ) = vertices (ip, jp, k ) + centerNele
           vertices (ip, j , kp) = vertices (ip, j , kp) + centerNele
           vertices (i , jp, kp) = vertices (i , jp, kp) + centerNele
           vertices (ip, jp, kp) = vertices (ip, jp, kp) + centerNele

           inBlock =       (i >= iminBlock) &
                     .and. (i <= imaxBlock) &
                     .and. (j >= jminBlock) &
                     .and. (j <= jmaxBlock) &
                     .and. (k >= kminBlock) &
                     .and. (k <= kmaxBlock)

           if (inBlock) then
               ed_cellDensity (i,j,k) = cellDensity
               ed_cellZbar    (i,j,k) = zbar
           end if

        enddo
     enddo
  enddo

  vertices = eighth * vertices
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
!                    kminDerv = kminBlock
!                    kmaxDerv = kmaxBlock + 1
!
!
  do k = kminDerv, kmaxDerv
     km = k - 1
     kp = k + 1
     do j = jminDerv, jmaxDerv
        jm = j - 1
        jp = j + 1
        do i = iminDerv, imaxDerv
           im = i - 1
           ip = i + 1

           v000 = vertices (i , j , k )
           vm00 = vertices (im, j , k )
           vp00 = vertices (ip, j , k )
           v0m0 = vertices (i , jm, k )
           vmm0 = vertices (im, jm, k )
           vpm0 = vertices (ip, jm, k )
           v0p0 = vertices (i , jp, k )
           vmp0 = vertices (im, jp, k )
           vpp0 = vertices (ip, jp, k )
           v00m = vertices (i , j , km)
           vm0m = vertices (im, j , km)
           vp0m = vertices (ip, j , km)
           v0mm = vertices (i , jm, km)
           vmmm = vertices (im, jm, km)
           vpmm = vertices (ip, jm, km)
           v0pm = vertices (i , jp, km)
           vmpm = vertices (im, jp, km)
           vppm = vertices (ip, jp, km)
           v00p = vertices (i , j , kp)
           vm0p = vertices (im, j , kp)
           vp0p = vertices (ip, j , kp)
           v0mp = vertices (i , jm, kp)
           vmmp = vertices (im, jm, kp)
           vpmp = vertices (ip, jm, kp)
           v0pp = vertices (i , jp, kp)
           vmpp = vertices (im, jp, kp)
           vppp = vertices (ip, jp, kp)

           dx = ed_mc (v000 - vm00 , vp00 - v000)                  ! rescaled derivatives need 1/2
           dy = ed_mc (v000 - v0m0 , v0p0 - v000)                  ! factor -> ed_mc has 1/2 factor
           dz = ed_mc (v000 - v00m , v00p - v000)                  ! -> overall factor = 1

           dxdy = (  ed_mc (v0p0 - vmp0 , vpp0 - v0p0) &           ! rescaled derivatives need 1/4
                   - ed_mc (v0m0 - vmm0 , vpm0 - v0m0) &           ! factor -> ed_mc has 1/2 factor
                   + ed_mc (vp00 - vpm0 , vpp0 - vp00) &           ! and relevant vertices like vpp0
                   - ed_mc (vm00 - vmm0 , vmp0 - vm00) ) * fourth  ! occur 2x -> overall factor = 1/4

           dxdz = (  ed_mc (v00p - vm0p , vp0p - v00p) &
                   - ed_mc (v00m - vm0m , vp0m - v00m) &
                   + ed_mc (vp00 - vp0m , vp0p - vp00) &
                   - ed_mc (vm00 - vm0m , vm0p - vm00) ) * fourth

           dydz = (  ed_mc (v00p - v0mp , v0pp - v00p) &
                   - ed_mc (v00m - v0mm , v0pm - v00m) &
                   + ed_mc (v0p0 - v0pm , v0pp - v0p0) &
                   - ed_mc (v0m0 - v0mm , v0mp - v0m0) ) * fourth

           dxdydz = (  ed_mc (v0pp - vmpp , vppp - v0pp) &         ! rescaled derivatives need 1/8
                     - ed_mc (v0pm - vmpm , vppm - v0pm) &         ! factor -> ed_mc has 1/2 factor
                     - ed_mc (v0mp - vmmp , vpmp - v0mp) &         ! and relevant vertices like vppp
                     + ed_mc (v0mm - vmmm , vpmm - v0mm) &         ! occur 3x -> overall factor = 1/12
                     + ed_mc (vp0p - vpmp , vppp - vp0p) &
                     - ed_mc (vp0m - vpmm , vppm - vp0m) &
                     - ed_mc (vm0p - vmmp , vmpp - vm0p) &
                     + ed_mc (vm0m - vmmm , vmpm - vm0m) &
                     + ed_mc (vpp0 - vppm , vppp - vpp0) &
                     - ed_mc (vpm0 - vpmm , vpmp - vpm0) &
                     - ed_mc (vmp0 - vmpm , vmpp - vmp0) &
                     + ed_mc (vmm0 - vmmm , vmmp - vmm0) ) * twelfth

           derivatives (1,i,j,k) = dx
           derivatives (2,i,j,k) = dy
           derivatives (3,i,j,k) = dz
           derivatives (4,i,j,k) = dxdy
           derivatives (5,i,j,k) = dxdz
           derivatives (6,i,j,k) = dydz
           derivatives (7,i,j,k) = dxdydz

        enddo
     enddo
  enddo
!
!
!     ...Scatter the vertices/derivative data into the appropriate cell places for obtaining the
!        tricubic expansion coefficients. The order of the vertices for each cell will be such
!        that it corresponds to the order needed for evaluating the tricubic expansion
!        coefficients:
!
!
!                   1)       0 0 0  (lower left vertex of cell)
!                   2)       1 0 0  (vertex on i-axis of cell)
!                   3)       0 1 0  (vertex on j-axis of cell)
!                   4)       0 0 1  (vertex on k-axis of cell)
!                   5)       1 1 0  (vertex on ij-plane of cell)
!                   6)       1 0 1  (vertex on ik-plane of cell)
!                   7)       0 1 1  (vertex on jk-plane of cell)
!                   8)       1 1 1  (upper right vertex of cell)
!
!
!                             k
!                             
!                             |
!                             | 7  -----------  8
!                             |  /|          /|
!                               / |         / |
!                              /  |        /  |
!                           4  -----------  6 |
!                             |   | j     |   |
!                             | 3  -----------  5
!                             |  /        |  /
!                             | /         | /
!                             |/          |/
!                              -----------    ---- i
!                           1              2
!
!
  do k = kminBlock, kmaxBlock
     kp = k + 1
     do j = jminBlock, jmaxBlock
        jp = j + 1
        do i = iminBlock, imaxBlock
           ip = i + 1

           ed_cellCubicNele ( 1,i,j,k) =    vertices (   i , j , k )   ! vertex values
           ed_cellCubicNele ( 2,i,j,k) =    vertices (   ip, j  ,k )   !
           ed_cellCubicNele ( 3,i,j,k) =    vertices (   i , jp ,k )   !
           ed_cellCubicNele ( 4,i,j,k) =    vertices (   i , j  ,kp)   !
           ed_cellCubicNele ( 5,i,j,k) =    vertices (   ip, jp ,k )   !
           ed_cellCubicNele ( 6,i,j,k) =    vertices (   ip, j  ,kp)   !
           ed_cellCubicNele ( 7,i,j,k) =    vertices (   i , jp ,kp)   !
           ed_cellCubicNele ( 8,i,j,k) =    vertices (   ip, jp ,kp)   !

           ed_cellCubicNele ( 9,i,j,k) = derivatives (1, i , j , k )   ! d/dx values
           ed_cellCubicNele (10,i,j,k) = derivatives (1, ip, j  ,k )   !
           ed_cellCubicNele (11,i,j,k) = derivatives (1, i , jp ,k )   !
           ed_cellCubicNele (12,i,j,k) = derivatives (1, i , j  ,kp)   !
           ed_cellCubicNele (13,i,j,k) = derivatives (1, ip, jp ,k )   !
           ed_cellCubicNele (14,i,j,k) = derivatives (1, ip, j  ,kp)   !
           ed_cellCubicNele (15,i,j,k) = derivatives (1, i , jp ,kp)   !
           ed_cellCubicNele (16,i,j,k) = derivatives (1, ip, jp ,kp)   !

           ed_cellCubicNele (17,i,j,k) = derivatives (2, i , j , k )   ! d/dy values
           ed_cellCubicNele (18,i,j,k) = derivatives (2, ip, j  ,k )   !
           ed_cellCubicNele (19,i,j,k) = derivatives (2, i , jp ,k )   !
           ed_cellCubicNele (20,i,j,k) = derivatives (2, i , j  ,kp)   !
           ed_cellCubicNele (21,i,j,k) = derivatives (2, ip, jp ,k )   !
           ed_cellCubicNele (22,i,j,k) = derivatives (2, ip, j  ,kp)   !
           ed_cellCubicNele (23,i,j,k) = derivatives (2, i , jp ,kp)   !
           ed_cellCubicNele (24,i,j,k) = derivatives (2, ip, jp ,kp)   !

           ed_cellCubicNele (25,i,j,k) = derivatives (3, i , j , k )   ! d/dz values
           ed_cellCubicNele (26,i,j,k) = derivatives (3, ip, j  ,k )   !
           ed_cellCubicNele (27,i,j,k) = derivatives (3, i , jp ,k )   !
           ed_cellCubicNele (28,i,j,k) = derivatives (3, i , j  ,kp)   !
           ed_cellCubicNele (29,i,j,k) = derivatives (3, ip, jp ,k )   !
           ed_cellCubicNele (30,i,j,k) = derivatives (3, ip, j  ,kp)   !
           ed_cellCubicNele (31,i,j,k) = derivatives (3, i , jp ,kp)   !
           ed_cellCubicNele (32,i,j,k) = derivatives (3, ip, jp ,kp)   !

           ed_cellCubicNele (33,i,j,k) = derivatives (4, i , j , k )   ! d2/dxdy values
           ed_cellCubicNele (34,i,j,k) = derivatives (4, ip, j  ,k )   !
           ed_cellCubicNele (35,i,j,k) = derivatives (4, i , jp ,k )   !
           ed_cellCubicNele (36,i,j,k) = derivatives (4, i , j  ,kp)   !
           ed_cellCubicNele (37,i,j,k) = derivatives (4, ip, jp ,k )   !
           ed_cellCubicNele (38,i,j,k) = derivatives (4, ip, j  ,kp)   !
           ed_cellCubicNele (39,i,j,k) = derivatives (4, i , jp ,kp)   !
           ed_cellCubicNele (40,i,j,k) = derivatives (4, ip, jp ,kp)   !

           ed_cellCubicNele (41,i,j,k) = derivatives (5, i , j , k )   ! d2/dxdz values
           ed_cellCubicNele (42,i,j,k) = derivatives (5, ip, j  ,k )   !
           ed_cellCubicNele (43,i,j,k) = derivatives (5, i , jp ,k )   !
           ed_cellCubicNele (44,i,j,k) = derivatives (5, i , j  ,kp)   !
           ed_cellCubicNele (45,i,j,k) = derivatives (5, ip, jp ,k )   !
           ed_cellCubicNele (46,i,j,k) = derivatives (5, ip, j  ,kp)   !
           ed_cellCubicNele (47,i,j,k) = derivatives (5, i , jp ,kp)   !
           ed_cellCubicNele (48,i,j,k) = derivatives (5, ip, jp ,kp)   !

           ed_cellCubicNele (49,i,j,k) = derivatives (6, i , j , k )   ! d2/dydz values
           ed_cellCubicNele (50,i,j,k) = derivatives (6, ip, j  ,k )   !
           ed_cellCubicNele (51,i,j,k) = derivatives (6, i , jp ,k )   !
           ed_cellCubicNele (52,i,j,k) = derivatives (6, i , j  ,kp)   !
           ed_cellCubicNele (53,i,j,k) = derivatives (6, ip, jp ,k )   !
           ed_cellCubicNele (54,i,j,k) = derivatives (6, ip, j  ,kp)   !
           ed_cellCubicNele (55,i,j,k) = derivatives (6, i , jp ,kp)   !
           ed_cellCubicNele (56,i,j,k) = derivatives (6, ip, jp ,kp)   !

           ed_cellCubicNele (57,i,j,k) = derivatives (7, i , j , k )   ! d3/dxdydz values
           ed_cellCubicNele (58,i,j,k) = derivatives (7, ip, j  ,k )   !
           ed_cellCubicNele (59,i,j,k) = derivatives (7, i , jp ,k )   !
           ed_cellCubicNele (60,i,j,k) = derivatives (7, i , j  ,kp)   !
           ed_cellCubicNele (61,i,j,k) = derivatives (7, ip, jp ,k )   !
           ed_cellCubicNele (62,i,j,k) = derivatives (7, ip, j  ,kp)   !
           ed_cellCubicNele (63,i,j,k) = derivatives (7, i , jp ,kp)   !
           ed_cellCubicNele (64,i,j,k) = derivatives (7, ip, jp ,kp)   !

        enddo
     enddo
  enddo
!
!
!     ...Calculate all Nele tricubic expansion coefficients in one shot.
!
!
  call  ut_triCubicCoefficients (numberOfCells, ed_cellCubicNele)
!
!
!     ...Repeat the same procedure for the Tele values.
!
!
  vertices = 0.0

  do k = kminData,kmaxData
     kp = k + 1
     do j = jminData,jmaxData
        jp = j + 1
        do i = iminData,imaxData
           ip = i + 1

           centerTele = blockData (TELE_VAR,i,j,k)

           vertices (i , j , k ) = vertices (i , j , k ) + centerTele
           vertices (ip, j , k ) = vertices (ip, j , k ) + centerTele
           vertices (i , jp, k ) = vertices (i , jp, k ) + centerTele
           vertices (i , j , kp) = vertices (i , j , kp) + centerTele
           vertices (ip, jp, k ) = vertices (ip, jp, k ) + centerTele
           vertices (ip, j , kp) = vertices (ip, j , kp) + centerTele
           vertices (i , jp, kp) = vertices (i , jp, kp) + centerTele
           vertices (ip, jp, kp) = vertices (ip, jp, kp) + centerTele

        enddo
     enddo
  enddo

  vertices = eighth * vertices

  do k = kminDerv, kmaxDerv
     km = k - 1
     kp = k + 1
     do j = jminDerv, jmaxDerv
        jm = j - 1
        jp = j + 1
        do i = iminDerv, imaxDerv
           im = i - 1
           ip = i + 1

           v000 = vertices (i , j , k )
           vm00 = vertices (im, j , k )
           vp00 = vertices (ip, j , k )
           v0m0 = vertices (i , jm, k )
           vmm0 = vertices (im, jm, k )
           vpm0 = vertices (ip, jm, k )
           v0p0 = vertices (i , jp, k )
           vmp0 = vertices (im, jp, k )
           vpp0 = vertices (ip, jp, k )
           v00m = vertices (i , j , km)
           vm0m = vertices (im, j , km)
           vp0m = vertices (ip, j , km)
           v0mm = vertices (i , jm, km)
           vmmm = vertices (im, jm, km)
           vpmm = vertices (ip, jm, km)
           v0pm = vertices (i , jp, km)
           vmpm = vertices (im, jp, km)
           vppm = vertices (ip, jp, km)
           v00p = vertices (i , j , kp)
           vm0p = vertices (im, j , kp)
           vp0p = vertices (ip, j , kp)
           v0mp = vertices (i , jm, kp)
           vmmp = vertices (im, jm, kp)
           vpmp = vertices (ip, jm, kp)
           v0pp = vertices (i , jp, kp)
           vmpp = vertices (im, jp, kp)
           vppp = vertices (ip, jp, kp)

           dx = ed_mc (v000 - vm00 , vp00 - v000)
           dy = ed_mc (v000 - v0m0 , v0p0 - v000)
           dz = ed_mc (v000 - v00m , v00p - v000)

           dxdy = (  ed_mc (v0p0 - vmp0 , vpp0 - v0p0) &
                   - ed_mc (v0m0 - vmm0 , vpm0 - v0m0) &
                   + ed_mc (vp00 - vpm0 , vpp0 - vp00) &
                   - ed_mc (vm00 - vmm0 , vmp0 - vm00) ) * fourth

           dxdz = (  ed_mc (v00p - vm0p , vp0p - v00p) &
                   - ed_mc (v00m - vm0m , vp0m - v00m) &
                   + ed_mc (vp00 - vp0m , vp0p - vp00) &
                   - ed_mc (vm00 - vm0m , vm0p - vm00) ) * fourth

           dydz = (  ed_mc (v00p - v0mp , v0pp - v00p) &
                   - ed_mc (v00m - v0mm , v0pm - v00m) &
                   + ed_mc (v0p0 - v0pm , v0pp - v0p0) &
                   - ed_mc (v0m0 - v0mm , v0mp - v0m0) ) * fourth

           dxdydz = (  ed_mc (v0pp - vmpp , vppp - v0pp) &
                     - ed_mc (v0pm - vmpm , vppm - v0pm) &
                     - ed_mc (v0mp - vmmp , vpmp - v0mp) &
                     + ed_mc (v0mm - vmmm , vpmm - v0mm) &
                     + ed_mc (vp0p - vpmp , vppp - vp0p) &
                     - ed_mc (vp0m - vpmm , vppm - vp0m) &
                     - ed_mc (vm0p - vmmp , vmpp - vm0p) &
                     + ed_mc (vm0m - vmmm , vmpm - vm0m) &
                     + ed_mc (vpp0 - vppm , vppp - vpp0) &
                     - ed_mc (vpm0 - vpmm , vpmp - vpm0) &
                     - ed_mc (vmp0 - vmpm , vmpp - vmp0) &
                     + ed_mc (vmm0 - vmmm , vmmp - vmm0) ) * twelfth

           derivatives (1,i,j,k) = dx
           derivatives (2,i,j,k) = dy
           derivatives (3,i,j,k) = dz
           derivatives (4,i,j,k) = dxdy
           derivatives (5,i,j,k) = dxdz
           derivatives (6,i,j,k) = dydz
           derivatives (7,i,j,k) = dxdydz

        enddo
     enddo
  enddo

  do k = kminBlock, kmaxBlock
     kp = k + 1
     do j = jminBlock, jmaxBlock
        jp = j + 1
        do i = iminBlock, imaxBlock
           ip = i + 1

           ed_cellCubicTele ( 1,i,j,k) =    vertices (   i , j , k )   ! vertex values
           ed_cellCubicTele ( 2,i,j,k) =    vertices (   ip, j  ,k )   !
           ed_cellCubicTele ( 3,i,j,k) =    vertices (   i , jp ,k )   !
           ed_cellCubicTele ( 4,i,j,k) =    vertices (   i , j  ,kp)   !
           ed_cellCubicTele ( 5,i,j,k) =    vertices (   ip, jp ,k )   !
           ed_cellCubicTele ( 6,i,j,k) =    vertices (   ip, j  ,kp)   !
           ed_cellCubicTele ( 7,i,j,k) =    vertices (   i , jp ,kp)   !
           ed_cellCubicTele ( 8,i,j,k) =    vertices (   ip, jp ,kp)   !

           ed_cellCubicTele ( 9,i,j,k) = derivatives (1, i , j , k )   ! d/dx values
           ed_cellCubicTele (10,i,j,k) = derivatives (1, ip, j  ,k )   !
           ed_cellCubicTele (11,i,j,k) = derivatives (1, i , jp ,k )   !
           ed_cellCubicTele (12,i,j,k) = derivatives (1, i , j  ,kp)   !
           ed_cellCubicTele (13,i,j,k) = derivatives (1, ip, jp ,k )   !
           ed_cellCubicTele (14,i,j,k) = derivatives (1, ip, j  ,kp)   !
           ed_cellCubicTele (15,i,j,k) = derivatives (1, i , jp ,kp)   !
           ed_cellCubicTele (16,i,j,k) = derivatives (1, ip, jp ,kp)   !

           ed_cellCubicTele (17,i,j,k) = derivatives (2, i , j , k )   ! d/dy values
           ed_cellCubicTele (18,i,j,k) = derivatives (2, ip, j  ,k )   !
           ed_cellCubicTele (19,i,j,k) = derivatives (2, i , jp ,k )   !
           ed_cellCubicTele (20,i,j,k) = derivatives (2, i , j  ,kp)   !
           ed_cellCubicTele (21,i,j,k) = derivatives (2, ip, jp ,k )   !
           ed_cellCubicTele (22,i,j,k) = derivatives (2, ip, j  ,kp)   !
           ed_cellCubicTele (23,i,j,k) = derivatives (2, i , jp ,kp)   !
           ed_cellCubicTele (24,i,j,k) = derivatives (2, ip, jp ,kp)   !

           ed_cellCubicTele (25,i,j,k) = derivatives (3, i , j , k )   ! d/dz values
           ed_cellCubicTele (26,i,j,k) = derivatives (3, ip, j  ,k )   !
           ed_cellCubicTele (27,i,j,k) = derivatives (3, i , jp ,k )   !
           ed_cellCubicTele (28,i,j,k) = derivatives (3, i , j  ,kp)   !
           ed_cellCubicTele (29,i,j,k) = derivatives (3, ip, jp ,k )   !
           ed_cellCubicTele (30,i,j,k) = derivatives (3, ip, j  ,kp)   !
           ed_cellCubicTele (31,i,j,k) = derivatives (3, i , jp ,kp)   !
           ed_cellCubicTele (32,i,j,k) = derivatives (3, ip, jp ,kp)   !

           ed_cellCubicTele (33,i,j,k) = derivatives (4, i , j , k )   ! d2/dxdy values
           ed_cellCubicTele (34,i,j,k) = derivatives (4, ip, j  ,k )   !
           ed_cellCubicTele (35,i,j,k) = derivatives (4, i , jp ,k )   !
           ed_cellCubicTele (36,i,j,k) = derivatives (4, i , j  ,kp)   !
           ed_cellCubicTele (37,i,j,k) = derivatives (4, ip, jp ,k )   !
           ed_cellCubicTele (38,i,j,k) = derivatives (4, ip, j  ,kp)   !
           ed_cellCubicTele (39,i,j,k) = derivatives (4, i , jp ,kp)   !
           ed_cellCubicTele (40,i,j,k) = derivatives (4, ip, jp ,kp)   !

           ed_cellCubicTele (41,i,j,k) = derivatives (5, i , j , k )   ! d2/dxdz values
           ed_cellCubicTele (42,i,j,k) = derivatives (5, ip, j  ,k )   !
           ed_cellCubicTele (43,i,j,k) = derivatives (5, i , jp ,k )   !
           ed_cellCubicTele (44,i,j,k) = derivatives (5, i , j  ,kp)   !
           ed_cellCubicTele (45,i,j,k) = derivatives (5, ip, jp ,k )   !
           ed_cellCubicTele (46,i,j,k) = derivatives (5, ip, j  ,kp)   !
           ed_cellCubicTele (47,i,j,k) = derivatives (5, i , jp ,kp)   !
           ed_cellCubicTele (48,i,j,k) = derivatives (5, ip, jp ,kp)   !

           ed_cellCubicTele (49,i,j,k) = derivatives (6, i , j , k )   ! d2/dydz values
           ed_cellCubicTele (50,i,j,k) = derivatives (6, ip, j  ,k )   !
           ed_cellCubicTele (51,i,j,k) = derivatives (6, i , jp ,k )   !
           ed_cellCubicTele (52,i,j,k) = derivatives (6, i , j  ,kp)   !
           ed_cellCubicTele (53,i,j,k) = derivatives (6, ip, jp ,k )   !
           ed_cellCubicTele (54,i,j,k) = derivatives (6, ip, j  ,kp)   !
           ed_cellCubicTele (55,i,j,k) = derivatives (6, i , jp ,kp)   !
           ed_cellCubicTele (56,i,j,k) = derivatives (6, ip, jp ,kp)   !

           ed_cellCubicTele (57,i,j,k) = derivatives (7, i , j , k )   ! d3/dxdydz values
           ed_cellCubicTele (58,i,j,k) = derivatives (7, ip, j  ,k )   !
           ed_cellCubicTele (59,i,j,k) = derivatives (7, i , jp ,k )   !
           ed_cellCubicTele (60,i,j,k) = derivatives (7, i , j  ,kp)   !
           ed_cellCubicTele (61,i,j,k) = derivatives (7, ip, jp ,k )   !
           ed_cellCubicTele (62,i,j,k) = derivatives (7, ip, j  ,kp)   !
           ed_cellCubicTele (63,i,j,k) = derivatives (7, i , jp ,kp)   !
           ed_cellCubicTele (64,i,j,k) = derivatives (7, ip, jp ,kp)   !

        enddo
     enddo
  enddo

  call  ut_triCubicCoefficients (numberOfCells, ed_cellCubicTele)
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
end subroutine ed_blockData3DRec
