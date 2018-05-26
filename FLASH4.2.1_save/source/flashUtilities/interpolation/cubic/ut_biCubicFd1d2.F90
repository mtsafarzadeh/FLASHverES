!!****if* source/flashUtilities/interpolation/cubic/ut_biCubicFd1d2
!!
!! NAME
!!
!!  ut_biCubicFd1d2
!!
!! SYNOPSIS
!!
!!  ut_biCubicFd1d2 (real, intent (in) :: a (1:16),
!!                   real, intent (in) :: x,
!!                   real, intent (in) :: y)
!!
!! DESCRIPTION
!!
!!  Calculates the function and rescaled pure 1st and 2nd derivative value for a pair
!!  [x,y] of rescaled [0,1] coordinates and the 16 bicubic expansion coefficients.
!!  The bicubic expansion reads, for one square, in terms of rescaled [0,1] x,y coordinates:
!!
!!                                   3   3            i j
!!                        F (x,y) = sum sum  a (i,j) x y
!!                                  i=0 j=0
!!
!!  The order of the supplied expansion coefficients a (i,j) must be such, that the
!!  j index has the highest ranking, followed by the i index. The overall location index
!!  of the a (i,j) inside the 16-dimensional vector is given by the following formula:
!!
!!                location index of (i,j)  =  1 + i + 4j
!!
!!  In the present function routine, the pair of rescaled [0,1] coordinates is first
!!  transformed into a 16-dimensional monomial vector, in which each mononial is placed
!!  according to the above stated location index. Then the function and derivative
!!  values are calculated. The rescaled derivatives are given by the formulae:
!!
!!
!!                         3   3                i-1 j
!!                 d/dx = sum sum  i * a (i,j) x   y
!!                        i=1 j=0
!!
!!                         3   3                i j-1
!!                 d/dy = sum sum  j * a (i,j) x y
!!                        i=0 j=1
!!
!!                         3   3                          i-2 j
!!               d2/dx2 = sum sum  i * (i - 1) * a (i,j) x   y
!!                        i=2 j=0
!!
!!                         3   3                          i j-2
!!               d2/dy2 = sum sum  j * (j - 1) * a (i,j) x y
!!                        i=0 j=2
!!
!!
!!  and are therefore also sums of appropriate expansion coefficients x monomial products.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th bicubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!  y     : rescaled [0,1] y coordinate
!!
!! NOTES
!!
!!  1) The function is defined as a real array of size 5:
!!
!!           ut_biCubicFd1d2 (1) = the function value
!!           ut_biCubicFd1d2 (2) = the rescaled d/dx value
!!           ut_biCubicFd1d2 (3) = the rescaled d/dy value
!!           ut_biCubicFd1d2 (4) = the rescaled d2/dx2 value
!!           ut_biCubicFd1d2 (5) = the rescaled d2/dy2 value
!!
!!  2) The code checks, if the supplied pair [x,y] is rescaled.
!!
!!***

function ut_biCubicFd1d2 (a,x,y)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:16)
  real, intent (in) :: x,y

  logical :: notRescaled

  real    :: p,q,r

  real    :: ut_biCubicFd1d2 (1:5)     ! declares the function as an array
  real    :: xy              (1:16)    ! will hold the monomials
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [ut_biCubicFd1d2]'    )
      call Logfile_stamp     (y, ' = rescaled y coordinate [ut_biCubicFd1d2]'    )
      call Driver_abortFlash ('[ut_biCubicFd1d2] ERROR: [x,y] pair not rescaled!')
  end if
!
!
!     ...Generate the monomial vector.
!
!
  xz ( 1)    = 1.0
  xz ( 2)    = x
  xz ( 3)    = x * xz (2)
  xz ( 4)    = x * xz (3)
  xz ( 5: 8) = z * xz (1:4)
  xz ( 9:12) = z * xz (5:8)
  xz (13:16) = z * xz (9:12)
!
!
!     ...Calculate the function value.
!
!
  ut_biCubicFd1d2 (1) = sum (xy (1:16) * a (1:16))
!
!
!     ...Calculate the d/dx and d/dy values. For the derivative evaluations, only
!        a subset of the monomials and the bicubic expansion coefficients are needed.
!        Below is a picture to show which sections of the two arrays are needed
!        for each individual type of derivative.
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
!                                 xy (i,j)
!
!
!         d/dx  -->  The picture shows a particular j-index section with i-index range
!                    from 0 to 3. Each box contains only 1 element.
!
!         d/dy  -->  The picture shows the entire vectors. Each box corresponds to
!                    a particular j-index and contains 4 elements with i-index
!                    range from 0 to 3.
!
!
  p = sum ( xy (1:13:4) * a (2:14:4) )
  q = sum ( xy (2:14:4) * a (3:15:4) )
  r = sum ( xy (3:15:4) * a (4:16:4) )

  ut_biCubicFd1d2 (2) = p + q + q + r + r + r    ! this is d/dx

  p = sum ( xy (1: 4) * a ( 5: 8) )
  q = sum ( xy (5: 8) * a ( 9:12) )
  r = sum ( xy (9:12) * a (13:16) )

  ut_biCubicFd1d2 (3) = p + q + q + r + r + r    ! this is d/dy
!
!
!     ...Calculate the d2/dx2 and d2/dy2 values. For the derivative evaluations, only
!        a subset of the monomials and the bicubic expansion coefficients are needed.
!        Below is a picture to show which sections of the two arrays are needed
!        for each individual type of derivative.
!
!
!                       a (i,j)
!
!
!                          --                  <--- shown is a section of 4 consecutive
!                         |XX|      --              boxes each containing an equal amount
!                         |XX|     |  |             of elements. The slope lines indicate
!                         |--|    /|  |             which boxes are multiplied together
!                         |  |  2/ |--|             (element by element). The numerical
!                         |  |  /  |  |             value on the slope lines indicate the
!                         |--| /  /|  |             multiplying factor. Boxes marked XX
!                         |  |/ 6/ |--|             are not needed.
!                         |  |  /  |XX|
!                         |--| /   |XX|
!                         |  |/    |--|
!                         |  |     |XX|
!                          --      |XX|
!                                   --
!
!                                xy (i,j)
!
!
!       d2/dx2  -->  The picture shows a particular j-index section with i-index range
!                    from 0 to 3. Each box contains only 1 element.
!
!       d2/dy2  -->  The picture shows the entire vectors. Each box corresponds to
!                    a particular j-index and contains 4 elements with i-index
!                    range from 0 to 3.
!
!
  p = sum ( xy (1:13:4) * a (3:15:4) )
  q = sum ( xy (2:14:4) * a (4:16:4) )

  ut_biCubicFd1d2 (4) = p + p + q + q + q + q + q + q    ! this is d2/dx2

  p = sum ( xy (1:4) * a ( 9:12) )
  q = sum ( xy (5:8) * a (13:16) )

  ut_biCubicFd1d2 (5) = p + p + q + q + q + q + q + q    ! this is d2/dz2
!
!
!     ...Ready!
!
!
  return
end function ut_biCubicFd1d2
