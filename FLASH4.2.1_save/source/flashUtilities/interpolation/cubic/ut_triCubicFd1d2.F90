!!****if* source/flashUtilities/interpolation/cubic/ut_triCubicFd1d2
!!
!! NAME
!!
!!  ut_triCubicFd1d2
!!
!! SYNOPSIS
!!
!!  ut_triCubicFd1d2 (real, intent (in) :: a (1:64),
!!                    real, intent (in) :: x,
!!                    real, intent (in) :: y,
!!                    real, intent (in) :: z)
!!
!! DESCRIPTION
!!
!!  Calculates the function and rescaled pure 1st and 2nd derivative value for a triple
!!  [x,y,z] of rescaled [0,1] coordinates and the 64 tricubic expansion coefficients.
!!  The tricubic expansion reads, for one cube, in terms of rescaled [0,1] x,y,z coordinates:
!!
!!                                   3   3   3              i j k
!!                      F (x,y,z) = sum sum sum  a (i,j,k) x y z
!!                                  i=0 j=0 k=0
!!
!!  The order of the supplied expansion coefficients a (i,j,k) must be such, that the
!!  k index has the highest ranking, followed by the j index and the i index. The overall
!!  location index of the a (i,j,k) inside the 64-dimensional vector is given by the
!!  following formula:
!!
!!                location index of (i,j,k)  =  1 + i + 4j + 16k
!!
!!  In the present function routine, the triple of rescaled [0,1] coordinates is first
!!  transformed into a 64-dimensional monomial vector, in which each mononial is placed
!!  according to the above stated location index. Then the function and derivative
!!  values are calculated. The rescaled derivatives are given by the formulae:
!!
!!
!!                         3   3   3                  i-1 j k
!!                 d/dx = sum sum sum  i * a (i,j,k) x   y z
!!                        i=1 j=0 k=0
!!
!!                         3   3   3                  i j-1 k
!!                 d/dy = sum sum sum  j * a (i,j,k) x y   z
!!                        i=0 j=1 k=0
!!
!!                         3   3   3                  i j k-1
!!                 d/dz = sum sum sum  k * a (i,j,k) x y z
!!                        i=0 j=0 k=1
!!
!!                         3   3   3                            i-2 j k
!!               d2/dx2 = sum sum sum  i * (i - 1) * a (i,j,k) x   y z
!!                        i=2 j=0 k=0
!!
!!                         3   3   3                            i j-2 k
!!               d2/dy2 = sum sum sum  j * (j - 1) * a (i,j,k) x y   z
!!                        i=0 j=2 k=0
!!
!!                         3   3   3                            i j k-2
!!               d2/dz2 = sum sum sum  k * (k - 1) * a (i,j,k) x y z
!!                        i=0 j=0 k=2
!!
!!
!!  and are therefore also sums of appropriate expansion coefficients x monomial products.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th tricubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!  y     : rescaled [0,1] y coordinate
!!  z     : rescaled [0,1] z coordinate
!!
!! NOTES
!!
!!  1) The function is defined as a real array of size 7:
!!
!!           ut_triCubicFd1d2 (1) = the function value
!!           ut_triCubicFd1d2 (2) = the rescaled d/dx value
!!           ut_triCubicFd1d2 (3) = the rescaled d/dy value
!!           ut_triCubicFd1d2 (4) = the rescaled d/dz value
!!           ut_triCubicFd1d2 (5) = the rescaled d2/dx2 value
!!           ut_triCubicFd1d2 (6) = the rescaled d2/dy2 value
!!           ut_triCubicFd1d2 (7) = the rescaled d2/dz2 value
!!
!!  2) The code checks, if the supplied triple [x,y,z] is rescaled.
!!
!!***

function ut_triCubicFd1d2 (a,x,y,z)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:64)
  real, intent (in) :: x,y,z

  logical :: notRescaled

  real    :: p,q,r

  real    :: ut_triCubicFd1d2 (1:7)     ! declares the function as an array
  real    :: xyz              (1:64)    ! will hold the monomials
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. (z < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0) .or. (z > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [ut_triCubicFd1d2]'        )
      call Logfile_stamp     (y, ' = rescaled y coordinate [ut_triCubicFd1d2]'        )
      call Logfile_stamp     (z, ' = rescaled z coordinate [ut_triCubicFd1d2]'        )
      call Driver_abortFlash ('[ut_triCubicFd1d2] ERROR: [x,y,z] triple not rescaled!')
  end if
!
!
!     ...Generate the monomial vector.
!
!
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
!
!
!     ...Calculate the function value.
!
!
  ut_triCubicFd1d2 (1) = sum (xyz (1:64) * a (1:64))
!
!
!     ...Calculate the d/dx, d/dy and d/dz values. For the derivative evaluations,
!        only a subset of the monomials and the tricubic expansion coefficients are
!        needed. Below is a picture to show which sections of the two arrays are needed
!        for each individual type of derivative.
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
!         d/dx  -->  The picture shows a particular jk-index pair section with
!                    i-index range from 0 to 3. Each box contains only 1 element.
!
!         d/dy  -->  The picture shows a particular k-index section with i- and
!                    j-index ranges from 0 to 3. Each box corresponds to a particular
!                    j-index and contains 4 element with i-index range from 0 to 3.
!
!         d/dz  -->  The picture shows the entire vectors. Each box corresponds to
!                    a particular k-index and contains 16 elements with i- and j-index
!                    ranges from 0 to 3.
!
!
  p = sum ( xyz (1:61:4) * a (2:62:4) )
  q = sum ( xyz (2:62:4) * a (3:63:4) )
  r = sum ( xyz (3:63:4) * a (4:64:4) )

  ut_triCubicFd1d2 (2) = p + q + q + r + r + r    ! this is d/dx

  p = sum (  xyz ( 1: 4) * a ( 5: 8) &
           + xyz (17:20) * a (21:24) &
           + xyz (33:36) * a (37:40) &
           + xyz (49:52) * a (53:56) )

  q = sum (  xyz ( 5: 8) * a ( 9:12) &
           + xyz (21:24) * a (25:28) &
           + xyz (37:40) * a (41:44) &
           + xyz (53:56) * a (57:60) )

  r = sum (  xyz ( 9:12) * a (13:16) &
           + xyz (25:28) * a (29:32) &
           + xyz (41:44) * a (45:48) &
           + xyz (57:60) * a (61:64) )

  ut_triCubicFd1d2 (3) = p + q + q + r + r + r    ! this is d/dy

  p = sum (  xyz ( 1:16) * a (17:32) )
  q = sum (  xyz (17:32) * a (33:48) )
  r = sum (  xyz (33:48) * a (49:64) )

  ut_triCubicFd1d2 (4) = p + q + q + r + r + r    ! this is d/dz
!
!
!     ...Calculate the d2/dx2, d2/dy2 and d2/dz2 values. For the derivative evaluations,
!        only a subset of the monomials and the tricubic expansion coefficients are
!        needed. Below is a picture to show which sections of the two arrays are needed
!        for each individual type of derivative.
!
!
!                       a (i,j,k)
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
!                                xyz (i,j,k)
!
!
!       d2/dx2  -->  The picture shows a particular jk-index pair section with
!                    i-index range from 0 to 3. Each box contains only 1 element.
!
!       d2/dy2  -->  The picture shows a particular k-index section with i- and
!                    j-index ranges from 0 to 3. Each box corresponds to a particular
!                    j-index and contains 4 element with i-index range from 0 to 3.
!
!       d2/dz2  -->  The picture shows the entire vectors. Each box corresponds to
!                    a particular k-index and contains 16 elements with i- and j-index
!                    ranges from 0 to 3.
!
!
  p = sum ( xyz (1:61:4) * a (3:63:4) )
  q = sum ( xyz (2:62:4) * a (4:64:4) )

  ut_triCubicFd1d2 (5) = p + p + q + q + q + q + q + q    ! this is d2/dx2

  p = sum (  xyz ( 1: 4) * a ( 9:12) &
           + xyz (17:20) * a (25:28) &
           + xyz (33:36) * a (41:44) &
           + xyz (49:52) * a (57:60) )

  q = sum (  xyz ( 5: 8) * a (13:16) &
           + xyz (21:24) * a (29:32) &
           + xyz (37:40) * a (45:48) &
           + xyz (53:56) * a (61:64) )

  ut_triCubicFd1d2 (6) = p + p + q + q + q + q + q + q    ! this is d2/dy2

  p = sum (  xyz ( 1:16) * a (33:48) )
  q = sum (  xyz (17:32) * a (49:64) )

  ut_triCubicFd1d2 (7) = p + p + q + q + q + q + q + q    ! this is d2/dz2
!
!
!     ...Ready!
!
!
  return
end function ut_triCubicFd1d2
