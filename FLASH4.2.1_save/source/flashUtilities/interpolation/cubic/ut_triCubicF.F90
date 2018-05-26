!!****if* source/flashUtilities/interpolation/cubic/ut_triCubicF
!!
!! NAME
!!
!!  ut_triCubicF
!!
!! SYNOPSIS
!!
!!  ut_triCubicF (real, intent (in) :: a (1:64),
!!                real, intent (in) :: x,
!!                real, intent (in) :: y,
!!                real, intent (in) :: z)
!!
!! DESCRIPTION
!!
!!  Calculates the function value for a triple [x,y,z] of rescaled [0,1] coordinates and
!!  the 64 tricubic expansion coefficients. The tricubic expansion reads, for one cube,
!!  in terms of rescaled [0,1] x,y,z coordinates:
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
!!  according to the above stated location index. Then the function is calculated.
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
!!  1) The code checks, if the supplied triple [x,y,z] is rescaled.
!!
!!***

real function ut_triCubicF (a,x,y,z)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:64)
  real, intent (in) :: x,y,z

  logical :: notRescaled

  real    :: xyz (1:64)    ! will hold the monomials
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. (z < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0) .or. (z > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [ut_triCubicF]'        )
      call Logfile_stamp     (y, ' = rescaled y coordinate [ut_triCubicF]'        )
      call Logfile_stamp     (z, ' = rescaled z coordinate [ut_triCubicF]'        )
      call Driver_abortFlash ('[ut_triCubicF] ERROR: [x,y,z] triple not rescaled!')
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
  ut_triCubicF = sum (xyz (1:64) * a (1:64))
!
!
!     ...Ready!
!
!
  return
end function ut_triCubicF
