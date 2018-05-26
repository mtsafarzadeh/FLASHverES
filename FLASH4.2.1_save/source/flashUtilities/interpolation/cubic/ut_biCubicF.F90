!!****if* source/flashUtilities/interpolation/cubic/ut_biCubicF
!!
!! NAME
!!
!!  ut_biCubicF
!!
!! SYNOPSIS
!!
!!  ut_biCubicF (real, intent (in) :: a (1:16),
!!               real, intent (in) :: x,
!!               real, intent (in) :: y)
!!
!! DESCRIPTION
!!
!!  Calculates the function value for a pair [x,y] of rescaled [0,1] coordinates and
!!  the 16 bicubic expansion coefficients. The bicubic expansion reads, for one square,
!!  in terms of rescaled [0,1] x,y coordinates:
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
!!  according to the above stated location index. Then the function value is calculated.
!!
!! ARGUMENTS
!!
!!  a (i) : the i-th bicubic expansion coefficient
!!  x     : rescaled [0,1] x coordinate
!!  y     : rescaled [0,1] y coordinate
!!
!! NOTES
!!
!!  1) The code checks, if the supplied pair [x,y] is rescaled.
!!
!!***

real function ut_biCubicF (a,x,y)

  use Driver_interface,   ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp

  implicit none

  real, intent (in) :: a (1:16)
  real, intent (in) :: x,y

  logical :: notRescaled

  real    :: xy (1:16)    ! will hold the monomials
!
!
!     ...Check status of transmitted rescaled coordinates.
!
!
  notRescaled = (x < 0.0) .or. (y < 0.0) .or. &
                (x > 1.0) .or. (y > 1.0)

  if (notRescaled) then
      call Logfile_stamp     (x, ' = rescaled x coordinate [ut_biCubicF]'    )
      call Logfile_stamp     (y, ' = rescaled y coordinate [ut_biCubicF]'    )
      call Driver_abortFlash ('[ut_biCubicF] ERROR: [x,y] pair not rescaled!')
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
  ut_biCubicF = sum (xy (1:16) * a (1:16))
!
!
!     ...Ready!
!
!
  return
end function ut_biCubicF
