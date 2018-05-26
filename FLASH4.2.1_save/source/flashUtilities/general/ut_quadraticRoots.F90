!!****if* source/flashUtilities/general/ut_quadraticRoots
!!
!! NAME
!!
!!  ut_quadraticRoots
!!
!! SYNOPSIS
!!
!!  call ut_quadraticRoots (real,    intent (in)  :: q1,
!!                          real,    intent (in)  :: q0,
!!                          logical, intent (in)  :: rescale,
!!                          integer, intent (out) :: nReal,
!!                          real,    intent (out) :: root (1:2,1:2))
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the quadratic polynomial:
!!
!!                 x^2 + q1 * x + q0
!!
!!  An option for rescaling the coefficients is provided for cases of very large
!!  initial coefficients. The code does not decide automatically, if rescaling is
!!  needed. The explicit rescaling option must be provided by the user. If some
!!  applications are known to lead to extremely large quadratic coefficients, then
!!  it is always safe to enforce rescaling. The cost of the rescaling is an extra
!!  square root and two divisions.
!!
!!  The code also deals with the case in which the discriminant is very close
!!  to zero, but the individual terms are large and opposite. In this case the
!!  code automatically sets the discriminant equal to zero and thus induces a
!!  degenerate real root.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) Since there can be only one complex conjugate pair root, no order
!!           is necessary.
!!
!! ARGUMENTS
!!
!!  q1         : coefficient of x term
!!  q0         : independent coefficient
!!  rescale    : rescaling indicator
!!  nReal      : number of real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!
!! NOTES
!!
!!***

subroutine ut_quadraticRoots (q1, q0,         &
                              rescale,        &
                                       nReal, &
                                       root   )
  implicit none

  real,    intent (in)  :: q1, q0
  logical, intent (in)  :: rescale
  integer, intent (out) :: nReal
  real,    intent (out) :: root (1:2,1:2)

  real    :: a0, a1
  real    :: k, x, y, z

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real,    parameter :: accuracy = 1.d-12
  real,    parameter :: half     = 5.d-1
  real,    parameter :: one      = 1.d0
  real,    parameter :: zero     = 0.d0
!
!
!     ...Handle special cases.
!
!
  if (q0 == zero .and. q1 == zero) then

      nReal = 2

      root (:,Re) = zero
      root (:,Im) = zero

  else if (q0 == zero) then

      nReal = 2

      root (1,Re) = max (zero, - q1)
      root (2,Re) = min (zero, - q1)
      root (:,Im) = zero

  else if (q1 == zero) then

      x = sqrt (abs (q0))

      if (q0 < zero) then

          nReal = 2

          root (1,Re) = x
          root (2,Re) = - x
          root (:,Im) = zero

      else

          nReal = 0

          root (:,Re) = zero
          root (1,Im) = x
          root (2,Im) = - x

      end if

  else
!
!
!     ...The general case. Do rescaling (if requested).
!
!
      if (rescale) then

          x = abs (q1)
          y = sqrt (abs (q0))

          if (x > y) then
              k  = x
              z  = one / x
              a1 = sign (one , q1)
              a0 = (q0 * z) * z
          else
              k  = y
              a1 = q1 / y
              a0 = sign (one , q0)
          end if

      else
          a1 = q1
          a0 = q0
      end if
!
!
!     ...Determine the roots of the quadratic. Note, that either a1 or a0 might
!        have become equal to zero due to underflow. But both cannot be zero.
!
!
      if (a0 /= zero) then

          x = a1 * a1
          y = a0 + a0
          y = y + y
          z = max (abs (x) , abs (y))
          y = x - y
          z = abs (y / z)

          if (z < accuracy) then                    ! this catches the cases where the discriminant
              y = zero                              ! is considered equal to zero, but the a1^2 and the 4a0
          end if                                    ! are large and opposite in magnitude.

          if (y > zero) then

              y = sqrt (y)
              y = - half * (a1 + sign (one,a1) * y)
              z = a0 / y

              nReal = 2

              root (1,Re) = max (y,z)                ! 1st real root from x^2 + a1 * x + a0
              root (2,Re) = min (y,z)                ! 2nd real root from x^2 + a1 * x + a0
              root (:,Im) = zero

          else if (y == zero) then

              z = - half * a1

              nReal = 2

              root (1,Re) = z                        ! degenerate real root from x^2 + a1 * x + a0
              root (2,Re) = z                        ! degenerate real root from x^2 + a1 * x + a0
              root (:,Im) = zero

          else

              y = half * sqrt (abs (y))
              z = - half * a1

              nReal = 0

              root (1,Re) = z
              root (2,Re) = z
              root (1,Im) = y                        ! complex conjugate pair of roots
              root (2,Im) = - y                      ! from x^2 + a1 * x + a0

          end if

      else

          nReal = 2

          root (1,Re) = max (zero, - a1)             ! nonzero root from x^2 + a1 * x
          root (2,Re) = min (zero, - a1)             ! and zero root from x^2 + a1 * x
          root (:,Im) = zero

      end if
!
!
!     ...Rescale the roots (if needed).
!
!
      if (rescale) then
          root = root * k
      end if

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ut_quadraticRoots
