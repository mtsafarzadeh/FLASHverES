!!****if* source/flashUtilities/general/ut_cubicRoots
!!
!! NAME
!!
!!  ut_cubicRoots
!!
!! SYNOPSIS
!!
!!  call ut_cubicRoots (real,    intent (in)  :: c2,
!!                      real,    intent (in)  :: c1,
!!                      real,    intent (in)  :: c0,
!!                      integer, intent (out) :: nReal,
!!                      real,    intent (out) :: root (1:3,1:2))
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the cubic polynomial:
!!
!!                 x^3 + c2 * x^2 + c1 * x + c0
!!
!!  The first real root (which always exists) is obtained using an optimized
!!  Newton-Raphson scheme, except for near triply degenerate roots, where
!!  the analytical formula is being used. The other remaining roots are obtained
!!  through deflation into a quadratic.
!!
!!  The cubic root solver can handle any size of cubic coefficients and there is
!!  no danger of overflow due to proper rescaling of the cubic polynomial.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) Since there can be only one complex conjugate pair root, no order
!!           is necessary.
!!
!!        3) All real roots preceede the complex ones.
!!
!! ARGUMENTS
!!
!!  c2         : coefficient of x^2 term
!!  c1         : coefficient of x term
!!  c0         : independent coefficient
!!  nReal      : number of different real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!
!! NOTES
!!
!!***

subroutine ut_cubicRoots (c2, c1, c0,         &
                                       nReal, &
                                       root   )
  implicit none

  real,    intent (in)  :: c2, c1, c0
  integer, intent (out) :: nReal
  real,    intent (out) :: root (1:3,1:2)

  character (len = 12) :: cubicType

  logical :: converged
  logical :: rescaleQuadratic

  integer :: n
  integer :: nReal2

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real    :: a0, a1, a2
  real    :: k, s, t, u, x, y, z
  real    :: tolerance
  real    :: xPrev, xShift

  real    :: root2 (1:2,1:2)

  real, parameter :: accuracy    = 1.d-14

  real, parameter :: zero        = 0.d0
  real, parameter :: one100th    = 1.d-2
  real, parameter :: one27th     = 1.d0 / 27.d0
  real, parameter :: two27th     = 2.d0 / 27.d0
  real, parameter :: eleven27th  = 11.d0 / 27.d0
  real, parameter :: sixteen27th = 16.d0 / 27.d0
  real, parameter :: one9th      = 1.d0 / 9.d0
  real, parameter :: third       = 1.d0 / 3.d0
  real, parameter :: half        = 5.d-1
  real, parameter :: one         = 1.d0
  real, parameter :: two         = 2.d0

  real, parameter :: pi          = 3.1415926535897932384d0

  real, parameter :: p1          = 1.09574d0           !
  real, parameter :: q1          = 3.23900d-1          ! Newton-Raphson coeffs for class 1 and 2
  real, parameter :: r1          = 3.23900d-1          !
  real, parameter :: s1          = 9.57439d-2          !

  real, parameter :: p3          = 1.14413d0           !
  real, parameter :: q3          = 2.75509d-1          ! Newton-Raphson coeffs for class 3
  real, parameter :: r3          = 4.45578d-1          !
  real, parameter :: s3          = 2.59342d-2          !

  real, parameter :: q4          = 7.71845d-1          ! Newton-Raphson coeffs for class 4
  real, parameter :: s4          = 2.28155d-1          !

  real, parameter :: p51         = 8.78558d-1          !
  real, parameter :: p52         = 1.92823d-1          !
  real, parameter :: p53         = 1.19748d0           !
  real, parameter :: p54         = 3.45219d-1          !
  real, parameter :: q51         = 5.71888d-1          !
  real, parameter :: q52         = 5.66324d-1          !
  real, parameter :: q53         = 2.83772d-1          ! Newton-Raphson coeffs for class 5 and 6
  real, parameter :: q54         = 4.01231d-1          !
  real, parameter :: r51         = 7.11154d-1          !
  real, parameter :: r52         = 5.05734d-1          !
  real, parameter :: r53         = 8.37476d-1          !
  real, parameter :: r54         = 2.07216d-1          !
  real, parameter :: s51         = 3.22313d-1          !
  real, parameter :: s52         = 2.64881d-1          !
  real, parameter :: s53         = 3.56228d-1          !
  real, parameter :: s54         = 4.45532d-3          !
!
!
!     ...Handle special cases.
!
!            1) all terms zero
!            2) only quadratic term is nonzero -> linear equation.
!            3) only independent term is zero -> quadratic equation.
!
!
  if (c0 == zero .and. c1 == zero .and. c2 == zero) then

      cubicType = 'allzero'

  else if (c0 == zero .and. c1 == zero) then

      k  = one
      a2 = c2

      cubicType = 'linear'

  else if (c0 == zero) then

      k  = one
      a2 = c2
      a1 = c1

      cubicType = 'quadratic'
      rescaleQuadratic = .true.

  else
!
!
!     ...The general case. Rescale cubic polynomial, such that largest absolute coefficient
!        is (exactly!) equal to 1. Honor the presence of a special cubic case that might have
!        been obtained (due to underflow in the coefficients).
!
!
      x = abs (c2)
      y = sqrt (abs (c1))
      z = abs (c0) ** third

      u = max (x,y,z)

      if (u == x) then

          k  = one / x
          t  = k * k
          a2 = sign (one , c2)
          a1 = (c1 * k) * k
          a0 = ((c0 * k) * k) * k

      else if (u == y) then

          k  = one / y
          a2 = c2 * k
          a1 = sign (one , c1)
          a0 = ((c0 * k) * k) * k

      else

          k  = one / z
          a2 = c2 * k
          a1 = (c1 * k) * k
          a0 = sign (one , c0)

      end if

      k = one / k

      if (a0 == zero .and. a1 == zero .and. a2 == zero) then
          cubicType = 'allzero'
      else if (a0 == zero .and. a1 == zero) then
          cubicType = 'linear'
      else if (a0 == zero) then
          cubicType = 'quadratic'
          rescaleQuadratic = .false.
      else
          cubicType = 'general'
      end if

  end if
!
!
!     ...Select the case.
!
!        1) Only zero roots.
!
!
  select case (cubicType)

    case ('allzero')

      nReal = 3

      root (:,Re) = zero
      root (:,Im) = zero
!
!
!     ...2) The linear equation case -> additional 2 zeros.
!
!
    case ('linear')

      x = - a2 * k

      nReal = 3

      root (1,Re) = max (zero, x)
      root (2,Re) = zero
      root (3,Re) = min (zero, x)
      root (:,Im) = zero
!
!
!     ...3) The quadratic equation case -> additional 1 zero.
!
!
    case ('quadratic')

      call ut_quadraticRoots (a2, a1,                   &
                              rescaleQuadratic,         &
                                                nReal2, &
                                                root2   )

      if (nReal2 == 2) then

          x = root2 (1,1) * k         ! real roots of quadratic are ordered x >= y
          y = root2 (2,1) * k

          nReal = 3

          root (1,Re) = max (x, zero)
          root (2,Re) = max (y, min (x, zero))
          root (3,Re) = min (y, zero)
          root (:,Im) = zero

      else

          nReal = 1

          root (1,Re) = zero
          root (2,Re) = root2 (1,Re) * k
          root (3,Re) = root2 (2,Re) * k
          root (1,Im) = zero
          root (2,Im) = root2 (1,Im) * k
          root (3,Im) = root2 (2,Im) * k

      end if
!
!
!     ...3) The general cubic case. Set the best Newton-Raphson root estimates for the cubic.
!           The easiest and most robust conditions are checked first. The most complicated
!           ones are last and only done when absolutely necessary.
!
!
    case ('general')

      if (a0 == one) then

          x = - p1 + q1 * a1 - a2 * (r1 - s1 * a1)

          s = a2
          t = a1
          u = a0
          xShift = zero

      else if (a0 == - one) then

          x = p1 - q1 * a1 - a2 * (r1 - s1 * a1)

          s = a2
          t = a1
          u = a0
          xShift = zero

      else if (a1 == one) then

          if (a0 > zero) then
              x = a0 * (- q4 - s4 * a2)
          else
              x = a0 * (- q4 + s4 * a2)
          end if

          s = a2
          t = a1
          u = a0
          xShift = zero

      else if (a1 == - one) then

          y = - two27th
          y = y * a2
          y = y * a2 - third
          y = y * a2

          if (a0 < y) then
              x = + p3 - q3 * a0 - a2 * (r3 + s3 * a0)       ! + guess
          else
              x = - p3 - q3 * a0 - a2 * (r3 - s3 * a0)       ! - guess
          end if

          s = a2
          t = a1
          u = a0
          xShift = zero

      else if (a2 == one) then

          s = a1 - third
          t = a0 - one27th

          if (abs (s) < accuracy .and. abs (t) < accuracy) then     ! triple -1/3 root

              x = - third * k

              nReal = 3

              root (1:3,Re) = x
              root (1:3,Im) = zero

              return

          else

              y = third * a1 - two27th

              if (a1 <= third) then
                  if (a0 > y) then
                      x = - p51 - q51 * a0 + a1 * (r51 - s51 * a0)   ! - guess
                  else
                      x = + p52 - q52 * a0 - a1 * (r52 + s52 * a0)   ! + guess
                  end if
              else
                  if (a0 > y) then
                      x = - p53 - q53 * a0 + a1 * (r53 - s53 * a0)   ! <-1/3 guess
                  else
                      x = + p54 - q54 * a0 - a1 * (r54 + s54 * a0)   ! >-1/3 guess
                  end if
              end if

              if (abs (s) < one100th .and. abs (t) < one100th) then  ! use shifted root
                  u = - third * s + t
                  t = s
                  s = zero
                  xShift = third
                  x = x + xShift
              else
                  s = a2
                  t = a1
                  u = a0
                  xShift = zero
              end if

          end if

      else if (a2 == - one) then

          s = a1 - third
          t = a0 + one27th

          if (abs (s) < accuracy .and. abs (t) < accuracy) then     ! triple 1/3 root

              x = third * k

              nReal = 3

              root (1:3,Re) = x
              root (1:3,Im) = zero

              return

          else

              y = two27th - third * a1

              if (a1 <= third) then
                  if (a0 < y) then
                      x = + p51 - q51 * a0 - a1 * (r51 + s51 * a0)   ! +1 guess
                  else
                      x = - p52 - q52 * a0 + a1 * (r52 - s52 * a0)   ! -1 guess
                  end if
              else
                  if (a0 < y) then
                      x = + p53 - q53 * a0 - a1 * (r53 + s53 * a0)   ! >1/3 guess
                  else
                      x = - p54 - q54 * a0 + a1 * (r54 - s54 * a0)   ! <1/3 guess
                  end if
              end if

              if (abs (s) < one100th .and. abs (t) < one100th) then  ! use shifted root
                  u = third * s + t
                  t = s
                  s = zero
                  xShift = - third
                  x = x + xShift
              else
                  s = a2
                  t = a1
                  u = a0
                  xShift = zero
              end if

          end if

      end if
!
!
!     ...Perform Newton-Raphson iterations.
!
!
      converged = .false.

      do while (.not.converged)

         tolerance = abs (x) * accuracy

         z = x + s
         y = x + z
         z = z * x + t
         y = y * x + z
         z = z * x + u
         xPrev = x
         x = x - z / y

         converged = abs (xPrev - x) <= tolerance

      end do

      x = x - xShift                  ! unshift root
      y = x * k                       ! save original cubic real root for deflation analysis
!
!
!     ...Forward / backward deflate rescaled cubic (if needed) to check for other real roots.
!        The deflation analysis is performed on the rescaled cubic. The actual deflation must
!        be performed on the original cubic, not the rescaled one. Otherwise deflation errors
!        will be enhanced when undoing the rescaling on the extra roots.
!
!
      t = x * x                       ! rescaled root squared -> x^2
      s = max (t * x , abs (t * a2))  ! maximum between |x^3| and |a2*x^2|
      t = abs (x * a1)                ! value of |a1*x|
      u = abs (a0)                    ! value of |a0|
      z = max (s,t,u)                 ! maximum of all terms

      if (z == s) then
          x = one / y
          a0 = - c0 * x               ! a0 -> backward deflation on unscaled cubic
          a1 = (a0 - c1) * x          ! a1 -> backward deflation on unscaled cubic
      else if (z == t) then
          a1 = c2 + y                 ! a1 ->  forward deflation on unscaled cubic
          a0 = - c0 / y               ! a0 -> backward deflation on unscaled cubic
      else
          a1 = c2 + y                 ! a1 ->  forward deflation on unscaled cubic
          a0 = c1 + a1 * y            ! a0 ->  forward deflation on unscaled cubic
      end if

      rescaleQuadratic = .true.

      call ut_quadraticRoots (a1, a0,                   &
                              rescaleQuadratic,         &
                                                nReal2, &
                                                root2   )

      if (nReal2 == 2) then

          x = root2 (1,Re)        ! real roots of quadratic are ordered x >= z
          z = root2 (2,Re)        ! use 'z', because 'y' is original cubic real root

          nReal = 3

          root (1,Re) = max (x, y)
          root (2,Re) = max (z, min (x, y))
          root (3,Re) = min (z, y)
          root (:,Im) = zero

      else

          nReal = 1

          root (1,Re) = y
          root (2,Re) = root2 (1,Re)
          root (3,Re) = root2 (2,Re)
          root (1,Im) = zero
          root (2,Im) = root2 (1,Im)
          root (3,Im) = root2 (2,Im)

      end if

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine ut_cubicRoots
