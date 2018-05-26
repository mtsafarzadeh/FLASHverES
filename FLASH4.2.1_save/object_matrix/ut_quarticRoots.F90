!!****if* source/flashUtilities/general/ut_quarticRoots
!!
!! NAME
!!
!!  ut_quarticRoots
!!
!! SYNOPSIS
!!
!!  call ut_quarticRoots (real,    intent (in)  :: q3,
!!                        real,    intent (in)  :: q2,
!!                        real,    intent (in)  :: q1,
!!                        real,    intent (in)  :: q0,
!!                        logical, intent (in)  :: printInfo,
!!                        integer, intent (out) :: nReal,
!!                        real,    intent (out) :: root (1:4,1:2))
!!
!! DESCRIPTION
!!
!!  Calculates all real + complex roots of the quartic polynomial:
!!
!!                 x^4 + q3 * x^3 + q2 * x^2 + q1 * x + q0
!!
!!  An option for printing a detailed info about the intermediate stages in solving
!!  the quartic is available. Since the code has not yet been extensively tested,
!!  this enables a detailed check in case something went wrong and the roots obtained
!!  are not proper.
!!
!!  The quartic root solver can handle any size of quartic coefficients and there is
!!  no danger of overflow, due to proper rescaling of the quartic polynomial.
!!
!!  The order of the roots is as follows:
!!
!!        1) For real roots, the order is according to their algebraic value
!!           on the number scale (largest positive first, largest negative last).
!!
!!        2) For complex conjugate pair roots, the order is according to the
!!           algebraic value of their real parts (largest positive first). If
!!           the real parts are equal, the order is according to the algebraic
!!           value of their imaginary parts (largest first).
!!
!!        3) All real roots preceede the complex ones.
!!
!! ARGUMENTS
!!
!!  q3         : coefficient of x^3 term
!!  q2         : coefficient of x^2 term
!!  q1         : coefficient of x term
!!  q0         : independent coefficient
!!  printInfo  : if true, detailed info will be printed about intermediate stages
!!  nReal      : number of different real roots found
!!  root (n,1) : real part of n-th root
!!  root (n,2) : imaginary part of n-th root
!!
!! NOTES
!!
!!***

subroutine ut_quarticRoots (q3, q2, q1, q0,        &
                            printInfo,             &
                                            nReal, &
                                            root   )
  
  implicit none

  real,    intent (in)  :: q3, q2, q1, q0
  logical, intent (in)  :: printInfo
  integer, intent (out) :: nReal
  real,    intent (out) :: root (1:4,1:2)

  character (len = 12) :: quarticType

  logical :: bisection
  logical :: minimum
  logical :: notZero
  logical :: rescaleQuadratic
  logical :: rootQs, rootQt, rootQu

  integer :: n
  integer :: nReal2, nReal3

  integer, parameter :: Re = 1
  integer, parameter :: Im = 2

  real    :: a0, a1, a2, a3
  real    :: a, b, c, d, k, s, t, u, v, w, x, y, z
  real    :: maxQs, maxQt, maxQu
  real    :: signQs, signQt, signQu
  real    :: tolQs, tolQt, tolQu

  real    :: interval (1:4)
  real    :: root2    (1:2,1:2)
  real    :: root3    (1:3,1:2)

  real, parameter :: accuracy    = 1.d-16
  real, parameter :: tolerance   = 1.d+14

  real, parameter :: zero        = 0.d0
  real, parameter :: one256th    = 3.90625d-3
  real, parameter :: one100th    = 1.d-2
  real, parameter :: one16th     = 6.25d-2
  real, parameter :: three8th    = 3.75d-1
  real, parameter :: three4th    = 7.5d-1
  real, parameter :: fourth      = 2.5d-1
  real, parameter :: third       = 1.d0 / 3.d0
  real, parameter :: half        = 5.d-1
  real, parameter :: one         = 1.d0
  real, parameter :: two         = 2.d0
  real, parameter :: three       = 3.d0
  real, parameter :: four        = 4.d0
!
!
!     ...Start.
!
!
  if (printInfo) then
      write (*,'(a)'        ) ' ------------------------------------------------'
      write (*,'(a,es24.16)') ' initial quartic q3    = ',q3
      write (*,'(a,es24.16)') ' initial quartic q2    = ',q2
      write (*,'(a,es24.16)') ' initial quartic q1    = ',q1
      write (*,'(a,es24.16)') ' initial quartic q0    = ',q0
      write (*,'(a)'        ) ' ------------------------------------------------'
  end if
!
!
!     ...Handle special cases. Since the cubic solver handles all its
!        special cases by itself, we need to check only for two cases:
!
!            1) independent term is zero -> solve cubic and include
!               the zero root
!
!            2) the biquadratic case.
!
!
  if (q0 == zero) then

      k  = one
      a3 = q3
      a2 = q2
      a1 = q1

      quarticType = 'cubic'

  else if (q3 == zero .and. q1 == zero) then

      k  = one
      a2 = q2
      a0 = q0

      quarticType      = 'biquadratic'
      rescaleQuadratic = .true.

  else
!
!
!     ...The general case. Rescale quartic polynomial, such that largest absolute coefficient
!        is (exactly!) equal to 1. Honor the presence of a special quartic case that might have
!        been obtained (due to underflow in the coefficients).
!
!
      s = abs (q3)
      t = sqrt (abs (q2))
      u = abs (q1) ** third
      x = abs (q0) ** fourth
      y = max (s,t,u,x)

      if (y == s) then

          k  = one / s
          a3 = sign (one , q3)
          a2 = (q2 * k) * k
          a1 = ((q1 * k) * k) * k
          a0 = (((q0 * k) * k) * k) * k

      else if (y == t) then

          k  = one / t
          a3 = q3 * k
          a2 = sign (one , q2)
          a1 = ((q1 * k) * k) * k
          a0 = (((q0 * k) * k) * k) * k

      else if (y == u) then

          k  = one / u
          a3 = q3 * k
          a2 = (q2 * k) * k
          a1 = sign (one , q1)
          a0 = (((q0 * k) * k) * k) * k

      else

          k  = one / x
          a3 = q3 * k
          a2 = (q2 * k) * k
          a1 = ((q1 * k) * k) * k
          a0 = sign (one , q0)

      end if

      k = one / k

      if (printInfo) then
          write (*,'(a,es24.16)') ' rescaling factor      = ',k
          write (*,'(a)'        ) ' ------------------------------------------------'
          write (*,'(a,es24.16)') ' rescaled quartic q3   = ',a3
          write (*,'(a,es24.16)') ' rescaled quartic q2   = ',a2
          write (*,'(a,es24.16)') ' rescaled quartic q1   = ',a1
          write (*,'(a,es24.16)') ' rescaled quartic q0   = ',a0
          write (*,'(a)'        ) ' ------------------------------------------------'
      end if

      if (a0 == zero) then
          quarticType = 'cubic'
      else if (a3 == zero .and. a1 == zero) then
          quarticType      = 'biquadratic'
          rescaleQuadratic = .false.
      else
          quarticType = 'general'
      end if

  end if
!
!
!     ...Select the case.
!
!        1) The quartic with independent term = 0 -> solve cubic and add a zero root.
!
!
  select case (quarticType)

    case ('cubic')

      call ut_cubicRoots (a3, a2, a1,         &
                                      nReal3, &
                                      root3   )
      if (nReal3 == 3) then

          x = root3 (1,Re) * k       ! real roots of cubic are ordered x >= y >= z
          y = root3 (2,Re) * k
          z = root3 (3,Re) * k

          nReal = 4

          root (1,Re) = max (x, zero)
          root (2,Re) = max (y, min (x, zero))
          root (3,Re) = max (z, min (y, zero))
          root (4,Re) = min (z, zero)
          root (:,Im) = zero

      else                           ! there is only one real cubic root here

          x = root3 (1,Re) * k

          nReal = 2

          root (1,Re) = max (x, zero)
          root (2,Re) = min (x, zero)
          root (3,Re) = root3 (2,Re) * k
          root (4,Re) = root3 (3,Re) * k
          root (1,Im) = zero
          root (2,Im) = zero
          root (3,Im) = root3 (2,Im) * k
          root (4,Im) = root3 (3,Im) * k

      end if
!
!
!     ...2) The quartic with x^3 and x terms = 0 -> solve biquadratic.
!
!
    case ('biquadratic')

      call ut_quadraticRoots (q2, q0,                   &
                              rescaleQuadratic,         &
                                                nReal2, &
                                                root2   )
      if (nReal2 == 2) then

          x = root2 (1,Re)         ! real roots of quadratic are ordered x >= y
          y = root2 (2,Re)

          if (y >= zero) then

              x = sqrt (x) * k
              y = sqrt (y) * k

              nReal = 4

              root (1,Re) = x
              root (2,Re) = y
              root (3,Re) = - y
              root (4,Re) = - x
              root (:,Im) = zero

          else if (x >= zero .and. y < zero) then

              x = sqrt (x)       * k
              y = sqrt (abs (y)) * k

              nReal = 2

              root (1,Re) = x
              root (2,Re) = - x
              root (3,Re) = zero
              root (4,Re) = zero
              root (1,Im) = zero
              root (2,Im) = zero
              root (3,Im) = y
              root (4,Im) = - y

          else if (x < zero) then

              x = sqrt (abs (x)) * k
              y = sqrt (abs (y)) * k

              nReal = 0

              root (:,Re) = zero
              root (1,Im) = y
              root (2,Im) = x
              root (3,Im) = - x
              root (4,Im) = - y

          end if

      else                             ! complex conjugate pair biquadratic roots x +/- iy.
              
          x = root2 (1,Re)
          y = root2 (1,Im)
          z = sqrt (x * x + y * y)
          y = sqrt (half * (z - x)) * k
          x = sqrt (half * (z + x)) * k

          nReal = 0

          root (1,Re) = x
          root (2,Re) = x
          root (3,Re) = - x
          root (4,Re) = - x
          root (1,Im) = y
          root (2,Im) = - y
          root (3,Im) = y
          root (4,Im) = - y

      end if
!
!
!     ...3) The general quartic case. Search for stationary points. Set the first
!           derivative polynomial (cubic) equal to zero and find its roots.
!
!
    case ('general')

      x = three4th * a3
      y = half     * a2
      z = fourth   * a1

      if (printInfo) then
          write (*,'(a,es24.16)') ' dQ(x)/dx cubic c2     = ',x
          write (*,'(a,es24.16)') ' dQ(x)/dx cubic c1     = ',y
          write (*,'(a,es24.16)') ' dQ(x)/dx cubic c0     = ',z
          write (*,'(a)'        ) ' ------------------------------------------------'
      end if

      call ut_cubicRoots (x, y, z,         &
                                   nReal3, &
                                   root3   )

      if (nReal3 == 1) then

          s = root3 (1,Re)                                     ! Q'(x) root s
          t = - two                                            ! dummy root value, outside the -ve range
          c = s * s                                            ! Q'(x) root s squared
          b = c * s                                            ! Q'(x) root s cubed
          a = b * s                                            ! Q'(x) root s 4-th power (always +ve)
          b = abs (b * a3)                                     ! absolute value of a3 * s^3 in Q(s)
          c = abs (c * a2)                                     ! absolute value of a2 * s^2 in Q(s)
          d = abs (s * a1)                                     ! absolute value of a1 * s   in Q(s)
          x = s + a3                                           ! start Horner scheme for Q(s)
          x = x * s + a2                                       ! Horner step for Q(s)
          x = x * s + a1                                       ! Horner step for Q(s)
          x = x * s + a0                                       ! Q(s), minimum, possibly < 0
          maxQs  = max (a, b, c, d, abs (a0))                  ! maximum value of all absolute Q(s) terms
          signQs = sign (one, x)                               ! sign of Q(s)
          signQt = one                                         ! activate bisection between [s,t=-2]
          signQu = zero                                        ! don't trigger bisections in [t,u] and [u,-2]
          rootQs = abs (x) * tolerance < maxQs                 ! Q'(x) root s is also root of Q(x)
          rootQt = .false.                                     ! don't trigger exact root test at t
          rootQu = .false.                                     ! don't trigger exact root test at u

          if (printInfo) then
              write (*,'(a,es24.16)'  ) ' dQ(x)/dx root s       = ',s
              write (*,'(a)'          ) ' ------------------------------------------------'
              write (*,'(a,2es24.16)' ) ' Q(s), max abs term    = ',x,maxQs
              write (*,'(a)'          ) ' ------------------------------------------------'
              write (*,'(a,es24.16)'  ) ' Sign of Q(s)          = ',signQs
              write (*,'(a,es24.16,a)') ' Sign of Q(t)          = ',signQt,' (activate bisection [s,t])'
              write (*,'(a,es24.16,a)') ' Sign of Q(u)          = ',signQu,' (no bisections [t,u] and [u,-2])'
              write (*,'(a)'          ) ' ------------------------------------------------'
              write (*,'(a,L1)'       ) ' Exact root at s       = ',rootQs
              write (*,'(a,L1,a)'     ) ' Exact root at t       = ',rootQt,' (no exact root test at t)'
              write (*,'(a,L1,a)'     ) ' Exact root at u       = ',rootQu,' (no exact root test at u)'
              write (*,'(a)'          ) ' ------------------------------------------------'
          end if

      else if (nReal3 == 3) then

          s = root3 (1,Re)                                     ! real root order is s >= t >= u
          t = root3 (2,Re)
          u = root3 (3,Re)

          c = s * s                                            ! Q'(x) root s squared
          b = c * s                                            ! Q'(x) root s cubed
          a = b * s                                            ! Q'(x) root s 4-th power (always +ve)
          b = abs (b * a3)                                     ! absolute value of a3 * s^3 in Q(s)
          c = abs (c * a2)                                     ! absolute value of a2 * s^2 in Q(s)
          d = abs (s * a1)                                     ! absolute value of a1 * s   in Q(s)
          x = s + a3                                           ! start Horner scheme for Q(s)
          x = x * s + a2                                       ! Horner step for Q(s)
          x = x * s + a1                                       ! Horner step for Q(s)
          x = x * s + a0                                       ! Q(s), minimum, possibly < 0
          maxQs  = max (a, b, c, d, abs (a0))                  ! maximum value of all absolute Q(s) terms
          tolQs  = abs (x) * tolerance                         ! value of Q(s) * tolerance -> Q(x) root condition
          signQs = sign (one, x)                               ! sign of Q(s)
          rootQs = tolQs < maxQs                               ! Q'(x) root s is also root of Q(x)

          c = t * t                                            ! Q'(x) root t squared
          b = c * t                                            ! Q'(x) root t cubed
          a = b * t                                            ! Q'(x) root t 4-th power (always +ve)
          b = abs (b * a3)                                     ! absolute value of a3 * t^3 in Q(t)
          c = abs (c * a2)                                     ! absolute value of a2 * t^2 in Q(t)
          d = abs (t * a1)                                     ! absolute value of a1 * t   in Q(t)
          y = t + a3                                           ! start Horner scheme for Q(t)
          y = y * t + a2                                       ! Horner step for Q(t)
          y = y * t + a1                                       ! Horner step for Q(t)
          y = y * t + a0                                       ! Q(t), maximum, possibly > 0
          maxQt  = max (a, b, c, d, abs (a0))                  ! maximum value of all absolute Q(t) terms
          tolQt  = abs (y) * tolerance                         ! value of Q(t) * tolerance -> Q(x) root condition
          signQt = sign (one, y)                               ! sign of Q(t)
          rootQt = tolQt < maxQt                               ! Q'(x) root t is also root of Q(x)

          c = u * u                                            ! Q'(x) root u squared
          b = c * u                                            ! Q'(x) root u cubed
          a = b * u                                            ! Q'(x) root u 4-th power (always +ve)
          b = abs (b * a3)                                     ! absolute value of a3 * u^3 in Q(u)
          c = abs (c * a2)                                     ! absolute value of a2 * u^2 in Q(u)
          d = abs (u * a1)                                     ! absolute value of a1 * u   in Q(u)
          z = u + a3                                           ! start Horner scheme for Q(u)
          z = z * u + a2                                       ! Horner step for Q(u)
          z = z * u + a1                                       ! Horner step for Q(u)
          z = z * u + a0                                       ! Q(u), minimum, possibly < 0
          maxQu  = max (a, b, c, d, abs (a0))                  ! maximum value of all absolute Q(u) terms
          tolQu  = abs (z) * tolerance                         ! value of Q(u) * tolerance -> Q(x) root condition
          signQu = sign (one, z)                               ! sign of Q(u)
          rootQu = tolQu < maxQu                               ! Q'(x) root u is also root of Q(x)

          if (printInfo) then
              write (*,'(a,es24.16)'  ) ' dQ(x)/dx root s       = ',s
              write (*,'(a,es24.16)'  ) ' dQ(x)/dx root t       = ',t
              write (*,'(a,es24.16)'  ) ' dQ(x)/dx root u       = ',u
              write (*,'(a)'          ) ' ------------------------------------------------'
              write (*,'(a,2es24.16)' ) ' Q(s), max abs term    = ',x,maxQs
              write (*,'(a,2es24.16)' ) ' Q(t), max abs term    = ',y,maxQt
              write (*,'(a,2es24.16)' ) ' Q(u), max abs term    = ',z,maxQu
              write (*,'(a)'          ) ' ------------------------------------------------'
              write (*,'(a,es24.16)'  ) ' Sign of Q(s)          = ',signQs
              write (*,'(a,es24.16)'  ) ' Sign of Q(t)          = ',signQt
              write (*,'(a,es24.16)'  ) ' Sign of Q(u)          = ',signQu
              write (*,'(a)'          ) ' ------------------------------------------------'
              write (*,'(a,L1)'       ) ' Exact root at s       = ',rootQs
              write (*,'(a,L1,a)'     ) ' Exact root at t       = ',rootQt
              write (*,'(a,L1,a)'     ) ' Exact root at u       = ',rootQu
              write (*,'(a)'          ) ' ------------------------------------------------'
          end if

      end if
!
!
!     ...Collect bisection ranges and initial bisection or exact roots. Work from the right
!        edge (+2) of the +ve x-axis to the left edge (-2) of the x-axis. Any real roots found
!        (either as exact roots of Q(x) or during bisection) will thus be in proper order.
!
!        Deal with the following Q(x) graph section:
!
!
!                                            Q(x)
!                                             |
!                                     /       |
!                                    /        |
!                                   /         |
!                       ------------------|---|---> x-axis
!                                 /      +2   |
!                            \   /            |
!                              s              | <- Q(s)
!                                             |
!
!
!
      nReal = 0                                         ! counter of real roots in Q(x)

      if (rootQs) then                                  ! exact root of Q(x) at s ?

          nReal = 1

          root     (1,Re) = s                           ! exact root at s
          interval (1   ) = zero                        ! bisection interval zero -> exact root

          signQs = zero                                 ! don't trigger a bisection in [t,s]

      else if (signQs < zero) then                      ! root between s and +2 ?

          nReal = nReal + 1

          if (s < zero) then                            ! is the Q(x) axis between s and +2 ?
              if (a0 > zero) then
                  root     (nReal,Re) = s               ! initial root (search towards more +ve values)
                  interval (nReal   ) = - s             ! bisection interval [s,0]    s < 0
              else
                  root     (nReal,Re) = zero            ! initial root (search towards more +ve values)
                  interval (nReal   ) = two             ! bisection interval [0,2]    s < 0 
              end if
          else
              Root     (nReal,Re) = s                   ! initial root (search towards more +ve values)
              interval (nReal   ) = two - s             ! bisection interval [s,2]   s > 0
          end if

      end if
!
!
!     ...Deal with the following Q(x) graph section:
!
!
!                                             Q(x)
!                                              |
!                             t                | <- Q(t)
!                           /   \              |
!                                \             |
!                       -----------------------|---> x-axis
!                                  \           |
!                                   \   /      |
!                                     s        | <- Q(s)
!                                              |
!
!
      if (rootQt) then                                  ! exact root of Q(x) at t ?

          nReal = nReal + 1

          root     (nReal,Re) = t                       ! exact root at t
          interval (nReal   ) = zero                    ! bisection interval zero -> exact root

          signQt = zero                                 ! don't trigger a bisection in [u,t]

      else if (signQt * signQs < zero) then             ! crossing the x-axis -> root between t and s ?

          nReal = nReal + 1

          if (t < zero .and. s > zero) then             ! is the x = 0 axis between t and s ?
              if (a0 > zero) then
                  root     (nReal,Re) = s               ! initial root (search towards more -ve values)
                  interval (nReal   ) = - s             ! bisection interval [0,s]    t < 0 and s > 0
               else
                  root     (nReal,Re) = zero            ! initial root (search towards more -ve values)
                  interval (nReal   ) = t               ! bisection interval [t,0]    t < 0 and s > 0
              end if
          else
              root     (nReal,Re) = s                   ! initial root (search towards more -ve values)
              interval (nReal   ) = t - s               ! bisection interval [t,s]    t,s < 0 or t,s > 0
          end if

      end if
!
!
!     ...Deal with the following Q(x) graph section:
!
!
!                                             Q(x)
!                                              |
!                                      t       | <- Q(t)
!                                    /   \     |
!                                   /          |
!                       -----------------------|---> x-axis
!                                 /            |
!                            \   /             |
!                              u               | <- Q(u)
!                                              |
!
!
      if (rootQu) then                                  ! exact root of Q(x) at u ?

          nReal = nReal + 1

          root     (nReal,Re) = u                       ! exact root at u
          interval (nReal   ) = zero                    ! bisection interval zero -> exact root

          signQu = zero                                 ! don't trigger a bisection in [-2,u]

      else if (signQu * signQt < zero) then             ! crossing the x-axis -> root between u and t ?

          nReal = nReal + 1

          if (u < zero .and. t > zero) then             ! is the x = 0 axis between u and t ?
              if (a0 > zero) then
                  root     (nReal,Re) = u               ! initial root (search towards more +ve values)
                  interval (nReal   ) = - u             ! bisection interval [u,0]    u < 0 and t > 0
              else
                  root     (nReal,Re) = zero            ! initial root (search towards more +ve values)
                  interval (nReal   ) = t               ! bisection interval [0,t]    u < 0 and t > 0
              end if
          else
              root     (nReal,Re) = u                   ! initial root (search towards more +ve values)
              interval (nReal   ) = t - u               ! bisection interval [u,t]    u,t < 0 or u,t > 0
          end if

      end if
!
!
!     ...Deal with the following Q(x) graph section:
!
!
!                                             Q(x)
!                                              |
!                              \               |
!                               \              |
!                                \             |
!                       ---|-------------------|---> x-axis
!                         -2       \           |
!                                   \   /      |
!                                     u        | <- Q(u)
!                                              |
!
!
      if (signQu < zero) then                           ! crossing the x-axis -> root between -2 and u ?

          nReal = nReal + 1

          if (u > zero) then                            ! is the x = 0 axis between -2 and u ?
              if (a0 > zero) then
                  root     (nReal,Re) = u               ! initial root (search towards more -ve values)
                  interval (nReal   ) = - u             ! bisection interval [0,u]    u > 0
              else
                  root     (nReal,Re) = zero            ! initial root (search towards more -ve values)
                  interval (nReal   ) = - two           ! bisection interval [-2,0]   u > 0 
              end if
          else
              root     (nReal,Re) = u                   ! initial root (search towards more -ve values)
              interval (nReal   ) = - two - u           ! bisection interval [-2,u]   u < 0
          end if

      end if
!
!
!     ...Do all necessary bisections (if any).
!
!
      if (nReal > 0) then

          if (nReal > 4) then
              write (*,*) ' Quartic polynomial has > 4 real roots! WRONG! '
              stop
          end if

          do n = 1,nReal
             t = interval (n)                                     ! initial interval (if 0 -> exact root)
             if (t /= zero) then                                  ! do bisection, if interval length > 0
                 s = root (n,Re)                                  ! initial root
                 do while (abs (t) > abs (s * accuracy))          ! bisection iterates
                    t = t * half                                  ! new interval
                    u = s                                         ! save copy of root (avoid roundoff errors)
                    s = s + t                                     ! new root trial s
                    x = s + a3                                    ! start Horner scheme for Q(s)
                    x = x * s + a2                                ! Horner step for Q(s)
                    x = x * s + a1                                ! Horner step for Q(s)
                    x = x * s + a0                                ! quartic Q(x) value at new root trial
                    if (x >= zero) s = u                          ! retain former root if step was bad
                    if (printInfo) write (*,'(a,es24.16)'  )      ' Bisection root        = ',s
                 end do
                 root (n,Re) = s                                  ! store the bisected real root
                 if (printInfo) write (*,'(a)') ' ------------------------------------------------'
             end if
          end do
!
!
!     ...Find remaining roots (if any).
!
!
          if (nReal == 4) then            ! all 4 real roots found

              root (1,Re) = root (1,Re) * k
              root (2,Re) = root (2,Re) * k
              root (3,Re) = root (3,Re) * k
              root (4,Re) = root (4,Re) * k
              root (1,Im) = zero
              root (2,Im) = zero
              root (3,Im) = zero
              root (4,Im) = zero

              return

          else if (nReal < 4) then            ! deflate to cubic, if 1,2 or 3 real roots

              w = root (1,Re)                 ! 1st rescaled root of Q(x) -> w
              c = w * w                       ! 1st rescaled root squared -> w^2
              b = c * w                       ! 1st rescaled root cubed   -> w^3
              a = max (b * w , abs (b * a3))  ! maximum between |w^4| and |a3*w^3|
              b = abs (c * a2)                ! value of |a2*w^2|
              c = abs (w * a1)                ! value of |a1*w|
              d = abs (a0)                    ! value of |a0|
              x = max (a,b,c,d)               ! maximum of all terms
              w = w * k                       ! form 1st original Q(x) root

              if (x == a) then
                  z = one / w
                  u = - q0 * z                ! u -> backward deflation on original Q(x)
                  t = (u - q1) * z            ! t -> backward deflation on original Q(x)
                  s = (t - q2) * z            ! s -> backward deflation on original Q(x)
              else if (x == b) then
                  z = one / w
                  u = - q0 * z                ! u -> backward deflation on original Q(x)
                  t = (u - q1) * z            ! t -> backward deflation on original Q(x)
                  s = q3 + w                  ! s ->  forward deflation on original Q(x)
              else if (x == c) then
                  s = q3 + w                  ! s ->  forward deflation on original Q(x)
                  t = q2 + s * w              ! t ->  forward deflation on original Q(x)
                  u = - q0 / w                ! u -> backward deflation on original Q(x)
              else
                  s = q3 + w                  ! s ->  forward deflation on original Q(x)
                  t = q2 + s * w              ! t ->  forward deflation on original Q(x)
                  u = q1 + t * w              ! u ->  forward deflation on original Q(x)
              end if

              root (1,Re) = w                 ! store 1st original Q(x) root
              root (1,Im) = zero              !

              if (printInfo) then
                  write (*,'(a,es24.16)') ' Residual cubic c2     = ',s
                  write (*,'(a,es24.16)') ' Residual cubic c1     = ',t
                  write (*,'(a,es24.16)') ' Residual cubic c0     = ',u
                  write (*,'(a)'        ) ' ------------------------------------------------'
              end if

          end if
!
!
!     ...At this stage, two situations can happern:
!
!             1) only 1 real root of Q(x) was found -> solve deflated unscaled cubic x^3 + sx^2 + tx + u
!             2) more than 1 real root of Q(x) is present -> deflate unscaled cubic to unscaled quadratic.
!
!
          if (nReal > 1) then                 ! deflate to quadratic, if 2 or 3 real roots

              w = root (2,Re) * k             ! 2nd original root of Q(x) -> w
              b = w * w                       ! 2nd original root squared -> w^2
              a = max (b * w , abs (b * s))   ! maximum between |w^3| and |s*w^2|
              b = abs (w * t)                 ! value of |t*w|
              c = abs (u)                     ! value of |u|
              x = max (a,b,c)                 ! maximum of all terms

              if (x == a) then
                  z = one / w
                  y = t                       ! save cubic 't'
                  t = - u * z                 ! t -> backward deflation on unscaled cubic
                  s = (t - y) * z             ! s -> backward deflation on unscaled cubic
              else if (x == b) then
                  s = s + w                   ! s ->  forward deflation on unscaled cubic
                  t = - u / w                 ! t -> backward deflation on unscaled cubic
              else
                  s = s + w                   ! s ->  forward deflation on unscaled cubic
                  t = t + s * w               ! t ->  forward deflation on unscaled cubic
              end if

              root (2,Re) = w             ! store 2nd original Q(x) root
              root (2,Im) = zero

              if (printInfo) then
                  write (*,'(a,es24.16)') ' Residual quadratic q1 = ',s
                  write (*,'(a,es24.16)') ' Residual quadratic q0 = ',t
                  write (*,'(a)'        ) ' ------------------------------------------------'
              end if

          else                            ! only 1 real root was found -> solve deflated cubic

              call ut_cubicRoots (s,t,u,         &
                                         nReal3, &
                                         root3   )
              if (nReal3 == 3) then

                  x = root3 (1,Re)        ! real roots of cubic are ordered x >= y >= z
                  y = root3 (2,Re)
                  z = root3 (3,Re)
                  w = root  (1,Re)        ! root of Q(x) found so far

                  nReal = 4

                  root (1,Re) = max (x, w)
                  root (2,Re) = max (y, min (x, w))
                  root (3,Re) = max (z, min (y, w))
                  root (4,Re) = min (z, w)
                  root (2,Im) = zero
                  root (3,Im) = zero
                  root (4,Im) = zero

              else                        ! there is only one real cubic root here

                  x = root3 (1,Re)
                  w = root  (1,Re)        ! root of Q(x) found so far

                  nReal = 2

                  root (1,Re) = max (x, w)
                  root (2,Re) = min (x, w)
                  root (3,Re) = root3 (2,Re)
                  root (4,Re) = root3 (3,Re)
                  root (2,Im) = zero
                  root (3,Im) = root3 (2,Im)
                  root (4,Im) = root3 (3,Im)

              end if

              return

          end if
!
!
!     ...At this stage, again two situations can happern:
!
!             1) only 2 real root of Q(x) were found -> solve deflated unscaled quadratic x^2 +sx + t
!             2) 3 real roots of Q(x) are present -> deflate unscaled quadratic to unscaled linear.
!
!
          if (nReal > 2) then                 ! deflate to linear, if 3 real roots

              z = root (3,Re) * k             ! 3rd original root of Q(x) -> z
              a = max (z * z , abs (z * s))   ! maximum between |z^2| and |s*z|
              b = abs (t)                     ! value of |t|

              if (a >= b) then
                  w = t / z                   ! backward deflation to 4th original Q(x) root
              else
                  w = - (s + z)               !  forward deflation to 4th original Q(x) root
              end if

              x = root (1,Re)                 ! roots ordered x > y > z -> root w unknown
              y = root (2,Re)

              nReal = 4

              root (1,Re) = max (x, w)
              root (2,Re) = max (y, min (x, w))
              root (3,Re) = max (z, min (y, w))
              root (4,Re) = min (z, w)
              root (3,Im) = zero
              root (4,Im) = zero

          else                            ! only 2 real roots were found -> solve deflated quadratic

              rescaleQuadratic = .true.   ! coefficients of quadratic might be large

              call ut_quadraticRoots (s, t,                     &
                                      rescaleQuadratic,         &
                                                        nReal2, &
                                                        root2   )

              if (nReal2 == 2) then

                  x = root2 (1,Re)        ! real roots of quadratic are ordered x >= y
                  y = root2 (2,Re)
                  s = root  (1,Re)        ! real roots found so far of Q(x), ordered s >= t
                  t = root  (2,Re)

                  root (1,Re) = max (x,s)
                  root (2,Re) = min (max (x,t) , max (y,s))
                  root (3,Re) = max (min (x,t) , min (y,s))
                  root (4,Re) = min (y,t)
                  root (3,Im) = zero
                  root (4,Im) = zero

                  nReal = 4

              else                        ! the case of two extra complex roots

                  root (3,Re) = root2 (1,Re)
                  root (4,Re) = root2 (2,Re)
                  root (3,Im) = root2 (1,Im)
                  root (4,Im) = root2 (2,Im)

              end if                      

          end if

      else             ! # of real roots found is 0
!
!
!     ...If no real roots have been found by now, only complex roots are possible.
!        Find real parts of roots first.
!
!
          s = a3 * half
          t = s * s - a2
          u = s * t + a1                                      ! value of sqrt [H(-a3/4)] at stationary point -a3/4

          if (printInfo) then
              write (*,'(a,es24.16)') ' sqrt [H(-a3/4)] value = ',u
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if

          if (q3 /= zero) then
              s = q1 / q3
              minimum = (q0 > s * s)
          else
              minimum   = (4 * q0 > q2 * q2)                  ! H''(-a3/4) > 0 -> minimum
          end if

          notZero   = (abs (u) >= accuracy)                   ! H(-a3/4) is considered > 0 at stationary point
          bisection = notZero .or. (.not.notZero .and. minimum)

          v = a2 + a2                                         ! 2a2 needed for calculating derivatives of Q(x)
          w = a3 + a3 + a3                                    ! 3a3 needed for calculating derivatives of Q(x)

          if (bisection) then                                 ! bisection applied to smallest real component

              x = - fourth * a3                               ! local minimum of H(x)

              if (a3 > zero) then                             ! bisect [x,+2]  (x = negative)
                  s = two                                     ! initial root (search towards more -ve values)
                  t = x - two                                 ! negative bisection interval
              else                                            ! bisect [-2,x]  (x = positive)
                  s = - two                                   ! initial root (search towards more +ve values)
                  t = x + two                                 ! positive bisection interval
              end if

              do while (abs (t) > abs (s * accuracy))         ! bisection iterates on H(x)
                 t = t * half                                 ! new interval
                 u = s                                        ! save copy of root (avoid roundoff errors)
                 s = s + t                                    ! new root trial s
                 a = s + a3                                   ! start Horner scheme for Q(s)
                 a = a * s + a2                               ! Horner step for Q(s)
                 a = a * s + a1                               ! Horner step for Q(s)
                 a = a * s + a0                               ! quartic Q(s) at new root trial -> 1st term of H(s)
                 y = s + s                                    ! 2s
                 z = y + y                                    ! 4s
                 x = z + w                                    ! start Horner scheme for Q'(s)
                 x = x * s + v                                ! Horner step for Q'(s)
                 x = x * s + a1                               ! Q'(s)
                 b = - x                                      ! initialize 2nd term of H(s)
                 c = x * x                                    ! finished 3rd term of H(s)
                 x = z + y + w                                ! start Horner scheme for Q''(s)
                 x = x * s + a2                               ! Q''(s)
                 b = b * x                                    ! update 2nd term of H(s)
                 x = z + a3                                   ! Q'''(s)
                 a = a * x * x                                ! finished 1st term of H(s)
                 b = b * x                                    ! finished 2nd term of H(s)
                 x = a + b + c                                ! H(s) value at new root trial
                 if (x >= zero) s = u                         ! undo root step if bad
                 if (printInfo) write (*,'(a,es24.16)'  )     ' Bisection H(x) root   = ',s
              end do
              if (printInfo) write (*,'(a)') ' ------------------------------------------------'

              a = s * k                                       ! 1st real component -> a
              b = - half * q3 - a                             ! 2nd real component -> b

              z = a + a                                       ! 2a
              z = z + z                                       ! 4a
              z = z + q3                                      ! Q'''(a)
              x = one / z                                     ! 1 / Q'''(a)
              z = z + q3 + q3                                 ! start Horner scheme for Q'(a)
              z = z * a + q2 + q2                             ! Horner step for Q'(a)
              z = z * a + q1                                  ! Q'(a)
              x = x * z                                       ! Q'(a) / Q'''(a)
              x = max (x,zero)                                ! ensure >= 0 value
              c = sqrt (x)                                    ! 1st imaginary component -> c

              z = b + b                                       ! 2b
              z = z + z                                       ! 4b
              z = z + q3                                      ! Q'''(b)
              x = one / z                                     ! 1 / Q'''(b)
              z = z + q3 + q3                                 ! start Horner scheme for Q'(b)
              z = z * b + q2 + q2                             ! Horner step for Q'(b)
              z = z * b + q1                                  ! Q'(b)
              x = x * z                                       ! Q'(b) / Q'''(b)
              x = max (x,zero)                                ! ensure >= 0 value
              d = sqrt (x)                                    ! 2nd imaginary component -> d

          else                                                ! no bisection -> real components equal

              a = - fourth * q3                               ! 1st real component -> a
              b = a                                           ! 2nd real component -> b = a

              u = a + a
              s = u * (a + u) - q2                            ! - Q''(a) -> linear term of quadratic
              s = min (s,zero)                                ! force s to be =< 0
              y = - q0 - q0
              y = y + y
              y = y + q1 * u + q2 * (q2 - u * u)              ! discriminant
              t = fourth * (s * s - y)                        ! Q(a) -> independent term of quadratic
              y = max (y,zero)                                ! force discriminant to be >= 0
              y = sqrt (y)                                    ! square root of discriminant
              y = - half * (s - y)                            ! 1st root^2 of biquadratic
              z = t / y                                       ! 2nd root^2 of biquadratic

              if (y < zero .or. z < zero) then                ! safety net against -ve roots
                  c = zero
                  d = zero
              else
                  c = max (y,z)                               ! to ensure c > d
                  d = min (y,z)                               ! to ensure c > d
                  c = sqrt (c)                                ! 1st positive root of biquadratic
                  d = sqrt (d)                                ! 2nd positive root of biquadratic
              end if

          end if

          if (a > b) then
              root (1,Re) = a
              root (2,Re) = a
              root (3,Re) = b
              root (4,Re) = b
              root (1,Im) = c
              root (2,Im) = - c
              root (3,Im) = d
              root (4,Im) = - d
          else if (a < b) then
              root (1,Re) = b
              root (2,Re) = b
              root (3,Re) = a
              root (4,Re) = a
              root (1,Im) = d
              root (2,Im) = - d
              root (3,Im) = c
              root (4,Im) = - c
          else
              root (1,Re) = a
              root (2,Re) = a
              root (3,Re) = a
              root (4,Re) = a
              root (1,Im) = c
              root (2,Im) = -c
              root (3,Im) = d
              root (4,Im) = - d
          end if

      end if    ! # of real roots 'if'

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine ut_quarticRoots
