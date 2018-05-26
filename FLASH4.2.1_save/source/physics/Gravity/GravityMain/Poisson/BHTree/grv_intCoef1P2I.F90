!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_intCoef1P2I
!!
!! NAME
!!
!!  grv_intCoef1P2I
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!   grv_intCoef1P2I integrates by Simpson method coefficients
!!   in the Fourrier series of the long distance interactions in
!!   the Ewald formula for the gravitational potential with periodic 
!!   boundaries in 1 direction and isolated boundaries in the other two
!!   directions. The coefficients are calculated by function grv_coef1P2I.
!!   
!!
!! ARGUMENTS
!!
!!  l            - integer
!!  eta, dzeta   - real value
!!
!! RESULT
!!
!!  Value of the ...
!!
!! NOTES
!!
!!***

#include "Flash.h"

real function grv_intCoef1P2I(l,eta,dzeta)

#if defined(__INTEL_COMPILER)
  use IFPORT
#endif

  implicit none

  real,intent(in) :: eta,dzeta
  integer,intent(in) :: l
  real xmin,xmax,x,dx
  real s,f0,f1,f2,grv_coef1P2I
  integer n,k

! number of integration steps and range of integration
  n=300
  xmin=1.0e-8
  xmax=5.0

  s = 0.0
  x=xmin
  dx=0.5*(xmax-xmin)/real(n)

! integration by Simpson method
  do k=0,n
    f0=grv_coef1P2I(x,l,eta,dzeta)
    f1=grv_coef1P2I(x+dx,l,eta,dzeta)
    f2=grv_coef1P2I(x+2.0*dx,l,eta,dzeta)
    s=s+(f0+f2+4.0*f1)
    x=x+2.0*dx
  enddo

  grv_intCoef1P2I=s*dx/3.0

end

! integrated function
real function grv_coef1P2I(x,l,eta,dzeta)
#if defined FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_besj0
#elif defined(__INTEL_COMPILER)
  use IFPORT
#endif

  implicit none

  real,intent(in) :: eta,dzeta,x
  integer,intent(in) :: l

#if defined USER_BESSEL_J0
#define BESSEL_J0 USER_BESSEL_J0
  real, external :: USER_BESSEL_J0
#elif defined(FLASH_USE_SPECFUN)
#define BESSEL_J0 grv_besj0
#else
#define BESSEL_J0 bessel_j0
  intrinsic bessel_j0
#endif

  if (l.eq.0) then
    grv_coef1P2I=(BESSEL_J0(x*eta)-1.0)*exp(-dzeta*x*x)/x
  else
    grv_coef1P2I=x*BESSEL_J0(x*eta)*exp(-dzeta*x*x)/(x**2+l**2)
  endif

end
