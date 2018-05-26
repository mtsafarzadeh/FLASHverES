!!****if* source/physics/sourceTerms/Chemistry/ChemistryIntegrate3/chem_netIntRosen
!!
!! NAME
!!     chem_netIntRosen  
!!
!! SYNOPSIS
!!  subroutine bn_rosenMa28(real(IN)       ::y(:),
!!                          real(IN)       ::dydx(:),
!!                          integer(IN)    ::n,
!!                          real(IN)       ::x,
!!                          real(IN)       ::htry,
!!                          real(IN)       ::eps,
!!                          real(IN)       ::yscal(:),
!!                          real(OUT)      ::hdid,
!!                          real(OUT)      ::hnext, 
!!                          procedure(IN)  :: derivs,
!!                          procedure(IN)  :: jakob,
!!                          real(OUT)	   :: jcounts) 
!!
!! DESCRIPTION
!!
!!  fourth order rosenbrock step for integrating stiff ode's with monitoring 
!!  of local truncation error to adjust the stepsize. input are the dependent  
!!  variable y(1:n) and its derivative dydx(1:n), at the starting value of the  
!!  independent variable x. also input are the stepsize to be attempted htry,  
!!  the desired fractional accuracy eps, and the vector yscal(1:n) against  
!!  which the error is scaled. derivs is a routine that computes the righthand  
!!  side of the system of ode's. on output the y and x are replaced by their  
!!  new values, hdid is the stepsize that was actually accomplished and hnext  
!!  is the estimate for the next stepsize. 
!!   
!!  nmax is the maximum number of ode's.  
!!  grow and shrnk are the extremes of the time step growth  
!!  errcon = (grow/safety)**(1/pgrow) handles the case when errmax = 0.0 
!!  
!!  routine modified to assume dfdx = 0.0 for reaction networks
!!  
!!  sparse ma28 algebra version
!!
!! ARGUMENTS
!!
!!   y       - dependent variable, array of size y(1:n)
!!   dydx    - derivative of dependent variable, array of size dydx(1:n)
!!   n       - number of dependent variables
!!   x       - independent variable
!!   htry    - attempted stepsize 
!!   eps     - desired fractional accuracy
!!   yscal   - vector of size yscal(1:n) for scaling error
!!   hdid    - stepsize used
!!   hnext   - estimate for next stepsize
!!   derivs  - procedure(IN) name of the routine that contains the odes
!!   jakob   - procedure(IN) name of the routine that contains the jacobian of the odes
!!   bjakob  - procedure(IN) name of the routine that sets the pointers of the sparse jacobian
!!
!! NOTES
!!  this file contains 2 routines to drive the integration of 
!!  nuclear reaction networks with rosenbrock time stepper
!!  and 2 choices for the linear algebra.
!!  routine bn_rosenMa28 drives a rosenbrock step with ma28 algebra
!!  routine bn_rosenGift drives a rosenbrock step with gift algebra
!!
!!***
!!---------------------------------------------------------------------------------

subroutine Chemistry_netIntRosen(y,dydx,n,x,htry,eps,yscal,hdid,hnext,  & 
     &                      derivs,jakob,jcounts) 

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stampMessage
  ! can't use jakob interface; see notes in bnNetwork_interface for the mystery.
!!  use bnIntegrate_interface, ONLY: derivs, bjakob
  
  use Chemistry_interface, ONLY : Chemistry_sparsePointers
  use Chemistry_data
  use Chemistry_dataNetworkSize
#include "constants.h"

  implicit none

!! argument declarations
  external               jakob, derivs
  integer, intent(IN) :: n
  real, intent(IN)    :: yscal(n), dydx(n), htry, eps
  real, intent(INOUT) :: x, y(n)
  real, intent(OUT)   :: hdid, hnext, jcounts


!! local variables-----------------------------------------------
  integer, save       :: i,jtry
  integer, parameter  :: nmax=30, maxtry=400
  integer, save       :: nDummy ! to make the arguments consistent of jakob
  real, save          :: errmax,h,xsav,dysav(nmax),err(nmax),g1(nmax),  & 
       &                 g2(nmax),g3(nmax),g4(nmax),ysav(nmax),xx 

  real, parameter     ::  safety=0.9e0, grow=1.5e0, pgrow=-0.25e0,  & 
       &                  shrnk=0.5e0,  pshrnk=-1.0e0/3.0e0,errcon=0.1296e0

!!  shampine parameter set 
  real, parameter     ::  gamf =  1.0e0/2.0e0,    a21 =  2.0e0,  & 
       &                  a31 =  48.0e0/25.0e0,  a32 =  6.0e0/25.0e0,  & 
       &                  c21 = -8.0e0,          c31 =  372.0e0/25.0e0,  & 
       &                  c32 =  12.0e0/5.0e0,   c41 = -112.0e0/125.0e0,  & 
       &                  c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0,  & 
       &                  b1  =  19.0e0/9.0e0,   b2  =  1.0e0/2.0e0,  & 
       &                  b3  =  25.0e0/108.0e0, b4  =  125.0e0/108.0e0,  & 
       &                  e1  =  17.0e0/54.0e0,  e2  =  7.0e0/36.0e0,  & 
       &                  e3  =  0.0e0,          e4  =  125.0e0/108.0e0,  & 
       &                  c1x =  1.0e0/2.0e0,    c2x = -3.0e0/2.0e0,  & 
       &                  c3x =  121.0e0/50.0e0, c4x =  29.0e0/250.0e0,  & 
       &                  a2x =  1.0e0,          a3x = 3.0e0/5.0e0


!!  for the ma28 package
  integer, parameter  :: naij=200, n5 = 5*nmax, n8=8*nmax
  logical, save       :: firstCall = .false.
  integer, save       :: iloc(naij),jloc(naij), flag, nzo,  &
       &                 ivect(naij),jvect(naij),ikeep(n5),iw(n8)
  real, save          :: dfdy(naij),amat(naij),w(nmax),u

!!----------------------------------------------------------------------


!!  get and copy the nonzero locations 
  if (.NOT. firstCall) then 
     firstCall = .true.

     nzo     = 0
     do i=1,naij
        iloc(i) = 0
        jloc(i) = 0
     end do

   !!  get the sparse pattern
     call Chemistry_sparePointers(iloc,jloc,nzo,naij)
   !!  print *,'nzo: ', nzo
   !!  copy the location
     do i=1,nzo
        ivect(i) = iloc(i)
        jvect(i) = jloc(i)
     enddo

   !!  force the diagonal to be the pivot elements 
     do i=1,nzo 
        amat(i) = 1.0e-10 
        if (ivect(i) .eq. jvect(i)) amat(i) = 1.0e0 
     enddo

     u  = 0.1e0
     call ma28ad(n,nzo,amat,naij,iloc,naij,jloc,u,ikeep,iw,w,flag)

     if (flag .lt. 0) then
        call Logfile_stampMessage( '[chem_rosenMa28] negative flag returned from ma28ad')
        call Logfile_stampMessage( '[chem_rosenMa28] more than stpmax steps required')
        call Driver_abortFlash('ERROR: more than stpmax steps required in bn_rosenMa28')
     end if

  end if  ! of firstCall =.false.
  !!nDummy = n

!!  store the initial values 
  xsav = x 
  do i=1,n 
     ysav(i)  = y(i) 
     dysav(i) = dydx(i) 
  enddo

!!  get the sparse jacobian in sparse_dfdy 
  call jakob(xsav,ysav,dfdy,nzo,nDummy) 

!!  main loop 
  h = htry 
  do jtry = 1,maxtry 

   !!   form the a matrix and decompose it 
     xx = 1.0e0/(gamf*h)
     do i=1,nzo 
        amat(i) = -dfdy(i) 
        if (ivect(i) .eq. jvect(i)) amat(i) = xx + amat(i) 
     enddo

   !!  numeric decomp
   !!  print *,'n: ', n, ', nzo: ', nzo
     call ma28bd(n,nzo,amat,naij,ivect,jvect,jloc,ikeep,iw,w,flag) 

     if (flag .lt. 0) then 
        !!if (bn_meshMe .EQ. MASTER_PE) print *, 'ERROR negative flag in ma28bd flag',flag
        call Logfile_stampMessage( '[chem_rosenMa28] negative flag returned from ma28bd')
        call Driver_abortFlash('ERROR: negative return flag in ma28bd')
     end if


   !!  set up and solve the right hand side for g1 
     do i=1,n 
        g1(i) = dysav(i) 
     enddo
     call ma28cd(n,amat,naij,jloc,ikeep,g1,w,1) 
     g1(iELEC) = g1(iHP) + g1(iHEP) + 2.0*g1(iHEPP)
     !g1(iELEC) = g1(iHP) + g1(iHEP) + 2.0*g1(iHE2P)

   !!  compute intermediate values of y,x and dydx 
     do i=1,n 
        y(i) = ysav(i) + a21 * g1(i) 
     enddo
     x = xsav + a2x * h 
     call derivs(x,y,dydx) 


   !!  set up and solve the right hand side for g2 
     do i=1,n 
        g2(i) = dydx(i) + c21*g1(i)/h 
     enddo
     call ma28cd(n,amat,naij,jloc,ikeep,g2,w,1) 
     g2(iELEC) = g2(iHP) + g2(iHEP) + 2.0*g2(iHEPP)
     !g2(iELEC) = g2(iHP) + g2(iHEP) + 2.0*g2(iHE2P)

   !!  compute intermediate values of y,x and dydx 
     do i=1,n 
        y(i) = ysav(i) + a31*g1(i) + a32*g2(i) 
     enddo
     x = xsav + a3x*h 
     call derivs(x,y,dydx) 

   !!  set up and solve the right hand side for g3 
     do i=1,n 
        g3(i)  = dydx(i) + (c31*g1(i) + c32*g2(i))/h 
     enddo
     call ma28cd(n,amat,naij,jloc,ikeep,g3,w,1) 
     g3(iELEC) = g3(iHP) + g3(iHEP) + 2.0*g3(iHEPP)
     !g3(iELEC) = g3(iHP) + g3(iHEP) + 2.0*g3(iHE2P)

   !!  set up and solve the right hand side for g4 
     do i=1,n 
        g4(i)  = dydx(i) + (c41*g1(i) + c42*g2(i) + c43*g3(i))/h 
     end do
     call ma28cd(n,amat,naij,jloc,ikeep,g4,w,1) 
     g4(iELEC) = g4(iHP) + g4(iHEP) + 2.0*g4(iHEPP)
     !g4(iELEC) = g4(iHP) + g4(iHEP) + 2.0*g4(iHE2P)

   !!  compute the third and fourth order estimates of y 
     do i=1,n 
        y(i)   = ysav(i) + b1*g1(i) + b2*g2(i) + b3*g3(i) + b4*g4(i) 
        err(i) = e1*g1(i) + e2*g2(i) + e3*g3(i) + e4*g4(i) 
     enddo
     y(iELEC) = y(iHP) + y(iHEP) + 2.0*y(iHEPP)
     !y(iELEC) = y(iHP) + y(iHEP) + 2.0*y(iHE2P)
     x = xsav + h 

     if (x .eq. xsav) then 
        !!if (bn_meshMe .EQ. MASTER_PE) print *, 'step size not significant in bn_rosenMa28'
        call Logfile_stampMessage( '[chem_rosenMa28] Step size not significant!')
        call Driver_abortFlash('ERROR: step size not significant in bn_rosenMa28')
     end if


   !!  determine the scaled accuracy 
     errmax = 0.0e0 
     do i=1,n 
        errmax = max(errmax,abs(err(i)/yscal(i))) 
     enddo
     errmax = errmax/eps 

   !!  if the step succeded, compute the size of the next step and return 
     if (errmax .le. 1.0) then 
        hdid = h 
        if (errmax .gt. errcon) then 
           hnext = safety * h * errmax**pgrow 
        else 
           hnext = grow * h 
        end if
	jcounts = jcounts + jtry
        return 

      !!   if the step did not succeed cut the stepsize and try again 
     else 
        hnext = safety * h * errmax**pshrnk 
        h     = sign(max(abs(hnext),shrnk*abs(h)),h) 
     end if
  enddo

!!  too many tries
  !!if (bn_meshMe .EQ. MASTER_PE) print *, 'ERROR exceeded maxtry in routine bn_rosenMa28' 
  call Logfile_stampMessage('[chem_rosenMa28] ERROR exceeded maxtry')
  call Driver_abortFlash('ERROR: exceeded maxtry in bn_rosenMa28')

end subroutine Chemistry_netIntRosen

