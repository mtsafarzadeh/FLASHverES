!!****Chemisty_netIntRosen
!!
!!
!!
!! NAME
!!  Chemistry_netIntRosen
!!
!! SYNOPSIS
!!	Integrates our ODEs forward in time
!!
!! ARGUMENTS
!!	
!!  subroutine Chemistry_netIntRosen( real(INOUT)   :: y(:),
!!				      real(INOUT)   :: dydx(:),
!!				      int(IN)	    :: n,
!!				      real(IN)      :: htry, 
!!				      real(IN)      :: eps,
!!				      real(IN)      :: yscal(:),
!!				      real(OUT)     :: hdid,
!!			  	      real(OUT)     :: hnext,
!!				      PRO(IN)	    :: derivs,
!!				      PRO(IN)	    :: jakob,
!!				      real(OUT)	    :: jcounts)
!!  
!!
!!*****

subroutine Chemistry_netIntRosen(y,dydx,n,x,htry,eps,yscal,hdid,hnext, &
                                derivs,jakob,jcounts,den,temp)
 
use Chemistry_data
use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
use Logfile_interface, ONLY : Logfile_stampMessage
implicit none

#include "constants.h"
#include "Flash.h"

external derivs, jakob
integer, intent(IN)    :: n
real, intent(IN)      :: yscal(n),dydx(n),htry,eps, den, temp
real, intent(INOUT)   :: x, y(n)
real, intent(OUT)   :: hdid, hnext, jcounts

integer, parameter          :: nmax = 40, maxtry=1000
real, save      :: errmax, h, xsav, dysav(nmax), err(nmax), g1(nmax), &
                             &  g2(nmax),g3(nmax),g4(nmax),ysav(nmax),xx
real, parameter :: safety = 0.9e0, grow=1.5e0, pgrow=-0.25e0, &
                             &  shrnk=0.5e0,pshrnk=-1.0e0/3.0e0,errcon=0.1296e0

integer, save               :: i,j,jtry


!! shampine parameter set

real, parameter :: gamf = 1.0e0/2.0e0,  a21 = 2.0e0, &
                             & a31 = 48.0e0/25.0e0, a32 = 6.0e0/25.0e0, &
                             &  c21 = -8.0e0, c31 = 372.0e0/25.0e0, &
                             &  c32 = 12.e0/5.0e0, c41 = -112.0e0/125.0e0, &
                             &  c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0, &
                             &  b1 = 19.0e0/9.0e0, b2 = 1.0e0/2.0e0, &
                             &  b3 = 25.0e0/108.0e0, b4 = 125.0e0/108.0e0, &
                             &  e1 = 17.0e0/54.0e0, e2 = 7.0e0/36.0e0, &
                             &  e3 = 0.0e0, e4 = 125.0e0/108.0e0, &
                             &  c1x = 1.0e0/2.0e0, c2x = -3.0e0/2.0e0, &
                             &  c3x = 121.0e0/50.0e0, c4x = 29.0e0/250.0e0, &
                             &  a2x = 1.0e0, a3x = 3.0e0/5.0e0

integer, parameter       :: nmaxp1 = nmax + 1
real, save   :: dfdy(nmax,nmax), dmat(nmax,nmax)
real, save   :: av(nmax,nmaxp1), bb(nmax), yysav(nmax)
real         :: d
integer, save :: indx(nmax)


!!Here we go

!!Store initial values


!print *, 'IN CHEM_NETINTROSEN'
!print *, 'htry: ', htry


xsav = x
do i=1,n
   if(y(i) .lt. 0.0e0) then
      ysav(i) = 1.0e-30
   else
    ysav(i) = y(i)
   endif
    yysav(i) = y(i)
    dysav(i) = dydx(i)
enddo
    !! get the dense jacobian

    call jakob(xsav,ysav,dfdy,n,nmax)

    h = htry

    do jtry=1,maxtry

    !! Form the matrix

       xx = 1.0e0/(gamf*h)

       do j=1,n
         do i=1,n
            dmat(i,j) = -dfdy(i,j)
         enddo
       enddo

       do i=1,n
           dmat(i,i) = xx + dmat(i,i)
       enddo

    !! setup and solve the right hand side for g1
       do i=1,n
          g1(i) = dysav(i)
       enddo

       do j=1,n
          do i=1,n
            av(i,j) = dmat(i,j)
          enddo
       enddo

       do i=1,n
          bb(i) = g1(i)
          av(i,n+1) = g1(i)
       enddo
       !call ludcmp(av,n,nmax,indx,d)
       !call lubksb(av,n,nmax,indx,bb)
       call leqs(av,bb,n-1,nmax)

       do i=1,n
          g1(i) = bb(i)
       enddo
       g1(iELEC) = g1(iHP) + g1(iHEP) +g1(iCP) + g1(iOP) + g1(iO2P) + g1(iOHP) + g1(iCOP) + g1(iCHP) &
       		 & + g1(iCH2P) + g1(iHCOP) + g1(iHOCP) + g1(iH2OP) + g1(iH3OP) + g1(iCH3P)  + g1(iH2P) &
		 & - g1(iHM) - g1(iCM) - g1(iOM)

       do i=1,n
          y(i) = ysav(i) + a21 * g1(i)
       enddo
       x = xsav + a2x * h

       call derivs(x,y,dydx,den)

!! setup and solve the right hand side for g2
       do i=1,n
          g2(i) = dydx(i) + c21 * g1(i)/h
       enddo

       do j=1,n
          do i=1,n
             av(i,j) = dmat(i,j)
          enddo
       enddo

       do i=1,n
          bb(i) = g2(i)
          av(i,n+1) = g2(i)
       enddo
       !call ludcmp(av,n,nmax,indx,d)
       !call lubksb(av,n,nmax,indx,bb)
       call leqs(av,bb,n-1,nmax)

       do i=1,n
          g2(i) = bb(i)
       enddo
         g2(iELEC) = g2(iHP) + g2(iHEP) +g2(iCP) + g2(iOP) + g2(iO2P) + g2(iOHP) + g2(iCOP) + g2(iCHP) &
       		 & + g2(iCH2P) + g2(iHCOP) + g2(iHOCP) + g2(iH2OP) + g2(iH3OP) + g2(iCH3P)  + g2(iH2P) &
		 & - g2(iHM) - g2(iCM) - g2(iOM)

!! Compute intermediate values
      do i=1,n
         y(i) = ysav(i) + a31*g1(i) + a32*g2(i)
      enddo

      x = xsav + a3x*h

      call derivs(x,y,dydx,den)

!! setup and solve the 3rd RHS
      do i=1,n
         g3(i) = dydx(i) + (c31*g1(i) + c32*g2(i))/h
      enddo

      do j=1,n
         do i=1,n
            av(i,j) = dmat(i,j)
         enddo
      enddo

      do i=1,n 
          bb(i) = g3(i)
          av(i,n+1) = g3(i)
      enddo
      !call ludcmp(av,n,nmax,indx,d)
      !call lubksb(av,n,nmax,indx,bb)
      call leqs(av,bb,n-1,nmax)

      do i=1,n
         g3(i) = bb(i)
      enddo
         g3(iELEC) = g3(iHP) + g3(iHEP) +g3(iCP) + g3(iOP) + g3(iO2P) + g3(iOHP) + g3(iCOP) + g3(iCHP) &
       		 & + g3(iCH2P) + g3(iHCOP) + g3(iHOCP) + g3(iH2OP) + g3(iH3OP) + g3(iCH3P)  + g3(iH2P) &
		 & - g3(iHM) - g3(iCM) - g3(iOM)

!! setup and solve 4th RHS
      do i=1,n
         g4(i) = dydx(i) + (c41*g1(i) + c42*g2(i) + c43*g3(i))/h
      enddo

      do j=1,n
         do i=1,n
            av(i,j) = dmat(i,j)
         enddo
      enddo

      do i=1,n
         bb(i) = g4(i)
         av(i,n+1) = g4(i)
      enddo
      !call ludcmp(av,n,nmax,indx,d)
      !call lubksb(av,n,nmax,indx,bb)
      call leqs(av,bb,n-1,nmax)

      do i=1,n
         g4(i) = bb(i)
      enddo
         g4(iELEC) = g4(iHP) + g4(iHEP) +g4(iCP) + g4(iOP) + g4(iO2P) + g4(iOHP) + g4(iCOP) + g4(iCHP) &
       		 & + g4(iCH2P) + g4(iHCOP) + g4(iHOCP) + g4(iH2OP) + g4(iH3OP) + g4(iCH3P)  + g4(iH2P) &
		 & - g4(iHM) - g4(iCM) - g4(iOM)


!! compute the 3rd and 4th order estimates for y
      do i=1,n
         y(i) = ysav(i) + b1*g1(i) + b2*g2(i) + b3*g3(i) + b4*g4(i)
	 !write(*,'(6(A,1pe10.2))') 'ysav: ', ysav(i) , ' g1: ', g1(i) , ' g2: ', g2(i), ' g3: ', g3(i), ' g4: ', g4(i), ' yn: ', y(i)
         err(i) = e1*g1(i) + e2*g2(i) + e3*g3(i) + e4*g4(i)
      enddo
	 y(iELEC) = y(iHP) + y(iHEP) +y(iCP) + y(iOP) + y(iO2P) + y(iOHP) + y(iCOP) + y(iCHP) &
       		 & + y(iCH2P) + y(iHCOP) + y(iHOCP) + y(iH2OP) + y(iH3OP) + y(iCH3P)  + y(iH2P) &
		 & - y(iHM) - y(iCM) - y(iOM)


      x = xsav + h

     if (x .eq. xsav) then
         write(*,*) 'PROBLEM IN CHEMISTRY_NETINTROSEN'
         write(*,*) 'STEP SIZE NOT SIGNIFICANT IN CHEMISTRY_NETINTROSEN'
         call abort()
     endif

     if(y(1) .ne. y(1)) then  !!This should be a NAN
       write(*,*), ' ERROR IN CHEM_NETINTROSEN'
       do i=1,n
          write(*,'(A,I,8(A,1pe15.6))'), 'I: ', i, ' y(i): ', y(i), ' ys(i): ', ysav(i), ' den: ', den, ' temp: ', temp, ' g1: ', g1(i), ' g2(i): ', g2(i), ' g3(i): ', g3(i), ' g4(i): ', g4(i)
       enddo
       write(*,*), ' PRINTING JACOB '
       do i=1,n
         do j=1,n
           write(*,'(I,I,1pe15.6)'), i,j,dfdy(i,j)
         enddo
       enddo
       write(*,*), ' PRINTING RATES '
       do i=1,ireaction
         write(*,'(I,1pe15.6)'), i, ratraw(i)
       enddo
     endif


!! determine the scaled accuracy
     errmax = 0.0e0
     do i=1,n
        errmax = max(errmax,abs(err(i)/yscal(i)))
     !   write(*,'(a,i4)') 'i: ', i
     !   write(*,'(A,1pe10.2,A,1pe10.2,A,1pe10.2)') ' errmax: ', errmax, '  err_i', err(i), '  yscal_i: ', yscal(i)
     !   write(*,'(A,1pe10.2)') 'err: ', abs(err(i)/yscal(i))
     enddo
     errmax = errmax/eps
!     write(*,'(A,1pe10.2)') 'errmax: ', errmax

!! if the step succeded, compute the size of the next step and return
     if(errmax .le. 1.0) then
        
        hdid = h
        if(errmax .gt. errcon) then
           hnext = safety * h * errmax**pgrow
        else
           hnext = grow * h
        endif
        jcounts = jcounts + jtry
!        print *, 'STEP OK'
!        print *, 'h: ', h, ' hnext: ', hnext
       !if(jcounts .gt. 100000e0) then
       !    write(*,*) 'ERROR IN NETINTROSEN'
       !    write(*,*) 'JCOUNTS > 1000000'
       !    write(*,*) 'htry',hnext
          ! do i=1,NSPECIES
            !  write(*,*) "y: ", y(i)
          ! enddo
       !    call abort()
       ! endif

        return
   else !!Step did not work, cut stepsize and try again

     !if(jtry .gt. 0.99*maxtry) then
     !  write(*,*) 'Getting to be a problem'
     !  do i=1,NSPECIES
     !     write(*,*) "Y: ", y(i) ,' ysav:' , ysav(i)
     !  enddo
     !endif

     hnext = safety * h * errmax**pshrnk
     h = sign(max(abs(hnext),shrnk*abs(h)),h)
!     print *, 'STEP FAIL'
!     print *, 'hnext: ', h
   endif

enddo

!!If we get here, we took too many tries

 ! write(*,*) 'ERROR IN CHEMISTRY_NETINTROSEN'
 ! write(*,*) 'EXCEEDED MAXTRY ! '
  call Driver_abortFLASH('Exceeded Maxtry')

end subroutine Chemistry_netIntRosen

  subroutine leqs(a,b,n,np)
  implicit none
  save


  integer  n,np,n1,i,j,k,l,imax,jj
  real     a(np,np),b(np),r,c
 
!  print *, 'inside leqs'  

  n1=n-1
  do i=1,n
     r = abs(a(i,1))
     do j=2,n
        c = abs(a(i,j))
        if (r .lt. c) r=c
     enddo

     do j=1,n
        a(i,j) = a(i,j)/r
     enddo
     b(i) = b(i)/r
  enddo


  do j=1,n1
     l = j+1
     do i=l,n
        r = -a(i,j)
        if (r .eq. 0.0) goto 50
        r = r/a(j,j)
        do k=l,n
           a(i,k) = a(i,k) + r*a(j,k)
        enddo
        b(i) = b(i) + r*b(j)
50      continue
     enddo
  enddo


  b(n) = b(n)/a(n,n)
  do l=1,n1
     i = n-l
     r = 0.0e0
     imax = i+1
     do j = imax, n
        jj = i+n+1-j
         r = r+a(i,jj)*b(jj)
     enddo
     b(i) = (b(i)-r)/a(i,i)
  enddo
  return
  end
  
  
  
    subroutine ludcmp(a,n,np,indx,d)
    INTEGER n,np,indx(n),NMAX
    REAL d,a(np,np),TINY
    PARAMETER (NMAX=500, TINY=1.0e-20)

    INTEGER i,imax,j,k
    REAL aamax,dum,sum,vv(NMAX)
    d=1.
    do i=1,n
       aamax = 0.
       do j=1,n
          if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
       enddo
       if (aamax .eq. 0.) stop 'singular matrix in ludcmp'
       vv(i) = 1./aamax
    enddo 

    do j=1,n
         do i=1,j-1
            sum = a(i,j)
            do k=1,i-1
               sum = sum-a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
         enddo
         aamax = 0.
         do i=j,n
             sum=a(i,j)
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j) = sum
             dum = vv(i)*abs(sum)
             if (dum .ge. aamax) then
                imax = i
                aamax = dum
             endif
         enddo
         if (j .ne. imax) then
             do k=1,n
                  dum = a(imax,k)
                  a(imax,k) = a(j,k)
                  a(j,k) = dum
             enddo
             d=-d
             vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j) .eq. 0.) a(j,j) = TINY
         if (j .ne. n) then
            dum = 1./a(j,j)
            do i=j+1,n
                 a(i,j) = a(i,j)*dum
            enddo
         endif
     enddo
     return 
     END






  subroutine lubksb(a,n,np,indx,b)
  
  INTEGER n,np,indx(n)
  REAL a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum

  ii=0

  do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if (ii .ne. 0) then
          do j = ii,i-1
               sum = sum-a(i,j)*b(j)
          enddo
       else if (sum .ne. 0.) then
           ii=i
       endif
       b(i) = sum
  enddo
 

  do i=n,1,-1
       sum = b(i)
          do j=i+1,n
             sum = sum - a(i,j)*b(j)
          enddo
       b(i) = sum/a(i,i)
  enddo
  return
  END  
