    subroutine st_getsold(xoverexp,overexp,sold)
    implicit none
    real, intent(in)  :: overexp
    real, intent(in)  :: xoverexp
    real, intent(out) :: sold
    real    :: x0,x1
    real    :: numerator,denominator
    real    :: deltax
    integer :: count


!   Bypass everything an use 1/rho
    if(1.gt.0) then
      x0= -alog(overexp)
    else
!   if xoverexp is to big then return 1
    if (xoverexp .gt. 0.367879) then
     print *,'xoverexp is too big',xoverexp
     x0 = 1.
    else
!     otherwise use overexp to get a first guess
      x0 = -alog(overexp)
      deltax = 1
      count = 0
!     use newtons method to get an answer.  But theres a slope of 0 near 1
      do while(abs(deltax).gt.0.0002) 
!       solution for x0 minus the desired answer
        numerator   =  x0*exp(-x0)-xoverexp
!       if you are less than 1 increasing x0 increases the solution
        denominator = (1-x0)*exp(-x0)
        if(abs(denominator).lt.0.0001) denominator = 0.0001
        deltax = numerator/denominator

!       if you are close to x0=1 and might overshoot then tamp it down 
!       don't go more than 50% of the way to 1
        if(abs(deltax).gt.(0.8*abs(1-x0))) deltax = 0.8 * abs(1.-x0) * deltax/abs(deltax)
        x0 = x0-deltax
!        print *,'x0,answer,deltax',x0,x0*exp(-x0),deltax
        count = count+1
        if (count.gt.50) then
          print *,'There is a problem',xoverexp,x0
          deltax = 0
        endif
      enddo
    endif
    endif
    sold = x0
    return
    end 


