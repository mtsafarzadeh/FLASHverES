!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/hy_uhd_eigenVector
!!
!! NAME
!!
!!  hy_uhd_eigenVector
!!
!! SYNOPSIS
!!
!!  hy_uhd_eigenVector( real (OUT)         :: LeftEigvec(HY_WAVENUM,HY_VARINUM),
!!                      real (OUT)         :: RightEigvec(HY_VARINUM,HY_WAVENUM)
!!                      real (IN)          :: V(HY_VARINUM2),
!!                      integer(IN)        :: dir,
!!                      logical(IN)        :: cons,
!!                      real (IN)          :: C_fast,
!!                      real (IN),optional :: C_alfn,
!!                      real (IN),optional :: C_slow,
!!                      real (IN),optional :: A_f,
!!                      real (IN),optional :: A_s,
!!                      real (IN),optional :: B_beta(MDIM) )
!!
!! DESCRIPTION
!!
!!  This routine calculates MHD/Hydro eigenvectors in either primitive form (used in
!!  Riemann solver) or conservative form (used in conservative updates).
!!
!!
!! ARGUMENTS
!!
!!  LeftEigvec  - Left eigenvectors
!!  RightEigvec - Right eigenvectors
!!  V           - Primitive variables + gammas:
!!                (dens,velx,vely,velz,pres,(magx,magy,magz),gamc,game)
!!  dir         - x,y,z direction
!!  cons        - A logical switch to choose either primitive or conservative eigenvector
!!  C_fast      - Fast magnetoacoustic speed for MHD/Sound speed for Hydro
!!  C_alfn      - Alfven speed (needed for MHD only)
!!  C_slow      - Slow magnetoacoustic speed (needed for MHD only)
!!  A_f         - Normalization coefficient (needed for MHD only)
!!  A_s         - Normalization coefficient (needed for MHD only)
!!  B_beta      - Alfven velcoities in transversal direction (needed for MHD only)
!!
!!***

!#define DEBUG_HY_EIGEN

Subroutine hy_uhd_eigenVector&
     (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta)

  use Hydro_data,           ONLY : hy_meshMe
  use Logfile_interface,    ONLY : Logfile_open,Logfile_close

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! Arguments type declaration ---------------------------------------
  real, dimension(HY_WAVENUM,HY_VARINUM), intent(OUT) :: LeftEigvec
  real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
  real, dimension(HY_VARINUM2), intent(IN) :: V
  integer, intent(IN) :: dir
  logical, intent(IN) :: cons
  real, intent(IN) :: C_fast
  real, intent(IN),optional :: C_alfn,C_slow,A_f,A_s
  real, dimension(MDIM), intent(IN),optional  :: B_beta
  !! ------------------------------------------------------------------

  integer :: ii,jj
  real :: sqrtd,k,dinv,a2inv,a,a2,cf2,u2
  real :: test
  real, dimension(HY_VARINUM,HY_VARINUM) :: JacQ,JacQi
  real, dimension(HY_WAVENUM,HY_VARINUM) :: LeftEigvecTemp
  real, dimension(HY_VARINUM,HY_WAVENUM) :: RightEigvecTemp
  integer :: logUnit
  logical :: logUnitLocal=.true.

  ! parameters
  sqrtd =sqrt(V(HY_DENS))
  dinv  =1./V(HY_DENS)
  k=1.-V(HY_GAME) !! k=1.-game

  u2=dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))

  a = C_fast
  a2=a*a
  a2inv=1./a2

  ! initialize with zeros
  LeftEigvec = 0.
  RightEigvec = 0.

  !! -----------------------------------------------!
  !! Construct left eigenvectors                    !
  !! -----------------------------------------------!

  ! left/right - goint fast/slow waves
  select case(dir)
  case(DIR_X)
     LeftEigvec(HY_FASTLEFT,HY_VELX)= -C_fast
     LeftEigvec(HY_FASTRGHT,HY_VELX)=  C_fast

     LeftEigvec(HY_SLOWLEFT,HY_VELY)= -V(HY_DENS)
     LeftEigvec(HY_SLOWRGHT,HY_VELZ)=  V(HY_DENS)

  case(DIR_Y)
     LeftEigvec(HY_FASTLEFT,HY_VELY)= -C_fast
     LeftEigvec(HY_FASTRGHT,HY_VELY)=  C_fast

     LeftEigvec(HY_SLOWLEFT,HY_VELX)=  V(HY_DENS)
     LeftEigvec(HY_SLOWRGHT,HY_VELZ)= -V(HY_DENS)

  case(DIR_Z)
     LeftEigvec(HY_FASTLEFT,HY_VELZ)= -C_fast
     LeftEigvec(HY_FASTRGHT,HY_VELZ)=  C_fast

     LeftEigvec(HY_SLOWLEFT,HY_VELX)= -V(HY_DENS)
     LeftEigvec(HY_SLOWRGHT,HY_VELY)=  V(HY_DENS)
  end select
  
  LeftEigvec(HY_FASTLEFT, HY_PRES) =  dinv
  LeftEigvec(HY_FASTRGHT, HY_PRES) =  dinv


  ! scale fast waves with 1/2*a^2
  LeftEigvec(HY_FASTLEFT,:) = 0.5*a2inv*LeftEigvec(HY_FASTLEFT,:)
  LeftEigvec(HY_FASTRGHT,:) = 0.5*a2inv*LeftEigvec(HY_FASTRGHT,:)

  ! entropy wave
  LeftEigvec(HY_ENTROPY,HY_DENS) = 1.
  LeftEigvec(HY_ENTROPY,HY_PRES) =-a2inv

  !! -----------------------------------------------!
  !! Construct right eigenvectors                   !
  !! -----------------------------------------------!

  ! left/right - goint fast waves
  RightEigvec(HY_DENS,HY_FASTLEFT) = V(HY_DENS)
  RightEigvec(HY_DENS,HY_FASTRGHT) = V(HY_DENS)

  select case(dir)
  case(DIR_X)
     RightEigvec(HY_VELX,HY_FASTLEFT)= -C_fast
     RightEigvec(HY_VELX,HY_FASTRGHT)=  C_fast

     RightEigvec(HY_VELY,HY_SLOWLEFT)= -dinv
     RightEigvec(HY_VELZ,HY_SLOWRGHT)=  dinv

  case(DIR_Y)
     RightEigvec(HY_VELY,HY_FASTLEFT)= -C_fast
     RightEigvec(HY_VELY,HY_FASTRGHT)=  C_fast

     RightEigvec(HY_VELX,HY_SLOWLEFT)=  dinv
     RightEigvec(HY_VELZ,HY_SLOWRGHT)= -dinv

  case(DIR_Z)
     RightEigvec(HY_VELZ,HY_FASTLEFT)= -C_fast
     RightEigvec(HY_VELZ,HY_FASTRGHT)=  C_fast

     RightEigvec(HY_VELX,HY_SLOWLEFT)= -dinv
     RightEigvec(HY_VELY,HY_SLOWRGHT)=  dinv
  end select
  RightEigvec(HY_PRES,HY_FASTLEFT) =  V(HY_DENS)*a2
  RightEigvec(HY_PRES,HY_FASTRGHT) =  V(HY_DENS)*a2


  ! entropy wave
  RightEigvec(HY_DENS,HY_ENTROPY) = 1.


#ifdef DEBUG_HY_EIGEN
  do ii = 1,HY_WAVENUM
     test=dot_product(LeftEigvec(ii,:), RightEigvec(:,ii))
     if (abs(test-1.) > 1.e-4) then
        call Logfile_open(logUnit,logUnitLocal)
        write(logUnit,*)'dot product is not unity: test=',test,ii,dir   !DEBUG
        call Logfile_close(logUnitLocal)
     endif
  enddo

  do ii = 1,HY_WAVENUM
     if (ii < HY_WAVENUM) then
        test=dot_product(LeftEigvec(ii,:), RightEigvec(:,ii+1))
        if (abs(test) > 1.e-4) then
           call Logfile_open(logUnit,logUnitLocal)
           write(logUnit,*)'dot product is not zero: test=',test,ii,dir   !DEBUG
           call Logfile_close(logUnitLocal)
        endif
     else
        test=dot_product(LeftEigvec(ii,:), RightEigvec(:,ii-1))
        if (abs(test) > 1.e-4) then
           call Logfile_open(logUnit,logUnitLocal)
           write(logUnit,*)'dot product is not zero: test=',test,ii,dir   !DEBUG
           call Logfile_close(logUnitLocal)
        endif
     endif
  enddo
#endif


  !! Convert to eigenvectors in conservative forms
  if (cons) then
     ! Initialize Jacobian matrix Q and its inverse Q^(-1)
     JacQ = 0.
     JacQi= 0.

     ! Define Q
     JacQ (HY_DENS,HY_DENS) = 1.
     JacQ (HY_VELX,HY_DENS) = V(HY_VELX)
     JacQ (HY_VELX,HY_VELX) = V(HY_DENS)
     JacQ (HY_VELY,HY_DENS) = V(HY_VELY)
     JacQ (HY_VELY,HY_VELY) = V(HY_DENS)
     JacQ (HY_VELZ,HY_DENS) = V(HY_VELZ)
     JacQ (HY_VELZ,HY_VELZ) = V(HY_DENS)
     JacQ (HY_PRES,HY_DENS:HY_PRES)=(/0.5*u2,V(HY_DENS)*V(HY_VELX),&
                                             V(HY_DENS)*V(HY_VELY),&
                                             V(HY_DENS)*V(HY_VELZ),&
                                             -1./k/)

     ! Define Q^(-1)
     JacQi(HY_DENS,HY_DENS) = 1.
     JacQi(HY_VELX,HY_DENS) = -V(HY_VELX)*dinv
     JacQi(HY_VELX,HY_VELX) =  dinv
     JacQi(HY_VELY,HY_DENS) = -V(HY_VELY)*dinv
     JacQi(HY_VELY,HY_VELY) =  dinv
     JacQi(HY_VELZ,HY_DENS) = -V(HY_VELZ)*dinv
     JacQi(HY_VELZ,HY_VELZ) =  dinv
     JacQi(HY_PRES,HY_DENS:HY_PRES)=(/-.5*u2,V(HY_VELX),&
                                             V(HY_VELY),&
                                             V(HY_VELZ),&
                                             -1./)*k

     !! left eigenvectors
     do ii=1,HY_WAVENUM
        do jj=1,HY_VARINUM
           LeftEigvecTemp(ii,jj) = dot_product(LeftEigvec(ii,1:HY_VARINUM),JacQi(1:HY_VARINUM,jj))
        enddo
     enddo
     LeftEigvec = LeftEigvecTemp

     !! right eigenvectors
     do ii=1,HY_VARINUM
        do jj=1,HY_WAVENUM
           RightEigvecTemp(ii,jj) = dot_product(JacQ(ii,1:HY_VARINUM),RightEigvec(1:HY_VARINUM,jj))
        enddo
     enddo
     RightEigvec = RightEigvecTemp

#ifdef DEBUG_HY_EIGEN
  do ii = 1,HY_WAVENUM
     test=dot_product(LeftEigvec(ii,:), RightEigvec(:,ii))
     if (abs(test-1.) > 1.e-4) then
        call Logfile_open(logUnit,logUnitLocal)
        write(logUnit,*)'conserve dot product is not unity: test=',test,ii,dir   !DEBUG
        call Logfile_close(logUnitLocal)
     endif
  enddo

  do ii = 1,HY_WAVENUM
     if (ii < HY_WAVENUM) then
        test=dot_product(LeftEigvec(ii,:), RightEigvec(:,ii+1))
        if (abs(test) > 1.e-4) then
           call Logfile_open(logUnit,logUnitLocal)
           write(logUnit,*)'conserve dot product is not zero: test=',test,ii,dir   !DEBUG
           call Logfile_close(logUnitLocal)
        endif
     else
        test=dot_product(LeftEigvec(ii,:), RightEigvec(:,ii-1))
        if (abs(test) > 1.e-4) then
           call Logfile_open(logUnit,logUnitLocal)
           write(logUnit,*)'conserve dot product is not zero: test=',test,ii,dir   !DEBUG
           call Logfile_close(logUnitLocal)
        endif
     endif
  enddo
#endif



  endif


End Subroutine hy_uhd_eigenVector
