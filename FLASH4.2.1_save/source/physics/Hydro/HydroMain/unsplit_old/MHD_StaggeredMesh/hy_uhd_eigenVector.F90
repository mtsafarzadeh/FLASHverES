!!****if* source/physics/Hydro/HydroMain/unsplit_old/MHD_StaggeredMesh/hy_uhd_eigenVector
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
  use Driver_interface,     ONLY : Driver_abortFlash

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
  real :: sqrtd,sqrtdinv,signbN,k,dinv,a2inv,a,a2,u2
  real :: a_f1,a_f2,a_f3,a_s1,a_s2,a_s3,test
  real, dimension(HY_VARINUM,HY_VARINUM) :: JacQ,JacQi
  real, dimension(HY_WAVENUM,HY_VARINUM) :: LeftEigvecTemp
  real, dimension(HY_VARINUM,HY_WAVENUM) :: RightEigvecTemp
  integer :: logUnit
  logical :: logUnitLocal=.true.

  ! parameters
  sqrtd =sqrt(V(HY_DENS))
  sqrtdinv=1./sqrtd
  signbN = sign(1.,C_alfn)
  dinv  =1./V(HY_DENS)
  k=1.-V(HY_GAME) !! k=1.-game

  a2=V(HY_GAMC)*V(HY_PRES)*dinv

  if (a2 .le. 0.) then
     call Driver_abortFlash&
          ("[hy_uhd_eigenVector]: Imaginary sound speed has obtained! "//&
           "Please try other (more diffusive) slope limiter, flux, order, cfl, etc.")
  else
     a=sqrt(a2)
  endif


  a2inv=1./a2
  u2=dot_product(V(HY_VELX:HY_VELZ),V(HY_VELX:HY_VELZ))

  a_f1 = A_f*C_fast*signbN
  a_f2 = A_f*a/sqrtd
  a_f3 = A_f*a*sqrtd

  a_s1 = A_s*C_slow*signbN
  a_s2 = A_s*a/sqrtd
  a_s3 = A_s*a*sqrtd

  ! initialize with zeros
  LeftEigvec  = 0.
  RightEigvec = 0.


  !! -----------------------------------------------!
  !! Construct left eigenvectors                    !
  !! -----------------------------------------------!

  ! left/right - goint fast waves
  LeftEigvec(HY_FASTLEFT,HY_VELX:HY_VELZ) =  a_s1*B_beta
  LeftEigvec(HY_FASTRGHT,HY_VELX:HY_VELZ) = -LeftEigvec(HY_FASTLEFT,HY_VELX:HY_VELZ)

  select case(dir)
  case(DIR_X)
     LeftEigvec(HY_FASTLEFT,HY_VELX)= -A_f*C_fast
     LeftEigvec(HY_FASTRGHT,HY_VELX)= -LeftEigvec(HY_FASTLEFT,HY_VELX)
  case(DIR_Y)
     LeftEigvec(HY_FASTLEFT,HY_VELY)= -A_f*C_fast
     LeftEigvec(HY_FASTRGHT,HY_VELY)= -LeftEigvec(HY_FASTLEFT,HY_VELY)
  case(DIR_Z)
     LeftEigvec(HY_FASTLEFT,HY_VELZ)= -A_f*C_fast
     LeftEigvec(HY_FASTRGHT,HY_VELZ)= -LeftEigvec(HY_FASTLEFT,HY_VELZ)
  end select
  
  LeftEigvec(HY_FASTLEFT,HY_MAGX:HY_MAGZ) =  a_s2*B_beta  
  LeftEigvec(HY_FASTRGHT,HY_MAGX:HY_MAGZ) =  LeftEigvec(HY_FASTLEFT,HY_MAGX:HY_MAGZ)
  LeftEigvec(HY_FASTLEFT,HY_PRES) =  A_f*dinv
  LeftEigvec(HY_FASTRGHT,HY_PRES) =  LeftEigvec(HY_FASTLEFT,HY_PRES)

  ! left/right - goint slow waves
  LeftEigvec(HY_SLOWLEFT,HY_VELX:HY_VELZ) = -a_f1*B_beta
  LeftEigvec(HY_SLOWRGHT,HY_VELX:HY_VELZ) = -LeftEigvec(HY_SLOWLEFT,HY_VELX:HY_VELZ)
  select case(dir)
  case(DIR_X)
     LeftEigvec(HY_SLOWLEFT,HY_VELX)= -A_s*C_slow
     LeftEigvec(HY_SLOWRGHT,HY_VELX)= -LeftEigvec(HY_SLOWLEFT,HY_VELX)
  case(DIR_Y)
     LeftEigvec(HY_SLOWLEFT,HY_VELY)= -A_s*C_slow
     LeftEigvec(HY_SLOWRGHT,HY_VELY)= -LeftEigvec(HY_SLOWLEFT,HY_VELY)
  case(DIR_Z)
     LeftEigvec(HY_SLOWLEFT,HY_VELZ)= -A_s*C_slow
     LeftEigvec(HY_SLOWRGHT,HY_VELZ)= -LeftEigvec(HY_SLOWLEFT,HY_VELZ)
  end select

  LeftEigvec(HY_SLOWLEFT,HY_MAGX:HY_MAGZ) = -a_f2*B_beta
  LeftEigvec(HY_SLOWRGHT,HY_MAGX:HY_MAGZ) =  LeftEigvec(HY_SLOWLEFT,HY_MAGX:HY_MAGZ)
  LeftEigvec(HY_SLOWLEFT,HY_PRES) =  A_s*dinv
  LeftEigvec(HY_SLOWRGHT,HY_PRES) =  LeftEigvec(HY_SLOWLEFT,HY_PRES)

  ! scale fast and slow waves with 0.5/a^2
  LeftEigvec = 0.5*a2inv*LeftEigvec


  ! left/right -going alfven waves
  select case(dir)
  case(DIR_X)
     LeftEigvec(HY_ALFNLEFT,HY_VELY)=-B_beta(DIR_Z)
     LeftEigvec(HY_ALFNLEFT,HY_VELZ)= B_beta(DIR_Y)
     LeftEigvec(HY_ALFNRGHT,HY_VELY)= B_beta(DIR_Z)
     LeftEigvec(HY_ALFNRGHT,HY_VELZ)=-B_beta(DIR_Y)

     LeftEigvec(HY_ALFNLEFT,HY_MAGY)=-B_beta(DIR_Z)*sqrtdinv
     LeftEigvec(HY_ALFNLEFT,HY_MAGZ)= B_beta(DIR_Y)*sqrtdinv
     LeftEigvec(HY_ALFNRGHT,HY_MAGY)=-B_beta(DIR_Z)*sqrtdinv
     LeftEigvec(HY_ALFNRGHT,HY_MAGZ)= B_beta(DIR_Y)*sqrtdinv
  case(DIR_Y)
     LeftEigvec(HY_ALFNLEFT,HY_VELX)= B_beta(DIR_Z)
     LeftEigvec(HY_ALFNLEFT,HY_VELZ)=-B_beta(DIR_X)
     LeftEigvec(HY_ALFNRGHT,HY_VELX)=-B_beta(DIR_Z)
     LeftEigvec(HY_ALFNRGHT,HY_VELZ)= B_beta(DIR_X)

     LeftEigvec(HY_ALFNLEFT,HY_MAGX)= B_beta(DIR_Z)*sqrtdinv
     LeftEigvec(HY_ALFNLEFT,HY_MAGZ)=-B_beta(DIR_X)*sqrtdinv
     LeftEigvec(HY_ALFNRGHT,HY_MAGX)= B_beta(DIR_Z)*sqrtdinv
     LeftEigvec(HY_ALFNRGHT,HY_MAGZ)=-B_beta(DIR_X)*sqrtdinv
  case(DIR_Z)
     LeftEigvec(HY_ALFNLEFT,HY_VELX)=-B_beta(DIR_Y)
     LeftEigvec(HY_ALFNLEFT,HY_VELY)= B_beta(DIR_X)
     LeftEigvec(HY_ALFNRGHT,HY_VELX)= B_beta(DIR_Y)
     LeftEigvec(HY_ALFNRGHT,HY_VELY)=-B_beta(DIR_X)

     LeftEigvec(HY_ALFNLEFT,HY_MAGX)=-B_beta(DIR_Y)*sqrtdinv
     LeftEigvec(HY_ALFNLEFT,HY_MAGY)= B_beta(DIR_X)*sqrtdinv
     LeftEigvec(HY_ALFNRGHT,HY_MAGX)=-B_beta(DIR_Y)*sqrtdinv
     LeftEigvec(HY_ALFNRGHT,HY_MAGY)= B_beta(DIR_X)*sqrtdinv
  end select

  ! scale alfven waves with 0.5
  LeftEigvec(HY_ALFNLEFT,:) = 0.5*LeftEigvec(HY_ALFNLEFT,:)
  LeftEigvec(HY_ALFNRGHT,:) = 0.5*LeftEigvec(HY_ALFNRGHT,:)


  ! entropy wave
  LeftEigvec(HY_ENTROPY,HY_DENS) = 1.
  LeftEigvec(HY_ENTROPY,HY_PRES) = -a2inv


  !! -----------------------------------------------!
  !! Construct right eigenvectors                   !
  !! -----------------------------------------------!

  ! left/right - goint fast waves
  RightEigvec(HY_DENS,HY_FASTLEFT)   =  A_f*V(HY_DENS)
  RightEigvec(HY_DENS,HY_FASTRGHT)   =  RightEigvec(HY_DENS,HY_FASTLEFT)
  RightEigvec(HY_VELX:HY_VELZ,HY_FASTLEFT) =  a_s1*B_beta
  RightEigvec(HY_VELX:HY_VELZ,HY_FASTRGHT) = -a_s1*B_beta
  select case(dir)
  case(DIR_X)
     RightEigvec(HY_VELX,HY_FASTLEFT)= -A_f*C_fast
     RightEigvec(HY_VELX,HY_FASTRGHT)= -RightEigvec(HY_VELX,HY_FASTLEFT)
  case(DIR_Y)
     RightEigvec(HY_VELY,HY_FASTLEFT)= -A_f*C_fast
     RightEigvec(HY_VELY,HY_FASTRGHT)= -RightEigvec(HY_VELY,HY_FASTLEFT)
  case(DIR_Z)
     RightEigvec(HY_VELZ,HY_FASTLEFT)= -A_f*C_fast
     RightEigvec(HY_VELZ,HY_FASTRGHT)= -RightEigvec(HY_VELZ,HY_FASTLEFT)
  end select
  RightEigvec(HY_MAGX:HY_MAGZ,HY_FASTLEFT) =  a_s3*B_beta  
  RightEigvec(HY_MAGX:HY_MAGZ,HY_FASTRGHT) =  RightEigvec(HY_MAGX:HY_MAGZ,HY_FASTLEFT)
  RightEigvec(HY_PRES,HY_FASTLEFT) =  A_f*V(HY_DENS)*a2
  RightEigvec(HY_PRES,HY_FASTRGHT) =  RightEigvec(HY_PRES,HY_FASTLEFT)


  ! left/right - goint slow waves
  RightEigvec(HY_DENS,HY_SLOWLEFT)   =  A_s*V(HY_DENS)
  RightEigvec(HY_DENS,HY_SLOWRGHT)   =  RightEigvec(HY_DENS,HY_SLOWLEFT)
  RightEigvec(HY_VELX:HY_VELZ,HY_SLOWLEFT) = -a_f1*B_beta
  RightEigvec(HY_VELX:HY_VELZ,HY_SLOWRGHT) = -RightEigvec(HY_VELX:HY_VELZ,HY_SLOWLEFT)
  select case(dir)
  case(DIR_X)
     RightEigvec(HY_VELX,HY_SLOWLEFT)= -A_s*C_slow
     RightEigvec(HY_VELX,HY_SLOWRGHT)= -RightEigvec(HY_VELX,HY_SLOWLEFT)
  case(DIR_Y)
     RightEigvec(HY_VELY,HY_SLOWLEFT)= -A_s*C_slow
     RightEigvec(HY_VELY,HY_SLOWRGHT)= -RightEigvec(HY_VELY,HY_SLOWLEFT)
  case(DIR_Z)
     RightEigvec(HY_VELZ,HY_SLOWLEFT)= -A_s*C_slow
     RightEigvec(HY_VELZ,HY_SLOWRGHT)= -RightEigvec(HY_VELZ,HY_SLOWLEFT)
  end select
  RightEigvec(HY_MAGX:HY_MAGZ,HY_SLOWLEFT) = -a_f3*B_beta
  RightEigvec(HY_MAGX:HY_MAGZ,HY_SLOWRGHT) =  RightEigvec(HY_MAGX:HY_MAGZ,HY_SLOWLEFT)
  RightEigvec(HY_PRES,HY_SLOWLEFT) =  A_s*V(HY_DENS)*a2
  RightEigvec(HY_PRES,HY_SLOWRGHT) =  RightEigvec(HY_PRES,HY_SLOWLEFT)


  ! left/right -going alfven waves
  select case(dir)
  case(DIR_X)
     RightEigvec(HY_VELY,HY_ALFNLEFT)=-B_beta(DIR_Z)
     RightEigvec(HY_VELZ,HY_ALFNLEFT)= B_beta(DIR_Y)
     RightEigvec(HY_VELY,HY_ALFNRGHT)= B_beta(DIR_Z)
     RightEigvec(HY_VELZ,HY_ALFNRGHT)=-B_beta(DIR_Y)

     RightEigvec(HY_MAGY,HY_ALFNLEFT)=-B_beta(DIR_Z)*sqrtd
     RightEigvec(HY_MAGZ,HY_ALFNLEFT)= B_beta(DIR_Y)*sqrtd
     RightEigvec(HY_MAGY,HY_ALFNRGHT)=-B_beta(DIR_Z)*sqrtd
     RightEigvec(HY_MAGZ,HY_ALFNRGHT)= B_beta(DIR_Y)*sqrtd
  case(DIR_Y)
     RightEigvec(HY_VELX,HY_ALFNLEFT)= B_beta(DIR_Z)
     RightEigvec(HY_VELZ,HY_ALFNLEFT)=-B_beta(DIR_X)
     RightEigvec(HY_VELX,HY_ALFNRGHT)=-B_beta(DIR_Z)
     RightEigvec(HY_VELZ,HY_ALFNRGHT)= B_beta(DIR_X)

     RightEigvec(HY_MAGX,HY_ALFNLEFT)= B_beta(DIR_Z)*sqrtd
     RightEigvec(HY_MAGZ,HY_ALFNLEFT)=-B_beta(DIR_X)*sqrtd
     RightEigvec(HY_MAGX,HY_ALFNRGHT)= B_beta(DIR_Z)*sqrtd
     RightEigvec(HY_MAGZ,HY_ALFNRGHT)=-B_beta(DIR_X)*sqrtd
  case(DIR_Z)
     RightEigvec(HY_VELX,HY_ALFNLEFT)=-B_beta(DIR_Y)
     RightEigvec(HY_VELY,HY_ALFNLEFT)= B_beta(DIR_X)
     RightEigvec(HY_VELX,HY_ALFNRGHT)= B_beta(DIR_Y)
     RightEigvec(HY_VELY,HY_ALFNRGHT)=-B_beta(DIR_X)

     RightEigvec(HY_MAGX,HY_ALFNLEFT)=-B_beta(DIR_Y)*sqrtd
     RightEigvec(HY_MAGY,HY_ALFNLEFT)= B_beta(DIR_X)*sqrtd
     RightEigvec(HY_MAGX,HY_ALFNRGHT)=-B_beta(DIR_Y)*sqrtd
     RightEigvec(HY_MAGY,HY_ALFNRGHT)= B_beta(DIR_X)*sqrtd
  end select

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
     JacQ (HY_MAGX,HY_MAGX) = 1.
     JacQ (HY_MAGY,HY_MAGY) = 1.
     JacQ (HY_MAGZ,HY_MAGZ) = 1.
     JacQ (HY_PRES,HY_DENS:HY_MAGZ)=(/0.5*u2,V(HY_DENS)*V(HY_VELX),&
                                             V(HY_DENS)*V(HY_VELY),&
                                             V(HY_DENS)*V(HY_VELZ),&
                                             -1./k,     &
                                             V(HY_MAGX),&
                                             V(HY_MAGY),&
                                             V(HY_MAGZ)/)

     ! Define Q^(-1)
     JacQi(HY_DENS,HY_DENS) = 1.
     JacQi(HY_VELX,HY_DENS) = -V(HY_VELX)*dinv
     JacQi(HY_VELX,HY_VELX) =  dinv
     JacQi(HY_VELY,HY_DENS) = -V(HY_VELY)*dinv
     JacQi(HY_VELY,HY_VELY) =  dinv
     JacQi(HY_VELZ,HY_DENS) = -V(HY_VELZ)*dinv
     JacQi(HY_VELZ,HY_VELZ) =  dinv
     JacQi(HY_MAGX,HY_MAGX) = 1.
     JacQi(HY_MAGY,HY_MAGY) = 1.
     JacQi(HY_MAGZ,HY_MAGZ) = 1.
     JacQi(HY_PRES,HY_DENS:HY_MAGZ)=(/-.5*u2,V(HY_VELX),&
                                             V(HY_VELY),&
                                             V(HY_VELZ),&
                                             -1.,       &
                                             V(HY_MAGX),&
                                             V(HY_MAGY),&
                                             V(HY_MAGZ)/)*k

     !! left eigenvectors
     do ii=1,HY_WAVENUM
        do jj=1,HY_VARINUM
           LeftEigvecTemp(ii,jj) = dot_product(LeftEigvec(ii,:),JacQi(:,jj))
        enddo
     enddo
     LeftEigvec = LeftEigvecTemp

     !! right eigenvectors
     do ii=1,HY_VARINUM
        do jj=1,HY_WAVENUM
           RightEigvecTemp(ii,jj) = dot_product(JacQ(ii,:),RightEigvec(:,jj))
        enddo
     enddo
     RightEigvec = RightEigvecTemp

  endif


End Subroutine hy_uhd_eigenVector
