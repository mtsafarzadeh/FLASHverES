!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhGenEwaldFieldLevel
!!
!! NAME
!!
!!  Gravity_bhGenEwaldFieldLevel
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhGenEwaldFieldLevel(
!!                      integer(in) :: ewald_periodicity,
!!                      integer(in) :: nx,
!!                      integer(in) :: ny,
!!                      integer(in) :: nz,
!!                      real(in)    :: ewald_xmax,
!!                      real(in)    :: ewald_ymax,
!!                      real(in)    :: ewald_zmax,
!!                      real(inout) :: field(-1:,-1:,-1:)
!!        )
!!
!! DESCRIPTION
!!
!!   Generates one level of the nested grid of the Ewald field. Called by
!!   grv_bhGenerateEwaldField.
!!
!! ARGUMENTS
!!
!!   ewald_periodicity : bit array denoting which direction has periodic (1)
!!                       or isolated (0) boundaries
!!   nx : number of points in ewald field in direction x
!!   ny : number of points in ewald field in direction y
!!   nz : number of points in ewald field in direction z
!!   ewald_xmax : size of ewald field in direction x
!!   ewald_ymax : size of ewald field in direction y
!!   ewald_zmax : size of ewald field in direction z
!!   field - ewald field array
!!
!!***

#include "Flash.h"

subroutine Gravity_bhGenEwaldFieldLevel(nx, ny, nz, ewald_xmax, ewald_ymax, ewald_zmax, &
           & ewald_periodicity, field)
  use Logfile_interface, ONLY : Logfile_stamp
  use Gravity_data, ONLY : grv_bhEwaldSeriesN, &
    grv_meshNumProcs, grv_meshMe, grv_meshComm, grv_bhEwaldNRef, &
    grv_bhLx, grv_bhLy, grv_bhLz, grv_bhDirectionQuad
  use grv_bhInterface, ONLY : grv_intCoef1P2I !, DERFCX, DERFC
#ifdef FLASH_USE_SPECFUN
  use grv_bhInterface, ONLY : grv_derfc, grv_derfcx, grv_derf
#endif

  implicit none
#include "Flash_mpi.h"
#include "constants.h"

  integer, intent(IN) :: ewald_periodicity, nx, ny, nz
  real, intent(IN) ::  ewald_xmax, ewald_ymax, ewald_zmax
  real, intent(INOUT) :: field(-1:,-1:,-1:)

  real, parameter :: pi = PI
  integer :: chunk, i1d, i, j, k, ni, nj, nk
  integer :: hi, hj, hk
  integer :: es_nrx, es_nry, es_nrz, es_nfx, es_nfy, es_nfz, es_radius2
  real :: ewald_alpha, ewald_dzeta, ewald_gamma, ewald_beta
  real :: ewald_eta
  real :: ewc1, ewc2, ewc3, ewc4, ewc6, ewc7
  real :: ratio_p, ratio_pinv, ratio_p1, ratio_p2, ratio_pinv1, ratio_pinv2
  real :: Linv, Linv2, ierr
  real :: x, y, z, xni, yni, zni, rni, grav_integral
  real :: hir, hjr, hkr, kx
  real :: cr1, cr2, cr3, cf1, cf2, cf3
  real :: loc_ewald(-1:nx,-1:ny,-1:nz)

! If compilation fails in one of the next lines, or there are other problems
! related to ERFC, see file README.erfc in the BHTree/Wunsch implementation
! directory!

#ifdef USER_ERFC
#define ERFC USER_ERFC
  real, external :: USER_ERFC
#elif defined(FLASH_USE_SPECFUN)
#define ERFC grv_derfc
#else
#define ERFC erfc
  intrinsic erfc
#endif

#ifdef USER_ERFC_SCALED
#define ERFC_SCALED USER_ERFC_SCALED
  real, external :: USER_ERFC_SCALED
#elif defined(FLASH_USE_SPECFUN)
#define ERFC_SCALED grv_derfcx
#else
#define ERFC_SCALED erfc_scaled
  intrinsic erfc_scaled
#endif

#ifdef USER_ERF
#define ERF USER_ERF
  real, external :: USER_ERF
#elif defined(FLASH_USE_SPECFUN)
#define ERF grv_derf
#else
#define ERF erf
  intrinsic erf
#endif

! If compilation fails in one of the preceding lines, or there are other problems
! related to ERFC, see file README.erfc in the BHTree/Wunsch implementation
! directory!

! prepare constants for generating Ewald field
! default values regardless of problem orientation
! range of coefficients in sums in real and Fourier space
  es_nrx = grv_bhEwaldSeriesN
  es_nry = grv_bhEwaldSeriesN
  es_nrz = grv_bhEwaldSeriesN
  es_nfx = grv_bhEwaldSeriesN
  es_nfy = grv_bhEwaldSeriesN
  es_nfz = grv_bhEwaldSeriesN
  es_radius2 = grv_bhEwaldSeriesN*grv_bhEwaldSeriesN

! axis ratio of ellipsis
  cr1 = 1.0
  cr2 = 1.0
  cr3 = 1.0
  cf1 = 1.0
  cf2 = 1.0
  cf3 = 1.0

! particular values for problem with given kind of boundary conditions and orientation
! note ewald_periodicity = 1,2,4 means 1 direction periodic and 2 isolated,
! ewald_periodicity = 3,5,6 means 2 directions periodic and 1 isolated,
! ewald_periodicity = 7 means boundary conditions isolated in tree directions 
  if (ewald_periodicity == 1) then
    Linv = 1.0/grv_bhLx 
    es_nry = 0
    es_nrz = 0
    es_nfy = 0
    es_nfz = 0
  else if (ewald_periodicity == 2) then
    Linv = 1.0/grv_bhLy 
    es_nrx = 0
    es_nrz = 0
    es_nfx = 0
    es_nfz = 0
  else if (ewald_periodicity == 4) then
    Linv = 1.0/grv_bhLz
    es_nrx = 0
    es_nry = 0
    es_nfx = 0
    es_nfy = 0
  else if (ewald_periodicity == 6) then
    Linv = 1.0/grv_bhLy
    ratio_p = grv_bhLz/grv_bhLy
    es_nrx = 0
    es_nrz = ceiling(grv_bhEwaldSeriesN/ratio_p)
    es_nfx = 0
    es_nfz = ceiling(grv_bhEwaldSeriesN*ratio_p)
    cr3 = ratio_p**2
    cf3 = 1.0/(ratio_p**2)
  else if (ewald_periodicity == 5) then
    Linv = 1.0/grv_bhLz
    ratio_p = grv_bhLx/grv_bhLz
    es_nry = 0
    es_nrx = ceiling(grv_bhEwaldSeriesN/ratio_p)
    es_nfy = 0
    es_nfx = ceiling(grv_bhEwaldSeriesN*ratio_p)
    cr1 = ratio_p**2
    cf1 = 1.0/(ratio_p**2)
  else if (ewald_periodicity == 3) then
    Linv = 1.0/grv_bhLx
    ratio_p = grv_bhLy/grv_bhLx
    es_nry = ceiling(grv_bhEwaldSeriesN/ratio_p)
    es_nrz = 0
    es_nfy = ceiling(grv_bhEwaldSeriesN*ratio_p)
    es_nfz = 0
    cr2 = ratio_p**2
    cf2 = 1.0/(ratio_p**2)
  else if (ewald_periodicity == 7) then
! this case we haven't finished yet 
    Linv = 1.0/grv_bhLx 
    ratio_p1 = grv_bhLy/grv_bhLx
    ratio_p2 = grv_bhLz/grv_bhLx
    es_nry = ceiling(grv_bhEwaldSeriesN/ratio_p1)
    es_nrz = ceiling(grv_bhEwaldSeriesN/ratio_p2)
    es_nfy = ceiling(grv_bhEwaldSeriesN*ratio_p1)
    es_nfz = ceiling(grv_bhEwaldSeriesN*ratio_p2)
    cr2 = ratio_p1**2
    cr3 = ratio_p2**2
    cf2 = 1.0/(ratio_p1**2)
    cf3 = 1.0/(ratio_p2**2)
    ratio_pinv1 = 1.0/ratio_p1
    ratio_pinv2 = 1.0/ratio_p2
  endif

! set following constants according to orientation and BCs
  if (ewald_periodicity == 7) then
    Linv2 = Linv/(pi*ratio_p1*ratio_p2)
  else if ((ewald_periodicity == 1).OR.(ewald_periodicity == 2).OR.(ewald_periodicity == 4)) then 
    Linv2 = 2.0*Linv
  else if ((ewald_periodicity == 3).OR.(ewald_periodicity == 5).OR.(ewald_periodicity == 6)) then
    Linv2 = 2.0*Linv/(pi*ratio_p)
  endif

  ratio_pinv = 1.0/ratio_p
  !ewald_alpha = 2.0*Linv*ratio_pinv
  ewald_alpha = 2.0*Linv
  ewald_dzeta = pi*pi*Linv*Linv/(ewald_alpha*ewald_alpha)
 
  ! constants which simplify equations below
  ewc1 = 0.25*pi
  ewc2 = 1.0/sqrt(ewald_dzeta)
  ewc3 = 0.5*ewc2
  ewc4 = 2.0*sqrt(ewald_dzeta/pi)
  ewc6 = 0.25/ewald_dzeta

! compute Ewald field
  ! size of a chunk of data for a given processor
  chunk = 1 + ((nz+2)*(ny+2)*(nx+2) / grv_meshNumProcs)

  do k = -1,nz
    do j = -1,ny
      do i = -1,nx

        ! on all CPUs set each element to zero at first
        field(i,j,k) = 0.0
        loc_ewald(i,j,k) = 0.0

        ! calculate the 1D index: 0..ewald_field_z*ny*ewald_field_x-1
        i1d = (k+1)*(ny+2)*(nx+2) + (j+1)*(nx+2) + (i+1)
        ! check if this point should be calculated on this processor
        if ((i1d >= grv_meshMe*chunk) .and. (i1d < (grv_meshMe+1)*chunk)) then
          ! coordinates of the point
          x = i * ewald_xmax / (nx-1)
          y = j * ewald_ymax / (ny-1)
          z = k * ewald_zmax / (nz-1)



          ! subtract term 1/r - this enables more accurate interpolation in
          ! function grv_bhEwald()
          ! this subtraction is compensated in routines Gravity_bhNodeContrib 
          ! and Gravity_bhBotNodeContrib
          !rni = sqrt(x*x + y*y + z*z)
          !loc_ewald(i,j,k) = loc_ewald(i,j,k) - 1.0/(rni + 1.0d-99)

          ! first term - short term interactions
          do ni = -es_nrx,es_nrx
            do nj = -es_nry,es_nry
              do nk = -es_nrz,es_nrz
              ! terms with non-negligible contributions must lie inside ellipse
                if ((cr1*ni*ni+cr2*nj*nj+cr3*nk*nk) <= es_radius2) then
                  xni = x + ni*grv_bhLx
                  yni = y + nj*grv_bhLy
                  zni = z + nk*grv_bhLz
                  rni = sqrt(xni*xni + yni*yni + zni*zni)

                  loc_ewald(i,j,k) = loc_ewald(i,j,k) + ERFC(ewald_alpha*rni) / (rni + 1.0d-99)

                endif
              enddo
            enddo
          enddo

          ! second term - long term interactions
          do hi = -es_nfx,es_nfx
            do hj = -es_nfy,es_nfy
              do hk = -es_nfz,es_nfz
                ! terms with non-negligible contributions must lie inside ellipse 
                ! (perpendicular to the previous ellipse)
               if ((cf1*hi*hi+cf2*hj*hj+cf3*hk*hk) <= es_radius2) then

!               periodic boundary conditions in one direction
                if ((ewald_periodicity == 1).OR.(ewald_periodicity == 2).OR. & 
                &  (ewald_periodicity == 4)) then
                  if (ewald_periodicity == 1) then
                    ewald_eta = 2*pi*Linv*sqrt(y**2+z**2)
                  else if (ewald_periodicity == 2) then
                    ewald_eta = 2*pi*Linv*sqrt(x**2+z**2)
                  else if (ewald_periodicity == 4) then
                    ewald_eta = 2*pi*Linv*sqrt(x**2+y**2)
                  endif

                  loc_ewald(i,j,k) = loc_ewald(i,j,k) + Linv2*grv_intCoef1P2I(max(abs(hi),abs(hj),abs(hk)),ewald_eta,ewald_dzeta) &  
                  &                * exp(-ewald_dzeta*(hi**2+hj**2+hk**2))*cos((2*pi)*(hi*x+hj*y+hk*z)*Linv)


!               periodic boundary conditions in two directions
                else if ((ewald_periodicity == 3).or.(ewald_periodicity == 5).or. &
                &   (ewald_periodicity == 6)) then
                  ! hr's (needed for ewald_beta) and ewald_gamma according to orientation
                  hir = real(hi)
                  hjr = real(hj)
                  hkr = real(hk)
                  if (ewald_periodicity == 6) then
                    ewald_gamma = 2*pi*x*Linv
                    hkr = ratio_pinv*real(hk)
                  else if (ewald_periodicity == 5) then
                    ewald_gamma = 2*pi*y*Linv
                    hir = ratio_pinv*real(hi)
                  else if (ewald_periodicity == 3) then
                    ewald_gamma = 2*pi*z*Linv
                    hjr = ratio_pinv*real(hj)
                  endif

                  if ((hi == 0).and.(hj == 0).and.(hk == 0)) then
                    ! analytical evaluation of appropriate integral
                    grav_integral = 2.0*(ewc1*(ewald_gamma*ERF(ewc3*ewald_gamma)+ &
                        & ewc4*exp(-ewc6*ewald_gamma*ewald_gamma)))
                    loc_ewald(i,j,k)=loc_ewald(i,j,k)-grav_integral*Linv2
                  else
                    ewald_beta = cos((2*pi)*(hir*x+hjr*y+hkr*z)*Linv)
! rename h2=ewc7 ?
                    ewc7 = sqrt(hir*hir+hjr*hjr+hkr*hkr)
                    ! analytical evaluation of appropriate integral 
                    ! (in order to avoid behaviour like infty*zero we performed 
                    ! decomposition into Chebyschev polynomials)
                    grav_integral = ewc1*(exp(-ewc7*ewald_gamma+ewald_dzeta*ewc7**2) &
                    &             * ERFC(ewc2*(ewc7*ewald_dzeta-0.5*ewald_gamma)) &
                    &             + exp(-ewc3*ewc3*ewald_gamma*ewald_gamma) &
                    &             * ERFC_SCALED(ewc2*(ewc7*ewald_dzeta+0.5*ewald_gamma)))/ewc7
                    ! alternative evaluation:
                    !grav_integral = ewc1*exp(-ewc7*ewald_gamma+ewald_dzeta*ewc7*ewc7)* &
                    !              & (ERFC(ewc2*(ewc7*ewald_dzeta-0.5*ewald_gamma)) &
                    !              & +exp(2.0*ewc7*ewald_gamma)* &
                    !              & ERFC(ewc2*(ewc7*ewald_dzeta+0.5*ewald_gamma)))/ewc7
                    ! alternative evaluation2:
                    !grav_integral = ewc1*exp(-ewc3*ewc3*ewald_gamma*ewald_gamma) &
                    !&             * (grv_erfcapp(ewc2*(ewc7*ewald_dzeta-0.5*ewald_gamma)) & 
                    !&             +  grv_erfcapp(ewc2*(ewc7*ewald_dzeta+0.5*ewald_gamma)))/ewc7
                    !
                    loc_ewald(i,j,k) = loc_ewald(i,j,k) &
                    &                + exp(-ewald_dzeta*(hir*hir+hjr*hjr+hkr*hkr)) &
                    &                * grav_integral*ewald_beta*Linv2
                  endif
                  
                  ! periodical BCs in three directions
                else if (ewald_periodicity == 7) then
                  if ((hi**2+hj**2+hk**2).gt.0) then
                    ewc7=hi**2+(hj*ratio_pinv1)**2+(hk*ratio_pinv2)**2
                    kx = 2*pi*(hi*x + hj*y*ratio_pinv1 + hk*z*ratio_pinv2)*Linv
                    loc_ewald(i,j,k) = loc_ewald(i,j,k) + Linv2*exp(-ewald_dzeta*ewc7) &
                    &                * cos(kx) / ewc7
                  endif

                endif

               endif
              enddo
            enddo
          enddo
        endif ! end my chunk

      enddo
    enddo
  enddo

  call MPI_AllReduce(loc_ewald,field,(nx+2)*(ny+2)*(nz+2),FLASH_REAL,FLASH_SUM,grv_meshComm,ierr)  
  field(0,0,0) = field(1,1,1)

  return
end subroutine Gravity_bhGenEwaldFieldLevel

