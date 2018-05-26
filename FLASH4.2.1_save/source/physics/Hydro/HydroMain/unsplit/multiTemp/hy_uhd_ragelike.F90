!!****if* source/physics/Hydro/HydroMain/unsplit/multiTemp/hy_uhd_ragelike
!!
!! NAME
!!
!!  hy_uhd_ragelike
!!
!!
!!  SYNOPSIS
!!
!!  hy_uhd_ragelike(real, intent(in)  :: uextra,
!!                  real, intent(in)  :: soln(NUNK_VARS),
!!                  real, intent(in)  :: dens_old,
!!                  real, intent(in)  :: duele_adv,
!!                  real, intent(in)  :: duion_adv,
!!                  real, intent(in)  :: durad_adv,
!!                  real, intent(in)  :: xc,
!!                  real, intent(in)  :: yc, 
!!                  real, intent(in)  :: zc,
!!                  real, intent(out) :: uele_new,
!!                  real, intent(out) :: uion_new,
!!                  real, intent(out) :: urad_new)
!!  
!!
!! DESCRIPTION
!!
!! Ragelike EOS
!!  
!!
!! ARGUMENTS
!!   uextra           - u to be divided 
!!   soln(NUNK_VARS)  - solution array
!!   dens_old         - old density 
!!   duele_adv        - advected du_ele
!!   duion_adv        - advected du_ion
!!   durad_adv        - advected du_rad
!!   xc               - x coordinate
!!   yc               - y coordinate   
!!   zc               - z coordinate
!!   uele_new         - output uele
!!   uion_new         - output uion
!!   urad_new         - output urad
!!
!! DESCRIPTION
!!
!!   This routine divides uextra in uele,ion and rad given their _adv counterparts
!!
!! SIDE EFFECTS
!!
!!
!!***


subroutine hy_uhd_ragelike(uextra, soln, dens_old, &
     duele_adv, duion_adv, durad_adv, &
     xc, yc, zc, &
     uele_new, uion_new, urad_new)

#include "Flash.h"
#include "constants.h"

  use Logfile_interface,    ONLY: Logfile_stampMessage
  use Hydro_interface,      ONLY: Hydro_recalibrateEints
  use hy_uhd_MultiTempData, ONLY: hy_3Ttry_D

  implicit none

  ! Arguments:
  real, intent(in)  :: uextra
  real, intent(in)  :: soln(NUNK_VARS)
  real, intent(in)  :: dens_old
  real, intent(in)  :: duele_adv
  real, intent(in)  :: duion_adv
  real, intent(in)  :: durad_adv
  real, intent(in)  :: xc 
  real, intent(in)  :: yc 
  real, intent(in)  :: zc
  real, intent(out) :: uele_new
  real, intent(out) :: uion_new
  real, intent(out) :: urad_new

  ! Local variables:
  real :: pele_adv
  real :: pion_adv
  real :: prad_adv
  real :: PeP
  real :: PiP
  real :: PrP
  real :: dens_new
  real :: utot_new
  real :: uele_old
  real :: uion_old
  real :: urad_old
  real :: uele_adv
  real :: uion_adv
  real :: urad_adv

  character(len=MAX_STRING_LENGTH) :: errmsg

  ! Subroutine body:
  dens_new = soln(DENS_VAR)
  uele_old = soln(EELE_VAR) * dens_old
  uion_old = soln(EION_VAR) * dens_old
  urad_old = soln(ERAD_VAR) * dens_old
  uele_adv = uele_old + duele_adv
  uion_adv = uion_old + duion_adv
  urad_adv = urad_old + durad_adv

  ! This will call EOS to compute the most up-to-date
  ! pressures. It is important to use up-to-date values when
  ! computing the advected quantities:

  call hy_uhd_getPressure( &
       uele_adv / dens_new, &
       uion_adv / dens_new, &
       urad_adv / dens_new, &
       soln, &
       pele_adv, &
       pion_adv, &
       prad_adv)

  ! Try using pressure ratios to update internal energies:
  if (hy_3Ttry_D==3.0) then
     PrP = 0.0
     PeP = pele_adv / (pele_adv + pion_adv)
     PiP = pion_adv / (pele_adv + pion_adv)
  else
     PrP = prad_adv / (pele_adv + pion_adv + prad_adv)
     PeP = pele_adv / (pele_adv + pion_adv + prad_adv)
     PiP = pion_adv / (pele_adv + pion_adv + prad_adv)
  end if

  call Hydro_recalibrateEints(uextra,PiP,PeP,PrP)
  uele_new = uele_old + (duele_adv + PeP)
  uion_new = uion_old + (duion_adv + PiP)
  urad_new = urad_old + (durad_adv + PrP)

  ! Now, check for negative internal energies:
  if(uele_new <= 0.0 .or. uion_new <= 0.0 .or. urad_new < 0.0 ) then

     call Logfile_stampMessage(&
          "[hy_uhd_unsplitUpdateMultiTemp] Negative internal energies in 3T update", &
          .true.)
     call Logfile_stampMessage(&
          "  Falling back to using internal energy ratios instead of pressure ratios", &               
          .true.)
     call Logfile_stampMessage("", .true.)

     ! Negative internal energy detected, try to update
     ! using internal energies instead...

     if (hy_3Ttry_D==3.0) then
        PrP = 0.0
        PeP = uele_adv / (uele_adv + uion_adv)
        PiP = uion_adv / (uele_adv + uion_adv)
     else
        ! PrP = U(PRAD_VAR,i,j,k) / U(PRES_VAR,i,j,k)
        PrP = urad_adv / (uele_adv + uion_adv + urad_adv)
        PeP = uele_adv / (uele_adv + uion_adv + urad_adv)
        PiP = uion_adv / (uele_adv + uion_adv + urad_adv)
     end if

     call Hydro_recalibrateEints(uextra,PiP,PeP,PrP)
     uele_new = uele_old + (duele_adv + PeP)
     uion_new = uion_old + (duion_adv + PiP)
     urad_new = urad_old + (durad_adv + PrP)

  end if

  ! Perform error check for negative internal energies:
  if( uele_new <= 0.0 .or. uion_new <= 0.0 .or. urad_new < 0.0 ) then

     utot_new = uele_adv + uion_adv + urad_adv + uextra

     call Logfile_stampMessage("[hy_uhd_ragelike] ERROR DETECTED", .true.)
     call Logfile_stampMessage("  Zero or negative advected specific internal energy detected", .true.)

     write(errmsg, "(a,1pe20.13)") "  DENS_OLD = ", dens_old
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  DENS_NEW = ", dens_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  UTOT_NEW = ", utot_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  UEXTRA = ", uextra
     call Logfile_stampMessage(trim(errmsg), .true.)

     call Logfile_stampMessage("  NEW INTERNAL ENERGY DENSITIES:")              
     write(errmsg, "(a,1pe20.13)") "    UELE_NEW = ", uele_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    UION_NEW = ", uion_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    URAD_NEW = ", urad_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     call Logfile_stampMessage("  OLD SPECIFIC INTERNAL ENERGIES:")
     write(errmsg, "(a,1pe20.13)") "    EELE_VAR = ", soln(EELE_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    EION_VAR = ", soln(EION_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    ERAD_VAR = ", soln(ERAD_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     call Logfile_stampMessage("  NEW SPECIFIC INTERNAL ENERGIES:")
     write(errmsg, "(a,1pe20.13)") "    EELE     = ", uele_new/dens_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    EION     = ", uion_new/dens_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    ERAD     = ", urad_new/dens_new
     call Logfile_stampMessage(trim(errmsg), .true.)

     call Logfile_stampMessage("  ADVECTED INTERNAL ENERGY DENSITIES:")
     write(errmsg, "(a,1pe20.13)") "    uele_adv  = ", uele_adv
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    uion_adv  = ", uion_adv
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "    urad_adv  = ", urad_adv
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  PELE_VAR = ", soln(PELE_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  PION_VAR = ", soln(PION_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  PRAD_VAR = ", soln(PRAD_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  PELE_ADV = ", pele_adv
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  PION_ADV = ", pion_adv
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  PRAD_ADV = ", prad_adv
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  TELE_VAR = ", soln(TELE_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  TION_VAR = ", soln(TION_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  TRAD_VAR = ", soln(TRAD_VAR)
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  CELL X = ", xc
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  CELL Y = ", yc
     call Logfile_stampMessage(trim(errmsg), .true.)

     write(errmsg, "(a,1pe20.13)") "  CELL Z = ", zc
     call Logfile_stampMessage(trim(errmsg), .true.)

     call Logfile_stampMessage("  You might be able to get things to keep running by trying", .true.)
     call Logfile_stampMessage("  more diffusive hydro options (reduce CFL, order, etc...)", .true.)

     call Driver_abortFlash("[hy_uhd_ragelike] Negative 3T internal energy, CHECK LOG")

  end if



end subroutine hy_uhd_ragelike
