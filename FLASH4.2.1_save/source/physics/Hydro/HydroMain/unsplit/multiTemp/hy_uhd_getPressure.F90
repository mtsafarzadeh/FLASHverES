!!****if* source/physics/Hydro/HydroMain/unsplit/multiTemp/hy_uhd_getPressure
!!
!! NAME
!!
!!  hy_uhd_getPressure
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_getPressure(real, intent(in)  :: eele,
!!                          real, intent(in)  :: eion,
!!                          real, intent(in)  :: erad,
!!                          real, intent(in)  :: soln(NUNK_VARS),
!!                          real, intent(out) :: pele,
!!                          real, intent(out) :: pion,
!!                          real, intent(out) :: prad) 
!!
!! DESCRIPTION
!! 
!! Calculates pele, ion and rad given energies and solution vector
!!
!! ARGUMENTS
!!
!!  eele        - electron energy 
!!  eion        - ion energy 
!!  erad        - radiation energy 
!!  soln        - solution vector     
!!  pele        - electron pressure
!!  pion        - ion pressure
!!  prad        - radiation pressure
!!
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine hy_uhd_getPressure( &
     eele, eion, erad, &
     soln, &
     pele, pion, prad)

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  use Eos_interface, ONLY: Eos

  implicit none

  ! Arguments:
  real, intent(in)  :: eele
  real, intent(in)  :: eion
  real, intent(in)  :: erad
  real, intent(in)  :: soln(NUNK_VARS)
  real, intent(out) :: pele
  real, intent(out) :: pion
  real, intent(out) :: prad
  
  ! Local variables:
  integer :: n

  integer :: mode, vecLen
  real    :: massfrac(NSPECIES)
  real    :: eos_arr(EOS_NUM)
  logical :: mask(EOS_VARS+1:EOS_NUM)

  ! Subroutine body:
  vecLen = 1

  ! Load inputs into EOS array. Since we are calling with
  ! MODE_DENS_EI_GATHER
  eos_arr = 0.0
  eos_arr(EOS_DENS) = soln(DENS_VAR)
  eos_arr(EOS_EINTELE) = eele
  eos_arr(EOS_EINTION) = eion
  eos_arr(EOS_EINTRAD) = erad
  eos_arr(EOS_TEMPELE) = soln(TELE_VAR)
  eos_arr(EOS_TEMPION) = soln(TION_VAR)
  eos_arr(EOS_TEMPRAD) = soln(TRAD_VAR)

  ! Use the mask array to ask EOS for pressures:
  mask = .false.
  mask(EOS_PRESELE) = .true.
  mask(EOS_PRESION) = .true.
  mask(EOS_PRESRAD) = .true.

  ! Load species, if defined. Checking avoids annoying compiler
  ! warnings:
#if (NSPECIES > 0)
  do n = 1, NSPECIES
     massfrac(n) = soln(SPECIES_BEGIN-1+n)
  enddo
#endif

  call Eos(MODE_DENS_EI_GATHER,vecLen,eos_arr,massfrac,mask)

  ! Recover the output:
  pele = eos_arr(EOS_PRESELE)
  pion = eos_arr(EOS_PRESION)
  prad = eos_arr(EOS_PRESRAD)

end subroutine hy_uhd_getPressure
