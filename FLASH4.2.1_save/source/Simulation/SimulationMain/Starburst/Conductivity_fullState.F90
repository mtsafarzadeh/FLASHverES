!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity_fullState
!!
!! NAME
!!  Conductivity_fullState
!!
!! SYNOPSIS
!!  call Conductivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                     OPTIONAL,real(out)   :: isochoricCond,
!!                     OPTIONAL,real(out)   :: diffCoeff,
!!                     OPTIONAL,integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron conductivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!!  Returns thermal conductivity and/or diffusivity coefficients.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   isochoricCond  :   isochoric conductivity
!!   diffCoeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!***

#include "Flash.h"
#include "constants.h"  
#include "Eos.h"


subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_interface, ONLY: Conductivity
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele, cond_navo
  use Eos_interface, ONLY: Eos, Eos_getAbarZbar, Eos_getTempData

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff
  real,    OPTIONAL, intent(OUT)  :: isochoricCond
  integer, OPTIONAL, intent(IN) :: component

  real, pointer :: massfrac(:)

  real :: isochoricCondLoc, diffCoeffLoc
  integer :: componentLoc, tempToUse

  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen

  real :: xtemp, xden
  real :: nele, ll
  real :: abar, zbar

  real, parameter :: cexp = 2.5

  isochoricCondLoc = 0.0
  diffCoeffLoc = 0.0
  
#if defined(DENS_VAR) && defined(TEMP_VAR)
  if (present(component)) then
     componentLoc = component
  else
     componentLoc = 0
  end if

  tempToUse = TEMP_VAR

  select case (componentLoc)
  case(1)
#ifdef TION_VAR
     tempToUse = TION_VAR
#endif
  case(2)
#ifdef TELE_VAR
     tempToUse = TELE_VAR
#endif
  case(3)
#ifdef TRAD_VAR
     tempToUse = TRAD_VAR
#endif
  end select

    if (cond_useConductivity) then
       call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)

       xden = solnVec(DENS_VAR)
       xtemp = solnVec(tempToUse)
       nele = zbar * xden * cond_navo / abar
       call cond_logLambda(xtemp, nele, zbar, ll)

       isochoricCondLoc = (8.0/PI)**1.5*cond_boltz**3.5 / (sqrt(cond_mele)*cond_qele**4) * &
            xtemp**cexp / (ll * (zbar + 3.3))

       if (present(diffCoeff)) then

          vecLen = 1
          !mode = MODE_DENS_TEMP_GATHER
          mode = MODE_DENS_TEMP
          call Eos_getTempData(solnVec,eos_arr,mode)
          eos_arr(EOS_DENS) = xden
          mask = .false.
          mask(EOS_CV) = .true.
          !mask(EOS_CVELE)  = .true.
          mask(EOS_DET) = .true.
          if (NSPECIES > 0) then
             massfrac => solnVec(SPECIES_BEGIN:SPECIES_END)
             call Eos(mode,vecLen,eos_arr,massfrac,mask)
          else
             call Eos(mode,vecLen,eos_arr,mask=mask)
          end if
          diffCoeffLoc = isochoricCondLoc/(xden*eos_arr(EOS_CV))
       end if

    end if

#endif

  if(present(isochoricCond)) isochoricCond = isochoricCondLoc
  if(present(diffCoeff)) diffCoeff = diffCoeffLoc
  
end subroutine Conductivity_fullState
