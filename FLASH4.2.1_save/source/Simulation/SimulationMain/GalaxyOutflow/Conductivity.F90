!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity
!!
!! NAME
!!
!!  Conductivity
!!
!! SYNOPSIS
!!
!!  call Conductivity(real, intent(IN)  :: xtemp,
!!                    real, intent(IN)  :: xden,
!!                    real, dimension(NSPECIES), intent(IN)  :: massfrac,
!!                    real, intent(OUT)  :: cond,
!!                    real, intent(OUT)  :: diff_coeff,
!!                    integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron conductivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   cond       :   conductivity
!!   diff_coeff :   diffusion coefficient ( = cond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are reqested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!***

#include "Flash.h"

subroutine Conductivity(xtemp,xden,massfrac,cond,diff_coeff, component)
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele, cond_navo
  use Eos_interface, ONLY:     Eos, Eos_getAbarZbar
  use Cosmology_data,    ONLY: csm_scaleFactor
  implicit none
  
#include "constants.h"  
#include "Eos.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, cond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen

  real :: nele, ll
  real :: abar, zbar
  real :: cond_phys
  real :: xden_phys
  real :: xtemp_phys

  real, parameter :: cexp = 2.5
  if (cond_useConductivity) then
     call Eos_getAbarZbar(abar=abar,zbar=zbar,massFrac=massfrac)

     ! use physical density here
     xden_phys = xden/(csm_scaleFactor)**3
     xtemp_phys = xtemp*(csm_scaleFactor)**2

     nele = zbar * xden_phys* cond_navo / abar
     call cond_logLambda(xtemp_phys, nele, zbar, ll)
    
     ! use physical temperature here 
     cond_phys = (8.0/PI)**1.5*cond_boltz**3.5 / (sqrt(cond_mele)*cond_qele**4) * &
          (xtemp_phys)**cexp / (ll * (zbar + 3.3))
!    physical cond is comoving cond /a
     cond = cond_phys*csm_scaleFactor
     vecLen = 1
     mode = MODE_DENS_TEMP
     eos_arr(EOS_TEMP) = xtemp_phys
     eos_arr(EOS_DENS) = xden_phys

     mask = .false.
     mask(EOS_CV)  = .true.
     mask(EOS_DET) = .true.

     call Eos(mode,vecLen,eos_arr,massfrac,mask)
     ! this is the physical diff_coefficient
     diff_coeff = cond_phys/((xden_phys)*eos_arr(EOS_CV))
     ! the physical diff_coeff is  a^2 comoving one
     ! so the comoving one diff/a*2
     !print *,'temp phys',xtemp_phys
     !print *,'dens phys',xden_phys
     !print *,'cond phys',cond_phys
     !print *,'diff coeff phys,',diff_coeff
     diff_coeff = diff_coeff/(csm_scaleFactor**2) 
  else
     cond = 0.0
     diff_coeff = 0.0
  end if
   !print *,'diff coeff comoving inside,',diff_coeff

end subroutine Conductivity





