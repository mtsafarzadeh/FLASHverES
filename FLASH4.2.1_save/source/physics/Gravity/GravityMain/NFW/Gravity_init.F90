!!****if* source/physics/Gravity/GravityMain/PointMass/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes the gravitational physics unit for Pointmass.
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Gravity_init()

  use Gravity_data
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

  real :: redshift,Omegam,hubble,rhocrit
#include "constants.h"

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshMe)

  call PhysicalConstants_get("newton", grv_newton)
  call RuntimeParameters_get("mvir",   grv_mvir)
  call RuntimeParameters_get("redshift",redshift)
  call RuntimeParameters_get("nfwc",   grv_nfwc)
  call RuntimeParameters_get("Omegam",Omegam)
  call RuntimeParameters_get("hubble",hubble)
  call RuntimeParameters_get("useGravity", useGravity)

  ! compute the critical density and check if the masses match
  rhocrit  = 1.878E-29*hubble*hubble*(Omegam*(1.+redshift)**3.+(1.-Omegam))
  grv_rvir = (3./(4.*3.1415)*grv_mvir*1.98E3/(180*rhocrit))**(1./3.)*1E10

  grv_rs   = grv_rvir/grv_nfwc
  grv_Fc   = alog(1+grv_nfwc)-grv_nfwc/(1+grv_nfwc)

  grv_factor = -grv_newton * grv_mvir/grv_Fc*1.988E33
  print *,'grv_newton',grv_newton
  print *,'grv_Fc',grv_Fc
  print *,'grv_mvir',grv_mvir*1.988E33
  print *,'grv_factor',grv_factor
  print *,'grv_rvir',grv_rvir
  print *,'sample acc',grv_factor/grv_rvir

!==============================================================================

!==============================================================================

  return
end subroutine Gravity_init
