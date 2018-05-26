!!****if* source/physics/Hydro/HydroMain/unsplit_old/multiTemp/hy_uhd_MultiTempData
!!
!! NAME
!!   
!!  hy_uhd_MultiTempData 
!!
!!
!! SYNOPSIS
!!
!!  use hy_uhd_MultiTempData
!!
!!
!! DESCRIPTION
!!
!!***


module hy_uhd_MultiTempData

#include "Flash.h"
#include "constants.h"

  !!*****Runtime parameters*****
  integer, save :: hy_eosMode
!!$  integer, save :: hy_ppmEintFluxConstructionMeth = 0
!!$  integer, save :: hy_ppmEnerFluxConstructionMeth = 0
!!$  integer, save :: hy_ppmEintCFluxConstructionMeth = 0
!!$  integer, save :: hy_ppmEnerCFluxConstructionMeth = 0

  ! DEV: The following should be in Hydro_data instead !? - KW
  logical,save :: hy_updateHydroFluxes

  !!*****End Runtime parameters*****


  !!*****End Directly Derived from Runtime parameters*****

  ! Maybe should become a runtime parameter.
  logical, parameter :: hy_alwaysCallDetectShock = .TRUE.

  logical, save :: hy_movingGrid
  
  !!*****Constants database
  real,    save :: hy_pi
  real,    save :: hy_eMass, hy_pMass, hy_eMassInUAmu


  integer, parameter :: hy_numXN     = NSPECIES+NMASS_SCALARS
  integer, parameter :: hy_numMS     = NMASS_SCALARS

  ! For testing ways to advect components and handle shock heating
  
  integer, save :: hy_3Ttry_B, &
                   hy_3Ttry_E, &
                   hy_3Ttry_F, &
                   hy_3Ttry_G

  real, save :: hy_3Ttry_D
  integer, save :: hy_3Ttry_B_rad

end module hy_uhd_MultiTempData
