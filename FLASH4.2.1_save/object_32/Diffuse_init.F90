!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_init
!!
!! NAME
!!
!!  Diffuse_init
!!
!!
!! SYNOPSIS
!!
!!  call Diffuse_init()
!!
!! Description
!!
!!  Initializes local data for Unit Diffuse defined in Module Diffuse_data.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse, and used by all three major diffuse subroutines,
!!  Diffuse_species(for mass diffusivity), useDiffusetherm (for 
!!  thermal Conductivity), and useDiffuseVisc (for viscosity).
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!    useDiffuse
!!        whether any method of the Diffuse unit should contribute to fluxes
!!    useDiffuseTherm
!!        whether Diffuse_therm should contribute to fluxes
!!    useDiffuseVisc
!!        whether Diffuse_visc should contribute to fluxes
!!    useDiffuseSpecies
!!        whether Diffuse_therm should contribute to fluxes [TO BE IMPLEMENTED]
!!    geometry [STRING]
!!        Grid geometry
!!    dt_diff_factor
!!        factor that scales the timestep returned by Diffuse_computeDt
!!    diffusion_cutoff_density
!!        density below which we no longer diffuse
!!    thermal_diff_method
!!***


subroutine Diffuse_init()

  use Diffuse_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getComm, &
       Driver_getNumProcs, Driver_getMype
  use Logfile_interface, ONLY : Logfile_stampMessage
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
   RuntimeParameters_mapStrToInt
  use Grid_interface, ONLY:   Grid_setFluxHandling
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"


  character(len=MAX_STRING_LENGTH),save :: str_geometry
  character(len=MAX_STRING_LENGTH) :: flmode_str
  real :: ssol

  ! Everybody should know this
  call Driver_getMype(MESH_COMM,diff_meshMe)
  call Driver_getNumProcs(MESH_COMM,diff_meshNumProcs)
  call Driver_getComm(MESH_COMM,diff_meshComm)

  call RuntimeParameters_get('dt_diff_factor', dt_diff_factor)
  call RuntimeParameters_get('useDiffuseTherm', useDiffuseTherm)
  call RuntimeParameters_get('useDiffuseVisc', useDiffuseVisc)
  call RuntimeParameters_get('useMagneticResistivity', useDiffuseMagneticResistivity)
  call RuntimeParameters_get('useDiffuseSpecies', useDiffuseSpecies)
  call RuntimeParameters_get("diff_useEleCond", diff_useEleCond)
  call RuntimeParameters_get("diff_useIonCond", diff_useIonCond)
  call RuntimeParameters_get('diffusion_cutoff_density', diffusion_cutoff_density)
  call RuntimeParameters_get('useDiffuse',useDiffuse)
  if (.not. useDiffuse) then
     write(6,*)'WARNING:  You have included the Diffuse physics unit '
     write(6,*)'   but have set the runtime parameter useDiffuse to FALSE'
     write(6,*)'   No Diffusion will occur but Diffuse_init will continue.'
  end if

  call RuntimeParameters_get ("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, diff_geometry)

!! Determine the geometries of the individual dimensions

  if (diff_geometry == CARTESIAN)then
     dirGeom(IAXIS) = XYZ
     dirGeom(JAXIS) = XYZ
     dirGeom(KAXIS) = XYZ
  elseif(diff_geometry == POLAR)then
     dirGeom(IAXIS) = RAD_CYL
     dirGeom(JAXIS) = PHI_CYL
     dirGeom(KAXIS) = XYZ
  elseif(diff_geometry == CYLINDRICAL) then
     dirGeom(IAXIS) = RAD_CYL
     dirGeom(JAXIS) = XYZ
     dirGeom(KAXIS) = PHI_CYL
  elseif(diff_geometry == SPHERICAL) then
     dirGeom(IAXIS) = RAD_SPH
     dirGeom(JAXIS) = THETA
     dirGeom(KAXIS) = PHI_SPH
  else
     call Driver_abortFlash("unsupported geometry ")
  end if

#ifdef FLASH_MULTISPECIES
  diff_singleSpeciesZ = 0.0     !value should not be used!
  diff_singleSpeciesA = 0.0     !value should not be used!
#else
  call RuntimeParameters_get("eos_singleSpeciesZ", diff_singleSpeciesZ)
  call RuntimeParameters_get("eos_singleSpeciesA", diff_singleSpeciesA)
#endif

  ! Store physical constants:
  call PhysicalConstants_get("speed of light",diff_speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",ssol)
  diff_asol    = 4.0 * ssol / diff_speedlt
  call PhysicalConstants_get("electron mass", diff_mele)
  call PhysicalConstants_get("Boltzmann", diff_boltz)
  call PhysicalConstants_get("Avogadro", diff_avo)

  call diff_saInit

  call RuntimeParameters_get('diff_eleFlMode', flmode_str)
  call makeLowercase(flmode_str)
  call RuntimeParameters_mapStrToInt(flmode_str, diff_eleFlMode)

  call RuntimeParameters_get('diff_ionFlMode', flmode_str)
  call makeLowercase(flmode_str)
  call RuntimeParameters_mapStrToInt(flmode_str, diff_ionFlMode)

  call RuntimeParameters_get('diff_eleFlCoef', diff_eleFlCoef)
  call RuntimeParameters_get('diff_ionFlCoef', diff_ionFlCoef)

  call diff_fbInit()

  return
end subroutine Diffuse_init
