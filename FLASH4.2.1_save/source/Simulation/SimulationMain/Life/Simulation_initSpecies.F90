!!****if* source/Simulation/SimulationMain/Life/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for
!!  the TwoGamma setup, which advects two fluids with different Gamma values
!!
!!***

subroutine Simulation_initSpecies()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use Multispecies_interface, ONLY : Multispecies_setProperty

  implicit none
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
#include "constants.h"

  character(len=MAX_STRING_LENGTH) :: str

  ! ! Load all runtime parameters that are needed for setting up the
  ! ! multispecies database
  ! call RuntimeParameters_get('sim_abarTarg', sim_abarTarg)
  ! call RuntimeParameters_get('sim_zbarTarg', sim_zbarTarg)
  ! call RuntimeParameters_get('sim_abarCham', sim_abarCham)
  ! call RuntimeParameters_get('sim_zbarCham', sim_zbarCham)
  ! call RuntimeParameters_get('sim_eosTarg', str)
  ! call RuntimeParameters_mapStrToInt(str, sim_eosTarg)
  ! call RuntimeParameters_get('sim_eosCham', str)
  ! call RuntimeParameters_mapStrToInt(str, sim_eosCham)
  ! call RuntimeParameters_get('sim_zminTarg', sim_zminTarg)

  ! call Multispecies_setProperty(TARG_SPEC, A, sim_abarTarg)
  ! call Multispecies_setProperty(TARG_SPEC, Z, sim_zbarTarg)
  ! call Multispecies_setProperty(TARG_SPEC, GAMMA, 1.6666666666666666667e0)
  ! call Multispecies_setProperty(TARG_SPEC, MS_EOSTYPE, sim_eosTarg)
  ! call Multispecies_setProperty(TARG_SPEC, MS_ZMIN, sim_zminTarg)

  ! call Multispecies_setProperty(CHAM_SPEC, A, sim_abarCham)
  ! call Multispecies_setProperty(CHAM_SPEC, Z, sim_zbarCham)
  ! call Multispecies_setProperty(CHAM_SPEC, GAMMA, 1.6666666666666666667e0)
  ! call Multispecies_setProperty(CHAM_SPEC, MS_EOSTYPE, sim_eosCham)


end subroutine Simulation_initSpecies

