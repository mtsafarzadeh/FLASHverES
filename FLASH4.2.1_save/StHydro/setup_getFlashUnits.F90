!!****f* object/setup_flashUnits
!!
!! NAME
!!
!!  setup_getFlashUnits
!!
!!
!! SYNOPSIS
!!
!!
!!  setup_getFlashUnits(unit_names)
!!
!!  setup_getFlashUnits(character())
!!
!!
!! DESCRIPTION
!!
!!  Return a character array of size NUM_UNITS containing
!!  the names of all of the FLASH units used to assemble
!!  the current executable
!!
!!  The unit_names variable should be declared as
!!
!!    use flashUnits
!!
!!  
!!    character (len=MAX_STRING_LENGTH) :: flash_units(NUM_UNITS) 
!!
!!
!!  The length of each character string is set to MAX_STRING_LENGTH,
!!  which is defined in the automatically generated flash_defines.fh
!!
!!***

  subroutine setup_getFlashUnits(unit_names)

#include "constants.h"
    implicit none

    integer, PARAMETER :: NUM_UNITS = 21
    character (len=MAX_STRING_LENGTH) :: unit_names(NUM_UNITS)
    integer :: i

    i = 0

    i = i + 1; unit_names(i) = &
"Driver/DriverMain/Unsplit"
    i = i + 1; unit_names(i) = &
"Grid/GridBoundaryConditions"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/UG"
    i = i + 1; unit_names(i) = &
"Grid/localAPI"
    i = i + 1; unit_names(i) = &
"IO/IOMain/hdf5/serial/UG"
    i = i + 1; unit_names(i) = &
"IO/localAPI"
    i = i + 1; unit_names(i) = &
"PhysicalConstants/PhysicalConstantsMain"
    i = i + 1; unit_names(i) = &
"RuntimeParameters/RuntimeParametersMain"
    i = i + 1; unit_names(i) = &
"Simulation/SimulationMain/StirTurbHydro"
    i = i + 1; unit_names(i) = &
"flashUtilities/contiguousConversion"
    i = i + 1; unit_names(i) = &
"flashUtilities/general"
    i = i + 1; unit_names(i) = &
"flashUtilities/interpolation/oneDim"
    i = i + 1; unit_names(i) = &
"flashUtilities/nameValueLL"
    i = i + 1; unit_names(i) = &
"flashUtilities/rng"
    i = i + 1; unit_names(i) = &
"flashUtilities/system/memoryUsage/legacy"
    i = i + 1; unit_names(i) = &
"monitors/Logfile/LogfileMain"
    i = i + 1; unit_names(i) = &
"monitors/Timers/TimersMain/MPINative"
    i = i + 1; unit_names(i) = &
"physics/Eos/EosMain/Gamma"
    i = i + 1; unit_names(i) = &
"physics/Eos/localAPI"
    i = i + 1; unit_names(i) = &
"physics/Hydro/HydroMain/unsplit/Hydro_Unsplit"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Stir/StirMain/StirOrig"


    return

  end subroutine setup_getFlashUnits

  subroutine setup_getNumFlashUnits(numUnits)

#include "constants.h"
    implicit none

    integer, intent(out) :: numUnits
    integer, PARAMETER :: NUM_UNITS = 21

    numUnits = NUM_UNITS

    return

  end subroutine setup_getNumFlashUnits

