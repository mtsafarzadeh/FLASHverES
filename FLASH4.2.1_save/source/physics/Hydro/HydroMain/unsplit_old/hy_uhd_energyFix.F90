!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_energyFix
!!
!! NAME
!!
!!  hy_uhd_energyFix
!!
!! SYNOPSIS
!!
!!  hy_uhd_energyFix( integer (IN) :: blockID,
!!                    integer (IN) :: blkLimits(2,MDIM),
!!                    real(IN)     :: dt,
!!                    real(IN)     :: del(MDIM),
!!                    integer(IN)  :: eosMode)
!!
!! DESCRIPTION
!!
!!  This routine corrects energy in two different ways:
!!  The first choice is to fix the energy due to the differences of
!!  the magnetic pressures using the cell-centered magnetic 
!!  fields and the divergence-free cell face-centered magnetic fields.
!!  This correction is optional, but may be useful for low beta
!!  plasma flows. To enable this first correction during the simulation,
!!  two runtime parameters "hy_killdivb" and "hy_energyFixSwitch" 
!!  should be both turned on in flash.par file.
!!  The second correction is to use the internal energy evolution
!!  to avoid any negativity states of pressure in the Eos routines.
!!
!!  Note that this routine does more than just correcting energy, and
!!  also computes several quantities such as divergence of magnetic
!!  fields, total pressure, current density, and electric fields, etc.
!!  This routine also takes care of abundances.
!!
!!  This implementation is a stub implementation at this level.
!!
!! ARGUMENTS
!!
!!  blockID   - a local block ID
!!  blkLimits - an array that holds the lower and upper indices of the section
!!              of block without the guard cells
!!  dt        - time step
!!  del       - grid deltas in each direction
!!  eosMode   - a mode used in a call to Eos
!!
!!***


Subroutine hy_uhd_energyFix(blockID,blkLimits,dt,del,eosMode)

  implicit none

#include "Flash.h"
#include "constants.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  real, intent(IN) :: dt
  real, dimension(MDIM), intent(IN) :: del
  integer, intent(IN) :: eosMode
  !! -----------------------------------------------------

End Subroutine hy_uhd_energyFix
