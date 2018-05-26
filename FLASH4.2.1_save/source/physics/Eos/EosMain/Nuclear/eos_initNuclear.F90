!!****if* source/physics/Eos/EosMain/Nuclear/eos_initNuclear
!!
!! NAME
!!
!!  eos_initNuclear
!!
!! SYNOPSIS
!!  
!!  subroutine eos_initNuclear()
!!                 
!!
!! DESCRIPTION
!!
!!  Initialization for the Nuclear EOS appropriate for 
!!  core-collapse supernova simulations.  
!!
!! ARGUMENTS
!!
!! NOTES
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in the kernel 
!!      directory are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  
!!      stellarcollapse.org/equationofstate.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M. 2013, ApJ, 765, 29
!!
!!
!!***

#include "Eos.h"
subroutine eos_initNuclear()

  use Eos_data, ONLY : eos_type
  use eosmodule
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  eos_type=EOS_NUC

  call PhysicalConstants_get("Avogadro", avo)

  call RuntimeParameters_get('eos_file', eos_file)

  call readtable(eos_file)
  e_zeroPoint = energy_shift

  return
end subroutine eos_initNuclear
