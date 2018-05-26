!!****f* source/physics/Hydro/Hydro_recalibrateEints
!! NAME
!!
!!  Hydro_recalibrateEints
!! 
!! SYNOPSIS
!!
!!  call Hydro_recalibrateEints(
!!                     integer(IN) :: range(HIGH, MDIM),
!!                     integer(IN) :: blockID)
!!
!! DESCRIPTION
!!
!! This function recalibrates multiTemp component internal energies and energies
!! so that their sums agree with the overall values for eint.
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which Eos is to be applies. The example shows how
!!          the array describes the block section.
!!
!!   blockID: current block number
!!
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then recalibration is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then recalibration is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface should be defined in Fortran Module 
!!      Hydro_interface. All functions calling this routine should include
!!      a statement like
!!      use Hydro_interface, ONLY : Hydro_recalibrateEints
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices cannot
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos_wrapped
!!     Eos
!!     Eos.h
!!
!!***

subroutine Hydro_recalibrateEints(range,blockID)

  implicit none

#include "constants.h"

  integer, dimension(2,MDIM), intent(in) :: range
  integer,intent(in) :: blockID

  return
end subroutine Hydro_recalibrateEints

subroutine Hydro_recalibrateEintsForCell(eint,eion,eele,erad,e1,e2,e3)
  implicit none
  real,intent(in)    :: eint
  real,intent(INOUT) :: eion,eele
  real,intent(INOUT),OPTIONAL :: erad
  real,intent(INOUT),OPTIONAL :: e1,e2,e3

end subroutine Hydro_recalibrateEintsForCell
