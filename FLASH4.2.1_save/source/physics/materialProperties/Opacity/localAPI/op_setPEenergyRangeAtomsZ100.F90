!!****if* source/physics/materialProperties/Opacity/localAPI/op_setPEenergyRangeAtomsZ100
!!
!! NAME
!!
!!  op_setPEenergyRangeAtomsZ100
!!
!! SYNOPSIS
!!
!!  call op_setPEenergyRangeAtomsZ100 ()
!!
!! DESCRIPTION
!!
!!  This routine sets the energy range (boundaries) for determining the photoelectric
!!  cross sections according to the F.Biggs and R.Lighthill report:
!!
!!       Analytical Approximations for X-Ray Cross Sections II
!!       Frank Biggs and Ruth Lighthill
!!       Weapons Effects Research Department
!!       Sandia Laboratories, December 1971
!!
!!  The atomic elements set in this routine are in the range: 74 =< Z =< 100.
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEenergyRangeAtomsZ100 ()

  implicit none

  return
end subroutine op_setPEenergyRangeAtomsZ100
