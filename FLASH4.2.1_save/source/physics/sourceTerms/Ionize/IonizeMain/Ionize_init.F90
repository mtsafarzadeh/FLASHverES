!!****if* source/physics/sourceTerms/Ionize/IonizeMain/Ionize_init
!!
!! NAME
!!  
!!  Ionize_init
!!
!!
!! SYNOPSIS
!! 
!!  call Ionize_init()
!!
!!  
!! DESCRIPTION
!!
!!  Perform various initializations (apart from the problem-dependent ones)
!!  for the Ionize unit.
!!
!!
!! ARGUMENTS
!!
!!   
!!
!!***

subroutine Ionize_init()
  use ion_interface, ONLY : ion_readTable
  use Ionize_data,ONLY : ion_xfrac,ion_tneimin, ion_tneimax, ion_idx,&
       ion_dneimin, ion_dneimax, ion_smallx, useIonize, &
       ion_symbols, ion_nelect

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
#include "constants.h"
  
  implicit none
   
  

  

  call RuntimeParameters_get("smallx", ion_smallx)  
  call RuntimeParameters_get("tneimin", ion_tneimin)
  call RuntimeParameters_get("tneimax", ion_tneimax)
  call RuntimeParameters_get("dneimin", ion_dneimin)
  call RuntimeParameters_get("dneimax", ion_dneimax)
  call RuntimeParameters_get("useIonize",useIonize)

  !Read ionize coefficients and related data into the module level arrays: 
  !ion_nion12, ion_tp, ion_cfinz, ion_cfric.
  call ion_readTable() 

end subroutine Ionize_init
