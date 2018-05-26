!!****f* source/source_terms/Cool/CoolMain/Cool_data
!!
!! NAME
!!  
!!  Cool_data
!!
!!
!! SYNOPSIS
!! 
!!  use Cool_data
!!  
!! DESCRIPTION
!!
!! Store the data for the Cooling unit
!!
!! ARUMENTS
!! tradmin
!! tradmax
!! dradmin
!! dradmax
!!***

module Cool_data
      implicit none
      real, save :: amu
      real, save, dimension(353) :: metal_cooling,hhe_cooling
      real, save :: tradmin, tradmax, dradmin, dradmax, massfracH,nocool
end module Cool_data	
