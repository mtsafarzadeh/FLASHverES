!!****f* source/source_terms/Cool/CoolMain/Cool_init
!!
!! NAME
!!  
!!  Cool_init
!!
!!
!! SYNOPSIS
!! 
!!  call Cool_init(integer(IN) :: myPE)
!!  
!! DESCRIPTION
!!
!!  Apply the radloss source term operator to a block of zones.
!!  The radiative losses rate is used to update the
!!  internal energy in the zone. The radiative losses from an 
!!  optically thin plasma are from Sarazin (1986) Rev.Mod.Phys.
!!  and from Raymond (1976).
!!  (see also Peres et al. 1982, ApJ 252, 791)
!!
!!  After we call radloss, call the eos to update the pressure 
!!  and temperature based on the radiative losses.
!!
!!
!!***

subroutine Cool_init (myPE)
      use Cool_data
      use RuntimeParameters_interface, ONLY : RuntimeParameters_get
      use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

#include "Flash.h"

      implicit none
      integer, intent(IN) :: myPE
      integer :: i

!
!==============================================================================
!
      call PhysicalConstants_get("proton mass", amu)
      call RuntimeParameters_get("tradmin", tradmin)      
      call RuntimeParameters_get("tradmax", tradmax)
      call RuntimeParameters_get("dradmin", dradmin)
      call RuntimeParameters_get("dradmax", dradmax)
      call RuntimeParameters_get("massfracH", massfracH)
      call RuntimeParameters_get("nocool", nocool)
	
      open(unit=2,file="./Total_metals_cooling.dat")
      open(unit=3,file="./Metal_free_cooling.dat")
      do i= 1,352
       read(2,*) metal_cooling(i)
       read(3,*) hhe_cooling(i)
      end do 
      close(2)
      close(3)
      return
end subroutine Cool_init
