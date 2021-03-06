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
      real, save, dimension(241,161) :: hhe_cooling,  hhe_heating
      real, save, dimension(241,161) :: total_cooling,total_heating
      real, save :: tradmin, tradmax, dradmin, dradmax, massfracH
!     set to 0 to turn on cooling but not use it
      real, save :: nocool
!     Initial dust fraction in solar units
      real, save :: dustinit
!     Metallicity in solar units
      real, save ::  metallicity
!     UV flux as in the paper
      real, save :: G0
!     Initial dust grain size in microns
      real, save :: A0
!     Initial dust grain density in cgs units
      real, save :: graindensity
!     Initial dust grain mass in grams 
      real, save :: mdust0
!     Coefficients for dust destruction by thermal sputtering 
      real, save :: Adust,Bdust
!     Number of ionizing photons per second generated by the source
!     in units of 1E50 per second
      real, save :: Ngamma
!     Distance to the source
      real, save :: sourcedist 
!     Ngamma/4 pi / distance to source ^2 * 6.3E-18 cm^2 cross section
!     for hydrogen photionization 
      real, save :: Gammapi
!     if this is nonzero then we use selfshielding 
      logical, save :: selfshield
!     The size of a self-shielding blob in cm 
      real, save :: shieldlength
!     Compton luminosity in units of solar
      real, save :: Comptonlum
!     Compton flux in units ergs/s/cm2
      real, save :: Comptonflux
!     Compton temp in kev
      real, save :: Comptontemp
!     Compton temp in Kelvin  4 k T = <hv>
      real, save :: ComptontempK
end module Cool_data
