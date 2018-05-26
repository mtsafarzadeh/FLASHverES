!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!
!! 
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Eos_data, ONLY: eos_singleSpeciesA
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none
  integer :: i 
  real :: redshift
  real :: Omegab, Omegam, hubble
  real :: rhocrit,mcheck
  real :: integral,integral2,t,dt

#include "Flash.h"
#include "constants.h"



  character(len=MAX_STRING_LENGTH) :: str_geometry
  integer :: meshGeom

  ! this is virial mass in solar masses
  call RuntimeParameters_get("mvir",   sim_mvir)
  ! this is the nfw concentration parameter
  call RuntimeParameters_get("nfwc",   sim_nfwc)
  call RuntimeParameters_get("redshift",redshift)
  ! this is the energy of the blast in ergs
  call RuntimeParameters_get("Eblast",sim_Eblast)
  ! this is radius of the blast in cm
  call RuntimeParameters_get("rblast",sim_rblast)
  call RuntimeParameters_get("geometry", str_geometry)
  call RuntimeParameters_mapStrToInt(str_geometry, meshGeom)
  ! cosmological parameters to figure out the virial density
  call RuntimeParameters_get("Omegam",Omegam)
  call RuntimeParameters_get("Omegab",Omegab)
  call RuntimeParameters_get("hubble",hubble)
  if (meshGeom /= SPHERICAL .AND. NDIM /= 2) then
     call Driver_abortFlash("ERROR: invalid geometry for SodSpherical")
  endif

  ! compute the critical density and check if the masses match
  rhocrit  = 1.878E-29*hubble*hubble*(Omegam*(1.+redshift)**3.+(1.-Omegam))
  sim_rvir = (3./(4.*3.1415)*sim_mvir*1.98E3/(180*rhocrit))**(1./3.)*1E10

  ! Tvir =  G \mu mp/r_vir/2/kb
  sim_tvir = 6.67E-8*(1.6726219E-24*eos_singleSpeciesA)*sim_mvir*1.98E33&
                  /sim_rvir/2./1.38E-16

  ! this is F(c)
  sim_Fc   = alog(1+sim_nfwc)-sim_nfwc/(1+sim_nfwc)
  sim_Ac = 2.*sim_nfwc/sim_Fc

  
  ! compute the critical density and check if the masses match
  mcheck   = rhocrit*(sim_rvir/1E10)**3.*180*4.*3.1415/3./1.98E3
  print *,'Mass ',sim_mvir,' = ',mcheck
  print *,'Tvir ',sim_tvir
  print *,'rvir ',sim_rvir/3.086E21
  print *,'c',    sim_nfwc
  print *,'Fc',   sim_Fc
  print *,'Ac',   sim_Ac

  ! compute the central baryonic density
  integral = 0.
  dt = sim_nfwc/5000.
  do i = 0,5000
   t = dt*(i+0.5)
   integral = integral+(1+t)**(sim_Ac/t)*t*t*dt
  enddo
  print *,'integral', integral
  sim_rho0 = rhocrit*(180/3.)*Omegab/Omegam*sim_nfwc**3*exp(sim_Ac)/integral
  print *,'rho0',sim_rho0

  ! we put the explosion at rblast the gas mass there is
  integral2 = 0.
  dt = sim_nfwc/5000*(sim_rblast/sim_rvir)
  do i = 0,5000
   t = dt*(i+0.5)
   integral2 = integral2+(1+t)**(sim_Ac/t)*t*t*dt
  enddo
  sim_mblast = Omegab/Omegam*sim_mvir*integral2/integral 
  print *,'rblast',sim_rblast
  print *,'gas fraction',integral2/integral,sim_mblast
  
 

  call Logfile_stamp( "initializing SZSpherical problem",  &
       "[Simulation_init]")
     
end subroutine Simulation_init







