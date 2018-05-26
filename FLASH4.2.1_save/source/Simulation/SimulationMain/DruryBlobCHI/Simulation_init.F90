!!****if* source/Simulation/SimulationMain/magnetoHD/BrioWu/Simulation_init
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
!! ARGUMENTS
!!
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for BrioWu problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

  real :: kx, ky, kz, twopi, L, A1, phi, phi2, phi3
  integer :: ikxmin, ikymin, ikzmin, ikxmax, ikymax,ikzmax, ikx, iky, ikz
  integer :: i, j, k


#include "constants.h"
#include "Flash.h"
  
!  call PhysicalConstants_get("Boltzmann", sim_boltz)
sim_boltz=1.38e-16
!  call PhysicalConstants_get("amu", sim_amu)
sim_amu= 1.67e-24
  call RuntimeParameters_get('smallx',    sim_smallx)
  call RuntimeParameters_get('smallp',    sim_smallP)
  call RuntimeParameters_get('posn',      sim_posn)
  call RuntimeParameters_get('gamma',     sim_gamma)
  call RuntimeParameters_get('rho_cloud',  sim_rhocloud)
  call RuntimeParameters_get('fieldbeta',      sim_beta)
  call RuntimeParameters_get('chi',      sim_chi)
  call RuntimeParameters_get('killdivb',  sim_killdivb)
  call RuntimeParameters_get('xmin',      sim_xmin)
  call RuntimeParameters_get('xmax',      sim_xmax)
  call RuntimeParameters_get('ymin',      sim_ymin)
  call RuntimeParameters_get('ymax',      sim_ymax)
  call RuntimeParameters_get('zmin',      sim_zmin)
  call RuntimeParameters_get('zmax',      sim_zmax)
  call RuntimeParameters_get('diff_time', sim_diff_time)

  call RuntimeParameters_get('kxmin', ikxmin)
  call RuntimeParameters_get('kxmax', ikxmax)
  call RuntimeParameters_get('kymin', ikymin)
  call RuntimeParameters_get('kymax', ikymax)
  call RuntimeParameters_get('kzmin', ikzmin)
  call RuntimeParameters_get('kzmax', ikzmax)

  call RuntimeParameters_get('igrid', sim_igrid)
  call RuntimeParameters_get('jgrid', sim_jgrid)
  call RuntimeParameters_get('kgrid', sim_kgrid)

  sim_pamb = sim_rhocloud*sim_boltz*1.e4/0.6/sim_amu
  sim_rhoamb = sim_rhocloud/sim_chi

print *, sim_amu, sim_pamb,sim_rhoamb

 sim_gCell = .true.

  if (NDIM == 1) then
     sim_killdivb = .false.
  endif

     if (NDIM == 2) then
        ikzmin=1
        ikzmax=1
     endif

     L = (sim_xmax - sim_xmin)

     twopi = 2.*PI

     randomgrid = 0.
     if (NDIM .eq. 2) then
     do i=1,sim_igrid
     do j=1,sim_jgrid
      posx(i,j)=(sim_posn-sim_xmin)/float(sim_igrid)*real(i)
      posy(i,j)=sim_ymin+(sim_ymax-sim_ymin)/float(sim_jgrid)*real(j)
     enddo
     enddo
     endif

     do ikx = ikxmin, ikxmax,2
        kx = twopi*ikx / L

        do iky = ikymin, ikymax,2
           ky = twopi*iky / L

              call random_number(A1)
              call random_number(phi)
              call random_number(phi2)

     do i=1,sim_igrid
     do j=1,sim_jgrid
              randomgrid(i,j)=randomgrid(i,j)+A1*sin(kx*posx(i,j)+twopi*phi)*sin(ky*posy(i,j)+twopi*phi2)
        enddo
     enddo

     enddo
     enddo

     randomgrid=exp(2.*randomgrid/maxval(randomgrid))

end subroutine Simulation_init
