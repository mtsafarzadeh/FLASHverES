!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkFindList
!!
!! NAME
!!
!!  pt_sinkFindList
!!
!! SYNOPSIS
!!
!!  call pt_sinkFindList(real, intent(IN) :: x,
!!                       real, intent(IN) :: y,
!!                       real, intent(IN) :: z,
!!                       real, intent(IN) :: rad,
!!                       logical, intent(IN) :: create_part,
!!                       integer, dimension(maxsinks), intent(OUT) :: pindex_found,
!!                       integer, intent(OUT) :: np_found)
!!
!! DESCRIPTION
!!
!!  Searches the global and local sink particle lists to find any sink particle within
!!  a given search radius around any given position (x,y,z) and returns the list indexes
!!  of the found particels in pindex_found and the number found in np_found.
!!
!! ARGUMENTS
!!
!!   x - x center of search sphere
!!
!!   y - y center of search sphere
!!
!!   z - z center of search sphere
!!
!!   rad - radius of search sphere
!!
!!   create_part - logical flag indicating whether to create a local dummy
!!                 particle or not, in case the found particle did not
!!                 reside on the local processor and was thus only found in
!!                 the global sink particle list.
!!
!!   pindex_found - list of found particle indexes to be returned by this routine
!!
!!   np_found - number of found particles
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Martin Schroen, 2011
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine pt_sinkFindList(x, y, z, rad, create_part, pindex_found, np_found)

  use Particles_sinkData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

  real, intent(IN)    :: x, y, z, rad
  logical, intent(IN) :: create_part

  integer, dimension(maxsinks), intent(OUT) :: pindex_found
  integer, intent(OUT)                      :: np_found
  character(len=80), save :: grav_boundary_type
  logical, save           :: first_call = .true.
  integer                 :: found, found2, pno, lp
  real                    :: dist, redshift, onePlusRedshift, dx, dy, dz
  real, save              :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz

  if (first_call) then

     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
     if ((grav_boundary_type.ne."isolated").and.(grav_boundary_type.ne."periodic")) then
        call Driver_abortFlash("Sink particles can only be used with perioidic of isolated gravity type!")
     end if
     call RuntimeParameters_get('xmin', xmin)
     call RuntimeParameters_get('xmax', xmax)
     call RuntimeParameters_get('ymin', ymin)
     call RuntimeParameters_get('ymax', ymax)
     call RuntimeParameters_get('zmin', zmin)
     call RuntimeParameters_get('zmax', zmax)
     Lx = xmax-xmin
     Ly = ymax-ymin
     Lz = zmax-zmin
     first_call = .FALSE.

  endif

  call Cosmology_getRedshift(redshift)
  onePlusRedshift = redshift + 1.0

  np_found = 0

  do pno = 1, localnpf

     found = 0

     dx = x - particles_global(ipx,pno)
     dy = y - particles_global(ipy,pno)
     dz = z - particles_global(ipz,pno)
     if (grav_boundary_type.eq."periodic") then !PBC
       if (dx .lt. -0.5*Lx) dx = dx+Lx
       if (dx .gt. +0.5*Lx) dx = dx-Lx
       if (dy .lt. -0.5*Ly) dy = dy+Ly
       if (dy .gt. +0.5*Ly) dy = dy-Ly
       if (dz .lt. -0.5*Lz) dz = dz+Lz
       if (dz .gt. +0.5*Lz) dz = dz-Lz
     endif
     dist = sqrt(dx**2 + dy**2 + dz**2)

     if (dist .le. rad) found = pno

     if (found .gt. 0) then

        np_found = np_found + 1

        found2 = 0
        do lp = 1, localnp
          if (int(particles_local(iptag,lp)) .eq. int(particles_global(iptag,pno))) found2 = lp
        enddo

        if (found2 .gt. 0) then
          pindex_found(np_found) = found2
        else if (create_part) then
          ! it was a particle in the global list, so create a dummy particle if desired
          localnp = localnp + 1
          if (localnp .gt. sinks_maxSinks) &
            call Driver_abortFlash('sink_findList: Sink particle number exceeds sinks_maxSinks. Increase.')
          particles_local(:,localnp) = particles_global(:,found)
          particles_local(ipblk,localnp) = NONEXISTENT
          pindex_found(np_found) = localnp
        else
          ! do this even if found <= localnp (means, particle is already created)
          pindex_found(np_found) = found ! index in global list
        endif

     endif ! found .gt. 0

  end do ! loop over all particles in global list

  return

end subroutine pt_sinkFindList
