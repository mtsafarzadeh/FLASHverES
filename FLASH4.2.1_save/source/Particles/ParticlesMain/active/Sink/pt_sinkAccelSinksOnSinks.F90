!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkAccelSinksOnSinks
!!
!! NAME
!!
!!  pt_sinkAccelSinksOnSinks
!!
!! SYNOPSIS
!!
!!  call pt_sinkAccelSinksOnSinks(real(out) :: local_min_radius,
!!                                real(out) :: local_max_accel)
!!
!! DESCRIPTION
!!
!!  Computes sinks <-> sinks gravitational accelerations.
!!  Compute the acceleration of the sink particles from other sink particles
!!  and add it to preexisting acceleration array. Send out local_min_radius
!!  and local_max_accel to constrain subcycle time step in Particles_sinkAdvanceParticles.
!!
!! ARGUMENTS
!!
!!   local_min_radius - minimum radius between sinks
!!
!!   local_max_accel - maximum acceleration between sinks
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged by Christoph Federrath, 2013
!!
!!***

subroutine pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)

    use Particles_sinkData
    use pt_sinkInterface, only: pt_sinkGatherGlobal
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Driver_data, ONLY : dr_globalMe
    use Cosmology_interface, ONLY : Cosmology_getRedshift
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

    real, intent(out) :: local_min_radius, local_max_accel

    character(len=80), save :: softening_type_sinks, grav_boundary_type
    logical, save :: first_call = .true.

    real, save    :: newton, maxradius_pbc, softening_radius
    real          :: softening_radius_comoving, slope, hinv, h2inv
    integer, save :: softeningtype
    real, save    :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz, local_min_radius_init

    integer       :: p, pf, n, nx, ny, nz
    real          :: xp, yp, zp, dx, dy, dz, dx_inside, dy_inside, dz_inside
    real          :: radius, masspf, q, kernelvalue, r3, paccx, paccy, paccz
    real, save    :: LxPBC(-nrep_pbc:nrep_pbc), LyPBC(-nrep_pbc:nrep_pbc), LzPBC(-nrep_pbc:nrep_pbc)

    real :: redshift, oneplusz3

    if (first_call) then

       call RuntimeParameters_get("sink_softening_radius", softening_radius)
       call RuntimeParameters_get("sink_softening_type_sinks", softening_type_sinks)
       select case (softening_type_sinks)
       case ("spline")
          softeningtype=1
       case ("linear")
          softeningtype=2
       case default
          softening_type_sinks = "spline"
          softeningtype = 2
          if(dr_globalMe .eq. MASTER_PE) print*, "pt_sinkAccelSinksOnSinks: invalid grav softening type specified"
       end select
       if(dr_globalMe .eq. MASTER_PE) print*, "pt_sinkAccelSinksOnSinks: grav softening type = ", &
            & trim(softening_type_sinks)

       call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

       if ((grav_boundary_type.ne."isolated").and.(grav_boundary_type.ne."periodic")) then
          call Driver_abortFlash("Sink particles can only be used with perioidic of isolated gravity type!")
       end if

       call RuntimeParameters_get("xmin", xmin)
       call RuntimeParameters_get("xmax", xmax)
       call RuntimeParameters_get("ymin", ymin)
       call RuntimeParameters_get("ymax", ymax)
       call RuntimeParameters_get("zmin", zmin)
       call RuntimeParameters_get("zmax", zmax)

       Lx = xmax-xmin
       Ly = ymax-ymin
       Lz = zmax-zmin

       if (grav_boundary_type.eq."periodic") then
          maxradius_pbc = real(nrep_pbc)*min(Lx,Ly,Lz)
          ! prepare container for faster access in case of periodic BCs
          do n = -nrep_pbc, nrep_pbc
             LxPBC(n) = n*Lx
             LyPBC(n) = n*Ly
             LzPBC(n) = n*Lz
          enddo
       endif

       local_min_radius_init = sqrt(Lx**2+Ly**2+Lz**2)

       call PhysicalConstants_get("Newton", newton)

       first_call = .false.

    end if

    local_min_radius = local_min_radius_init
    local_max_accel = 0.0

    if (localnpf .eq. 0) return

    call Cosmology_getRedshift(redshift)
    softening_radius_comoving = softening_radius * (1.0 + redshift)

    hinv = 2.0 / softening_radius_comoving
    h2inv = hinv**2
    slope = 1.0 / softening_radius_comoving**3

    oneplusz3 = (1.0 + redshift) ** 3.0

    call pt_sinkGatherGlobal()

    if (grav_boundary_type .eq. "isolated") then

       do p = 1, localnp

          paccx = 0.0
          paccy = 0.0
          paccz = 0.0

          xp = particles_global(POSX_PART_PROP,p)
          yp = particles_global(POSY_PART_PROP,p)
          zp = particles_global(POSZ_PART_PROP,p)

          do pf = 1, localnpf

             if (pf .ne. p) then

                masspf = particles_global(MASS_PART_PROP,pf)
                dx = xp - particles_global(POSX_PART_PROP,pf)
                dy = yp - particles_global(POSY_PART_PROP,pf)
                dz = zp - particles_global(POSZ_PART_PROP,pf)
                radius = sqrt(dx*dx+dy*dy+dz*dz)
                local_min_radius = min(local_min_radius, radius)

                if (radius .lt. softening_radius_comoving) then
                   if (softeningtype .eq. 1) then   ! spline softening
                      q = radius*hinv
                      if ((q.gt.1.e-5).and.(q.lt.1.0)) &
                         & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                      if ((q.ge.1.0)  .and.(q.lt.2.0)) &
                         & kernelvalue = h2inv * &
                         & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                      paccx = paccx + masspf*dx*kernelvalue
                      paccy = paccy + masspf*dy*kernelvalue
                      paccz = paccz + masspf*dz*kernelvalue
                   end if
                   if (softeningtype.eq.2) then    ! linear kernel
                      paccx = paccx + masspf*dx*slope
                      paccy = paccy + masspf*dy*slope
                      paccz = paccz + masspf*dz*slope
                   end if
                else    ! Newtonian gravity of point mass
                   r3 = 1.0 / radius**3
                   paccx = paccx + masspf*dx*r3
                   paccy = paccy + masspf*dy*r3
                   paccz = paccz + masspf*dz*r3
                end if

             end if ! pf .ne. p

          end do   ! pf

          particles_local(ACCX_PART_PROP,p) = particles_local(ACCX_PART_PROP,p) - newton*paccx*oneplusz3
          particles_local(ACCY_PART_PROP,p) = particles_local(ACCY_PART_PROP,p) - newton*paccy*oneplusz3
          particles_local(ACCZ_PART_PROP,p) = particles_local(ACCZ_PART_PROP,p) - newton*paccz*oneplusz3

          local_max_accel = max(local_max_accel, abs(particles_local(ACCX_PART_PROP,p)))
          local_max_accel = max(local_max_accel, abs(particles_local(ACCY_PART_PROP,p)))
          local_max_accel = max(local_max_accel, abs(particles_local(ACCZ_PART_PROP,p)))

       end do   ! local particles p

    end if

    if (grav_boundary_type .eq. "periodic") then

       do p = 1, localnp

          paccx = 0.0
          paccy = 0.0
          paccz = 0.0

          xp = particles_global(POSX_PART_PROP,p)
          yp = particles_global(POSY_PART_PROP,p)
          zp = particles_global(POSZ_PART_PROP,p)

          do pf = 1, localnpf

            if (pf .ne. p) then

             masspf = particles_global(MASS_PART_PROP,pf)
             dx_inside = xp - particles_global(POSX_PART_PROP,pf)
             dy_inside = yp - particles_global(POSY_PART_PROP,pf)
             dz_inside = zp - particles_global(POSZ_PART_PROP,pf)

             do nx = -nrep_pbc, nrep_pbc
                dx = dx_inside + LxPBC(nx)
                do ny = -nrep_pbc, nrep_pbc
                   dy = dy_inside + LyPBC(ny)
                   do nz = -nrep_pbc, nrep_pbc
                      dz = dz_inside + LzPBC(nz)

                      radius = sqrt(dx*dx+dy*dy+dz*dz)

                      if ((radius .gt. 0) .and. (radius .lt. maxradius_pbc)) then

                         local_min_radius = min(local_min_radius, radius)
                         if (radius .lt. softening_radius) then
                            if (softeningtype .eq. 1) then    ! spline softening

                               q = radius*hinv
                               if ((q.gt.1.e-5).and.(q.lt.1.0)) &
                                 & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                               if ((q.ge.1.0)  .and.(q.lt.2.0)) &
                                 & kernelvalue = h2inv * &
                                 & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                               paccx = paccx + masspf*dx*kernelvalue
                               paccy = paccy + masspf*dy*kernelvalue
                               paccz = paccz + masspf*dz*kernelvalue    
                            end if
                            if (softeningtype.eq.2) then ! linear kernel inside smoothing radius
                               paccx = paccx + masspf*dx*slope
                               paccy = paccy + masspf*dy*slope
                               paccz = paccz + masspf*dz*slope
                            endif
                         else  ! Newtonian graity of point mass
                            r3 = 1.0/radius**3
                            paccx = paccx + masspf*dx*r3
                            paccy = paccy + masspf*dy*r3
                            paccz = paccz + masspf*dz*r3
                         end if

                      end if  ! within maxradius_pbc

                   end do  ! nz
                end do   ! ny
             end do   ! nx

            end if ! pf .ne. p

          end do  ! pf

          particles_local(ACCX_PART_PROP,p) = particles_local(ACCX_PART_PROP,p) - newton*paccx*oneplusz3
          particles_local(ACCY_PART_PROP,p) = particles_local(ACCY_PART_PROP,p) - newton*paccy*oneplusz3
          particles_local(ACCZ_PART_PROP,p) = particles_local(ACCZ_PART_PROP,p) - newton*paccz*oneplusz3

          local_max_accel = max(local_max_accel, abs(particles_local(ACCX_PART_PROP,p)))
          local_max_accel = max(local_max_accel, abs(particles_local(ACCY_PART_PROP,p)))
          local_max_accel = max(local_max_accel, abs(particles_local(ACCZ_PART_PROP,p)))

       end do ! loop over local particles p

    end if    ! gravity boundary

    return

end subroutine pt_sinkAccelSinksOnSinks
