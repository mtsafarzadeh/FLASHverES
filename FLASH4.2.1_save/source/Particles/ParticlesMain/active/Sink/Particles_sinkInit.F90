!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkInit
!!
!! NAME
!!
!!  Particles_sinkInit
!!
!! SYNOPSIS
!!
!!  call Particles_sinkInit(logical, INTENT(in) :: restart)
!!
!! DESCRIPTION
!!
!!  Initializes sink particle properties.
!!
!! ARGUMENTS
!!
!!   restart - logical flag indicating restart or not
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine Particles_sinkInit(restart)

   ! initialize the sink particles arrays

   use Driver_data, ONLY : dr_globalMe
   use Particles_sinkData
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
        RuntimeParameters_mapStrToInt
   use Driver_interface, ONLY : Driver_abortFlash

   implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

   logical, INTENT(in) :: restart

   integer :: ierr

   call RuntimeParameters_get ("useSinkParticles", useSinkParticles)
   call RuntimeParameters_get ("sinks_maxSinks", sinks_maxSinks)

   if (.not. useSinkParticles) return

   allocate (NumParticlesPerBlock(MAXBLOCKS), stat=ierr)
   NumParticlesPerBlock(:) = 0

   if (dr_globalMe .eq. MASTER_PE) print*, "Particles_sinkInit: Initializing Sink Particles."

   ! Sink particle Properties (see Flash.h)
   ipx = POSX_PART_PROP
   ipy = POSY_PART_PROP
   ipz = POSZ_PART_PROP
   ipvx = VELX_PART_PROP
   ipvy = VELY_PART_PROP
   ipvz = VELZ_PART_PROP
   ipm = MASS_PART_PROP

   ipblk = BLK_PART_PROP
   iptag = TAG_PART_PROP
   ipcpu = PROC_PART_PROP

   iplx = X_ANG_PART_PROP
   iply = Y_ANG_PART_PROP
   iplz = Z_ANG_PART_PROP
   iplx_old = X_ANG_OLD_PART_PROP
   iply_old = Y_ANG_OLD_PART_PROP
   iplz_old = Z_ANG_OLD_PART_PROP
   ipt = CREATION_TIME_PART_PROP
   ipmdot = ACCR_RATE_PART_PROP 
   ipraccr = ACCR_RADIUS_PART_PROP
   iold_pmass = OLD_PMASS_PART_PROP
   ipmgas = MGAS_PART_PROP
   ipdtold = DTOLD_PART_PROP
   ipbflx = X_BFLUX_PART_PROP
   ipbfly = Y_BFLUX_PART_PROP
   ipbflz = Z_BFLUX_PART_PROP

   n_empty = sinks_maxSinks
   RunningParticles = .true.
   if (.not. restart) then
      localnp = 0
      localnpf = 0
   end if

   if (sinks_maxSinks .gt. maxsinks) then
      call Driver_abortFlash("Particles_sinkInit: sinks_maxSinks > maxsinks. Must increase maxsinks in Particles_sinkData.")
   endif

   if (.not. restart) then !if we starting from scratch

      if (.not. allocated(particles_local)) then
         allocate (particles_local(pt_sinkParticleProps, sinks_maxSinks), stat=ierr)
      endif
      if (ierr /= 0) then
         call Driver_abortFlash("Particles_sinkInit:  could not allocate particles_local array")
      endif

      if (.not. allocated(particles_global)) &
           allocate (particles_global(pt_sinkParticleProps, sinks_maxSinks), stat=ierr)
      if (ierr /= 0) then
         call Driver_abortFlash("Particles_sinkInit:  could not allocate particles_global array for sink particles")
      endif

      if (ierr /= 0) then
         call Driver_abortFlash("Particles_sinkInit:  could not allocate particles_local array")
      endif

      particles_local = NONEXISTENT
      particles_global = NONEXISTENT

   end if  ! end of .not. restart

   if (allocated(particles_local)) particles_local(VELX_PART_PROP,:)=0.0
   if (allocated(particles_global)) particles_global(VELX_PART_PROP,:)=0.0

   if (allocated(particles_local)) particles_local(VELY_PART_PROP,:)=0.0
   if (allocated(particles_global)) particles_global(VELY_PART_PROP,:)=0.0

   if (allocated(particles_local)) particles_local(VELZ_PART_PROP,:)=0.0
   if (allocated(particles_global)) particles_global(VELZ_PART_PROP,:)=0.0

   return

end subroutine Particles_sinkInit
