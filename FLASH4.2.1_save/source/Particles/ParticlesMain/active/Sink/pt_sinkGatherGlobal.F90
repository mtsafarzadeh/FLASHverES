!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkGatherGlobal
!!
!! NAME
!!
!!  pt_sinkGatherGlobal
!!
!! SYNOPSIS
!!
!!  call pt_sinkGatherGlobal()
!!
!! DESCRIPTION
!!
!!  Fills the global sink particle list with a copy of all sink particles on all
!!  processors. Note that the global sink particle list is not ordered in the same
!!  way on all processors. It always contains the local sink particles in the first
!!  1:localnp indexes after the call is finished.
!!
!! ARGUMENTS
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

subroutine pt_sinkGatherGlobal()

   ! Get particle information from all processors
   ! i.e., update particles_global from particles_local

   use Particles_sinkData, ONLY : particles_local, particles_global, localnp, localnpf, &
                            & recv_buff, send_buff, MAX_MSGS, &
                            & pt_sinkParticleProps, sinks_maxSinks
   use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Driver_interface, ONLY : Driver_abortFlash

   implicit none

   include "Flash_mpi.h"

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

   integer       :: pcount_recv, pcount_send, pmsgcount_recv, pmsgcount_send
   integer       :: localnpt, jproc, i, k, ierr, np_offset, reqr
   integer       :: statr(MPI_STATUS_SIZE)
   logical       :: send_receive
   integer, save :: MyPE, NumPEs
   logical, save :: first_call = .TRUE.

   if (first_call) then

     MyPE = dr_globalMe
     NumPEs = dr_globalNumProcs
     first_call = .false.

   end if

   ! the global list always contains the local particles in the first indexes (1:localnp)
   if (localnp .gt. 0) particles_global(:,1:localnp) = particles_local(:,1:localnp)

   ! Get the total number of particles on each processor and the total number overall
   localnpf = 0

   ! Get the total number of sink particles on all processors
   call MPI_ALLREDUCE(localnp, localnpf, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   
   if (localnpf .gt. sinks_maxSinks) then
     call Driver_abortFlash("pt_sinkGatherGlobal: Number of all sink particles exceeds sinks_maxSinks. Increase.")
   end if

   np_offset = localnp

   !-----------------------------------------------------------------------------
   ! loop over all of the processors.  All the data is moved to local processor 
   ! using MPI sends and receives.
   !-----------------------------------------------------------------------------

   do jproc = 0, NumPEs-1

     if (jproc .ne. MyPE) then

        call MPI_IRECV(localnpt,1,FLASH_INTEGER,jproc,889,MPI_COMM_WORLD,reqr,ierr)
        call MPI_SSEND(localnp ,1,FLASH_INTEGER,jproc,889,MPI_COMM_WORLD,ierr)

        call MPI_WAIT(reqr,statr,ierr)

        ! if there are no particles on this processor
        ! and jproc, do not go through with send/receive
        if ((localnpt .eq. 0) .and. (localnp .eq. 0)) send_receive = .false.

        ! Let's send data; one particle property at a time.
        do k = 1, pt_sinkParticleProps

           pcount_recv = 0
           pcount_send = 0
           send_receive = .true.

           do while (send_receive)

              ! Do not receive more than 12 (=MAX_MSGS) at a time
              pmsgcount_recv = min(localnpt - pcount_recv, MAX_MSGS)
              if (pmsgcount_recv .gt. 0) then

                 call MPI_IRECV(recv_buff, pmsgcount_recv, FLASH_REAL, &
                      jproc, 4711, MPI_COMM_WORLD, reqr, ierr)

              end if

              pmsgcount_send = min(localnp - pcount_send, MAX_MSGS)
              if (pmsgcount_send .gt. 0) then
                 do i = 1, pmsgcount_send
                    send_buff(i) = particles_local(k,pcount_send+i)
                 end do

                 call MPI_SSEND(send_buff, pmsgcount_send, FLASH_REAL, &
                      jproc, 4711, MPI_COMM_WORLD, ierr)

                 pcount_send = pcount_send + pmsgcount_send
              end if

              if (pmsgcount_recv .gt. 0) then

                 call MPI_WAIT(reqr, statr, ierr)
                 do i = 1, pmsgcount_recv
                    particles_global(k,np_offset+pcount_recv+i) = recv_buff(i)
                 end do
                 pcount_recv = pcount_recv + pmsgcount_recv
              end if

              if ((pcount_recv .ge. localnpt) .and. (pcount_send .ge. localnp)) &
                   send_receive = .false.
           end do   ! while send_receive
        end do  ! particle properties

        np_offset = np_offset + localnpt

     end if
   end do

end subroutine pt_sinkGatherGlobal
