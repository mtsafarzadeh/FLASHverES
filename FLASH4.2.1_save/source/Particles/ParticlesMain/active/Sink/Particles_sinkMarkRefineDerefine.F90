!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkMarkRefineDerefine
!!
!! NAME
!!
!!  Particles_sinkMarkRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Particles_sinkMarkRefineDerefine()
!!
!! DESCRIPTION
!!
!!  This routine takes care of grid refinement based on Jeans length and sink particles.
!!  If the local density exceeds a given value that is computed based on Jeans analysis
!!  the block containing that cell is marked for refinement. Refinement and derefinement
!!  are triggered based on the number of cells per Jeans length, which the user must
!!  supply (jeans_ncells_ref and jeans_ncells_deref). Good values for these parameters are
!!  jeans_ncells_ref = 32 and jeans_ncells_deref = 64 (Federrath et al. 2011, ApJ 731, 62),
!!  but the user can choose any real number, where jeans_ncells_ref <= 2*jeans_ncells_deref.
!!  If sink particles are present, they must be at the highest level of AMR, so this routine
!!  also flags all cells within the sink particle accretion radius for refinement to the
!!  highest level.
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
!!   cleaned by Christoph Federrath, 2013
!!
!!***

subroutine Particles_sinkMarkRefineDerefine()

  use RuntimeParameters_interface, only: RuntimeParameters_get

  implicit none
  
  logical, save :: first_call = .true.
  integer, save :: gr_refine_var_thresh
  logical, save :: gr_refineOnSinkParticles
  logical, save :: gr_refineOnJeansLength

  if (first_call) then
     gr_refine_var_thresh = 1
     call RuntimeParameters_get("refineOnJeansLength", gr_refineOnJeansLength)
     call RuntimeParameters_get("refineOnSinkParticles", gr_refineOnSinkParticles)
     first_call = .false.
  end if
  
  !! Sink Particles:
  if (gr_refineOnSinkParticles) call mark_blocks(gr_refine_var_thresh, -1.0, -1, 4)
  !! Jeans Length:
  if (gr_refineOnJeansLength) call mark_blocks(gr_refine_var_thresh, -1.0, -1, 3)
  
  return
  
contains
  
  subroutine mark_blocks(Var, var_th, icmp, input)
    use tree
    use paramesh_dimensions
    use physicaldata, ONLY : unk
    use Grid_data, ONLY : gr_maxRefine
    use Cosmology_interface, ONLY : Cosmology_getRedshift
    use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPhysicalSize, & 
         Grid_getCellCoords, Grid_getBlkIndexLimits, Grid_getBlkCenterCoords
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Driver_data, ONLY : dr_globalMe
    use Grid_interface, ONLY : Grid_getBlkPhysicalSize
    use Particles_sinkData, ONLY : localnp, localnpf, particles_local, ipblk, maxsinks
    use pt_sinkInterface, only: pt_sinkGatherGlobal, pt_sinkFindList
    implicit none

    include "Flash_mpi.h"

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
    ! Arguments

    integer, intent(IN) :: Var
    real,    intent(IN) :: var_th
    integer, intent(IN) :: icmp
    integer, intent(IN) :: input

  ! Local data

    integer :: b, ii, jj, kk, p, p_blknum
    real :: density, redshift

    real, dimension(:), allocatable :: xc, yc, zc
    integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
    integer              :: size_x, size_y, size_z, kp, jp, ip, np_found

    real, dimension(MDIM) :: blockSize
    logical, save :: first_call = .true.
    real, save :: accretion_radius
    real       :: accretion_radius_comoving
    integer, dimension(maxsinks)  :: pindex_found

    real          :: jeans_min(MAXBLOCKS)
    real          :: jeans_min_par(MAXBLOCKS)
    real          :: dens_loc(MAXBLOCKS)
    real          :: dens_max(MAXBLOCKS)
    integer       :: nsend,nrecv, ierr, j, lb
    real          :: cs2, maxd, jeans_numb
    integer       :: reqr(2*MAXBLOCKS),reqs(2*MAXBLOCKS)
    integer       :: stats(MPI_STATUS_SIZE,MAXBLOCKS)
    integer       :: statr(MPI_STATUS_SIZE,maxblocks)
    real, save    :: jeans_ncells_ref, jeans_ncells_deref
    real, save    :: Newton
    real          :: comoving_density, oneplusz, oneplusz_neg2, oneplusz_cu

  !-------------------------------------------------------------------------------

    if(first_call) then

       call RuntimeParameters_get("sink_accretion_radius", accretion_radius)

       call RuntimeParameters_get("jeans_ncells_ref", jeans_ncells_ref)
       call RuntimeParameters_get("jeans_ncells_deref", jeans_ncells_deref)

       call PhysicalConstants_get("Newton", Newton)

       first_call = .false.
    end if


    call Cosmology_getRedshift(redshift)

    oneplusz = redshift + 1.0
    oneplusz_neg2 = oneplusz**(-2.0)
    oneplusz_cu = oneplusz**3.0

    accretion_radius_comoving = accretion_radius * oneplusz


    SELECT CASE(input)

    CASE(1)   !OVERDESNTIY REFINEMENT

    ! not implemented in this version


    CASE(2)    ! UNDERDENSITY REFINEMENT

    ! not implemented in this version


    CASE(3)   ! JEANS LENGTH REFINEMENT


       do lb = 1, lnblocks

          jeans_min(lb) = 1.0e99
          dens_loc(lb)  = 1.0e-50
          dens_max(lb)  = 1.0e-50

          if (nodetype(lb) .eq. 1 .or. nodetype(lb) .eq. 2) then

             call Grid_getBlkPhysicalSize(lb, blockSize)

             maxd = blockSize(1) / real(NXB)
             maxd = max(maxd, blockSize(2) / real(NYB))
             maxd = max(maxd, blockSize(3) / real(NZB))

             maxd = maxd / oneplusz

             do kk = NGUARD*K3D+1,NGUARD*K3D+NZB
                do jj = NGUARD*K2D+1,NGUARD*K2D+NYB
                   do ii = NGUARD+1,NGUARD+NXB

                      comoving_density = unk(DENS_VAR,ii,jj,kk, lb)

                      cs2 = unk(PRES_VAR,ii,jj,kk, lb) / comoving_density
                      cs2 = cs2 * oneplusz_neg2
                      density = comoving_density * oneplusz_cu

                      jeans_numb = sqrt(PI*cs2 / Newton / density) / maxd

                      if (jeans_numb .lt. jeans_min(lb)) then
                         jeans_min(lb) = jeans_numb
                         dens_loc(lb) = density
                      end if

                   enddo
                enddo
             enddo

          endif !type of block

       enddo ! blocks

       ! communicate error of parent to children

       jeans_min_par(1:lnblocks) = 0.
       nrecv = 0
       do lb = 1,lnblocks
          if (parent(1,lb).gt.-1) then
             if (parent(2,lb).ne.dr_globalMe) then
                nrecv = nrecv + 1
                call MPI_IRecv(jeans_min_par(lb),1, &
                     MPI_DOUBLE_PRECISION, &
                     parent(2,lb), &
                     lb, &
                     MPI_COMM_WORLD, &
                     reqr(nrecv), &
                     ierr)
             else
                jeans_min_par(lb) = jeans_min(parent(1,lb))
             end if
          end if
       end do


       ! parents send error to children

       nsend = 0
       do lb = 1,lnblocks
          do j = 1,nchild
             if (child(1,j,lb).gt.-1) then
                if (child(2,j,lb).ne.dr_globalMe) then
                   nsend = nsend + 1
                   call MPI_ISend(jeans_min(lb), &
                        1, &
                        MPI_DOUBLE_PRECISION, &
                        child(2,j,lb), &  ! PE TO SEND TO
                        child(1,j,lb), &  ! THIS IS THE TAG
                        MPI_COMM_WORLD, &
                        reqs(nsend), &
                        ierr)
                end if
             end if
          end do
       end do


       if (nsend.gt.0) then
          call MPI_Waitall (nsend, reqs, stats, ierr)
       end if
       if (nrecv.gt.0) then
          call MPI_Waitall (nrecv, reqr, statr, ierr)
       end if


       ! label blocks for refinement

       do lb = 1, lnblocks

          if (nodetype(lb) .eq. 1) then

                ! refinement

                if (jeans_min(lb) .lt. jeans_ncells_ref) then
                   derefine(lb) = .false.
                   refine(lb) = .true.
                end if

                if (lrefine(lb) .ge. lrefine_max) refine(lb) = .false.

          end if       ! leaf blocks
       end do        ! blocks


       ! label blocks for derefinement

       do lb = 1, lnblocks

          if (nodetype(lb) .eq. 1) then

                if (.not. refine(lb) .and. .not. stay(lb) & 
                     .and. jeans_min(lb) .gt. jeans_ncells_deref &
                     .and. jeans_min_par(lb) .gt. jeans_ncells_deref) then
                   derefine(lb) = .true.
                else
                   derefine(lb) = .false.
                end if

                if (lrefine(lb) .ge. lrefine_max) refine(lb) = .false.

          end if          ! leaf blocks
       end do            ! blocks


    CASE(4)   ! SINK PARTICLE REFINEMENT

       ! Any block with a sink particle in it should be at highest refinement level
       do p = 1, localnp
          p_blknum = particles_local(ipblk, p)
          if (lrefine(p_blknum) .lt. gr_maxRefine) then
             refine(p_blknum) = .true.
             derefine(p_blknum) = .false.
             stay(p_blknum) = .true.

          end if
          if (lrefine(p_blknum) .eq. gr_maxRefine) then
             derefine(p_blknum) = .false.
             stay(p_blknum) = .true.
          end if
       end do

       ! update particles_global array
       call pt_sinkGatherGlobal()

       ! Any cell within accretion_radius of sink particle should be at the
       ! highest refinement level (its block, to be precise)

       do b = 1, lnblocks

          if (nodetype(b).eq.1) then

             ! find cells (including GCs) within sink particle accretion radius

             call Grid_getBlkIndexLimits(b, blkLimits, blkLimitsGC)
             size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
             size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
             size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

             allocate(xc(size_x))
             allocate(yc(size_y))
             allocate(zc(size_z))

             call Grid_getCellCoords(IAXIS, b, CENTER, .true., xc, size_x)
             call Grid_getCellCoords(JAXIS, b, CENTER, .true., yc, size_y)
             call Grid_getCellCoords(KAXIS, b, CENTER, .true., zc, size_z)

             do kp = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                do jp = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
                   do ip = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

                      ! cell within accretion radius?

                      call pt_sinkFindList(xc(ip), yc(jp), zc(kp), accretion_radius_comoving, &
                          & .FALSE., pindex_found, np_found)
                      if (np_found .gt. 0) then
                         if (lrefine(b) .lt. gr_maxRefine) then
                            refine   (b) = .TRUE.
                            derefine (b) = .FALSE.
                            stay     (b) = .TRUE.

                         endif
                         if (lrefine(b) .eq. gr_maxRefine) then
                            derefine (b) = .FALSE.
                            stay     (b) = .TRUE.
                         endif
                      end if

                   end do
                end do
             end do

             deallocate(xc)
             deallocate(yc)
             deallocate(zc)

          end if      ! nodetype

       end do         ! loop over blocks

    END SELECT

    return

  end subroutine mark_blocks

end subroutine Particles_sinkMarkRefineDerefine
