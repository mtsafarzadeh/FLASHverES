!!****if* source/Simulation/SimulationMain/StirTurb/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities() 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4): solnData, scr

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use Simulation_data, ONLY: sim_velambient, sim_tcc, sim_shiftTime, & 
      sim_rhoblob,sim_blobradius,sim_blobx,sim_bloby, sim_velframe
  use IO_data, ONLY : io_restart, io_statsFileName, io_globalMe
  use Grid_interface, ONLY : Grid_computeUserVars, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getSingleCellVol, Grid_getSingleCellCoords

  implicit none

#include "mpif.h"
#include "constants.h"
#include "Flash.h"
  
  
  real, intent(in) :: simTime
  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer :: size(MDIM)

  integer, parameter ::  nGlobalSum = 28          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k
  real :: dvol, del(MDIM)
  real :: centerCoords(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, DIMENSION(:,:,:,:), POINTER :: scr
  real :: vyblob   ! This is the velocity of the blob

  integer :: point(MDIM)

  ! Make sure the vorticity is up-to-date
  call Grid_computeUserVars()

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
 
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)



     ! Sum contributions from the indicated range of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
              call Grid_getSingleCellCoords(point, blockList(lb), CENTER, EXTERIOR, centerCoords)
    
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
              lsum(8) = lsum(8) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
              lsum(9) = lsum(9) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
              lsum(10) = lsum(10) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &  
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
              lsum(11) = lsum(11) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &  
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol 
              lsum(12) = lsum(12) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &  
                   &                                solnData(ENER_VAR,i,j,k)*dvol
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           
              lsum(13) = lsum(13) + 0.5* solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           
#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
              lsum(14) = lsum(14) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &  
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR
          
              if ((solnData(DENS_VAR,i,j,k).ge.sim_rhoblob/10).and.(abs(centerCoords(IAXIS)).le.3*sim_blobradius))then 
                 lsum(17) = lsum(17) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k)*dvol
                 if (solnData(DENS_VAR,i,j,k).ge.sim_rhoblob/3.) then 
                    lsum(16) = lsum(16) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k)*dvol
                    lsum(18) = lsum(18) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &
                                                 solnData(VELY_VAR,i,j,k)*dvol 
                    lsum(19) = lsum(19) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &
                                                  abs(centerCoords(IAXIS))*dvol
                    lsum(20) = lsum(20) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &
                                                  centerCoords(JAXIS)*dvol
                    lsum(21) = lsum(21) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k) * &
                                                  centerCoords(JAXIS)*centerCoords(JAXIS)*dvol
                    if (solnData(DENS_VAR,i,j,k).ge.sim_rhoblob/1.) then 
                      lsum(15) = lsum(15) + solnData(BLOB_MSCALAR,i,j,k)*solnData(DENS_VAR,i,j,k)*dvol
                    endif
                 endif
              endif 
              lsum(22) = lsum(22) + solnData(MAGX_VAR,i,j,k)*dvol
              lsum(23) = lsum(23) + solnData(MAGY_VAR,i,j,k)*dvol
              lsum(24) = lsum(24) + solnData(MAGZ_VAR,i,j,k)*dvol
              lsum(25) = lsum(25) + solnData(MAGP_VAR,i,j,k)*dvol
              lsum(26) = lsum(26) + solnData(PRES_VAR,i,j,k)*dvol
         enddo
       enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)




  enddo

  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Allreduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MPI_Comm_World, error)



  gsum(18)=gsum(18)/gsum(16)
  vyblob   = gsum(18) 
  gsum(19)=gsum(19)/gsum(16)/sim_blobradius
  gsum(20)=gsum(20)/gsum(16)/sim_blobradius
  gsum(21)=gsum(21)/gsum(16)/sim_blobradius/sim_blobradius
  gsum(21)=sqrt(gsum(21)-gsum(20)*gsum(20))
  sim_blobx = gsum(19)*sim_blobradius
  sim_bloby = gsum(20)*sim_blobradius
  gsum(27) = sim_velframe !sim_velambient
  gsum(28) = simTime/sim_tcc


  ! This shifts the frame periodically to keep the blob on the grid
  if((simTime.ge.3.*sim_tcc).and.(simTime-sim_shiftTime.ge.sim_tcc)) then
  ! if the blob is shifted to the front by more than -0.5 blob radius
  ! and vyblob > 0 so that a velocity shift would move the blob closer
  ! to the front... then skip the shift for the moment
    if((sim_bloby.ge.-0.5).or.(vyblob.le.0)) then
    do lb = 1, count

      !get the index limits of the block
      call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

      ! get a pointer to the current block of data
      call Grid_getBlkPtr(blockList(lb), solnData)
      
      ! Shift the Frame of this block to the curret blob frame 
      do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
              solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) & 
                                      - solnData(VELY_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)/2. 
              solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k)-vyblob*1.2
              solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) &
                                      + solnData(VELY_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)/2. 
           enddo
        enddo
      enddo
      call Grid_releaseBlkPtr(blockList(lb), solnData)
    enddo
!    sim_velambient = sim_velambient-vyblob*1.2
    sim_velframe = sim_velframe+vyblob*1.2
    print *, sim_velambient, sim_velframe
    endif
    sim_shiftTime  = simTime
  endif
     
  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', &
                'mass blob                 ', &
                'x-momentum blob           ', &
                'y-momentum blob           ', & 
                'z-momentum blob           ', &
                'E_total    blob           ', &
                'E_kinetic  blob           ', &
                'E_internal blob           ', & 
                'mass blob (>rho0)         ', &
                'mass blob (>rho0/3.)      ', &
                'mass blob (>rho0/10.)     ', &
                'x-velocty blob (>rho0/10.)', &
                'x-center blob  (>rho0/10.)', &
                'y-center blob  (>rho0/10.)', &
                'y^2 of blob (>rho0/10.)   ', &
                'volume sum of magx        ', &
                'volume sum of magy        ', &
                'volume sum of magz        ', &
                'volume sum of magp        ', &
                'volume sum of pres        ', &
                'velframe                  ', &
                'time/tcc                  '
           
10         format (2x,50(a25, :, 1X))
           
        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



