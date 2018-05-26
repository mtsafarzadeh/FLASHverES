!!****if* source/Simulation/SimulationMain/StirTurbHydro/IO_writeIntegralQuantities
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

subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalMe
  use Grid_interface, ONLY : Grid_computeUserVars, Grid_getCellCoords, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getSingleCellVol, Grid_getDeltas, Grid_getBlkData
  use Runtimeparameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

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

  integer, parameter ::  nGlobalSum = 22         ! Number of globally-summed quantities
  real :: gsum(nGlobalSum)                        ! Global summed quantities
  real :: lsum(nGlobalSum)                        ! Local summed quantities

  integer :: i, j, k
  real :: dvol, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, DIMENSION(:,:,:,:), POINTER :: scr

  integer :: point(MDIM)

  real, save :: sim_pi, sim_grv_const, sim_gamma
  logical :: gcell = .true.

!! 
!! SS : variables whose time series we output 
!!
  real :: lmin_dens, lmax_dens, gmin_dens, gmax_dens 
  real :: rms_velh_dw, rms_velv_dw, rms_vel_dw, rms_mach, rms_force
  real :: mean_pres, mean_temp
  
  call RuntimeParameters_get('gamma', sim_gamma)
  call PhysicalConstants_get( 'Newton', sim_grv_const)               !! Gravitational constant G
  call PhysicalConstants_get( 'Pi', sim_pi)                          !! Value of pi

  ! Make sure the vorticity is up-to-date
  call Grid_computeUserVars()

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  gmin_dens = 1.0
  gmax_dens = 1.0
  lmin_dens = 1.0
  lmax_dens = 1.0

  call Grid_getListOfBlocks(LEAF, blockList, count)
  

!! SS : start looping over the blocks 

  do lb = 1, count

     ! get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)


     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
             
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 

              lmin_dens = min(lmin_dens, solnData(DENS_VAR,i,j,k))
              lmax_dens = max(lmin_dens, solnData(DENS_VAR,i,j,k))

#endif           
              ! Total volume 
              lsum(2) = lsum(2) + dvol 

              ! Pressure 
              lsum(3) = lsum(3) + solnData(PRES_VAR,i,j,k)*dvol 


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! x-momentum
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              ! y-momentum
              lsum(5) = lsum(5) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              ! z-momentum
              lsum(6) = lsum(6) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol

              !! horizontal velocity dispersion squared - density weighted
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k)*(solnData(VELX_VAR,i,j,k)**2 + solnData(VELY_VAR,i,j,k)**2)*dvol

              !! vertical velocity dispersion squared - density weighted 

               lsum(8) = lsum(8) + solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)**2*dvol 
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(9) = lsum(9) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol 
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy - density weighted 
              lsum(10) = lsum(10) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol     

              ! density weighted velocity squared 

              lsum(11) = lsum(11) + solnData(DENS_VAR,i,j,k)*(solnData(VELX_VAR,i,j,k)**2 + & 
                                                              solnData(VELY_VAR,i,j,k)**2+solnData(VELZ_VAR,i,j,k)**2)*dvol      

     
              ! Sound speed squared 

              lsum(12) = lsum(12) + sim_gamma*solnData(PRES_VAR,i,j,k)*dvol
#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(13) = lsum(13) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

               ! Quantities related to the artificial forcing

               ! Forcing squared 

               lsum(14) = lsum(14) + (solnData(ACCX_VAR,i,j,k)**2 + solnData(ACCY_VAR,i,j,k)**2 + solnData(ACCZ_VAR,i,j,k)**2)*dvol

               !  rho*accx
               lsum(15) = lsum(15) + solnData(DENS_VAR,i,j,k)* &
                      &                         solnData(ACCX_VAR,i,j,k)*dvol

              ! rho*accy
               lsum(16) = lsum(16) + solnData(DENS_VAR,i,j,k)* &
                        &                       solnData(ACCY_VAR,i,j,k)*dvol

              ! rho*accz
               lsum(17) = lsum(17) + solnData(DENS_VAR,i,j,k)* &
                         &                      solnData(ACCZ_VAR,i,j,k)*dvol

              ! rho * accx *velx
               lsum(18) = lsum(18) + 2.0*solnData(DENS_VAR,i,j,k)* &
                          solnData(VELX_VAR,i,j,k)*solnData(ACCX_VAR,i,j,k)*dvol

              ! rho * accy * vely
               lsum(19) = lsum(19) + 2.0*solnData(DENS_VAR,i,j,k)* &
                          solnData(VELY_VAR,i,j,k)*solnData(ACCY_VAR,i,j,k)*dvol 

              ! rho * accz * velz 
               lsum(20) = lsum(20) + 2.0*solnData(DENS_VAR,i,j,k)* &
                          solnData(VELZ_VAR,i,j,k)*solnData(ACCZ_VAR,i,j,k)*dvol 
 
              ! Temperature 
               lsum(21) = lsum(21) + solnData(TEMP_VAR,i,j,k)*dvol           


           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)


  enddo     !! SS : end of the loop over all the blocks 

  !! SS : Calculate the sum of vorticity over all cells in the different blocks

  do lb = 1, count
     call Grid_getBlkPtr(blockList(lb), scr, SCRATCH_CTR)


     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              
              lsum(22) = lsum(22) + (scr(MVRT_SCRATCH_CENTER_VAR,i,j,k) * dvol)/2
           end do
        end do
     end do

     call Grid_releaseBlkPtr(blockList(lb), scr,SCRATCH_CTR)

  end do  !! end of loop over all the blocks
  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)
  
!! SS : Reduce calls for the min and max density
  call MPI_Reduce (lmin_dens, gmin_dens, 1, MPI_Double_Precision, MPI_Min, &
       &                MASTER_PE, MPI_Comm_World, error)

  call MPI_Reduce (lmax_dens, gmax_dens, 1, MPI_Double_Precision, MPI_Max, &
       &                MASTER_PE, MPI_Comm_World, error)


  if (io_globalMe  == MASTER_PE) then

                 mean_pres     = gsum(3)/gsum(2)
                 mean_temp     = gsum(21)/gsum(2)
                 rms_velh_dw   = sqrt(gsum(7)/gsum(1))     
                 rms_velv_dw   = sqrt(gsum(8)/gsum(1))     
                 rms_vel_dw    = sqrt(gsum(11)/gsum(1))
                 rms_mach      = sqrt(gsum(11)/gsum(12))
                 rms_force     = sqrt(gsum(14)/gsum(2))

     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time(1)                    ', &
                'mass(2)                      ', &
                'x-momentum(3)                ', &
                'y-momentum(4)                ', & 
                'z-momentum(5)                ', &
                'E_total(6)                   ', &
                'E_kinetic(7)                 ', &
                'E_internal(8)                ', &
                'rms_vel_horz(9)              ', &
                'rms_vel_vert(10)              ', &
                'rms_vel_dw(11)                ', &
                'rms_mach(12)                  ', &
                'rms_force(13)                 ', &
                'mean_pres(14)                 ', &
                'mean_temp(15)                 ', &
                'Min. dens(16)                 ', &
                'Max. dens(17)                 ', &
                'Mag_vorticity(18)             '
           
10         format (2x,50(a25, :, 1X))
           
        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12)  & ! Write the global sums to the file.
              simtime,    &                  ! Time 
              gsum(1),    &                  ! Mass 
              gsum(4),    &                  ! x-momentum 
              gsum(5),    &                  ! y-momentum 
              gsum(6),    &                  ! z-momentum 
              gsum(9),    &                  ! E_total
              gsum(10),   &                  ! E_kin 
              gsum(13),   &                  ! E_int
              rms_velh_dw, &                 ! rms value of horz. vel. 
              rms_velv_dw, &                 ! rms value of vert. vel. 
              rms_vel_dw,  &                 ! rms velocity - 3d 
              rms_mach,    &                 ! rms Mach number 
              rms_force,   &                 ! rms value of the forcing 
              mean_pres,   &                 ! Mean pressure 
              mean_temp,   &                 ! Mean temperature
              gmin_dens,   &                 ! Min. density
              gmax_dens,   &                 ! Max. density
              gsum(21)                       ! Mag. of vorticity

12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



