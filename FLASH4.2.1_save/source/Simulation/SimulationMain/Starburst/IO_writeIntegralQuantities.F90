!!****if* source/Simulation/SimulationMain/StirTurbScalar/IO_writeIntegralQuantities
!!    SS : change the above directory StirTurb to StirTurbScalar
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
!!   SS : Modified to add diagonostics for the magnetic energy, divB, plasma 
!!        beta and other stuff. 
!!***

!!REORDER(4): solnData, scr

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalMe
  use Grid_interface, ONLY : Grid_computeUserVars, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getSingleCellVol, Grid_getGlobalIndexLimits
  use Logfile_interface, ONLY : Logfile_stamp
  use IO_interface, ONLY : IO_writeCheckpoint, IO_setScalar
  use Driver_interface, ONLY : Driver_getNStep, Driver_abortFlash, Driver_getSimTime
! SS : added a call to runtime parameters for Gamma, magnetic and metals keyword
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
   use Multispecies_interface, ONLY : Multispecies_getAvg
   use Eos_interface, ONLY : Eos, Eos_wrapped
   use Simulation_data, ONLY : sim_rhoAmbient, sim_xmax
   use Driver_data, ONLY: dr_tmax
  implicit none

!!#include "mpif.h"
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
#include "Flash_mpi.h"
  
  real, intent(in) :: simTime
  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: killunit = 69
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  character (len=MAX_STRING_LENGTH), save :: killname
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer :: size(MDIM)

  integer, parameter ::  nGlobalSum = 15   ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities
  real :: lmin_dens, lmax_dens, gmin_dens, gmax_dens                           ! min and max density variables
  real :: lmin_temp, lmax_temp, gmin_temp, gmax_temp
  real :: lmin_gam, lmax_gam, gmin_gam, gmax_gam
  real, save :: sim_pi, sim_gamma    ! constants
  !!WJG: MY STUFF
  real :: avg_gam, curr_ke, curr_sim_time, spec_sum
  !real , dimension(NSPECIES) ::  weights
  !integer, dimension(NSPECIES) :: speciesMask
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: start_avg
  real, dimension(NSPECIES) :: end_avg
  real, dimension(NSPECIES) :: dfoverf
  real			    :: spec_cutoff, nhtot
  integer		    :: raiseJ21, didraise, curr_step, kick_step
  integer, dimension(MDIM)  :: GlobalIndexLimits
  logical                   :: fileexist

  integer :: i, j, k
  real :: dvol, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, DIMENSION(:,:,:,:), POINTER :: scr

  integer :: point(MDIM)
  
  !! SS:  value of PI called from constants.h 
  sim_pi = PI

  ! Make sure the vorticity is up-to-date
  call Grid_computeUserVars()
  ! SS : Call runtime parameters if any needed to be called 
  call RuntimeParameters_get('gamma', sim_gamma)
  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.

! SS : added the following for monitoring min and max dens
  gmin_dens = 1.0
  gmax_dens = 1.0e-50
  lmin_dens = 1.0
  lmax_dens = 1.0e-50
  gmin_temp = 1.0e50
  gmax_temp = 1.0e-50
  lmin_temp = 1.0e50
  lmax_temp = 1.0e-50
  gmin_gam = 1.0e50
  gmax_gam = 1.0e-50
  lmin_gam = 1.0e50
  lmax_gam = 1.0e-50
  dfoverf = 0.0e0
  nhtot = 0.0e0
  curr_sim_time = 0.0e0

  call Grid_getGlobalIndexLimits(GlobalIndexLimits)
  call Grid_getListOfBlocks(LEAF, blockList, count)
  call Driver_getNStep(curr_step)
  call Driver_getSimTime(curr_sim_time)
  !!WJG: Set speciesMask here since I only need to do this once
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated range of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
	      !!WJG: Set weights and species to get the average gamma for sound speed stuff  
	      eosData(:) = 1.0e-20
	      eosData(EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
	      eosData(EOS_DENS) = solnData(DENS_VAR,i,j,k)
	      call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
	      avg_gam = eosData(EOS_GAMC)
	      !!WJG: Get average gamma
	      !!call MultiSpecies_getAvg(GAMMA,avg_gam,weights,speciesMask)
	      !!avg_gam = avg_gam * NSPECIES
              lmin_gam = min(lmin_gam, avg_gam)
              lmax_gam = max(lmax_gam, avg_gam)
	      !!write(*,'(A,1pe10.3)') 'gam: (should be 1.66667) ', avg_gam
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
                lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
!! SS : compute the min and max densities 
                lmin_dens = min(lmin_dens, solnData(DENS_VAR,i,j,k))
                lmax_dens = max(lmax_dens, solnData(DENS_VAR,i,j,k))
                lmin_temp = min(lmin_temp, solnData(TEMP_VAR,i,j,k))
                lmax_temp = max(lmax_temp, solnData(TEMP_VAR,i,j,k))
#endif           
! Total Volume 
                lsum(2) = lsum(2) + dvol
! Pressure
                lsum(3) = lsum(3) + solnData(PRES_VAR,i,j,k)*dvol
#ifdef DENS_VAR   
! beginning of DENS_VAR loop

#ifdef VELX_VAR      
              ! momentum
                lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * solnData(VELX_VAR,i,j,k)*dvol 
#endif ! end of VELX_VAR loop

#ifdef VELY_VAR      
                lsum(5) = lsum(5) + solnData(DENS_VAR,i,j,k) * solnData(VELY_VAR,i,j,k)*dvol
#endif ! end of VELY_VAR loop

#ifdef VELZ_VAR      
                lsum(6) = lsum(6) + solnData(DENS_VAR,i,j,k) * solnData(VELZ_VAR,i,j,k)*dvol
#endif  ! end of VELZ_VAR loop
              ! total energy
#ifdef ENER_VAR
                lsum(7) = lsum(7) + solnData(ENER_VAR,i,j,k) * solnData(DENS_VAR,i,j,k)*dvol 
#endif ! end of ENER_VAR loop

#ifdef EINT_VAR
              ! internal energy
                lsum(8) = lsum(8) + solnData(DENS_VAR,i,j,k) * solnData(EINT_VAR,i,j,k)*dvol
#endif  ! end of EINT_VAR loop

#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! density weighted kinetic energy
                lsum(9) = lsum(9) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                    &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                    &                              solnData(VELY_VAR,i,j,k)**2+ & 
                    &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

              ! density weighted velocity squared 
                 lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k)*(solnData(VELX_VAR,i,j,k)**2+solnData(VELY_VAR,i,j,k)**2+solnData(VELZ_VAR,i,j,k)**2)*dvol
#endif
#endif
#endif ! end of VELX_VAR, VELY_VAR and VELZ_VAR loop

! Sound speed squared  
		 lsum(11) = lsum(11) + avg_gam*solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)*dvol

#ifdef VELX_VAR
#ifdef VELY_VAR
#ifdef VELZ_VAR
                 !pdv - only includes hydrodynamic pressure; dvol=dx*dy*dz and all the 1/dx in the 
                 ! the differentiation gives a total of dv**(2/3) which is like a surface term
                  lsum(12) = lsum(12) + solnData(PRES_VAR, i, j, k) * &
                       &                0.5* (solnData(VELX_VAR,i+1,j,k) - solnData(VELX_VAR,i-1,j,k) + &
                       &                      solnData(VELY_VAR,i,j+1,k) - solnData(VELY_VAR,i,j-1,k) + &
                       &                      solnData(VELZ_VAR,i,j,k+1) - solnData(VELZ_VAR,i,j,k-1)) * dvol**(2./3.)

#endif
#endif
#endif

!! Volume Weighted EINT
		lsum(13) = lsum(13) + solnData(EINT_VAR,i,j,k) * dvol

!! Volume Weighted TEMP
		lsum(14) = lsum(14) + solnData(TEMP_VAR,i,j,k) * dvol

!! Density Weighted TEMP
		lsum(15) = lsum(15) + solnData(DENS_VAR,i,j,k) * solnData(TEMP_VAR,i,j,k) * dvol

!! WJG: Adding Species Stuff. I am adding density weighted values for each species
!!      Evan also wants temperature & KE weighted by the density of each species
#endif ! ifdef DENS_VAR
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)
  enddo

  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)
!! SS : added the following reduce calls for the min and max density
  call MPI_Reduce (lmin_dens, gmin_dens, 1, MPI_Double_Precision, MPI_Min, &
       &                MASTER_PE, MPI_Comm_World, error)
  call MPI_Reduce (lmax_dens, gmax_dens, 1, MPI_Double_Precision, MPI_Max, &
       &                MASTER_PE, MPI_Comm_World, error)
  call MPI_Reduce (lmin_temp, gmin_temp, 1, MPI_Double_Precision, MPI_Min, &
       &                MASTER_PE, MPI_Comm_World, error)
  call MPI_Reduce (lmax_temp, gmax_temp, 1, MPI_Double_Precision, MPI_Max, &
       &                MASTER_PE, MPI_Comm_World, error)
  call MPI_Barrier(MPI_COMM_WORLD,error)

  if (io_globalMe  == MASTER_PE) then
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
				'time(1)                      ', &
				'mass(2)                      ', &
				'x-momentum(3)                ', &
				'y-momentum(4)                ', &
				'z-momentum(5)                ', &
				'DW-E_total(6)                ', &
				'DW-E_kinetic(7)              ', &
				'DW-E_internal(8)             ', &
				'DW-Temp(9)                   ', &
				'VW-E_internal(10)            ', &
				'VW-Temp(11)                  ', &
                		'Min Density(12)              ', &
                     		'Max Density(13)              ', &
                     		'Min Temperature(14)          '

10         format (2x,235(a25, :, 1X))

        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11)
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) &  ! Write the global sums to the file
!! SS : all quantities to be written 
			simtime, &			!!time
			gsum(1) , &			!!mass
			gsum(4) , &			!!x-mom
			gsum(5) , &			!!y-mom
			gsum(6) , &			!!z-mom
			gsum(7)/gsum(1) , &		!!DW-E_tot
			gsum(9)/gsum(1) , &		!!DW-E_Kin
			gsum(8)/gsum(1) , &		!!DW-E_int
			gsum(15)/gsum(1), &		!!DW-Temp
			gsum(14)/gsum(2), &		!!VW-Eint
			gsum(13)/gsum(2), &		!!VW-Temp
		        gmin_dens        , &            !!Min Dens
            		gmax_dens         , &            !!Max Dens
            		gmin_temp         , &            !!Min Temp
            		gmax_temp                        !!Max Temp

12   format (1x, 235(es25.18, :, 1x))
     close (funit)          ! Close the file.
     
  endif !!ends master PE check
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



