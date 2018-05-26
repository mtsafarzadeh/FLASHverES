!!****if* source/Simulation/SimulationMain/Life/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklst - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
  use Simulation_data, ONLY : sim_targetMass, sim_meshComm, sim_inputEnergy, &
                              sim_driverType, sim_pulseLength
  use Grid_interface,  ONLY : Grid_getSingleCellVol, &
                              Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_fillGuardCells, &
                              Grid_releaseBlkPtr, Grid_getSingleCellVol
  use Eos_interface,   ONLY : Eos_wrapped, Eos_getAbarZbar
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: i,j,k, lb
  real, pointer   :: blkPtr(:,:,:,:)
  real :: vol
  real :: targ_mass, mass_frac
  real :: rho
  real :: eele
  real :: elas
  real :: mfrac
  real :: abar, zbar
  real :: gradpe(MDIM), gradne(MDIM), cross(MDIM), bier(MDIM)
  real :: delta(MDIM) ! The cell width

  real :: block_mass, target_mass

  integer :: ierr
  integer :: is, ie, js, je, ks, ke
  real, allocatable :: nele(:,:,:)

  ! *************************************
  ! *                                   *
  ! *     COMPUTE TOTAL TARGET MASS     *
  ! *                                   *
  ! *************************************
  sim_targetMass = 0.0
  target_mass = 0.0

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     block_mass = 0.0
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              call Grid_getSingleCellVol(blklst(lb), EXTERIOR, (/i,j,k/), vol)
              mfrac = blkPtr(TARG_SPEC,i,j,k)
              rho = blkPtr(DENS_VAR,i,j,k)              
              block_mass = block_mass + vol * rho * mfrac
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)

     target_mass = target_mass + block_mass
  end do  

  ! Sum the target mass across all processes:  
  call mpi_allreduce(target_mass, sim_targetMass, 1, FLASH_REAL, & 
       MPI_SUM, sim_meshComm, ierr )     

  ! Deposit energy:

  if(stime < sim_pulseLength .and. &
       (sim_driverType == "uniform" .or. sim_driverType == "unispec")) then

     ! Compute amount of energy to dump in on this time
     ! step:
     elas = sim_inputEnergy * dt/sim_pulseLength

     if (sim_driverType == "uniform") then
        elas = elas/sim_targetMass
     end if

     do lb = 1, blkcnt
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)
        
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 
                 
                 eele  = blkPtr(EELE_VAR, i,j,k) 
                 mfrac = blkPtr(TARG_SPEC,i,j,k)
                 
                 eele = eele + mfrac*elas
                 blkPtr(EELE_VAR,i,j,k) = eele
              enddo
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
        
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blklst(lb))
     end do
  end if

end subroutine Simulation_adjustEvolution
