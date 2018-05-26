!!****if* source/Simulation/SimulationMain/StirTurbScalar/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block. 
!!
!!  SS : This routine also includes magnetic fields to study dynamo 
!!       action in turbulent fluids with the new Unsplit Staggered 
!!       Mesh solver available with FLASH4. The routine is heavily 
!!       modified from the Simulation_initBlock.F90 routine provided 
!!       in the MHD Orszag Tang vortex problem in FLASH 4. 
!!  
!!  SS : Add scalars to this routine. We will add 3 forced scalars 
!!       CONC, CONB and CONE which will be forced at the driving scale 
!!       of turbulence. Optional are three decaying scalars, DECC, DECD, DECE
!!     
!! ARGUMENTS
!!
!!  blockID - the number of the block to update
!!
!! 
!!
!!***

!!REORDER(4): solnData, face[xy]Data

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Eos_data, ONLY: eos_singleSpeciesA
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos, Eos_wrapped
  use Cosmology_data, ONLY: csm_rho0, csm_sound0,csm_scalefactor
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone, rhoZone, tempZone, presZone, tempAmbient
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  !ES Variables for White Noise INitial Conditions 
  real, allocatable  :: rvec(:)         ! for the random number generator
  real               :: random
  integer            :: rvecSize=0      ! number of random numbers needed
  integer, parameter :: iseed = -8679
  integer            :: iseedUse
  integer            :: icount
  integer            :: i4,j4,k4,i2,j2,k2,icount4,icount2

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
  random   = (xCoord(6)+sqrt(2.)*yCoord(6)+sqrt(5.)*zCoord(6))/1E18*csm_scaleFactor
  iseedUse = iseed*blockID+random   !use a different random seed on each block
  print *,'randomseed',iseedUse,blockID
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)


  !..get a blocks worth of random numbers between 0.0 and 1.0
  rvecSize = blkLimitsGC(HIGH,IAXIS)*blkLimitsGC(HIGH,JAXIS)*blkLimitsGC(HIGH,KAXIS)*2.
  allocate(rvec(rvecSize))
  call sim_ranmar(iseedUse, rvec, rvecSize)

  ! Loop over cells in the block.
  icount = 0
  tempAmbient = csm_sound0**2*eos_singleSpeciesA*sim_amu/sim_boltz/sim_gamma
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     k4 = (k-blkLimitsGC(LOW,KAXIS))/4       
     k4 = k4*(blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS))*(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS))/4/4
     k2 = (k-blkLimitsGC(LOW,KAXIS))/2       
     k2 = k2*(blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS))*(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS))/2/2
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        j4 = (j-blkLimitsGC(LOW,JAXIS))/4       
        j4 = j4*(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS))/4
        j2 = (j-blkLimitsGC(LOW,JAXIS))/2       
        j2 = j2*(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS))/2
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           i4 = (i-blkLimitsGC(LOW,IAXIS))/4       
           i2 = (i-blkLimitsGC(LOW,IAXIS))/2       
           icount = icount+1
           icount2 = k2+j2+i2+1 
           icount4 = k4+j4+i4+1 
           rhoZone = csm_rho0*& 
                     (1.0+sim_noise_amplitude*(1.0-0.5/sqrt(8.)-0.25/8.)*(1.0-2.0*rvec(icount4))&
                         +sim_noise_amplitude*(0.5-0.25/sqrt(8.))*(1.0-2.0*rvec(icount2))&
                         +sim_noise_amplitude*0.25*(1.0-2.0*rvec(icount)))
           tempZone = tempAmbient/&
                     (1.0+(sim_noise_amplitude*(1.0-0.5/sqrt(8.)-0.1/8.)*(1.0-2.0*rvec(icount4))&
                          +sim_noise_amplitude*(0.5-0.25/sqrt(8.))*(1.0-2.0*rvec(icount2))&
                          +sim_noise_amplitude*0.25*(1.0-2.0*rvec(icount))) *(1.-sim_isothermal))
           presZone= rhoZone*tempZone*sim_boltz/sim_amu/eos_singleSpeciesA
           eintZone = presZone*sim_gamma/(sim_gamma-1.)/rhoZone
           enerZone = eintZone + ekinZone
           ekinZone = 0.

           solnData(DENS_VAR,i,j,k) = rhoZone*sim_scalefactor**3
           solnData(TEMP_VAR,i,j,k) = tempZone/sim_scalefactor**2
           solnData(VELX_VAR,i,j,k) = 0.
           solnData(VELY_VAR,i,j,k) = 0.
           solnData(VELZ_VAR,i,j,k) = 0.
           solnData(PRES_VAR,i,j,k) = presZone*sim_scalefactor 
           solnData(ENER_VAR,i,j,k) = enerZone/sim_scalefactor**2
           solnData(EINT_VAR,i,j,k) = eintZone/sim_scalefactor**2
           solnData(GAMC_VAR,i,j,k) = sim_gamma
           solnData(GAME_VAR,i,j,k) = sim_gamma

           solnData(DUST_MSCALAR,i,j,k) = sim_dustinit
 
        enddo
     enddo
  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock






