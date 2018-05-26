!****if* source/Simulation/SimulationMain/StirTurbScalar/Simulation_initBlock
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
  use Cool_data, ONLY: dustinit
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos, Eos_wrapped
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
  real :: enerZone, ekinZone, eintZone, rhoZone, tempZone, presZone
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  !ES Variables for White Noise INitial Conditions 
  real, allocatable  :: rvec(:)         ! for the random number generator
  real               :: random
  integer            :: rvecSize=0      ! number of random numbers needed
  integer, parameter :: iseed = -8679
  integer            :: iseedUse
  integer            :: icount,icount2,isize2
  integer            :: itwo,jtwo,ktwo
  real               :: norm,normrho,normtemp

  !ES CHemistry stuff
  real               :: smallx
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData

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
  call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeX)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeX)



!Setup initial composition of all species
!Question, Cellular had some small fraction for each species, can this not be
!zero?
smallx = 1.0e-30;  !Trying something small
massFraction(:) = smallx

if (H_SPEC > 0) massFraction(H_SPEC) = max(sim_xH, smallx)
if (HP_SPEC > 0) massFraction(HP_SPEC) = max(sim_xHP, smallx)
if (HM_SPEC > 0) massFraction(HM_SPEC) = max(sim_xHM, smallx)
if (H2_SPEC > 0) massFraction(H2_SPEC) = max(sim_xH2, smallx)
if (H2P_SPEC > 0) massFraction(H2P_SPEC) = max(sim_xH2P, smallx)
if (H3P_SPEC > 0) massFraction(H3P_SPEC) = max(sim_xH3P, smallx)
if (HE_SPEC > 0) massFraction(HE_SPEC) = max(sim_xHE, smallx)
if (HEP_SPEC > 0) massFraction(HEP_SPEC) = max(sim_xHEP, smallx)
if (C_SPEC > 0) massFraction(C_SPEC) = max(sim_xC, smallx)
if (CP_SPEC > 0) massFraction(CP_SPEC) = max(sim_xCP, smallx)
if (CM_SPEC > 0) massFraction(CM_SPEC) = max(sim_xCM, smallx)
if (O_SPEC > 0) massFraction(O_SPEC) = max(sim_xO, smallx)
if (OP_SPEC > 0) massFraction(OP_SPEC) = max(sim_xOP, smallx)
if (OM_SPEC > 0) massFraction(OM_SPEC) = max(sim_xOM, smallx)
if (C2_SPEC > 0) massFraction(C2_SPEC) = max(sim_xC2, smallx)
if (O2_SPEC > 0) massFraction(O2_SPEC) = max(sim_xO2, smallx)
if (O2P_SPEC > 0) massFraction(O2P_SPEC) = max(sim_xO2P, smallx)
if (OH_SPEC > 0) massFraction(OH_SPEC) = max(sim_xOH, smallx)
if (OHP_SPEC > 0) massFraction(OHP_SPEC) = max(sim_xOHP, smallx)
if (CO_SPEC > 0) massFraction(CO_SPEC) = max(sim_xCO, smallx)
if (COP_SPEC > 0) massFraction(COP_SPEC) = max(sim_xCOP, smallx)
if (CH_SPEC > 0) massFraction(CH_SPEC) = max(sim_xCH, smallx)
if (CHP_SPEC > 0) massFraction(CHP_SPEC) = max(sim_xCHP, smallx)
if (CH2_SPEC > 0) massFraction(CH2_SPEC) = max(sim_xCH2, smallx)
if (CH2P_SPEC > 0) massFraction(CH2P_SPEC) = max(sim_xCH2P, smallx)
if (HCOP_SPEC > 0) massFraction(HCOP_SPEC) = max(sim_xHCOP, smallx)
if (HOCP_SPEC > 0) massFraction(HOCP_SPEC) = max(sim_xHOCP, smallx)
if (H2O_SPEC > 0) massFraction(H2O_SPEC) = max(sim_xH2O, smallx)
if (H2OP_SPEC > 0) massFraction(H2OP_SPEC) = max(sim_xH2OP, smallx)
if (H3OP_SPEC > 0) massFraction(H3OP_SPEC) = max(sim_xH3OP, smallx)
if (ELEC_SPEC > 0) massFraction(ELEC_SPEC) = max(sim_xELEC, smallx)


! Now setup the other variables

  random   =  (xCoord(6)+sqrt(2.)*yCoord(6)+sqrt(5.)*zCoord(6))/1E18
  iseedUse = iseed*blockID+random   !use a different random seed on each block
  print *,'randomseed',iseedUse,blockID
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)


  !..get a blocks worth of random numbers between 0.0 and 1.0
  rvecSize = blkLimitsGC(HIGH,IAXIS)*blkLimitsGC(HIGH,JAXIS)*blkLimitsGC(HIGH,KAXIS)*4
  allocate(rvec(rvecSize))
  call sim_ranmar(iseedUse, rvec, rvecSize)
  do i = 1,rvecSize
!  the 3 here is to get the rms to come out right
!  sim_noise_amplitude is the rms if d is small
    rvec(i) = (2*rvec(i)-1.)*sim_noise_amplitude*sqrt(3.)
  enddo
  norm = (exp(sim_noise_amplitude*sqrt(3.))-exp(-sim_noise_amplitude*sqrt(3.))) & 
         /(2.*sim_noise_amplitude*sqrt(3.))
  normrho  = 1./norm
  normtemp = norm*(1.-sim_isothermal)+sim_isothermal

  ! Loop over cells in the block.
  icount = 0
  isize2 = (blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS))/4
! if sim_noise = 1 then we are doubling and halving
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     ktwo = int(k/4)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        jtwo = int(j/4)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
          itwo = int(i/4)
          icount   = icount+1
          icount2  = itwo+isize2*(jtwo+isize2*ktwo)+1
          rhoZone  = sim_rhoAmbient *normrho                      & 
                       *exp(rvec(icount2))*(1.+0.2*rvec(icount))
          tempZone = sim_tempAmbient*(normtemp)                   &
                       *(sim_isothermal+ (1-sim_isothermal)/      & 
                       (exp(rvec(icount2))*(1.+0.2*rvec(icount))))


          eosData(EOS_TEMP) = rhoZone
          eosData(EOS_DENS) = tempZone
          eosData(EOS_PRES) = 0.0e0
          do n=SPECIES_BEGIN,SPECIES_END
            solnData(n, i,j,k) =  massFraction(n)
          enddo
          call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
          sim_gamma = eosData(EOS_GAMC)

          presZone = eosData(EOS_PRES)
          eintZone = eosData(EOS_EINT)
          ekinZone = 0.
          enerZone = eintZone + ekinZone

          solnData(DENS_VAR,i,j,k) = rhoZone
          solnData(TEMP_VAR,i,j,k) = tempZone
          solnData(VELX_VAR,i,j,k) = 0.
          solnData(VELY_VAR,i,j,k) = 0.
          solnData(VELZ_VAR,i,j,k) = 0.
          solnData(PRES_VAR,i,j,k) = presZone 
          solnData(ENER_VAR,i,j,k) = enerZone
          solnData(EINT_VAR,i,j,k) = eintZone
          solnData(GAMC_VAR,i,j,k) = sim_gamma
          solnData(GAME_VAR,i,j,k) = sim_gamma
          solnData(DUST_MSCALAR,i,j,k) = dustinit

        enddo
     enddo
  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock


