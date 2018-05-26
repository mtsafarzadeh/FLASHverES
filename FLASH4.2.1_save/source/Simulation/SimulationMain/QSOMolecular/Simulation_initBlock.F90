!!****if* source/Simulation/SimulationMain/Marcus/Simulation_initBlock
!!
!! NAME
!!
!!  init_block
!!
!! 
!! SYNOPSIS
!!
!!  Simulation/init_Block(integer(IN) :: blockID,
!!                        integer(IN) :: myPE)	                   
!!
!!
!!
!! DESCRIPTION
!!
!!
!!   Takes the one-d profile defined in Simulation_init and applys it
!!  This version sets up a spherical cluster in hydrostatic equilibrium
!!  in a 2D/3D cartesian grid. 
!!  The density and temperature profiles are give analytically.
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!  myPE -             my processor number
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_initblock (blockID,myPE)



!!***used modules from FLASH.***

  use Simulation_data
  use Chemistry_data
  use Chemistry_coolData
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getBlkPtr, &
      Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_putPointData
  use Eos_interface, ONLY : Eos, Eos_wrapped
  use Cool_data, ONLY: dustinit
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

! HERE are the arguments
  integer, intent(in) :: blockID
  integer, intent(in) :: myPE

 

!!***block number and (cell) counters***
  integer  ::  i, j, k, n

!!***some dummy variables***
  real             :: dummy

!!** this is where we read and write the data
real, pointer, dimension(:,:,:,:)  :: solnData

!!***vectors that store cell dimensions **
  real, allocatable, dimension(:) :: x, y, z

!!***coordinates of grid cells***
  real :: xx,  yy, zz, radius

!!***variables that you set as you go 
  real   ::   vx, vy, vz, p, rho
  real   ::   e, ek, T

!!*** sizes in each of the three direction
  integer :: sizeX,sizeY,sizeZ
  integer :: istat

!! This says that we are grabing the guard cells when we grab the block
  logical :: gcell = .true.
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

! This is a part dedicated to the Chemistry. This is followed
! from the Cellular problem and the unitTest problem
  real :: abar, zbar
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real :: tempZone, rhoZone, velZone
  real :: gamma, smallx
  real :: ptot, eint, etot
  real, dimension(EOS_NUM) :: eosData


! DEDICATED TO TEST THE CHEMISTRY RATES / Cooling Rates
  integer	  :: specieMap
  real, dimension(NSPECIES)		   :: nnin
  real, dimension(7)			   :: clum
  real  :: ctemp, denrate, divv, lco, nnh, crout
  integer	 :: iNH, iTemp


! ES Variables for White Noise Initial Conditions 
  real, allocatable  :: rvec(:)         ! for the random number generator
  real               :: random
  integer            :: rvecSize=0      ! number of random numbers needed
  integer, parameter :: iseed = -8679
  integer            :: iseedUse
  integer            :: icount,icount2,isize2
  integer            :: itwo,jtwo,ktwo
  real               :: norm,normrho,normtemp





!  print *,'Simulation InitBlock'
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData)

!!***get coordinates of this block***

  rvecSize = blkLimitsGC(HIGH,IAXIS)*blkLimitsGC(HIGH,JAXIS)*blkLimitsGC(HIGH,KAXIS)*4
  allocate(rvec(rvecSize))
  call sim_ranmar(iseedUse, rvec, rvecSize)
  do i = 1,rvecSize
    rvec(i) = (2*rvec(i)-1.)*sim_noise_amplitude*sqrt(3.)
  enddo
  norm = (exp(sim_noise_amplitude*sqrt(3.))-exp(-sim_noise_amplitude*sqrt(3.))) &
         /(2.*sim_noise_amplitude*sqrt(3.))
  normrho  = 1./norm
  normtemp = norm*(1.-sim_isothermal)+sim_isothermal

  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(x(sizeX),stat=istat)
  allocate(y(sizeY),stat=istat)
  allocate(z(sizeZ),stat=istat)
  x = 0.0
  y = 0.0
  z = 0.0
 
  if (NDIM==3) call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,z,sizeZ)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,y,sizey)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,x,sizex)

!Setup initial composition of all species
!Question, Cellular had some small fraction for each species, can this not be zero?
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
  

!***************************************************
!!***Now loop over all cells in the current block***
!  
  print *,'Starting the InitBlock Loop'
  icount = 0
  isize2 = (blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS))/4
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
    ktwo = (k/4)
    do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
      jtwo = (j/4)
      do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
        itwo = (i/4)
        icount   = icount+1
        icount2  = itwo+isize2*(jtwo+isize2*ktwo)+1
        rhoZone  = sim_rhoAmbient *normrho                      &
                       *exp(rvec(icount2))*(1.+0.2*rvec(icount))
        tempZone = sim_tempAmbient*(normtemp)                   &
                       *(sim_isothermal+ (1-sim_isothermal)/      &
                       (exp(rvec(icount2))*(1.+0.2*rvec(icount))))      

        print *,'rho-T',rhoZone,tempZone
        do n=SPECIES_BEGIN,SPECIES_END
          solnData(n, i,j,k) =  massFraction(n)
        !  print *, 'Massfraction(',n,'):', massFraction(n)
        enddo

        eosData(EOS_DENS) = rhoZone
        eosData(EOS_TEMP) = tempZone
        eosData(EOS_PRES) = 0.0e0

        call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
        gamma = eosData(EOS_GAMC)
        sim_gamma = eosData(EOS_GAMC)
        ptot = eosData(EOS_PRES)

        eint = eosData(EOS_EINT)

        vx   = 0.
        vy   = 0.
        vz   = 0.
        ek   = 0.

        solnData(DENS_VAR,i,j,k) = eosData(EOS_DENS)
        solnData(PRES_VAR,i,j,k) = eosData(EOS_PRES)
        solnData(VELX_VAR,i,j,k) = vx
        solnData(VELY_VAR,i,j,k) = vy
        solnData(VELZ_VAR,i,j,k) = vz
        solnData(GAME_VAR,i,j,k) = eosData(EOS_GAMC)
        solnData(GAMC_VAR,i,j,k) = eosData(EOS_GAMC)
        solnData(EINT_VAR,i,j,k) = eosData(EOS_EINT)
        solnData(ENER_VAR,i,j,k) = eosData(EOS_EINT) + ek
        solnData(TEMP_VAR,i,j,k) = eosData(EOS_TEMP)
        solnData(DUST_MSCALAR,i,j,k) = dustinit
        solnData(METL_MSCALAR,i,j,k) = 0.

      enddo      
    enddo
  enddo


  !print *, ' MADE IT HERE'
  call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
  call Grid_releaseBlkPtr(blockID,solnData)

  deallocate(x)
  deallocate(y)
  deallocate(z)

 return

end subroutine Simulation_initBlock
