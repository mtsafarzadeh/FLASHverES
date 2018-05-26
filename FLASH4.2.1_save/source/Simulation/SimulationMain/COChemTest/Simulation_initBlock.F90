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
  use Cool_data
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getBlkPtr, &
      Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_putPointData
  use Eos_interface, ONLY : Eos, Eos_wrapped
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
  real   ::   vx, vy, vz, p, rho, metals
  real   ::   e, ek, T

!!*** sizes in each of the three direction
  integer :: sizeX,sizeY,sizeZ
  integer :: istat

!! This says that we are grabing the guard cells when we grab the block
  logical :: gcell = .true.
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

!This is a part dedicated to the Chemistry. This is followed
!from the Cellular problem and the unitTest problem
real :: abar, zbar
real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
real :: temp_zone, rho_zone, vel_zone
real :: gamma, smallx
real :: ptot, eint, etot
real, dimension(EOS_NUM) :: eosData


!!DEDICATED TO TEST THE CHEMISTRY RATES / Cooling Rates
integer	  :: specieMap
real, dimension(SPECIES_BEGIN:SPECIES_END) :: xxin
real, dimension(NSPECIES)		   :: nnin, ys2
real  :: ctemp, denrate, divv, lco, nnh, crout
integer	 :: iNH, iTemp



!  print *,'Simulation InitBlock'
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData)

!!***get coordinates of this block***

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
  zz=0.
  do k = 1, sizeZ
     if(NDIM==3) zz = z(k)

     do j = 1, sizeY
        yy = y(j)

        do i = 1, sizeX
           xx = x(i)
!!         ***distance to center.***           
           radius=sqrt(xx**2 + yy**2 + zz**2)


 
!           print *, 'Sim Temp: ', sim_c_temp
	   rho= sim_c_den  !Test Cloudy
	   T= sim_c_temp  !Test Cloud
	   p=rho*T*sim_gasConst !!This will stay the same for pressure
        
!Getting ready to find Gamma and some other Eos stuff
	eosData(EOS_TEMP) = sim_c_temp
	eosData(EOS_DENS) = sim_c_den
        eosData(EOS_PRES) = 0.0e0
!	print *, 'bad_P: ', p
	call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
	gamma = eosData(EOS_GAMC)
 	sim_gamma = eosData(EOS_GAMC)
	ptot = eosData(EOS_PRES)
!	print *, 'good_p: ', ptot
	eint = eosData(EOS_EINT)	

          metals = 0.0

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
	solnData(METL_MSCALAR,i,j,k) = sim_meta

	do n=SPECIES_BEGIN,SPECIES_END
 	  solnData(n, i,j,k) = 	massFraction(n)
	!  print *, 'Massfraction(',n,'):', massFraction(n)
	enddo
	
	enddo      
     enddo
  enddo


!!TEST COOL_ROVIB
!!make sure to uncomment write line in cool_rovib
!!run as: ./flash3 > data.txt
!! cat data.txt | tail -n 41001 | head -n 41000

!!  denrate = 1.0e-30
!!  xxin(:) = smallx
!!  xxin(H_SPEC)   = 0.76
!!  xxin(H2_SPEC)  = 1.0e-5
!!  xxin(CO_SPEC)  = 7.7e-6
!!  xxin(OH_SPEC)  = 5.0e-6
!!  xxin(H2O_SPEC) = 5.5e-6

!!  nnin(:) = smallx
!!  nnin(iH) = (denrate*6.022e23)/aion(iH) * xxin(H_SPEC)
!!  nnin(iH2) = (denrate*6.022e23)/aion(iH2) * xxin(H2_SPEC)
!!  nnin(iCO) = (denrate*6.022e23)/aion(iCO) * xxin(CO_SPEC)
!!  nnin(iH2O) = (denrate*6.022e23)/aion(iH2O) * xxin(H2O_SPEC)
!!  nnin(iOH) = (denrate*6.022e23)/aion(iOH) * xxin(OH_SPEC)
!!  lco = 0.0e0
!!  divv = 0.0e0
!!  crout = 0.0e0

!! do iTemp = 1,410
!!     do iNH = 1,100
!!       nnh = 9.0 + (iNH-1)*0.12
!!       ctemp = 5.0 + (iTemp-1)*10.
!!       divv = nnin(iH2O)/10**(nnh)
!!       call Cool_rovib(ctemp,denrate,nnin,crout,divv,lco)
!!     enddo
!!  enddo
!!  print *, 'CALLING ABORT'
!!  call abort()


!!Get Chemistry Rates
!!  denrate = 1.0e-30
!!  xxin(:) = smallx
!!  xxin(H_SPEC) = 0.76
!!  xxin(H2_SPEC) = 1.0e-5
!!  xxin(CO_SPEC) = 7.7e-6
!!  xxin(OH_SPEC) = 5.0e-6
!!  xxin(H2O_SPEC) = 5.5e-6
!!  xxin(HE_SPEC) = 0.24

!!  ys2(:) = smallx
!!  ys2(iH) = xxin(H_SPEC)/aion(iH)
!!  ys2(iH2) = xxin(H2_SPEC)/aion(iH2)
!!  ys2(iCO) = xxin(CO_SPEC)/aion(iCO)
!!  ys2(iOH) = xxin(OH_SPEC)/aion(iOH)
!!  ys2(iH2O) = xxin(H2O_SPEC)/aion(iH2O)
!!  ys2(iHE) = xxin(HE_SPEC)/aion(iHE)
  
!!  do iTemp = 1,501
!!      ctemp = 10**(1.0 + 0.012*(iTemp-1))
      !!print *, ' itemp: ', iTemp, ' ctemp: ', ctemp
!!      call Chemistry_networkRates(ctemp,denrate,ys2)
!!  enddo
!!  print *, ' CALLING ABORT '
!!  call abort()


  !print *, ' MADE IT HERE'
  call Eos_wrapped(MODE_DENS_TEMP,blkLimits,blockID)
  call Grid_releaseBlkPtr(blockID,solnData)

  deallocate(x)
  deallocate(y)
  deallocate(z)

 return

end subroutine Simulation_initBlock
