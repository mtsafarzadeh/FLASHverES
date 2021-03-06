!!****if* source/physics/Eos/EosMain/multiTemp/Eos_wrapped
!! NAME
!!
!!  Eos_wrapped
!! 
!! SYNOPSIS
!!
!!  call Eos_wrapped(  integer(IN) :: mode,
!!                     integer(IN) :: range(HIGH, MDIM),
!!                     integer(IN) :: blockID,
!!            optional,integer(IN) :: gridDataStruct )
!!
!! DESCRIPTION
!!
!! This function is provided for the user's convenience and acts as a simple
!! wrapper to the Eos interface. The Eos interface uses a single, flexible data
!! structure "eosData" to pass the thermodynamic quantities in and out of the
!! funtion (see Eos). The wrapper hides formation and use of eosData
!! from the users.
!!
!! While Eos does not know anything about blocks, Eos_wrapped takes its
!! input thermodynamic state variables from a given block's storage area.
!! It works by taking a selected section of a block described by array
!! "range" and translating it to eosData before calling the Eos routine.
!! Upon return from Eos, Eos_wrapper updates certain state variables in
!! the same section of the block's storage area. Which variables are taken
!! as input, and which are updated, depends on the "mode" argument.
!!
!! If you want to return the derived quantities defined from EOS_VAR+1:EOS_NUM
!! in Eos.h, then you must use the direct interface Eos().
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   mode : determines which variables are used as Eos input.
!!          The valid values are MODE_DENS_EI (where density and internal
!!          energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!          MODE_DENS_TEMP (density and temperature are inputs).
!!          These quantities are defined in constants.h, the argument is 
!!          forwarded unchanged to the Eos function call.
!!          Note that internal energy is grid variable EINT_VAR, not ENER_VAR.
!!
!! 
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which Eos is to be applies. The example shows how
!!          the array describes the block section.
!!
!!   blockID: current block number
!!
!!   gridDataStruct : the grid data structure on whose data Eos is to be applied
!!
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then Eos is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However, if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then Eos is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_wrapped
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices can't
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos.h
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Eos_wrapped(mode,range,blockID, gridDataStruct)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkType
  use Logfile_interface, ONLY: Logfile_stampMessage 
  use Eos_interface, ONLY : Eos, Eos_putData, Eos_getData

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: mode
  integer, dimension(2,MDIM), intent(in) :: range
  integer,intent(in) :: blockID
  integer, optional, intent(IN) :: gridDataStruct

  real, pointer:: solnData(:,:,:,:)

#ifndef FIXEDBLOCKSIZE
  real, allocatable :: eosData(:),massFraction(:)
#else
  real, dimension(NSPECIES*MAXCELLS) :: massFraction
  real, dimension(EOS_NUM*MAXCELLS) :: eosData
#endif

  logical,target,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: ierr, istat, dataStruct
  integer :: i,j,k, vecLen
  integer,dimension(MDIM) :: pos
  integer :: diagFlag, blkType

!! ---------------------------------------------------------------------------------
  ! Test calling arguments
#define DEBUG
#ifdef DEBUG
  ierr = 1
  select case (mode)
  case (MODE_DENS_PRES)
     ierr = 0
  case (MODE_DENS_TEMP)
     ierr = 0
  case (MODE_DENS_EI)
     ierr = 0
  case (MODE_EOS_NOP,MODE_EOS_WRAPPERONLY)
     ierr = 0
  case (MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_ALL,MODE_DENS_TEMP_EQUI, MODE_DENS_TEMP_GATHER)
     ierr = 0
  case (MODE_DENS_TEMP_ION,MODE_DENS_TEMP_ELE,MODE_DENS_TEMP_RAD)
     ierr = 0
  case (MODE_DENS_EI_COMP,MODE_DENS_EI_ALL,MODE_DENS_EI_SCATTER,MODE_DENS_EI_GATHER)
     ierr = 0
  case (MODE_DENS_EI_SELE_GATHER,MODE_DENS_EI_SELERAD_GATHER,MODE_DENS_EI_SHOCKSELE_GATHER)
     ierr = 0
  case (MODE_DENS_EI_RECAL_GATHER)
     ierr = 0
  end select

  if(ierr /= 0) then
     call Driver_abortFlash("[Eos_wrapped] "//&
          "invalid mode: must be MODE_DENS_PRES, MODE_DENS_TEMP, MODE_DENS_EI, or variants thereof, or MODE_EOS_NOP")
  end if
#endif

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *

  vecLen = range(HIGH,IAXIS)-range(LOW,IAXIS)+1
  if (vecLen==0) return ! * Return immediately for empty IAXIS range! (for efficiency and avoiding index range errors)

  ! Initializations:   grab the solution data from UNK (or other data structure)
  !   and determine the length of the data being operated upon

  if(present(gridDataStruct))then
     dataStruct=gridDataStruct
  else
     dataStruct=CENTER
  end if
  call Grid_getBlkPtr(blockID,solnData,dataStruct)

#ifndef FIXEDBLOCKSIZE
  allocate(massFraction(NSPECIES*vecLen))
  allocate(eosData(EOS_NUM*vecLen))
#endif

  eosMask = .FALSE.
  eosMask(EOS_EINTION) = .TRUE.
  eosMask(EOS_EINTELE) = .TRUE.
  eosMask(EOS_EINTRAD) = .TRUE.
  eosMask(EOS_PRESION) = .TRUE.
  eosMask(EOS_PRESELE) = .TRUE.
  eosMask(EOS_PRESRAD) = .TRUE.
#ifdef EOS_ENTRELE
  eosMask(EOS_ENTRELE) = .TRUE.
#endif
#ifdef EOS_ENTRRAD
  eosMask(EOS_ENTRRAD) = .TRUE.
#endif


  do k = range(LOW,KAXIS), range(HIGH,KAXIS)
     do j = range(LOW,JAXIS), range(HIGH,JAXIS)

        pos(JAXIS)=j
        pos(KAXIS)=k
        pos(IAXIS)=range(LOW,IAXIS)
        
        select case (mode)
        case(MODE_EOS_WRAPPERONLY)
           ! No Call At All!
        case(MODE_DENS_EI_SHOCKSELE_GATHER)
           do i = 1,vecLen
              call Eos_getData(IAXIS,pos,1,solnData,dataStruct,eosData,massFraction)
#ifdef SHKS_VAR
              if (solnData(SHKS_VAR,i,j,k) > 3.5) then
                 call Eos(MODE_DENS_EI_SELE_GATHER,1,eosData,massFraction,mask=eosMask)
              else
#endif
                 call Eos(MODE_DENS_EI_GATHER,1,eosData,massFraction,mask=eosMask)
#ifdef SHKS_VAR
              end if
#endif
              call Eos_putData(IAXIS,pos,1,solnData,dataStruct,eosData)
              pos(IAXIS) = pos(IAXIS) + 1 
           end do
        case default
           call Eos_getData(IAXIS,pos,vecLen,solnData,dataStruct,eosData,massFraction)
           call Eos(mode,vecLen,eosData,massFraction,mask=eosMask,diagFlag=diagFlag)
           if (diagFlag .NE. 0) then
              call Grid_getBlkType(blockID, blkType)
              print*,'diagFlag!',range(:,1:NDIM),' blk',blockID, blkType
           end if
           call Eos_putData(IAXIS,pos,vecLen,solnData,dataStruct,eosData)
        end select


        ! The following is now done in Eos_putData:
!!$        solnData(GAME_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(pres+1:pres+veclen)/&
!!$             (eosData(eint+1:eint+veclen) *eosData(dens+1:dens+veclen)) +1

     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,dataStruct)

#ifndef FIXEDBLOCKSIZE
  deallocate(eosData)
  deallocate(massFraction)
#endif
  return
end subroutine Eos_wrapped

!!$        do i = 1,vecLen
!!$           massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
!!$                solnData(SPECIES_BEGIN:SPECIES_END,range(LOW,IAXIS)+i-1,j,k)
!!$        end do

!!$#ifdef EINT_VAR
!!$        energyInternal(1:vecLen) = solnData(EINT_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
!!$        
!!$        do i = 1,vecLen
!!$           if (solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) > &
!!$                (1.+ eos_eintSwitch)*eosData(ekin+i)) then
!!$              energyInternal(i) = solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) - eosData(ekin+i)
!!$           end if
!!$           energyInternal(i) = max(energyInternal(i), eos_smalle)
!!$        end do
!!$#else
!!$        do i = 1,vecLen
!!$           energyInternal(i) = solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) - eosData(ekin+i)
!!$           energyInternal(i) = max(energyInternal(i), eos_smalle)
!!$        end do
!!$#endif
!!$
!!$
!!$        eosData(eint+1:eint+vecLen) = energyInternal(1:vecLen)


!!$        eosData(ekin+1:ekin+vecLen) = 0.5*(solnData(VELX_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2 + &
!!$             &                                  solnData(VELY_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2 + &
!!$             &                                  solnData(VELZ_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2)


!!$        eosData(pres+1:pres+vecLen) = &
!!$             solnData(PRES_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
!!$        eosData(dens+1:dens+vecLen) = &
!!$             solnData(DENS_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
!!$        eosData(temp+1:temp+vecLen) = &
!!$             solnData(TEMP_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
!!$        eosData(gamc+1:gamc+vecLen) = &
!!$             solnData(GAMC_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)


        ! check for zero values before calculating gamma
!!$        iFlag = 0
!!$        where ( (eosData(eint+1:eint+vecLen) .eq. 0.) .or. (eosData(dens+1:dens+vecLen) .eq. 0.))
!!$           iFlag(1:vecLen) = 1
!!$        end where
!!$
!!$        !maybe there was a wrong flag set
!!$        if (maxval(iFlag) .gt. 0) then
!!$           
!!$           if (eos_meshMe .EQ. MASTER_PE) then
!!$              write(*,*) "ERROR After calling Eos, eosData(EOS_EINT) or eosData(EOS_DENS) are zero"
!!$              write(*,*) "  Perhaps the initialization routine is wrong..... or"
!!$              write(*,*) "  perhaps the runtime parameter eosMode is wrong."
!!$              write(*,*) "  This routine Eos_wrapped was called with mode= ", mode
!!$              write(*,*) "     Check constants.h to determine value of MODE_DENS_??"
!!$           endif
!!$           call Logfile_stampMessage('[Eos_wrapped] ERROR Density or Internal Energy are zero after a call to EOS!')
!!$           call Driver_abortFlash('[Eos_wrapped] ERROR Density or Internal Energy are zero after a call to EOS!')
!!$        end if



!!$
!!$        solnData(PRES_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(pres+1:pres+vecLen)
!!$        solnData(TEMP_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(temp+1:temp+vecLen)
!!$        solnData(GAMC_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(gamc+1:gamc+vecLen)
!!$#ifdef EINT_VAR
!!$        solnData(EINT_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(eint+1:eint+veclen)
!!$#endif
!!$        solnData(ENER_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(eint+1:eint+veclen) + eosData(ekin+1:ekin+vecLen)
!!$#ifdef ENTR_VAR
!!$        solnData(ENTR_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
!!$             eosData(entr+1:entr+veclen)
!!$#endif

