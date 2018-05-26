subroutine eos_tabLogOutsideCounts(force)

  use eos_data,                   ONLY : eos_meshMe
  use eos_tabData,                ONLY : EOS_TAB_NALLTAB,           &
                                         eos_allTab, &
                                         eos_tabAllDiag

  use Logfile_interface

  implicit none

#include "constants.h"
#include "Flash.h"

  logical, intent(in) :: force

  integer :: species
  logical :: tempIsLog, densIsLog
  real    :: temp, dens

  real    :: physTemp, physDens
  character(len=16) :: a1,a2
  character(len=MAX_STRING_LENGTH) :: str1,str2,str3


  do species = 1,NSPECIES
     if (eos_tabAllDiag(species)%highTempCount > 0) then
400     format(I3)
401     format(I16)
        write(a1,400) species
        write(a2,401) eos_tabAllDiag(species)%highTempCount
501     format('Species',A,': Lookup for T > Tmax(table): ',A,' times, ') 
502     format('first for (',1P,G13.7,',',G13.7,'), ')
503     format('highest T was ',1P,G13.7)
        write(str1,501) trim(a1), trim(adjustl(a2))
        write(str2,502) eos_tabAllDiag(species)%firstHighTempEvent%temp, &
                        eos_tabAllDiag(species)%firstHighTempEvent%dens
        write(str3,503) eos_tabAllDiag(species)%highestTemp
        call Logfile_stampMessage(trim(str1)//trim(str2)//str3,force)
      end if
   end do
  do species = 1,NSPECIES
     if (eos_tabAllDiag(species)%highDensCount > 0) then
        write(a1,400) species
        write(a2,401) eos_tabAllDiag(species)%highDensCount

601     format('Species',A,': Lookup dens > densMax(table): ',A,' times, ') 
603     format('highest ion dens was ',1P,G13.7)
        write(str1,601) trim(a1), trim(adjustl(a2))
        write(str2,502) eos_tabAllDiag(species)%firstHighDensEvent%temp, &
                        eos_tabAllDiag(species)%firstHighDensEvent%dens
        write(str3,603) eos_tabAllDiag(species)%highestDens
        call Logfile_stampMessage(trim(str1)//trim(str2)//str3,force)
     end if
  end do
  

end subroutine eos_tabLogOutsideCounts
