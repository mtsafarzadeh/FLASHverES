!!****if* source/physics/sourceTerms/Stir/StirMain/StirScalar/st_calcPhases
!!    SS : added a directory StirScalar inside Stir 
!! NAME
!!
!!  st_calcPhases
!!
!! SYNOPSIS
!!
!!  st_calcPhases()
!!
!! DESCRIPTION
!!
!!     This routine updates the stirring phases from the OU phases.
!!     It copies them over, then subtracts out the divergence part.
!!
!!     SS : 27th Dec.12 : modified to include purely solenoidal forcing now
!!          In future, generalize it to include either mixed or purely compressive 
!!          and purely solenoidal forcings.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***
subroutine st_calcPhases()

  use Stir_data, ONLY : st_nmodes, st_mode, st_aka, st_akb, st_OUphases

#include "Flash.h"

  implicit none

!! SS : added one more variable kc
  real bjiR, bjiI, kb, kc, kk
  integer i,j

  do i = 1, st_nmodes 
     kb = 0.
     kc = 0. !! initialise kc
     kk = 0. 
     do j=1,NDIM
        kk = kk + st_mode(j,i)*st_mode(j,i) 
        kb = kb + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+0+1)
        kc = kc + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+1+1)
     enddo
     do j=1,NDIM
        bjiR = st_OUphases(6*(i-1)+2*(j-1) + 0 + 1)
        bjiI = st_OUphases(6*(i-1)+2*(j-1) + 1 + 1)

       !! Original FLASH4 expressions 
	!st_aka(j,i) = bjiR - st_mode(j,i)*kb/kk
        !st_akb(j,i) = bjiI 

        !! SS : Purely solenoidal forcing now. In future, generalize this to include mixed, 
        !!      purely compressive and solenoidal forcings
        st_aka(j,i) = bjiR - st_mode(j,i)*kb/kk
        st_akb(j,i) = bjiI - st_mode(j,i)*kc/kk
     enddo
  enddo

  return
end subroutine st_calcPhases
