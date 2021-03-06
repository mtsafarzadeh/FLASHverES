!!****if* source/Grid/GridSolvers/Multipole/gr_mpoleCenterOfMass
!!
!! NAME
!!
!!  gr_mpoleCenterOfMass
!!
!! SYNOPSIS
!!
!!  gr_mpoleCenterOfMass(integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes the centroid of the specified variable and returns
!!  its location in the mpole_common variables Xcm, Ycm, and Zcm.
!!  Also computes the total value of the quantity and leaves it
!!  in the variable Mtot.  If Mtot=0, the routine aborts.
!!
!! ARGUMENTS
!!
!!  idensvar -- the index of the density variable
!!
!!***
#ifdef DEBUG_ALL
#define DEBUG_GRAVITY
#endif


subroutine gr_mpoleCenterOfMass (idensvar)

  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use gr_mpoleData, ONLY : G_2DSPHERICAL, Xcm, Ycm,Zcm,mpole_geometry,Mtot

  implicit none
  
#include "constants.h"
#include "Flash.h"

  include "Flash_mpi.h"
  
  integer, intent(IN)  :: idensvar
  
  integer, parameter :: nsum = 4    ! Number of summed quantities;
  ! coordinate w/mpole_sum_local
  real    :: sum(nsum), lsum(nsum)
  
  integer :: blockCount, blockList(MAXBLOCKS)
  integer :: error, lb
  character(len=124) :: str_buffer
  
  !==========================================================================
  

  ! Sum quantities over all locally held leaf blocks.


  if(mpole_geometry /= G_2DSPHERICAL) then
     

     sum  = 0.
     lsum = 0.

     call Grid_getListOfBlocks(LEAF, blockList, blockCount)     

     do lb = 1, blockCount
        call gr_mpoleLocalSum (blockList(lb), nsum,  idensvar, lsum )
     enddo
     
       !Give all processors a copy of the global sums.
     call MPI_AllReduce (lsum, sum, nsum, FLASH_REAL, & 
          MPI_Sum, gr_meshComm, error)
     
  
     !Now normalize the center-of-mass coordinates.
     
     Mtot = sum(1)
     
     if (abs(Mtot) < TINY(Mtot)) &
          call Driver_abortFlash("FATAL ERROR: mpole_centerOfMass:  Mtot = 0")
     
     if (Mtot < 0.) then
        write (str_buffer,*) 'mpole_centerOfMass:  Mtot = ', Mtot
        call Logfile_stamp( str_buffer, 'warning')
     endif
     
     Xcm = sum(2) / Mtot
     Ycm = sum(3) / Mtot
     Zcm = sum(4) / Mtot
  end if
#ifdef DEBUG_GRAVITY
  print*,'center of mass found'
#endif
  return
end subroutine gr_mpoleCenterOfMass
