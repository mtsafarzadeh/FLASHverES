!!****if* source/Simulation/SimulationMain/StirTurb/IO_readUserArray
!!
!!  NAME
!!    IO_readUserArray
!!
!!  SYNOPSIS
!!    call IO_readUserArray()
!!
!!
!!  DESCRIPTION
!!
!!    This is the supplied interface for users to read in additional
!!    quantities to the checkpoint or plotfile.  This routine should be used
!!    for reading in various types of arrays.  If the array is a global quantity
!!    only the master processor needs to read in the data.  If it is a quantity 
!!    which is different on all processors then each processor must read in its 
!!    own section of the array. (For a serial IO implementation each processor would
!!    need to send its data to the master.)  The specific implementation is left up
!!    to the user.  
!!
!!    In each case the user should make a call to either 
!!    io_h5read_generic_int_arr (hdf5) or io_ncmpi_read_generic_iarr (pnetcdf)  or
!!    io_h5read_generic_real_arr (hdf5) or io_ncmpi_read_generic_rarr (pnetcdf)
!!    depending on the io implementation.
!!  
!!  ARGUMENTS
!!    
!!  SS : This routine is the one modified by Evan for the StirTurb problem in FLASH4.0  
!!
!!  NOTES 
!!
!!    This routine should NOT
!!    be used to read in grid scope data or to read in single scalar
!!    values.  To read in user defined grid scope variables the user should
!!    use the keyword 'GRIDVAR' to declare a grid scope variable in the Config
!!    files.  Then set the runtime parameters plot_grid_var_1, plot_grid_var_2,   
!!    to the name of the grid var to include them in the checkpoint files and
!!    plotfiles.
!!
!!    To read in single scalar quantities the use the IO_setScalar routine to
!!    add a scalar to the scalar output list.
!!
!!
!!  SEE ALSO
!!
!!    io_h5write_generic_int_arr
!!    io_h5read_generic_real_arr
!!    IO_setScalar
!!    
!!    For the pnetcdf implementation see
!!    io_ncmpi_write_generic_iarr
!!    io_ncmpi_read_generic_rarr
!!
!!***



subroutine IO_readUserArray ()
  
  use Stir_data, ONLY : st_mode, st_aka, st_akb, st_randseed, st_seedLen, &
       st_OUphases, st_maxmodes

  use IO_data, ONLY : io_chkptFileID, io_globalMe, io_globalComm
  use ut_randomInterface, ONLY : ut_randomSeed
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  

  integer :: offset, datasetNameLen
  integer :: ierr

  !!For the StirTurb problem we need to read in mode, aka, akb, OUPhases and randseed

  !!Each processor reads in data so we don't have to do a broadcast

   
  !this assumes you are restarting on the same machine!!
  
  call ut_randomSeed (ut_size = st_seedLen) 
  if (.not. allocated (st_randseed) ) then
     allocate (st_randseed (st_seedLen) )
  endif
  
  
  if (io_globalMe .eq. MASTER_PE) then
     print *, 'seed length = ', st_seedLen
!ES
endif   ! SS : kept this endif 
  
! SS : added this loop separately 

 if (io_globalMe .eq. MASTER_PE) then   
  
  offset = 0
  datasetNameLen = 8 !this makes things a lot easier because string handling from fortran to c is
  !a big pain
  call io_h5read_generic_int_arr( &
       io_chkptFileID, &
       st_randseed, &
       st_seedLen, &
       st_seedLen, &
       offset, &
       "randseed", &
       datasetNameLen)
  
  
  !for mulitdimensional arrays we just read it as a single dim
  datasetNameLen = 4
  call io_h5read_generic_real_arr( &
       io_chkptFileID, &
       st_mode, &
       3*st_maxmodes, &
       3*st_maxmodes, &
       offset, &
       "mode", &
       datasetNameLen)
  
  datasetNameLen = 3
  call io_h5read_generic_real_arr( &
       io_chkptFileID, &
       st_aka, &
       3*st_maxmodes, &
       3*st_maxmodes, &
       offset, &
       "aka", &
       datasetNameLen)
  
  datasetNameLen = 3
  call io_h5read_generic_real_arr( &
       io_chkptFileID, &
       st_akb, &
       3*st_maxmodes, &
       3*st_maxmodes, &
       offset, &
       "akb", &
       datasetNameLen)
  
  
  datasetNameLen = 8
  call io_h5read_generic_real_arr( &
       io_chkptFileID, &
       st_OUphases, &
       6*st_maxmodes, &
       6*st_maxmodes, &
       offset, &
       "ouphases", &
       datasetNameLen)

!ES I moved the endif and added the Bcast stuff here  
  endif
  call MPI_Barrier(io_globalComm, ierr)
  call MPI_Bcast(st_randseed,     st_seedLen, FLASH_INTEGER, MASTER_PE, io_globalComm, ierr)
  call MPI_Bcast(st_mode,      3*st_maxmodes, FLASH_REAL, MASTER_PE, io_globalComm, ierr)
  call MPI_Bcast(st_aka,       3*st_maxmodes, FLASH_REAL, MASTER_PE, io_globalComm, ierr)
  call MPI_Bcast(st_akb,       3*st_maxmodes, FLASH_REAL, MASTER_PE, io_globalComm, ierr)
  call MPI_Bcast(st_OUphases,  6*st_maxmodes, FLASH_REAL, MASTER_PE, io_globalComm, ierr)

! Print some communication for testing 

  !print *,'Processor No. : ',io_globalMe
  !print *,'Print the random seed',st_randseed(1),st_randseed(2)
  !print *,'Print values for any two modes',st_mode(3,2),st_mode(2,10)
  !print *,'Print values for aka',st_aka(1,1),st_aka(2,1)
  !print *,'Print values for akb',st_akb(1,3),st_akb(2,4)
  !print *,'Here are the first two OUphases',st_OUphases(1),st_OUphases(2)

end subroutine IO_readUserArray
