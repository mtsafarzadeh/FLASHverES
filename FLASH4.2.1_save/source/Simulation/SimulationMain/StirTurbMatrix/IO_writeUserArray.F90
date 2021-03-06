!!****if* source/Simulation/SimulationMain/StirTurb/IO_writeUserArray
!!
!!
!!  NAME
!!    IO_writeUserArray
!!
!!  SYNOPSIS
!!    call IO_writeUserArray()
!!
!!
!!  DESCRIPTION
!!
!!    This is the supplied interface for users to write out additional
!!    quantities to the checkpoint or plotfile.  This routine should be used
!!    for writing out various types of arrays.  If the array is a global quantity
!!    only the master processor needs to write out the data.  If it is a quantity 
!!    which is different on all processors then each processor must write out its 
!!    own section of the array. (For a serial IO implementation each processor would
!!    need to send its data to the master.)  The specific implementation is left up
!!    to the user.  
!!
!!    In each case the user should make a call to either 
!!    io_h5write_generic_int_arr (hdf5) or io_ncmpi_write_generic_iarr (pnetcdf) or
!!    io_h5write_generic_real_arr (hdf5) or io_ncmpi_write_generic_rarr (pnetcdf)
!!    depending on the io implementation.
!!  
!!  ARGUMENTS
!!    
!!
!!  NOTES 
!!
!!    This routine should NOT
!!    be used to write out grid scope data or to write out single scalar
!!    values.  To write out user defined grid scope variables the user should
!!    use the keyword 'GRIDVAR' to declare a grid scope variable in the Config
!!    files.  Then set the runtime parameters plot_grid_var_1, plot_grid_var_2,   
!!    to the name of the grid var to include them in the checkpoint files and
!!    plotfiles.
!!
!!    To write out single scalar quantities the use the IO_setScalar routine to
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

subroutine IO_writeUserArray ()
  
  use Stir_data, ONLY : st_mode, st_aka, st_akb, st_randseed, st_seedLen, &
       st_OUphases, st_maxmodes

  use IO_data, ONLY : io_chkptFileID, io_globalMe

  implicit none

#include "constants.h"


  integer :: offset, datasetNameLen


  !!For the StirTurb problem we need to write out mode, aka, akb, OUPhases and randseed
  !!In this case, these values are global quantities so only master processor needs to write
  !!out values.  This is handled within the h5 routines.  We send in the local size and total
  !!size of the datasets to be the same.
   
  offset = 0
  datasetNameLen = 8 !this makes things a lot easier because string handling from fortran to c is
  !a big pain
  call io_h5write_generic_int_arr(io_globalMe, &
       io_chkptFileID, &
       st_randseed, &
       st_seedLen, &
       st_seedLen, &
       offset, &
       "randseed", &
       datasetNameLen)
  
  !for mulitdimensional arrays we just write it as a single dim
  datasetNameLen = 4
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       st_mode, &
       3*st_maxmodes, &
       3*st_maxmodes, &
       offset, &
       "mode", &
       datasetNameLen)

  datasetNameLen = 3
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       st_aka, &
       3*st_maxmodes, &
       3*st_maxmodes, &
       offset, &
       "aka", &
       datasetNameLen)
  
  datasetNameLen = 3
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       st_akb, &
       3*st_maxmodes, &
       3*st_maxmodes, &
       offset, &
       "akb", &
       datasetNameLen)
  
  
  datasetNameLen = 8
  call io_h5write_generic_real_arr(io_globalMe, &
       io_chkptFileID, &
       st_OUphases, &
       6*st_maxmodes, &
       6*st_maxmodes, &
       offset, &
       "ouphases", &
       datasetNameLen)
  
  


end subroutine IO_writeUserArray
