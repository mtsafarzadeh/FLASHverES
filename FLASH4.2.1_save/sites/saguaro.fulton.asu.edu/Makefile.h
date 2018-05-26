# FLASH makefile definitions for the Intel ifc compilers on Linux

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#HDF4_PATH = /opt/HDF-4.1r4
#HDF5_PATH = $(CURDIR)/../hdf5-1.6.6/hdf5

ZLIB_PATH  =
#MPI_PATH   = /packages/intel-mvapich-1.0.1
PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH =
MPE_PATH   =

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  Do not call the MPI wrappers here. Use the native compiler and specify 
#  MPI_PATH. This allows us to replace mpi with an equivalent later.
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpicxx
LINK    = mpif90

#FCOMP   = ${MPI_PATH}/bin/mpif90
#CCOMP   = ${MPI_PATH}/bin/mpicc
#CPPCOMP = ${MPI_PATH}/bin/mpicc
#LINK    = ${MPI_PATH}/bin/mpif90

# pre-processor flag
PP     = -D

#-----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized code
#  code ("-opt"), one for debugging ("-debug"), and one for testing ("-test").
#  Passing these flags to the setup script will cause the value associated with
#  the corresponding keys (i.e. those ending in "_OPT", "_DEBUG", or "_TEST") to
#  be incorporated into the final Makefile. For example, passing "-opt" to the
#  setup script will cause the flags following "FFLAGS_OPT" to be assigned to
#  "FFLAGS" in the final Makefile. If none of these flags are passed, the default
#  behavior will match that of the "-opt" flag.
#  In general, "-opt" is meant to optimize compilation and linking. "-debug"
#  should enable runtime bounds checking, debugger symbols, and other compiler-
#  specific debugging options. "-test" is useful for testing different
#  combinations of compiler flags particular to your individual system.
#----------------------------------------------------------------------------

OPENMP = -openmp
FFLAGS_OPT   =  -c -r8 -i4 -O3 -real_size 64
FFLAGS_DEBUG =  -c -r8 -i4 -real_size 64
FFLAGS_TEST  =  -c -r8 -i4 -real_size 64

FFLAGS_MPI  = -I $(MPI_PATH)/include

F90FLAGS = -I ${HDF5_DIR}/include -DH5_USE_16_API

CFLAGS_OPT   = -c -O2 -D_LARGEFILE64_SOURCE -malign-double
CFLAGS_DEBUG = -c -g -D_LARGEFILE64_SOURCE -malign-double
CFLAGS_TEST  = -c -02 -D_LARGEFILE64_SOURCE -malign-double

CFLAGS_HDF5 = -I${HDF5_DIR}/include -DH5_USE_16_API
CFLAGS_MPI  = -I${MPICH_HOME}/include
CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

#LFLAGS_OPT   = -r8 -i4 -Vaxlib -o
#LFLAGS_DEBUG = -r8 -i4 -Vaxlib -o
#LFLAGS_TEST  = -r8 -i4 -Vaxlib -o
LFLAGS_OPT   = -r8 -i4 -o
LFLAGS_DEBUG = -r8 -i4 -g -o
LFLAGS_TEST  = -r8 -i4 -o


#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_OPT   = 
LIB_DEBUG = 
LIB_TEST  =

LIB_HDF4  = -L$(HDF4_PATH)/lib -lmfhdf -ldf -lz -ljpeg
LIB_HDF5  = -L${HDF5_DIR}/lib -lhdf5_fortran -lhdf5 -lz
LIB_MPI   = 
#LIB_MPI   = -I $(MPI_PATH)/include -L$(MPI_PATH)/lib -lmpich -L/usr/local/ofed/lib64/ -libverbs -lpthread -libumad -lrt

LIB_NCMPI =
LIB_MPE   =

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ = 


#----------------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------------- 

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo
