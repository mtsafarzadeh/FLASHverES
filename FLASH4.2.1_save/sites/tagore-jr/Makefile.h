# FLASH makefile definitions for ix86-64 Linux (gfortran compiler)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------
#HDF4_PATH =
HDF5_PATH = /opt/local/var/macports/software/hdf5/1.6.9_0/opt/local
#LIB_HDF5 = ${HDF5_PATH}/lib

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

LIB_NCMPI = /usr/local/pkgs/pnetcdf-1.1.0/gcc
MPE_PATH   =
MPI_PATH = /usr/local/pkgs/mpich2-1.1.1p1/gfortran/bin

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------
FCOMP   = $(MPI_PATH)/mpif90
CCOMP   = $(MPI_PATH)/mpicc
CPPCOMP = $(MPI_PATH)/mpiCC
LINK    = $(MPI_PATH)/mpif90

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for
#  flash_test, and is set for quick code generation, and (sometimes)
#  profiling.  The Makefile generated by setup will assign the generic token
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

FFLAGS_OPT = -c -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
FFLAGS_DEBUG = -g -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none
FFLAGS_TEST = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none

F90FLAGS =

CFLAGS_OPT = -O2  -c
CFLAGS_DEBUG = -g -c
CFLAGS_TEST = -c

# Platform symbol
CDEFINES += -DDarwin

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I ${HDF5_PATH}/include

CFLAGS_NCMPI = -I$(LIB_NCMPI)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -o
LFLAGS_TEST  = -o


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

#LIB_HDF4  = -lmfhdf -ldf -ljpeg -lz
#LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 /usr/lib64/libz.a
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5

LIB_PAPI  =
#LIB_MATH  = -ldfftw -ldrfftw

LIB_MPI   = 
#LIB_NCMPI = -L $(NCMPI_PATH)/lib -lpnetcdf
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

