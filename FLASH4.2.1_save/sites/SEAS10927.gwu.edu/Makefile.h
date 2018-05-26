# FLASH makefile definitions for x86-64 Linux (GNU compilers)
#
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH = /usr/local/mpich2
#MPI_PATH  = /usr/local/mpich2dbg
HDF4_PATH  =
#HDF5_PATH  = /Users/mvanella/Documents/Software/hdf5-1.8.6
HDF5_PATH  = /usr/local/hdf5-1.8.10

ZLIB_PATH  = /usr/local/lib

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH = 
MPE_PATH   =

BLAS_PATH  = /usr/local
HYPRE_PATH = /usr/local/Hypre2.9
SUPERLU_PATH = /Users/mvanella/Documents/Software/SuperLU_4.3
SUPERLUDIST_PATH = /usr/local/SuperLU_DIST_3.2
PARMETIS_PATH = /usr/local

export cur-dir := $(shell pwd)

# Set the location of top directory
export setup_dir = $(cur-dir)


#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------

FCOMP   = ${MPI_PATH}/bin/mpif90
CCOMP   = ${MPI_PATH}/bin/mpicc
CPPCOMP = ${MPI_PATH}/bin/mpicxx
LINK    = ${MPI_PATH}/bin/mpif90

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


FFLAGS_OPT =  -c -O2 -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -Wuninitialized

FFLAGS_DEBUG = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -pedantic -Wall -Wextra -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check \
-fbacktrace -fdump-core -finit-real=nan \
-finit-integer=-999999 -fimplicit-none

FFLAGS_TEST =  -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none

F90FLAGS =


CFLAGS_OPT =  -c -O2 -Wuninitialized

CFLAGS_DEBUG =  -ggdb -c -Wno-div-by-zero -Wundef  \
-Wconversion -Wstrict-prototypes -Wunreachable-code \
-pedantic -Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal \
-Wunsafe-loop-optimizations -Wpadded -fstack-check -fstack-protector-all


CFLAGS_TEST = -c

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
FFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI = -I${NCMPI_PATH}/include

CFLAGS_BLAS = -I${BLAS_PATH}/include
FFLAGS_BLAS = -I${BLAS_PATH}/include

CFLAGS_HYPRE = -I${HYPRE_PATH}/include
FFLAGS_HYPRE = -I${HYPRE_PATH}/include

CFLAGS_SUPERLU = -I${SUPERLU_PATH}/include
FFLAGS_SUPERLU = -I${SUPERLU_PATH}/include

CFLAGS_PARMETIS = -I${PARMETIS_PATH}/include
FFLAGS_PARMETIS = -I${PARMETIS_PATH}/include

CFLAGS_SUPERLUDIST = -I${SUPERLUDIST_PATH}/include
FFLAGS_SUPERLUDIST = -I${SUPERLUDIST_PATH}/include



#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -o
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

LIB_HDF4  = 
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5  -lz -lsz

LIB_PAPI  =
LIB_MATH  = 

LIB_MPI   = 
LIB_NCMPI = -L${NCMPI_PATH}/lib -lpnetcdf
LIB_MPE   =

LIB_BLAS         = -L${BLAS_PATH}/lib -lopenblas
LIB_HYPRE        = -L${HYPRE_PATH}/lib -lHYPRE  
LIB_SUPERLU      = -L${SUPERLU_PATH}/lib -lsuperlu_4.3
LIB_SUPERLUDIST  = -L${SUPERLUDIST_PATH}/lib -lsuperlu_dist_3.2 
LIB_PARMETIS     = -L${PARMETIS_PATH}/lib -lparmetis -lmetis
LIB_STDCXX = -lstdc++

#Specify TEC_PLOT=YES in order to link the tec plot library.
TEC_PLOT=YES
ifeq ($(TEC_PLOT), YES)
CONFIG_LIB = -I${setup_dir}/../source/Simulation/SimulationMain/INavierStokes -L${setup_dir}/../source/Simulation/SimulationMain/INavierStokes -ltecio
endif

#HY_PRE=YES
#ifeq ($(HY_PRE), YES)
#CONFIG_LIB = $CONFIG_LIB " -I/usr/local/Hypre2.7/include -L/usr/local/Hypre2.7/lib -HYPRE" 
#endif

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

ifeq ($(FLASHBINARY),true)

FFLAGS_WO_WARNALL = $(patsubst -pedantic,,$(FFLAGS))

#Files mix and match assumed shape arrays, assumed size arrays
#and scalars in function calls.  This is fine but it is viewed as
#a problem when using strict type checking compiler options.
fftpack.o : %.o : %.f90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)       $<
gr_pfftDcftForward.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)       $<
gr_pfftDcftInverse.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)       $<

endif


#Configure lines:
#
# mpich2-1.2.1p1:
#./configure --prefix=/opt/mpich2/1.2.1p1/gcc-4.4.3
#--enable-error-checking=all --enable-error-messages=all 
#--with-pm=gforker:mpd --enable-g=dbg,meminit --enable-fast=defopt 
#--enable-f77 --enable-f90 --enable-cxx --enable-romio --enable-sharedlibs=gcc 
#--with-mpe CC=gcc F77=gfortran F90=gfortran CXX=g++ 2>&1 | tee ../mpich2_1.2.1p1_gcc-4.4.3_build.out
#
# hdf5-1.8.4-patch1:
#./configure --prefix=/opt/hdf5/1.8.4-patch1/gcc-4.4.3 
#CC=mpicc FC=mpif90 CXX=mpicxx --enable-production --enable-debug=all --enable-shared
#--enable-parallel --enable-using-memchecker 2>&1 | tee ../hdf5_1.8.4-patch1_gcc-4.4.3_build.out
#
# parallel-netcdf-1.1.1:
#./configure --prefix=/opt/parallel-netcdf/1.1.1/gcc-4.4.3 
#--enable-fortran --with-mpi=/opt/mpich2/1.2.1p1/gcc-4.4.3 
#CFLAGS="${CFLAGS} -g" 2>&1 | tee ../parallel-netcdf-1.1.1_gcc-4.4.3_build.out
#
# valgrind-3.5.0: (Need to patch because we have a new version of glibc)
#patch -Np0 -i ../valgrind_glibc211.diff || return 1
#autoreconf
#./configure --prefix=/opt/valgrind/3.5.0/gcc-4.4.3
#--without-mpicc 2>&1 | tee ../valgrind-3.5.0_gcc-4.4.3_build.out