# FLASH makefile definitions for x86-64 Linux (GNU compilers)
#
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

LIB_BASE = /export/lustre_1/home/tfitz/compiledLibs/intel2013

#MPI_PATH = /usr/local
MPI_PATH = /cell_root/software/openmpi/1.6/intel/sys

#MPI_PATH  = /usr/local/mpich2dbg
HDF4_PATH  =
HDF5_PATH  = $(LIB_BASE)



ZLIB_PATH  = 


PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH = 
MPE_PATH   =

BLAS_PATH  = $(LIB_BASE)
HYPRE_PATH = $(LIB_BASE)/hypre2.9
SUPERLU_PATH = $(LIB_BASE)

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


FFLAGS_OPT   = -c -r8 -i4 -O3 -real_size 64 -diag-disable 10120
FFLAGS_DEBUG = -c -g -r8 -i4 -O0 -check bounds -check format \
-check output_conversion -warn all -warn error -real_size 64 -check uninit \
-traceback -fp-stack-check -diag-disable 10120 -fpe0 -check pointers
FFLAGS_TEST  = ${FFLAGS_OPT} -fp-model precise
FFLAGS_HYPRE = -I${HYPRE_PATH}/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include

F90FLAGS = -I ${HDF5_PATH}/include -DH5_USE_16_API \
-c -O3 -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 \
-diag-disable 10120
CFLAGS_OPT   = -c -O3 -diag-disable 10120 -D_LARGEFILE64_SOURCE 
CFLAGS_DEBUG = -c -O0 -g -traceback -debug all -debug extended \
-D_LARGEFILE64_SOURCE -diag-disable 10120 -ftrapuv -fp-stack-check
CFLAGS_TEST  = ${CFLAGS_OPT} -fp-model precise

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_16_API
CFLAGS_NCMPI = -I${NCMPI_PATH}/include

CFLAGS_BLAS = -I${BLAS_PATH}/include
FFLAGS_BLAS = -I${BLAS_PATH}/include

CFLAGS_HYPRE = -I${HYPRE_PATH}/include
FFLAGS_HYPRE = -I${HYPRE_PATH}/include

CFLAGS_SUPERLU = -I${SUPERLU_PATH}/include
FFLAGS_SUPERLU = -I${SUPERLU_PATH}/include

CFLAGS_MPI   = -I$(MPI_PATH)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -diag-disable 10120 -O3 -o
LFLAGS_DEBUG = -diag-disable 10120 -o
LFLAGS_TEST  = -diag-disable 10120 -O2 -o


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
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5 -lz -lsz

LIB_PAPI  =
LIB_MATH  = -limf -lm

LIB_MPI   = 
LIB_NCMPI = -L${NCMPI_PATH}/lib -lpnetcdf
LIB_MPE   =

LIB_BLAS  = -L${BLAS_PATH}/lib -lopenblas
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE  
LIB_SUPERLU = -L${SUPERLU_PATH}/lib -lsuperlu_4.3

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
AR = xiar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo

#----------------------------------------------------------------------------
# ifort hacks
#----------------------------------------------------------------------------

ifeq ($(FLASHBINARY),true)

FFLAGS_WO_WARNALL = $(patsubst -warn all,,$(FFLAGS))

#Turn off compiler error messages for paramesh files that use wrapper
#functions such as MPI_int_SSEND.
amr_migrate_tree_data.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_test_neigh_values.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_checkpoint_default.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_morton.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(F90FLAGS) $(FDEFINES)	$<

#Fortran 77 source lines exceed 72 characters
umap.o : %.o : %.F
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
fftsg.o : %.o : %.f
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
fftsg3d.o : %.o : %.f
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<

#Files mix and match assumed shape arrays, assumed size arrays
#and scalars in function calls.  This is fine but it is viewed as
#a problem when using strict type checking compiler options.
fftpack.o : %.o : %.f90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
gr_pfftDcftForward.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<
gr_pfftDcftInverse.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_WARNALL) $(FDEFINES)	$<

#ifort version 12.1.0 hangs during compilation of hy_ppm_sweep.F90
#unless we use -O0
hy_ppm_sweep.o : %.o : %.F90
	$(FCOMP) $(FFLAGS) -O0 $(F90FLAGS) $(FDEFINES)	$<

#ifort version 12.1.0 generates bad code for Grid_advanceDiffusion.F90
#from Grid/GridSolvers/HYPRE/Grid_advanceDiffusion.F90 when we use the
#-openmp option: this is very strange because this file contains no openmp.
FFLAGS_WO_OPENMP = $(patsubst $(OPENMP),-recursive -reentrancy threaded,$(FFLAGS))
Grid_advanceDiffusion.o : %.o : %.F90
	$(FCOMP) $(FFLAGS_WO_OPENMP) $(F90FLAGS) $(FDEFINES)	$<

endif

