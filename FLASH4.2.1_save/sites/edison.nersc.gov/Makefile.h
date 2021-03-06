# FLASH makefile definitions for Edison, a Cray XC30 machine at NERSC.
#
# The software stack is managed by the modules system.  The default
# programming environment is PrgEnv-intel/5.0.41.
#
# Load the following modules before compiling (names as of 6th Aug 13):
# * HDF5: module load cray-hdf5-parallel/1.8.11
# * PNETCDF: module load cray-parallel-netcdf/1.3.1.1

#----------------------------------------------------------------------------
# Compiler wrapper scripts
#----------------------------------------------------------------------------

FCOMP   =  ftn
CCOMP   =  cc
CPPCOMP =  CC
LINK    =  ftn


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

# Suppress the following compiler warning message:
# icc: command line warning #10120: overriding '-O2' with '-O3'

OPENMP = -openmp
FFLAGS_OPT   = -c -g -r8 -i4 -O3 -real_size 64 -diag-disable 10120
FFLAGS_DEBUG = -c -g -r8 -i4 -O0 -check bounds -check format \
-check output_conversion -warn all -warn error -real_size 64 -check uninit \
-traceback -fp-stack-check -diag-disable 10120 -fpe0 -check pointers
FFLAGS_TEST  = ${FFLAGS_OPT} -fp-model precise

CFLAGS_OPT   = -c -O3 -g -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 \
-diag-disable 10120
CFLAGS_DEBUG = -c -O0 -g -traceback -debug all -debug extended \
-D_LARGEFILE64_SOURCE -diag-disable 10120 -ftrapuv -fp-stack-check
CFLAGS_TEST  = ${CFLAGS_OPT} -fp-model precise

# No MPI path required because we are using compiler wrapper scripts
FFLAGS_MPI   =
CFLAGS_MPI   =

# No need to specify the path to the HDF5 header files because it is handled
# by the compiler wrapper scripts
CFLAGS_HDF5  = -DH5_USE_16_API

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -diag-disable 10120 -O3 -o
LFLAGS_DEBUG = -diag-disable 10120 -o
LFLAGS_TEST  = ${LFLAGS_OPT}

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

# No need to explicitly link against the HDF5 library because it is handled
# by the compiler wrapper scripts

LIB_HDF5    =
LIB_NCMPI   =
LIB_MPI     =

LIB_OPT     =
LIB_DEBUG   =
LIB_TEST    =
CONFIG_LIB  = ${IPM_INTEL}

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

MV    = mv -f
AR    = ar -r
RM    = rm -f
CD    = cd
RL    = ranlib
ECHO  = echo
AWK   = awk
CAT   = cat


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
