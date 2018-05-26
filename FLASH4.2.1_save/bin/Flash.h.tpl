## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for Generating the BuildStamp script
## For syntax of this file see "Readme.template"
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## ppDefines         -> Pre Processor Defines
##
## variableDefines   -> Defines related to variable names
## speciesDefines    -> Defines related to various species
## massscalarsDefines-> Defines related to various Mass Scalars
## nmassscalar_groups-> Number of groups of Mass scalars 
## mfrac_specDefine  -> 
## facevarDefines    -> Defines related to FACE Variables
## fluxDefines       -> Defines related to FLUXes
## partpropDefines   -> Defines related to Particle Properties
## scratchvarDefines   -> Defines related to scratch variables
## scratchcentervarDefines -> Defines related to scratch center variables
## scratchfacexvarDefines  -> Defines related to scratch facex variables
## scratchfaceyvarDefines  -> Defines related to scratch facey variables
## scratchfacezvarDefines  -> Defines related to scratch facez variables
## particletypeDefines -> Defines related to particle variables
##
## nguard            -> # of guard cells
## fixedblocksize    -> 0/1 (are we using fixed block size?)
## ndim              -> # of dimensions (user specified)
## nxb,nyb,nzb       -> zones per block in each dimension (user specified)
## maxblocks         -> User specified value for MAXBLOCKS
## nscratchvars      -> # of scratch variables
## nscratchcentervars -> # of scratch center variables
## nscratchfacexvars  -> # of scratch facex variables
## nscratchfaceyvars  -> # of scratch facey variables
## nscratchfacezvars  -> # of scratch facez variables
## nrealprops        -> # of real properties (of particles?)
## nparticletypes    -> # of particle types (passive, active, star, etc.)
## nflux             ->
## nfacevars         -> # of face variables
## nvars             -> # of Property variables
## nspecies          -> # of species
## nmass_scalars     -> # of mass scalars
## max_plot_vars     -> # of plot_var_N runtime parameters that can be used
##
## gridGeometry      -> grid geometry as specified at setup time, as one of the
##                       following unquoted strings:
##                       GRID_GEOM_CARTESIAN
##                       GRID_GEOM_SPHERICAL
##                       GRID_GEOM_CYLINDRICAL
##                       GRID_GEOM_POLAR
##                       GRID_GEOM_UNDEF
## curvilinear       -> whether curvilinear support was requested at setup time
##                       (0 or 1)
##
## NO other variables may be used with %(VARNAME)s syntax
##

#if 0
 *************************************************************************
  File:       Flash.h, generated by setup script
  Contains:   problem dependent parameters.

  Example:    When the setup script parses the Config files it
              counts up the number of declared grid variables.  It records 
              the final number in this file.  It also defines the indices 
              at which the vars are located.

              In other cases, some parameters need to be specified at 
              compile time.
              If a grid with a fixed block size is used (Paramesh or fixed 
              block size Uniform Grid) MAXBLOCKS, NXB, NYB, NZB, NGUARD etc 
              are defined in this Flash.h file so they can be used at 
              compile time.

  Note:       parameters which can be computed from other parameters. These 
              computations are now performed at the Flash.h level and not by 
              the setup script. The big advantage to this is that one can 
              tweak some entries of Flash.h and dependent quantities get 
              changed appropriately

              If you have some parameter computed from others enclose the 
              computation in parenthesis
              
              #define NEW (OLD + 5) and not #define NEW OLD + 5

              Even though both are valid syntax the second can produce 
              unexpected results. e.g. suppose OLD=2, then NEW*6 will 
              evaluate to 42 in the first case and 32 in the second
              
  Postscript: Since this header file is to be used with C, F77 and F90 
              there are lot of nice things we cannot do. 

              (1) No indendation since F77 will choke (so inner ifdef will 
                  not be indented)
              (2) All Comments will have to be in (#if 0,stuff, #endif ) 
                  construct. C, F77 and F90 do not agree on a common comment 
                  character
              (3) Hence this file is sprinkled with comments indicating 
                  start and finish of major ifdef statements

  !!!!Do not edit!  See problem config files or modify setup options to 
      change simulation parameters!!!!
  *********************************************************************
#endif

#ifndef __FLASH_HEADER_FILE__
#define __FLASH_HEADER_FILE__

#define CONSTANT_ZERO (0)
#define CONSTANT_ONE  (1)
#define CONSTANT_TWO  (2)

#define STRINGIFY_VAL(x) STRINGIFY_LIT(x)
#define STRINGIFY_LIT(x) #x
#define FILE_AT_LINE_C __FILE__ "@" STRINGIFY_VAL(__LINE__)

#ifdef __INTEL_COMPILER
#define FILE_AT_LINE (__FILE__ // "@" // STRINGIFY_VAL(__LINE__))
#else
#define FILE_AT_LINE __FILE__
#endif

#if 0
  **************************************************************************
  The #define definitions below indicate the index at which a grid scope 
  variable is located in the 'unk' (unknowns) data structure.  Property 
  variables are denoted with '_VAR' following the variable name.  Species 
  are also stored in the unk datastructure.  They are
  denoted with '_SPEC' following the species name.

  NPROP_VARS = the total number of grid property variables in the simulation
               (ie those ending in _VAR)
  NSPECIES = the number of species in the simulation 
             (ie those ending in _SPEC)
  NMASS_SCALARS = the number of mass scalars in a simulation
  NMASS_SCALAR_GROUPS = # of groups mass scalars to be renormed are split into
  NUNK_VARS = NPROP_VARS + NSPECIES + NMASS_SCALARS
  SPECIES_BEGIN = the index at which the first species is stored in the 
                  unk data structure
  SPECIES_BEGIN = NPROP_VARS + 1
  **************************************************************************
#endif

#define PROP_VARS_BEGIN CONSTANT_ONE
#define UNK_VARS_BEGIN (PROP_VARS_BEGIN)

%(variableDefines)s
#define NPROP_VARS %(nvars)s
#define PROP_VARS_END (PROP_VARS_BEGIN + NPROP_VARS - CONSTANT_ONE)

#define SPECIES_BEGIN (PROP_VARS_END + CONSTANT_ONE)
%(speciesDefines)s

%(mfrac_specDefine)s
#define NSPECIES %(nspecies)s
#define SPECIES_END (SPECIES_BEGIN + NSPECIES - CONSTANT_ONE)

#define MASS_SCALARS_BEGIN (SPECIES_END + CONSTANT_ONE)

%(massscalarsDefines)s
#define NMASS_SCALARS %(nmass_scalars)s
#define MASS_SCALARS_END (MASS_SCALARS_BEGIN+NMASS_SCALARS-CONSTANT_ONE)

#define NMASS_SCALAR_GROUPS %(nmassscalar_groups)s 

#define NUNK_VARS (NPROP_VARS + NSPECIES + NMASS_SCALARS)
#define UNK_VARS_END (UNK_VARS_BEGIN - CONSTANT_ONE + NUNK_VARS)

#if 0
  Macros for non-replicated variable arrays.
   - meshes are indexed beginning with 0.
   - local and global variables are indexed beginning at 1.
  
  NONREP_NLOCS -- gives the number of variables in the local subset of a mesh
  NONREP_LOC2GLOB -- gets the global index from a local index
  NONREP_GLOB2LOC -- gets the local index from a global index
  NONREP_MESHOFGLOB -- determines the mesh that owns a global index
#endif

#define NONREP_NLOCS(mesh,meshes,globs) (((globs)-(mesh)+(meshes)-1)/(meshes))
#define NONREP_LOC2GLOB(loc,mesh,meshes) (((loc)-1)*(meshes)+(mesh)+1)
#define NONREP_GLOB2LOC(glob,mesh,meshes) (((glob)-1)/(meshes)+1)
#define NONREP_MESHOFGLOB(glob,meshes) (mod((glob)-1,meshes))

#if 0
  For each non-replicated unk X, the following are defined:
    X_NONREP -- the integer id of this nonrep (starts at 1)
    X_NONREP_LOC2UNK(n) -- macro to convert a local index to its corresponding index in unk
    X_NONREP_MAXLOCS -- the maximum number of locals any processor can have (taken directly from CONFIG line)
    X_NONREP_RPCOUNT -- character literal of the runtime parameter indicating the number of variables in the nonrep array.
    
  The following are defined for code that has to deal with all nonrep variables in a general way (mostly IO).
    NONREP_COUNT -- the number of nonrep variable arrays in this simulation
    NONREP_NAMEF_FLAT_LWR -- a concatenation of all name formatters in lowercase for all nonreps
    NONREP_NAMEF_START -- integer array literal of length NONREP_COUNT+1 that stores for each nonrep id, where its name formatter begins in NONREP_NAMEF_FLAT_LWR.  The ending index of each name can be retrieved with NONREP_NAMEF_START(nonrep+1)-1
    NONREP_LOCUNK1 -- integer array literal of length NONREP_COUNT that foreach nonrep gives its starting unk index
    NONREP_MAXLOCS -- integer array literal of length NONREP_COUNT that gives the maximum number of locals foreach nonrep
    NONREP_UNKS_START -- integer array literal of length NONREP_COUNT+1 that gives that starting index in NONREP_UNKS_FLAT for a given nonrep.  The ending index can be computed as NONREP_UNKS_START(nonrep+1)-1
    NONREP_RPCOUNT_FLAT -- concatenation of each nonreps runtime parameter string designating the size of the variable array
    NONREP_RPCOUNT_START -- integer array literal that for each nonrep gives the starting position of its runtime parameter name in NONREP_RPCOUNT_FLAT.  The ending position is NONREP_RPCOUNT_START(nonrep+1)-1    
    * If you didnt pick up on this, conceptually anything named ***_FLAT and ***_START together form a jagged array of arrays
#endif
  
%(nonrepDefines)s

#if 0
  ************************************************************************
  The #define definitions below indicate the index at which a grid scope 
  FACE variable is located in the facevars data structure.  
  ************************************************************************
#endif

%(facevarDefines)s    
#define NFACE_VARS %(nfacevars)s

#if 0
  ************************************************************************
  The #define definitions below indicate the index at which a flux 
  variable is located in the flux data structure (data structure can vary 
  depending on grid).
  Property fluxes are denoted with '_FLUX' following the variable name.
  Flux Species are also stored in the flux datastructure.  

  NPROP_FLUX = the total number of grid property variables in the 
                 simulation (ie those ending in _FLUX)
  NSPECIES_FLUX = the number of flux species in the simulation
  NMASS_SCALARS_FLUX = flux mass scalars

  NFLUXES = NPROP_FLUX + NSPECIES_FLUX + NMASS_SCALARS_FLUX
  ************************************************************************
#endif

%(fluxDefines)s
#define NPROP_FLUX %(nflux)s
#define NSPECIES_FLUX %(nspecies)s
#define NMASS_SCALARS_FLUX %(nmass_scalars)s
#define NFLUXES (NPROP_FLUX + NSPECIES_FLUX + NMASS_SCALARS_FLUX)

#define PROP_FLUX_BEGIN CONSTANT_ONE
#define FLUXES_BEGIN (PROP_FLUX_BEGIN)

#define PROP_FLUX_END (PROP_FLUX_BEGIN + NPROP_FLUX - CONSTANT_ONE)
#define SPECIES_FLUX_BEGIN (PROP_FLUX_END + CONSTANT_ONE)
#define SPECIES_FLUX_END (SPECIES_FLUX_BEGIN + NSPECIES_FLUX - CONSTANT_ONE)
#define MASS_SCALARS_FLUX_BEGIN (SPECIES_FLUX_END + CONSTANT_ONE)
#define MASS_SCALARS_FLUX_END (MASS_SCALARS_FLUX_BEGIN + NMASS_SCALARS_FLUX - CONSTANT_ONE)
#define NVARS_TOTAL (NUNK_VARS + (NFACE_VARS * 3) + NFLUXES)

#if 0
  ************************************************************************
  The #define definitions below pertain to particles.
  The particles data structure  needs to be defined at compile time and thus
  the number of particle properties (such as, velocity, position etc)
  are recorded here.  For example, #define VELX_PART_PROP 2 indicates that
  x velocity of a particle is stored in the second position in the particles
  array.

  Since FLASH3 all particle properties are reals.
  If particles are not being used in a simulation then the number of
  particle properties = 0 (defined as NPART_PROPS)
  DEV: CD.  It is set to 1 and not 0 by the setup script.  I guess this is 
  to avoid zero sized arrays???
  ************************************************************************
#endif

%(partpropDefines)s

#define PART_PROPS_BEGIN CONSTANT_ONE
#define NPART_PROPS %(nrealprops)s
#define PART_PROPS_END (PART_PROPS_BEGIN + NPART_PROPS - CONSTANT_ONE)

#if 0
  ************************************************************************
  The #define definitions below pertain to particles types.
  Particles can be of type passive, active, star, etc.....

  PASSIVE, if it exists, should always be the first type
  ************************************************************************
#endif

%(particletypeDefines)s

#define PART_TYPES_BEGIN CONSTANT_ONE
#define NPART_TYPES %(nparticletypes)s
#define PART_TYPES_END (PART_TYPES_BEGIN + NPART_TYPES - CONSTANT_ONE)
    
#if 0
  ************************************************************************
  NDIM = number of dimensions in the simulation (or 2 if user did not 
  specify anything)
  MAXBLOCKS = maximum number of blocks allowed per processor

  ************************************************************************
#endif

#define NDIM %(ndim|2)s
#define MAXBLOCKS %(maxblocks|-1)s


#if 0
  ************************************************************************
  Definitions for scratch array.
  ************************************************************************
#endif
%(scratchvarDefines)s
#define NSCRATCH_GRID_VARS %(nscratchvars)s
#define SCRATCH_GRID_VARS_BEGIN (CONSTANT_ONE)
#define SCRATCH_GRID_VARS_END (SCRATCH_GRID_VARS_BEGIN - CONSTANT_ONE + NSCRATCH_GRID_VARS)


#if 0
  ************************************************************************
  Definitions for scratch_ctr array.
  ************************************************************************
#endif
%(scratchcentervarDefines)s
#define NSCRATCH_CENTER_VARS %(nscratchcentervars)s
#define SCRATCH_CENTER_VARS_BEGIN (CONSTANT_ONE)
#define SCRATCH_CENTER_VARS_END (SCRATCH_CENTER_VARS_BEGIN - CONSTANT_ONE + NSCRATCH_CENTER_VARS)


#if 0
  ************************************************************************
  Definitions for scratch_facevarx array.
  ************************************************************************
#endif
%(scratchfacexvarDefines)s
#define NSCRATCH_FACEX_VARS %(nscratchfacexvars)s
#define SCRATCH_FACEX_VARS_BEGIN (CONSTANT_ONE)
#define SCRATCH_FACEX_VARS_END (SCRATCH_FACEX_VARS_BEGIN - CONSTANT_ONE + NSCRATCH_FACEX_VARS)


#if 0
  ************************************************************************
  Definitions for scratch_facevary array.
  ************************************************************************
#endif
%(scratchfaceyvarDefines)s
#define NSCRATCH_FACEY_VARS %(nscratchfaceyvars)s
#define SCRATCH_FACEY_VARS_BEGIN (CONSTANT_ONE)
#define SCRATCH_FACEY_VARS_END (SCRATCH_FACEY_VARS_BEGIN - CONSTANT_ONE + NSCRATCH_FACEY_VARS)


#if 0
  ************************************************************************
  Definitions for scratch_facevarz array.
  ************************************************************************
#endif
%(scratchfacezvarDefines)s
#define NSCRATCH_FACEZ_VARS %(nscratchfacezvars)s
#define SCRATCH_FACEZ_VARS_BEGIN (CONSTANT_ONE)
#define SCRATCH_FACEZ_VARS_END (SCRATCH_FACEZ_VARS_BEGIN - CONSTANT_ONE + NSCRATCH_FACEZ_VARS)


#define MAX_PLOT_VARS %(max_plot_vars)s

#if 0
  ************************************************************************
  The #define definitions below pertain simulations where a fixed block 
  size. mesh is used. 'Fixed block size' means that the size of the block 
  is specified at compile time rather than runtime.  It also means that 
  all the blocks in a simulation have the same number of zones per block 
  and the same number of guardcells.

  Paramesh and a fixed block size uniform grid run in this mode. 

  NXB (NYB, NZB) = the number of zones in the x(y,z) direction respectively
  NGUARD = the number of guardcells in the block.  In fixed block size 
           mode the x,y and z dimensions of the block all must have the 
           same number of guardcells.

  GRID_ILO_GC = the index of the first zone in block in the i dimension 
                including guardcells.  (Always = 1)
  GRID_IHI_GC = the index of the last zone in a block in the i dimension 
                including guardcells.
  If, NGUARD = 4 and NXB = 8 then GRID_IHI_GC = 16 (fortran counting)  
  GRID_ILO    = the index of the first interior (non guardcell) zone in 
                a block in the i dimension.
  If, NGUARD = 4 then GRID_ILO = 5  (fortran counting)
  GRID_IHI    = the index of the last interior zone in a block in the i 
                dimension.              
  If, NGUARD = 4 and NXB = 8 then GRID_IHI = 12 (fortran counting)
  ************************************************************************
#endif


#define K1D CONSTANT_ONE

#if NDIM > CONSTANT_ONE
#define K2D CONSTANT_ONE
#else
#define K2D CONSTANT_ZERO
#endif

#if NDIM > CONSTANT_TWO
#define K3D CONSTANT_ONE
#else
#define K3D CONSTANT_ZERO
#endif

#define NGUARD %(nguard)s

#if 0
  *************
  Grid geometry
  *************
#endif
#define %(gridGeometry)s
#define GRID_CURVILINEAR %(curvilinear)s

#if GRID_CURVILINEAR == CONSTANT_ZERO
#undef GRID_CURVILINEAR
#endif

#if 0
  ************************************************************************
  FIXEDBLOCKSIZE is initially defined to be 0 or 1
  Soon FIXEDBLOCKSIZE is undefined if it was set to 0
  ************************************************************************
#endif

#define FIXEDBLOCKSIZE %(fixedblocksize)s

#if FIXEDBLOCKSIZE == CONSTANT_ZERO
#undef FIXEDBLOCKSIZE
#endif

#ifdef FIXEDBLOCKSIZE 

#define NXB %(nxb|-1)s
#define NYB %(nyb|-1)s
#define NZB %(nzb|-1)s

#define GRID_ILO_GC CONSTANT_ONE
#define GRID_JLO_GC CONSTANT_ONE
#define GRID_KLO_GC CONSTANT_ONE
        
#define GRID_IHI_GC (NXB + CONSTANT_TWO*NGUARD*K1D)
#define GRID_JHI_GC (NYB + CONSTANT_TWO*NGUARD*K2D)
#define GRID_KHI_GC (NZB + CONSTANT_TWO*NGUARD*K3D)
        
#define GRID_ILO (NGUARD*K1D + CONSTANT_ONE)
#define GRID_JLO (NGUARD*K2D + CONSTANT_ONE)
#define GRID_KLO (NGUARD*K3D + CONSTANT_ONE)

#define GRID_IHI (NGUARD*K1D + NXB)
#define GRID_JHI (NGUARD*K2D + NYB)
#define GRID_KHI (NGUARD*K3D + NZB)


#if 0
  ************************************************************************
  !!!DEV: MAXCELLS pertains to all simulations but someone else needs to 
     decide how this works with non fixed block size!

  It is the maxmimum of  GRID_{IJK}HI_GC. We compute the maximum using 
  preprocessor stuff
  ************************************************************************
#endif

#if GRID_IHI_GC > GRID_JHI_GC
#define FLASHPP_MAX_IJ_TEMP GRID_IHI_GC
#else
#define FLASHPP_MAX_IJ_TEMP GRID_JHI_GC
#endif

#if FLASHPP_MAX_IJ_TEMP > GRID_KHI_GC
#define MAXCELLS FLASHPP_MAX_IJ_TEMP
#else
#define MAXCELLS GRID_KHI_GC
#endif

#if 0
  ************************************************************************
  End of #ifdef FIXEDBLOCKSIZE follows
  ************************************************************************
#endif
#endif 

#if 0
  ************************************************************************
  FL_NON_PERMANENT_GUARDCELLS is initially defined to be 0 or 1
  Soon FL_NON_PERMANENT_GUARDCELLS is undefined if it was set to 0
  ************************************************************************
#endif

#define FL_NON_PERMANENT_GUARDCELLS %(npg)s

#if FL_NON_PERMANENT_GUARDCELLS == CONSTANT_ZERO
#undef FL_NON_PERMANENT_GUARDCELLS
#endif

#if 0
  ***************************************************************************
  The following set of #defines are specific to FLASH and picked up from 
  Config Files and default values in "globals.py". Some uses:
 
  * Allow code to determine which units have been compiled into the code. 
    In each Config file the unit can declare a Pre-processor symbol using 

    PPDEFINE fubar
    This gets translated to 
    #define fubar in Flash.h

  * Some configuration of FLASH specific to a unit or simulation. For example,

    PPDEFINE FLASH_MHD_DIM 1
    becomes
    #define FLASH_MHD_DIM 1 in Flash.h
    
  ***************************************************************************
#endif

%(ppDefines!\n)s


#if 0
  ************************************************************************
  If we are in reorder mode, this provides a clean way to transpose arrays
  if assigning slices.
  ************************************************************************
#endif

#ifdef INDEXREORDER
#define TRANSPOSE_IF_REORDER(arr) transpose(arr)
#else
#define TRANSPOSE_IF_REORDER(arr) arr
#endif

#define DRIFT_ENABLE 0

#if DRIFT_ENABLE
#define Grid_releaseBlkPtr Driver_driftSetSrcLoc(__FILE__,__LINE__); call Grid_releaseBlkPtr
#endif

#if 0
  ************************************************************************
  End of #ifdef ensuring that this header file is included only once
  ************************************************************************
#endif

#endif
