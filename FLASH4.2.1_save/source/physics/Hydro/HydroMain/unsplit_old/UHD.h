!!***h* source/physics/Hydro/HydroMain/unsplit/UHD.h
!!
!! This is the internal header file for the Unsplit MHD/Hydro units.
!!
!!***

#include "Flash.h"
#include "constants.h"

#define DIR_X 1
#define DIR_Y 2
#define DIR_Z 3


#define MINMOD  1
#define MC      2
#define HYBRID  3
#define VANLEER 4
#define LIMITED 5

#define HARTEN      1
#define HARTENHYMAN 2

#define UPDATE_INTERIOR 1
#define UPDATE_BOUND    2
#define UPDATE_ALL      3

#define FWDCONVERT 1
#define BWDCONVERT 2


!!--------------------------------------------------------------!!
!! [1] FOR UNSPLIT MHD IMPLEMENTATION --------------------------!!
!!--------------------------------------------------------------!!
#ifdef FLASH_USM_MHD

#define INJECTION_PROL 1
#define BALSARA_PROL   2

#define ROE   1
#define LLF   2
#define HLL   3
#define HLLC  4
#define HLLD  5
#define MARQ  6
#define HYBR  7

#ifdef CURX_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

#ifdef CURY_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

#ifdef CURZ_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

#ifdef DIVV_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

!! PRIMITIVE VARIABLES FOR MHD
#define HY_DENS 1
#define HY_VELX 2
#define HY_VELY 3
#define HY_VELZ 4
#define HY_PRES 5
#define HY_MAGX 6
#define HY_MAGY 7
#define HY_MAGZ 8
#define HY_GAMC 9
#define HY_GAME 10
#define HY_EINT 11

! MHD 3T 
#ifdef FLASH_UHD_3T
#define HY_EELE 12
#define HY_EION 13
#define HY_ERAD 14
#define HY_GRAV 15
#define HY_TEMP 16
#define HY_ABAR 17
#define HY_ZBAR 18
#else
#define HY_GRAV 12
#define HY_TEMP 13
#define HY_ABAR 14
#define HY_ZBAR 15
#endif


!! CONSERVATIVE VARIABLES FOR MHD 
#define HY_XMOM 2
#define HY_YMOM 3
#define HY_ZMOM 4
#define HY_ENER 5
!! NOTE THAT THE REST OF CONSERVATIVE VARIABLES 
!! ARE ALREADY DEFINED IN PRIMITIVE VARIABLES

!! FLUX VARIABLES FOR MHD
#define HY_DENS_FLUX 1
#define HY_XMOM_FLUX 2
#define HY_YMOM_FLUX 3
#define HY_ZMOM_FLUX 4
#define HY_ENER_FLUX 5
#define HY_MAGX_FLUX 6
#define HY_MAGY_FLUX 7
#define HY_MAGZ_FLUX 8
#define HY_EINT_FLUX 9
#define HY_PRES_FLUX 10

! MHD 3T
#ifdef FLASH_UHD_3T
#define HY_EELE_FLUX 11
#define HY_EION_FLUX 12
#define HY_ERAD_FLUX 13
#endif

!! WAVE STRUCTURES FOR MHD
#define HY_FASTLEFT 1
#define HY_ALFNLEFT 2
#define HY_SLOWLEFT 3
#define HY_ENTROPY  4
#define HY_SLOWRGHT 5
#define HY_ALFNRGHT 6
#define HY_FASTRGHT 7

!! EXTRA PARAMETERS FOR PURE MHD
#define HY_VARINUM HY_MAGZ

#ifdef FLASH_UHD_3T
#define HY_SCRATCH_NUM 14
#else
#define HY_SCRATCH_NUM 11
#endif



#else
!!--------------------------------------------------------------!!
!! [2] FOR UNSPLIT HYDRO IMPLEMENTATION ------------------------!!
!!--------------------------------------------------------------!!
#define ROE   1
#define LLF   2
#define HLL   3
#define HLLC  4
#define MARQ  5
#define MMAR  6
#define HYBR  7

!! PRIMITIVE VARIABLES FOR PURE HYDRO
#define HY_DENS 1
#define HY_VELX 2
#define HY_VELY 3
#define HY_VELZ 4
#define HY_PRES 5
#define HY_GAMC 6
#define HY_GAME 7
#define HY_EINT 8

#ifdef FLASH_UHD_3T
#define HY_EELE 9
#define HY_EION 10
#define HY_ERAD 11
#define HY_GRAV 12
#define HY_TEMP 13
#define HY_ABAR 14
#define HY_ZBAR 15
#else
#define HY_GRAV 9
#define HY_TEMP 10
#define HY_ABAR 11
#define HY_ZBAR 12
#endif

!! CONSERVATIVE VARIABLES FOR PURE HYDRO
#define HY_XMOM 2
#define HY_YMOM 3
#define HY_ZMOM 4
#define HY_ENER 5

!! FLUX VARIABLES FOR PURE HYDRO
#define HY_DENS_FLUX 1
#define HY_XMOM_FLUX 2
#define HY_YMOM_FLUX 3
#define HY_ZMOM_FLUX 4
#define HY_ENER_FLUX 5
#define HY_EINT_FLUX 6
#define HY_PRES_FLUX 7
#ifdef FLASH_UHD_3T
#define HY_EELE_FLUX 8
#define HY_EION_FLUX 9
#define HY_ERAD_FLUX 10
#endif

!! WAVE STRUCTURES FOR PURE HYDRO
#define HY_FASTLEFT 1
#define HY_SLOWLEFT 2
#define HY_ENTROPY  3
#define HY_SLOWRGHT 4
#define HY_FASTRGHT 5

!! EXTRA PARAMETERS FOR PURE HYDRO
#define HY_VARINUM HY_PRES

#ifdef FLASH_UHD_3T
#define HY_SCRATCH_NUM 11
#else
#define HY_SCRATCH_NUM 8
#endif

!! NOTE: the following scratch indices are the definitions of the LOCAL
!!       scratch arrays defined in Hydro_data.F90, NOT from Config file.
!! The following definitions are only applied to hydro case
!! because the USM solver SHOULD use global (not local as in the below) 
!! scratch arrays defined via Config file

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#ifndef FLASH_UHD_3T
!without 3T

#define XP01_SCRATCH_CENTER_VAR 1
#define XP02_SCRATCH_CENTER_VAR 2
#define XP03_SCRATCH_CENTER_VAR 3
#define XP04_SCRATCH_CENTER_VAR 4
#define XP05_SCRATCH_CENTER_VAR 5
#define XP06_SCRATCH_CENTER_VAR 6
#define XP07_SCRATCH_CENTER_VAR 7
#define XP08_SCRATCH_CENTER_VAR 8

#define XN01_SCRATCH_CENTER_VAR 9
#define XN02_SCRATCH_CENTER_VAR 10
#define XN03_SCRATCH_CENTER_VAR 11
#define XN04_SCRATCH_CENTER_VAR 12
#define XN05_SCRATCH_CENTER_VAR 13
#define XN06_SCRATCH_CENTER_VAR 14
#define XN07_SCRATCH_CENTER_VAR 15
#define XN08_SCRATCH_CENTER_VAR 16


#if NDIM >= 2
#define YP01_SCRATCH_CENTER_VAR 17
#define YP02_SCRATCH_CENTER_VAR 18
#define YP03_SCRATCH_CENTER_VAR 19
#define YP04_SCRATCH_CENTER_VAR 20
#define YP05_SCRATCH_CENTER_VAR 21
#define YP06_SCRATCH_CENTER_VAR 22
#define YP07_SCRATCH_CENTER_VAR 23
#define YP08_SCRATCH_CENTER_VAR 24

#define YN01_SCRATCH_CENTER_VAR 25
#define YN02_SCRATCH_CENTER_VAR 26
#define YN03_SCRATCH_CENTER_VAR 27
#define YN04_SCRATCH_CENTER_VAR 28
#define YN05_SCRATCH_CENTER_VAR 29
#define YN06_SCRATCH_CENTER_VAR 30
#define YN07_SCRATCH_CENTER_VAR 31
#define YN08_SCRATCH_CENTER_VAR 32


#if NDIM == 3
#define ZP01_SCRATCH_CENTER_VAR 33
#define ZP02_SCRATCH_CENTER_VAR 34
#define ZP03_SCRATCH_CENTER_VAR 35
#define ZP04_SCRATCH_CENTER_VAR 36
#define ZP05_SCRATCH_CENTER_VAR 37
#define ZP06_SCRATCH_CENTER_VAR 38
#define ZP07_SCRATCH_CENTER_VAR 39
#define ZP08_SCRATCH_CENTER_VAR 40

#define ZN01_SCRATCH_CENTER_VAR 41
#define ZN02_SCRATCH_CENTER_VAR 42
#define ZN03_SCRATCH_CENTER_VAR 43
#define ZN04_SCRATCH_CENTER_VAR 44
#define ZN05_SCRATCH_CENTER_VAR 45
#define ZN06_SCRATCH_CENTER_VAR 46
#define ZN07_SCRATCH_CENTER_VAR 47
#define ZN08_SCRATCH_CENTER_VAR 48

#endif
!end #if NDIM == 3
#endif
!end #if NDIM >= 2

#else
! with 3T
#define XP01_SCRATCH_CENTER_VAR 1
#define XP02_SCRATCH_CENTER_VAR 2
#define XP03_SCRATCH_CENTER_VAR 3
#define XP04_SCRATCH_CENTER_VAR 4
#define XP05_SCRATCH_CENTER_VAR 5
#define XP06_SCRATCH_CENTER_VAR 6
#define XP07_SCRATCH_CENTER_VAR 7
#define XP08_SCRATCH_CENTER_VAR 8
#define XP09_SCRATCH_CENTER_VAR 9
#define XP10_SCRATCH_CENTER_VAR 10
#define XP11_SCRATCH_CENTER_VAR 11

#define XN01_SCRATCH_CENTER_VAR 12
#define XN02_SCRATCH_CENTER_VAR 13
#define XN03_SCRATCH_CENTER_VAR 14
#define XN04_SCRATCH_CENTER_VAR 15
#define XN05_SCRATCH_CENTER_VAR 16
#define XN06_SCRATCH_CENTER_VAR 17
#define XN07_SCRATCH_CENTER_VAR 18
#define XN08_SCRATCH_CENTER_VAR 19
#define XN09_SCRATCH_CENTER_VAR 20
#define XN10_SCRATCH_CENTER_VAR 21
#define XN11_SCRATCH_CENTER_VAR 22


#if NDIM >= 2
#define YP01_SCRATCH_CENTER_VAR 23
#define YP02_SCRATCH_CENTER_VAR 24
#define YP03_SCRATCH_CENTER_VAR 25
#define YP04_SCRATCH_CENTER_VAR 26
#define YP05_SCRATCH_CENTER_VAR 27
#define YP06_SCRATCH_CENTER_VAR 28
#define YP07_SCRATCH_CENTER_VAR 29
#define YP08_SCRATCH_CENTER_VAR 30
#define YP09_SCRATCH_CENTER_VAR 31
#define YP10_SCRATCH_CENTER_VAR 32
#define YP11_SCRATCH_CENTER_VAR 33

#define YN01_SCRATCH_CENTER_VAR 34
#define YN02_SCRATCH_CENTER_VAR 35
#define YN03_SCRATCH_CENTER_VAR 36
#define YN04_SCRATCH_CENTER_VAR 37
#define YN05_SCRATCH_CENTER_VAR 38
#define YN06_SCRATCH_CENTER_VAR 39
#define YN07_SCRATCH_CENTER_VAR 40
#define YN08_SCRATCH_CENTER_VAR 41
#define YN09_SCRATCH_CENTER_VAR 42
#define YN10_SCRATCH_CENTER_VAR 43
#define YN11_SCRATCH_CENTER_VAR 44


#if NDIM == 3
#define ZP01_SCRATCH_CENTER_VAR 45
#define ZP02_SCRATCH_CENTER_VAR 46
#define ZP03_SCRATCH_CENTER_VAR 47
#define ZP04_SCRATCH_CENTER_VAR 48
#define ZP05_SCRATCH_CENTER_VAR 49
#define ZP06_SCRATCH_CENTER_VAR 50
#define ZP07_SCRATCH_CENTER_VAR 51
#define ZP08_SCRATCH_CENTER_VAR 52
#define ZP09_SCRATCH_CENTER_VAR 53
#define ZP10_SCRATCH_CENTER_VAR 54
#define ZP11_SCRATCH_CENTER_VAR 55


#define ZN01_SCRATCH_CENTER_VAR 56
#define ZN02_SCRATCH_CENTER_VAR 57
#define ZN03_SCRATCH_CENTER_VAR 58
#define ZN04_SCRATCH_CENTER_VAR 59
#define ZN05_SCRATCH_CENTER_VAR 60
#define ZN06_SCRATCH_CENTER_VAR 61
#define ZN07_SCRATCH_CENTER_VAR 62
#define ZN09_SCRATCH_CENTER_VAR 63
#define ZN10_SCRATCH_CENTER_VAR 64
#define ZN11_SCRATCH_CENTER_VAR 65
#define ZN12_SCRATCH_CENTER_VAR 66

#endif
!end #if NDIM == 3
#endif
!end #if NDIM >= 2

#endif
!end if ifndef FLASH_UHD_3T

#endif
!end if #ifndef FLASH_UHD_NEED_SCRATCHVARS

!! Memory Efficient Setup without using global SCRATCH arrays (SCRATCH_CTR)
#ifndef FLASH_UHD_NEED_SCRATCHVARS
#define HY_NSCRATCH_VARS (2*HY_SCRATCH_NUM*NDIM)
#endif


#endif
!!--------------------------------------------------------------!!
!! END OF MHD AND HYDRO DEFINITIONS ----------------------------!!
!!--------------------------------------------------------------!!

#define HY_NSPEC NSPECIES+NMASS_SCALARS


!! DEFINE TOTAL NUMBERS OF VARIABLES AND WAVES
#define HY_WAVENUM  HY_FASTRGHT
#define HY_VARINUM2 HY_VARINUM+2
#define HY_VARINUM3 HY_VARINUM+3
#define HY_VARINUM4 HY_VARINUM+4
#define HY_VARINUM7 HY_VARINUM+7

#if defined(FLASH_UHD_NEED_SCRATCHVARS)||defined(FLASH_USM_MHD)||(HY_NSPEC<=0)
#ifndef FLASH_UHD_3T
#define HY_VARINUMMAX HY_VARINUM4
#else
#define HY_VARINUMMAX HY_VARINUM7
#endif
#else
#ifndef FLASH_UHD_3T
#define HY_VARINUMMAX HY_VARINUM4+HY_NSPEC
#else
#define HY_VARINUMMAX HY_VARINUM7+HY_NSPEC
#endif
#endif

#ifndef FLASH_UHD_3T
#ifndef GRAVITY
#define HY_END_VARS HY_EINT
#else
#define HY_END_VARS HY_GRAV
#endif
#define HY_END_FLUX HY_PRES_FLUX
#else
#define HY_END_VARS HY_ERAD
#define HY_END_FLUX HY_ERAD_FLUX
#endif
