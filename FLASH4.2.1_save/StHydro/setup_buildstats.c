/* !!! DO NOT EDIT, FILES WRITTEN BY SETUP SCRIPT !!!
   
!!****f* object/setup_buildstats
!!
!! NAME
!!
!!  setup_buildstats
!!
!!
!! SYNOPSIS
!!
!!  call setup_buildstats(build_date, build_dir, build_machine, setup_call)
!!
!!  call setup_buildstats(character, character, character, character)
!!
!! 
!! DESCRIPTION
!!
!!  Simple subroutine generated at build time that returns the build date
!!  build directory, build machine, c and f flags and the full setup 
!!  command used to assemble the FLASH executable.
!!
!!
!!
!!***
*/

#include "Flash.h"
#include "constants.h"
#include "mangle_names.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void FTOC(setup_buildstats)(char* build_date, 
		    char* build_dir, 
		    char* build_machine, 
		    char* setup_call, 
		    char* c_flags, 
		    char* f_flags){



     strncpy(build_date, "Fri Dec 18 18:17:48 CST 2015",80);
     strncpy(build_dir, "/work/01734/evan1022/FLASH4.2.1_save/StHydro", 80);
     strncpy(build_machine, "Linux login1.stampede.tacc.utexas.edu 2.6.32-431.17.1.el6.x86_64 #1 SMP Wed May ", 80);
     strncpy(setup_call, "/work/01734/evan1022/FLASH4.2.1_save/bin/setup.py -auto -3d -nxb=32 -nyb=32 -nzb=32 +ug +uhd StirTurbHydro -objdir=StHydro ",400);
     strncpy(c_flags, "mpicc -I/opt/apps/intel15/hdf5/1.8.16/x86_64/include -DH5_USE_16_API -I /work/01734/evan1022/hypre-2.9.0b/src/install/include -I/opt/apps/intel15/mvapich2/2.1/include -c -O3 -g -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 -diag-disable 10120 -DMAXBLOCKS=1 -DNXB=32 -DNYB=32 -DNZB=32 -DN_DIM=3", 400);
     strncpy(f_flags, "mpif90 -c -g -r8 -i4 -O3 -real_size 64 -diag-disable 10120 -I /opt/apps/intel15/hdf5/1.8.16/x86_64/include -DH5_USE_16_API -DMAXBLOCKS=1 -DNXB=32 -DNYB=32 -DNZB=32 -DN_DIM=3", 400);


}

