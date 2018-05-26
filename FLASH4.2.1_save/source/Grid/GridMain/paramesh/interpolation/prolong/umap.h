! Any source file including this header file should include "Flash.h" before it. - KW

#ifdef FIXEDBLOCKSIZE
       PARAMETER ( mvx_m   = MAXCELLS)
#else
! This file should not be included by any source files when not using fixed block sizes,
! since this file should be only included when using a PARAMESH Grid implementation! - KW
#endif

       PARAMETER ( mui_m   = NUNK_VARS)

!       PARAMETER ( mvx_m   = max(nxb,nyb,nzb) + 2*nguard )
!       PARAMETER ( mui_m   = nvar )

! The following should work in most cases.

!       PARAMETER ( mvx_m   = 128 + 2 * 4 )
!       PARAMETER ( mui_m   = 35 )

!       PARAMETER ( mvm_m   = mvx_m ) ! not used for dimensioning anywhere - KW
        PARAMETER ( mvxu_m  = mvx_m * mui_m )
        PARAMETER ( mvxmu_m = mvx_m * mvxu_m )
