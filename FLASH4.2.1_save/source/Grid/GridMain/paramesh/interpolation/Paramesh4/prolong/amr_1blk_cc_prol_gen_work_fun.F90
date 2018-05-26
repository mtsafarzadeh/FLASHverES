!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_gen_work_fun(recvt, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)



!
!------------------------------------------------------------------------
!
! This routine is a wrapper routine which calls the functions
! which prolong data for WORK. The argument idest selects the layer
! within WORK on which prolongation is required.
! The argument interp should select the actual interpolation function to be
! called, but is ignored here.
!
!
! Original:     Peter MacNeice          January 2002
! Modified:     Klaus Weide             September 2006
!------------------------------------------------------------------------

      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      integer, intent(in) :: interp
      real,    intent(inout) :: recvt(:,:,:)

!------------------------------------

! Call the minimally changed subroutine for Paramesh2
      call amr_prolong_gen_work1_fun &
     &     (recvt,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
     &     mype,lb)

      return
      end subroutine amr_1blk_cc_prol_gen_work_fun
