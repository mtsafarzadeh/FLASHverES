!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_restrict_nc_user()




!------------------------------------------------------------------------
!
! This routine is a stub routine and is meant to be replaced by a
! user written routine if the default is not adequate
!
!------------------------------------------------------------------------

      implicit none

      return
      end subroutine amr_restrict_nc_user
