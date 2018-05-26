!!****ih* source/Grid/GridSolvers/HYPRE/gr_hypreLocalInterface
!!
!! This is a module file for the HYPRE solver in FLASH that defines
!! additional interfaces private to the GridSolvers/HYPRE implementations.
!!
!! NOTES
!!
!!  Adding explicit interfaces here is being tried as an alternative to
!!  writiting executable FORTRAN wrappers for each additional HYPRE routine
!!  we want to call from FLASH.
!!
!! SEE ALSO
!!
!!  gr_hypreF90CAdapters.c
!!***

#include "constants.h"

Module gr_hypreLocalInterface

#if 0
  ! Maybe one day... for now, not all compilers support this.
  interface
     subroutine hypre_sstructinnerprod(fx, fy, fresult, ierr) &
          bind(C,name="HYPRE_SStructInnerProd")
       implicit none
       integer*8,intent(IN)  :: fx, fy
       real,     intent(OUT) :: fresult
       integer,  intent(OUT) :: ierr
    end subroutine hypre_sstructinnerprod

     integer function hypre_pcggetconverged(solver, converged) &
          bind(C,name="HYPRE_PCGGetConverged")
       implicit none
       integer*8,VALUE :: solver ! intent(IN)
       integer,intent(OUT) :: converged
     end function hypre_pcggetconverged

     subroutine hypre_describeerror(hypreErr, description) &
          bind(C,name="HYPRE_DescribeError")
       implicit none
       integer,VALUE :: hypreErr ! intent(IN)
       character(len=1),dimension(MAX_STRING_LENGTH),intent(OUT) :: description
     end subroutine hypre_describeerror

  end interface

#else
  ! fallback... wrapper implementations are in gr_hypreF90CAdapters.c

  interface
     subroutine hypre_sstructinnerprod(fx, fy, fresult, ierr)
       implicit none
       integer*8,intent(IN)  :: fx, fy
       real,     intent(OUT) :: fresult
       integer,  intent(OUT) :: ierr
    end subroutine hypre_sstructinnerprod

     subroutine hypre_pcggetconverged(solver, converged, ierr)
       implicit none
       integer*8,intent(IN) :: solver
       integer,intent(OUT) :: converged
       integer,  intent(OUT) :: ierr
     end subroutine hypre_pcggetconverged

     subroutine hypre_describeerror(hypreErr, description)
       implicit none
       integer,intent(IN) :: hypreErr
       character(len=MAX_STRING_LENGTH),intent(OUT) :: description
     end subroutine hypre_describeerror

  end interface

#endif

end Module gr_hypreLocalInterface


