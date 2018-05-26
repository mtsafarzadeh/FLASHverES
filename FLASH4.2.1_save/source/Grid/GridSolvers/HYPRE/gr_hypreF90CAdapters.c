
/******************************************************************************
 *
 * HYPRE_PCG Fortran interface -
 *           KW based on hypre-2.9.0b/src/parcsr_ls/F90_HYPRE_parcsr_pcg.c
 *
 *****************************************************************************/

#include "_hypre_parcsr_ls.h"
#include "fortran.h"


/*--------------------------------------------------------------------------
 * HYPRE_PCGSetConvergenceFactorTol
 *--------------------------------------------------------------------------*/

void
hypre_F90_IFACE(hypre_pcgsetconvergencefactortol, HYPRE_PCGSETCONVERGENCEFACTORTOL)
   ( hypre_F90_Obj *solver,
     hypre_F90_Dbl *cfTol,
     hypre_F90_Int *ierr    )
{
   *ierr = (hypre_F90_Int)
      ( HYPRE_PCGSetConvergenceFactorTol(
           hypre_F90_PassObj (HYPRE_Solver, solver),
           hypre_F90_PassDbl (cfTol) ) );
}

/*--------------------------------------------------------------------------
 * HYPRE_PCGSetRecomputeResidual
 *--------------------------------------------------------------------------*/

void
hypre_F90_IFACE(hypre_pcgsetrecomputeresidual, HYPRE_PCGSETRECOMPUTERESIDUAL)
   ( hypre_F90_Obj *solver,
     hypre_F90_Int *recomputeResidual,
     hypre_F90_Int *ierr    )
{
   *ierr = (hypre_F90_Int)
      ( HYPRE_PCGSetRecomputeResidual(
           hypre_F90_PassObj (HYPRE_Solver, solver),
           hypre_F90_PassInt (recomputeResidual) ) );
}

/*--------------------------------------------------------------------------
 * HYPRE_PCGSetRecomputeResidualP
 *--------------------------------------------------------------------------*/

void
hypre_F90_IFACE(hypre_pcgsetrecomputeresidualp, HYPRE_PCGSETRECOMPUTERESIDUALP)
   ( hypre_F90_Obj *solver,
     hypre_F90_Int *recomputeResidualP,
     hypre_F90_Int *ierr    )
{
   *ierr = (hypre_F90_Int)
      ( HYPRE_PCGSetRecomputeResidualP(
           hypre_F90_PassObj (HYPRE_Solver, solver),
           hypre_F90_PassInt (recomputeResidualP) ) );
}


void
hypre_F90_IFACE(hypre_pcggetconverged, HYPRE_PCGGETCONVERGED)
   ( hypre_F90_Obj *solver,
     hypre_F90_Int *converged,
     hypre_F90_Int *ierr    )
{
   *ierr = (hypre_F90_Int)
      ( HYPRE_PCGGetConverged(
           hypre_F90_PassObj (HYPRE_Solver, solver),
           hypre_F90_PassIntRef (converged) ) );
}


#include "HYPRE_utilities.h"

void
hypre_F90_IFACE(hypre_describeerror, HYPRE_DESCRIBEERROR)
   ( hypre_F90_Int *hypre_ierr,
     char          *descr   )
{
  HYPRE_DescribeError(
           hypre_F90_PassInt (hypre_ierr),
           hypre_F90_PassObjRef (char, descr) );
}
