!!****f* source/Grid/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!!  This routine is normally called by the implementation of
!!  Grid_updateRefinement.
!!
!! ARGUMENTS
!!
!!  none
!! 
!! SIDE EFFECTS
!!  
!!  This routine works by modifying the global (processor-local) flag
!!  arrays
!!      newchild, refine, derefine, stay
!!  imported from the AMR implementation.
!!
!! NOTES
!!  
!!  The default implementation of this routine uses the Flash runtime
!!  parameters refine_var_{1,2,3,4}, refine_cutoff_{1,2,3,4},
!!  derefine_cutoff_{1,2,3,4}, and refine_filter_{1,2,3,4} to determine
!!  refinement.
!!  
!!  A non-directional guardcell fill for CENTER (and also EOS calls
!!  for all block cells, including guardcells, if any refinement
!!  variables refine_var_# require this to be current) must have been
!!  performed when this routine is invoked. Moreover, there must not
!!  be any intervening calls that would modify the solution data in
!!  unk (at least, for the variables to be used for refinement
!!  criteria).
!!  
!!  Users creating their own implementation of this interface or of
!!  Grid_updateRefinement should make sure that the above remains
!!  true; or should include the appropriate call(s) to
!!  Grid_fillGuardCells in their implementation.
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_fillGuardCells
!!***



subroutine Grid_markRefineDerefine()

implicit none
   
end subroutine Grid_markRefineDerefine














