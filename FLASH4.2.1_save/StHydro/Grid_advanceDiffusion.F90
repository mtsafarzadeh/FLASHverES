!!****f* source/Grid/Grid_advanceDiffusion
!!
!! NAME
!!  Grid_advanceDiffusion
!!
!! SYNOPSIS
!!
!!  call Grid_advanceDiffusion (integer, intent(IN) :: iVar,
!!                              integer(IN)         :: iSrc,
!!                              integer, intent(IN) :: iFactorB,
!!                              integer, intent(IN) :: iFactorA,
!!                              integer, intent(IN) :: bcTypes(6),
!!                              real,    intent(IN) :: bcValues(2,6),
!!                              real,    intent(IN) :: dt,
!!                              real,    intent(IN) :: chi,
!!                              real,    intent(IN) :: scaleFact,
!!                              real,    intent(IN) :: theta,
!!                              integer, intent(IN), OPTIONAL :: iFactorC,
!!                              integer, intent(IN), OPTIONAL :: iFactorD
!!                              integer, OPTIONAL, intent(IN) :: pass,
!!                              logical(IN) :: solnIsDelta)
!!
!! DESCRIPTION
!!
!!      This routine advances a generalized diffusion operator of the form
!!
!!         A*(df/dt) + C*f = div(B*grad(f)) + D ,
!!
!!      where
!!         f -> f(x,t) is the  Variable to be diffused (x=1D..3D position);
!!
!!         A,B,C,D are optional given scalar factors/terms that may depend
!!         on position; they are either physcially constant in time, or at
!!         least considered time-independent for the purpose of the operation
!!         implemented here (typically by computing their values from the
!!         solution state reached by the previous time step).
!!
!!      Presently it is used to do heat conduction and multigroup diffusion.
!!
!! ARGUMENTS
!!   iVar           : Variable on which the diffusion operatorion is performed (e.g., TEMP_VAR)
!!   iFactorA       :| Are factors in the equation with spatial variation.
!!   iFactorB       :| Factor C,D are optional and are generally used
!!   iFactorC       :| to represent emission/absorption in MGD.
!!   iFactorD       :| iFactorA is needed only for conduction.
!!   bcTypes        : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested.
!!   bcValues       : Values of iVar,iFactorB on boundary (DIRICHLET).                        
!!   dt             : The time step.
!!   scaleFact      : Factor by which the end solution is scaled (not used).
!!   chi            : useful for constant diffusion problems (not used).
!!   theta          : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   pass           : Ignored in unsplit solver.
!!                    pass=1 order of directional sweep X-Y-Z, 
!!                    pass=2 order of directional sweep Z-Y-X.
!!   iSrc           : Ignored.
!!   solnIsDelta    : Is the solution only a delta that the caller has to apply to the
!!                    temperature, rather than temperature itself (ignored).
!!
!!
!! NOTES
!!
!!  It is currently assumed in implementations
!!  that boundary condition types at all sides of the domain are the same. 
!!  However, there is no check to make sure that this is true.
!!  What this means to the user  is that only the first value in bcTypes is
!!  checked here, and may be assumed to give the type of boundary condition
!!  for all 2*NDIM directions.
!!
!! SEE ALSO
!! 
!!  Diffuse_advance1D
!!  
!!
!!***


subroutine Grid_advanceDiffusion (iVar, iSrc, iFactorB, iFactorA, bcTypes, bcValues, dt, chi, scaleFact, &
     theta, solnIsDelta, iFactorC, iFactorD, pass)       
  
  implicit none
  
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iSrc
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  real, intent(IN)    :: dt 
  real, intent(IN)    :: chi
  real, intent(IN)    :: scaleFact
  real, intent(IN)    :: theta
  logical, intent(IN) :: solnIsDelta
  integer, dimension(6),  intent(IN) :: bcTypes
  real   , dimension(2,6),intent(IN) :: bcValues
  integer, intent(IN), OPTIONAL :: pass
  integer, intent(IN), OPTIONAL :: iFactorC
  integer, intent(IN), OPTIONAL :: iFactorD   
  
end subroutine Grid_advanceDiffusion
