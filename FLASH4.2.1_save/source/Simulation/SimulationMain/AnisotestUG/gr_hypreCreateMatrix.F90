!!****if* source/Grid/GridSolvers/HYPRE/UG/gr_hypreCreateMatrix
!!
!!  NAME 
!!
!!   gr_hypreCreateMatrix
!!
!!  SYNOPSIS
!!
!!   call gr_hypreCreateMatrix(integer(IN)           :: iVar,
!!                             integer(IN)           :: iFactorB,
!!                             integer(IN)           :: iFactorA,
!!                             integer(IN)           :: bcTypes(6),
!!                             real(IN)              :: bcValues(2,6),
!!                             real(IN)              :: dt,
!!                             real(IN)              :: alpha,
!!                             integer(IN)           :: blockCount,
!!                             integer(IN)           :: blockList(blockCount),
!!                             logical(IN)           :: JacobiMatrix)
!!
!!
!!  DESCRIPTION 
!!      This routine computes one of the matrices A and B, depending on the
!!      logical input argument 'JacobiMatrix':
!!          Ax = b, where A is the matrix to be inverted
!!          B = MX, where M is a matrix whose product with iVar produces RHS B.
!!
!!      A*(df/dt) + C*f = div(B*grad(f)) + D
!!      f -> Variable to be diffused.
!!      C,D are optional factors (not implemented here, the caller should add them later.)
!!
!!
!! ARGUMENTS
!! 
!!   iVar         : Variable on which the diffusion operation is performed (e.g., TEMP_VAR)
!!   iFactorB     : a factors in the equation with spatial variation.
!!   iFactorA     : a factors in the equation with spatial variation.
!!   bcTypes      : Boundary condition types.  Should be chosen from the constants
!!                  GRID_PDE_BND_* defined in Grid_interface.F90, or the special value VACUUM
!!                  defined in constants.h.
!!                  Presently OUTFLOW and VACUUM are supported, DIRICHLET less well tested.
!!   bcValues     : Values of iVar,iFactorB (!DEV: ??) on boundary (currently used for DIRICHLET and GIVENGRAD BCs).                        
!!   dt           : The time step.
!!   alpha        : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank-Nicolson
!!   blockCount   : The number of blocks in the list.   
!!   blockList    : The list of blocks on which the solution must be updated.        
!!   JacobiMatrix : TRUE computes A; FALSE computes M.
!!
!! SIDE EFFECTS
!!
!!   On return, the elements of the HYPRE matrix A that represent the second-derivative
!!   operator (div B grad) term have been defined, and it is ready for use.
!!   The gr_hypreData module variable gr_hypreMatA holds the
!!   handle for the HYPRE Solver object.
!!
!! NOTES
!!
!!   This routine does not actually 'create' the matrix object in the sense of HYPRE.
!!   It expects a matrix object already created and initialized; that is done,
!!   together with initialization of the grid object, in gr_hypreSetupGrid.
!!
!!   Currently, gr_hypreCreateMatrixFcB is called from Grid_advanceDiffusion with
!!   JacobiMatrix==.FALSE. only when the implicitness parameter theta passed to
!!   Grid_advanceDiffusion is 0. (KW 2012-12-05, corrected 2014-12-05)
!!
!! SEE ALSO
!!
!!  Grid_interface
!!***

!!REORDER(4): solnVec

subroutine gr_hypreCreateMatrix(iVar, iFactorB, iFactorA, bcTypes, bcValues, dt, &
     alpha, blockCount, blockList, JacobiMatrix)
  
  use gr_hypreLocalInterface, ONLY: gr_hypreApplyBcToFace
  use gr_hypreData,     ONLY : gr_hypreLower, gr_hypreUpper, &
                               gr_hypreMatA, &
                               gr_hypreAnisoDiffusion
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_fillGuardCells, Grid_getBlkBC, &
    Grid_getBlkCornerID, Grid_getCellCoords, Grid_getBlkData, Grid_getDeltas
  use Timers_interface, ONLY : Timers_start, Timers_stop 
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET
  use Driver_interface, ONLY : Driver_abortFlash
  
  implicit none
#include "Flash.h"  
#include "constants.h"
#include "HYPREf.h"    
  
  !!-----------------------------------------------------------------------
  !!         ARGUMENTS
  !!-----------------------------------------------------------------------
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
  real,    intent(IN) :: dt
  real,    intent(IN) :: alpha
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  logical, intent(IN) :: JacobiMatrix
      
  !!-----------------------------------------------------------------------
  !!         LOCAL VARIABLES.
  !!-----------------------------------------------------------------------  
  integer :: ierr
  real :: bxim12,bxip12,bxjm12,bxjp12,bxkm12,bxkp12,byim12,byip12,byjm12,byjp12,bykm12
  real :: bykp12,bzim12,bzip12,bzjm12,bzjp12,bzkm12,bzkp12 
  real, dimension(MDIM)     :: del
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec, solnData !,facexdata,faceydata,facezdata
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits 
  integer :: datasize(MDIM), datasizeGC(MDIM)
  integer, parameter ::  mypart = 0  !! HYPRE part
  integer ::  var
  integer ::  blockID
  integer ::  nentries, stencil_indices(19)
  integer :: i, j, k,  lb
  integer, dimension(2,MDIM):: faces 
  real :: condimh, condiph
  real :: condjmh, condjph
  real :: condkmh, condkph
  real :: term1, term2, term3, term4, term5, term6
  real :: theta, eps
  real :: dirichlet_multiplier  
  integer :: dir, ii
  real, allocatable :: faceAreas  (:,:,:,:)   
  real, allocatable :: BoxVal(:)

  !!-----------------------------------------------------------------------
  !! Anisotropic diffuusion; intended for
  !! anisotropic heat conduction along magnetic fields with HYPRE.
  !! Clearly, this feature is intended for use with MHD.
  !! Pointers for face magnetic fields in each direction
!  real, pointer, dimension(:,:,:,:) :: Bfcx, Bfcy, Bfcz
!  real :: Bximh, Bxiph, Byimh, Byiph, Bzimh, Bziph
!  real :: Bxjmh, Bxjph, Byjmh, Byjph, Bzjmh, Bzjph
!  real :: Bxkmh, Bxkph, Bykmh, Bykph, Bzkmh, Bzkph
!!-----------------------------------------------------------------------


  call Timers_start("gr_hypreCreateMatrix") 
  
  if(JacobiMatrix) then
     theta = alpha
     dirichlet_multiplier = 1.0
  else
     theta = alpha - 1.0
     dirichlet_multiplier = 0.0     
  end if
  
!  if (.not. gr_hypreAnisoDiffusion) then
!     nentries = NDIM*2 + 1
!  else
     ! HYPRE for anisotropic diffusion only supports 2D and 3D.
     ! This is meant to be used with multidimensional MHD.
     if (NDIM == 1) then
        call Driver_abortFlash&
          ('The HYPRE GridSolver for anisotropic diffusion is only available for multidimensional setups')
     elseif (NDIM == 2) then
        nentries = 9
     elseif (NDIM == 3) then
        nentries = 19 !(27-8)
     endif
!  endif

  do i = 1, nentries
     stencil_indices(i) = i-1
  enddo
  
  var    = 0  !! var iterator.
  
  do lb = 1, blockCount 
     

     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID, solnVec)
!     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getDeltas(blockID, del)
     call Grid_getBlkBC (blockID, faces)


     !! For anisotropic thermal conduction along magnetic fields
!     if (gr_hypreAnisoDiffusion) then
!        if (NDIM > 1) then
!           ! Note that facevariables are only defined for multidimensional MHD
!           call Grid_getBlkPtr(blockID, facexdata, FACEX)
!           call Grid_getBlkPtr(blockID, faceydata, FACEY)
!           if (NDIM > 2) call Grid_getBlkPtr(blockID, facezdata, FACEZ)
!        endif
!     endif

if (lb .eq. 1) then
     datasize(1:MDIM)   = blkLimits(HIGH,1:MDIM)   - blkLimits(LOW,1:MDIM)  +1 
     datasizeGC(1:MDIM) = blkLimitsGC(HIGH,1:MDIM) - blkLimitsGC(LOW,1:MDIM)+1

     allocate(BoxVal(nentries*product(datasize(1:NDIM)))) !nentries * total number of grid cells per block 

!!$print*,product((/1,2,3/)),datasize(1:NDIM)
!!$print*,nentries,product(datasize(1:NDIM)),nentries*product(datasize(1:NDIM));stop
     !!-----------------------------------------------------------------------
     !!         COMPUTE CELL FACE AREAS
     !!-----------------------------------------------------------------------     

     allocate(faceAreas(NDIM, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))              
     
endif

     call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, & 
          blkLimitsGC(LOW,:), faceAreas(IAXIS,:,:,:), datasizeGC)        
     
#if NDIM >= 2
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), faceAreas(JAXIS,:,:,:), datasizeGC)       
     
#if NDIM == 3
     
     call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
          blkLimitsGC(LOW,:), faceAreas(KAXIS,:,:,:), datasizeGC)   
#endif
     
#endif

     eps = 0. !1.e-20
     ii = 1
     BoxVal = 0.0

     !! We construct the matrix A (of AX=b) here, without taking into account
     !! those values in GC regions. In other words, BoxVal at
     !! ii = blkLimits(LOW,IAXIS) or blkLimits(HIGH, IAXIS) does not
     !! include condimh or condiph in this do-loop.
     !! We include these in the below by calling gr_hypreApplyBcToFace separately.
     do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
              !print*,JacobiMatrix;stop
              condimh = 0.
              condiph = 0.
              
! Compute B field on faces

! First Bx
              bxim12 = 0.5!*(solnVec(MAGX_VAR, i-1, j, k)+solnVec(MAGX_VAR, i, j, k))
              bxip12 = 0.5!*(solnVec(MAGX_VAR, i+1, j, k)+solnVec(MAGX_VAR, i, j, k))
              bxjm12 = 0.5!*(solnVec(MAGX_VAR, i, j-1, k)+solnVec(MAGX_VAR, i, j, k))
              bxjp12 = 0.5!*(solnVec(MAGX_VAR, i, j+1, k)+solnVec(MAGX_VAR, i, j, k))
              bxkm12 = 0.5!*(solnVec(MAGX_VAR, i, j, k-1)+solnVec(MAGX_VAR, i, j, k))
              bxkp12 = 0.5!*(solnVec(MAGX_VAR, i, j, k+1)+solnVec(MAGX_VAR, i, j, k))

!By

              byim12 = 0.!5*(solnVec(MAGY_VAR, i-1, j, k)+solnVec(MAGY_VAR, i, j, k))
              byip12 = 0.!5*(solnVec(MAGY_VAR, i+1, j, k)+solnVec(MAGY_VAR, i, j, k))
              byjm12 = 0.!5*(solnVec(MAGY_VAR, i, j-1, k)+solnVec(MAGY_VAR, i, j, k))
              byjp12 = 0.!5*(solnVec(MAGY_VAR, i, j+1, k)+solnVec(MAGY_VAR, i, j, k))
              bykm12 = 0.!5*(solnVec(MAGY_VAR, i, j, k-1)+solnVec(MAGY_VAR, i, j, k))
              bykp12 = 0.!5*(solnVec(MAGY_VAR, i, j, k+1)+solnVec(MAGY_VAR, i, j, k))

!Bz

              bzim12 = 0.!5*(solnVec(MAGZ_VAR, i-1, j, k)+solnVec(MAGZ_VAR, i, j, k))
              bzip12 = 0.!5*(solnVec(MAGZ_VAR, i+1, j, k)+solnVec(MAGZ_VAR, i, j, k))
              bzjm12 = 0.!5*(solnVec(MAGZ_VAR, i, j-1, k)+solnVec(MAGZ_VAR, i, j, k))
              bzjp12 = 0.!5*(solnVec(MAGZ_VAR, i, j+1, k)+solnVec(MAGZ_VAR, i, j, k))
              bzkm12 = 0.!5*(solnVec(MAGZ_VAR, i, j, k-1)+solnVec(MAGZ_VAR, i, j, k))
              bzkp12 = 0.!5*(solnVec(MAGZ_VAR, i, j, k+1)+solnVec(MAGZ_VAR, i, j, k))

!print *, bzim12,bzip12,bzjm12,bzjp12,bzkm12,bzkp12

              !! i-1,j,k
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then                         

                term1 = - 0.5*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12**2/(bxim12**2 + byim12**2+bzim12**2+eps)

                term2 = 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*bxjp12*byjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term3 = -0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*bxjm12*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                term4 = 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bxkp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                term5 = -0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bxkm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condimh = (term1+term2+term3+term4+term5) *faceAreas(IAXIS,i,j,k)

                BoxVal(ii+1) = -condimh*theta*dt/del(IAXIS)

              end if
              
              !! i+1,j,k
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then

               term1 = - 0.5*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12**2/(bxip12**2+byip12**2+bzip12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*bxjp12*byjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term3 = 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*bxjm12*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                term4 = - 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bxkp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                term5 = 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bxkm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condiph = (term1+term2+term3+term4+term5) *faceAreas(IAXIS,i+1,j,k)

                BoxVal(ii+2) =  -condiph*theta*dt/(del(IAXIS))

              end if
              
!              BoxVal(ii) =  BoxVal(ii) + (theta*((Condimh/del(IAXIS))+(Condiph/del(IAXIS)))*dt)    
 
#if NDIM >= 2
              condjmh = 0.
              condjph = 0.

              !! i,j-1,k
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then                                    
!                 condjmh = 0.5*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*faceAreas(JAXIS,i,j,k)

                term1 = - 0.5*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                term2 =  0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*bxip12*byip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term3 = - 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*bxjm12*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                term4 = 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bykp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                term5 = - 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bykm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condjmh = (term1+term2+term3+term4+term5) *faceAreas(JAXIS,i,j,k)

                BoxVal(ii+3) = -condjmh*theta*dt/(del(JAXIS))

              end if
              
              !! i,j+1,k
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then                      
!                 condjph = 0.5*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*faceAreas(JAXIS,i,j+1,k)

                 term1 = - 0.5*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*byjp12**2/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                 term2 = - 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*bxip12*byip12/(bxip12**2+byip12**2+bzip12**2+eps)

                 term3 = 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*bxjm12*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                 term4 = - 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bykp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                 term5 = 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bykm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                 condjph = (term1+term2+term3+term4+term5) *faceAreas(JAXIS,i,j+1,k)

                 BoxVal(ii+4) = -condjph*theta*dt/(del(JAXIS))

              endif
              
!              BoxVal(ii) =  BoxVal(ii) +  theta*(Condjmh/del(JAXIS)+Condjph/del(JAXIS))*dt !! diag           
              


                    !! i-1,j-1,k                                                                                                                          
               if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then
               if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then

!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                            

                  term1 = - 0.125*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12*byim12/(bxim12**2+byim12**2+bzim12**2+eps)

                  term2 = - 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*bxjm12*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                  condjmh = (term1+term2) *faceAreas(JAXIS,i,j,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                                      
                  BoxVal(ii+18) = -condjmh*theta*dt/(del(JAXIS))

               endif
               endif


                     !! i+1,j+1,k                                                                                                                                                                                                                                        
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                             
                term1 = - 0.125*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12*byip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*bxjp12*byjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                condiph = (term1+term2) *faceAreas(JAXIS,i,j,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                                     

                BoxVal(ii+12) = -condiph*theta*dt/(del(JAXIS))

              endif
              endif

                    !! i-1,j+1,k                                                                                                                                                                                             
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                            
                    term1 = 0.125*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12*byim12/(bxim12**2 + byim12**2+bzim12**2+eps)

                    term2 = 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*bxjp12*byjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                    condjmh = (term1+term2) *faceAreas(JAXIS,i,j,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                       

                    BoxVal(ii+16) = -condjmh*theta*dt/(del(JAXIS))

               endif
               endif

              !! i+1,j-1,k                                                                                                                                                                                                                                        
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
               if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                            

                term1 = 0.125*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12*byip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term2 = 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*bxjm12*byjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                condjmh = (term1+term2) *faceAreas(JAXIS,i,j,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                                     
                BoxVal(ii+14) = -condjmh*theta*dt/(del(JAXIS))

              endif
              endif

              
#if NDIM == 3
              condkmh = 0.
              condkph = 0.
              
              !! i,j,k-1
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then                                    

                term1 = - 0.5*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bzkm12**2/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                term2 = 0.125*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12*bzip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term3 = - 0.125*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12*bzim12/(bxim12**2 + byim12**2+bzim12**2+eps)

                term4 = 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*byjp12*bzjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term5 = - 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*byjm12*bzjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                condkmh = (term1+term2+term3+term4+term5) *faceAreas(KAXIS,i,j,k)

                BoxVal(ii+5) = -condkmh*theta*dt/(del(KAXIS))               

              endif
              
              !! i,j,k+1
              if ((k /= blkLimits(HIGH, KAXIS)) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) then                      

!                 condkph = 0.5*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*faceAreas(KAXIS,i,j,k+1)
                term1 = - 0.5*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bzkp12**2/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12*bzip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term3 = 0.125*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12*bzim12/(bxim12**2 + byim12**2+bzim12**2+eps)

                term4 = - 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*byjp12*bzjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term5 = 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*byjm12*bzjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                condkph = (term1+term2+term3+term4+term5) *faceAreas(KAXIS,i,j,k+1)

                 BoxVal(ii+6) = -condkph*theta*dt/(del(KAXIS))                                        

              endif
              
!              BoxVal(ii) =  BoxVal(ii) +  theta*(Condkmh/del(KAXIS)+Condkph/del(KAXIS))*dt !! diag  


              !! i,j+1,k-1                                                                                                                                                                                                                                        
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then

!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                            
                term1 = 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*byjp12*bzjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term2 = 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bykm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                                    
                 BoxVal(ii+9) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i+1,j,k-1                                                                                                                                                                                                                                        
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then

!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                            
                term1 = 0.125*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12*bzip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term2 = 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bxkm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j,k) ! WHICH AREA DO WE USE HERE?

                BoxVal(ii+13) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i,j-1,k-1                                                                                                                                                                                                                                        
               if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
 
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction
                term1 = - 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*byjm12*bzjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bykm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?

                BoxVal(ii+10) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i-1,j,k-1                                       
              if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                            
                term1 = - 0.125*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12*bzim12/(bxim12**2 + byim12**2+bzim12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bxkm12*bzkm12/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                                   
                BoxVal(ii+17) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i,j+1,k+1                                                                                                                                                                                                                                
              if ((j /= blkLimits(HIGH, JAXIS)) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then

!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                       
                term1 = - 0.125*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*byjp12*bzjp12/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bykp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                            
                BoxVal(ii+7) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i+1,j,k+1                                                                                                                                                                                                                                
              if (i /= blkLimits(HIGH, IAXIS) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                           
                term1 = - 0.125*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12*bzip12/(bxip12**2+byip12**2+bzip12**2+eps)

                term2 = - 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bxkp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                  
                BoxVal(ii+11) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i,j-1,k+1                                                                                                                                                                                                                            
              if ((j /= blkLimits(LOW, JAXIS)) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                  
                term1 = 0.125*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*byjm12*bzjm12/(bxjm12**2+byjm12**2+bzjm12**2+eps)

                term2 = 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bykp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                            
                BoxVal(ii+8) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

              !! i-1,j,k+1                                                                                                                                                                                                                            
               if ((i /= blkLimits(LOW, IAXIS)) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) then
              if ((k /= blkLimits(LOW, KAXIS)) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) then
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                       
                    term1 = 0.125*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12*bzim12/(bxim12**2 + byim12**2+bzim12**2+eps)

                    term2 = 0.125*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bxkp12*bzkp12/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                    condkmh = (term1+term2) *faceAreas(KAXIS,i,j+1,k) ! WHICH AREA DO WE USE HERE?                                                                                                                                                   
                    BoxVal(ii+15) = -condkmh*theta*dt/(del(KAXIS))

              endif
              endif

                    !! i,j,k                                                                                                                                                                                                             
!MB Here is a first attempt at filling the matrix elements for anisotropic conduction                                                                                                                                                                      
                term1 = 0.5*(solnVec(iFactorB,i+1,j,k)+ solnVec(iFactorB,i,j,k))*bxip12**2/(bxip12**2+byip12**2+bzip12**2+eps)

                term2 = 0.5*(solnVec(iFactorB,i-1,j,k)+ solnVec(iFactorB,i,j,k))*bxim12**2/(bxim12**2 + byim12**2+bzim12**2+eps)

                term3 = 0.5*(solnVec(iFactorB,i,j+1,k)+ solnVec(iFactorB,i,j,k))*byjp12**2/(bxjp12**2+byjp12**2+bzjp12**2+eps)

                term4 = 0.5*(solnVec(iFactorB,i,j-1,k)+ solnVec(iFactorB,i,j,k))*byjm12**2/(bxjm12**2+byjm12**2+bzjm12**2+eps)


                term5 = 0.5*(solnVec(iFactorB,i,j,k+1)+ solnVec(iFactorB,i,j,k))*bzkp12**2/(bxkp12**2+bykp12**2+bzkp12**2+eps)

                term6 = 0.5*(solnVec(iFactorB,i,j,k-1)+ solnVec(iFactorB,i,j,k))*bzkm12**2/(bxkm12**2+bykm12**2+bzkm12**2+eps)

                condkmh = (term1+term2+term3+term4+term5+term6) *faceAreas(KAXIS,i,j,k)

                BoxVal(ii) = -condkmh*theta*dt/(del(KAXIS))

#endif             
              
#endif               
              ii = ii + nentries


              
           end do
        end do
     end do

     call HYPRE_SStructMatrixSetBoxValues(gr_hypreMatA, mypart, gr_hypreLower(lb,1:NDIM), & 
          gr_hypreUpper(lb,1:NDIM), var, nentries, stencil_indices(1:nentries), BoxVal(:), ierr)
     

     !! We fix BoxVal at the block boundaries.
     dir = ILO_FACE
     do i = IAXIS, NDIM
        do j = LOW, HIGH
           if (faces(j,i) /= NOT_BOUNDARY) then               
              call gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,mypart,var,iFactorB,bcTypes(dir),dir, &
                   bcValues(:,dir), dt, theta, del(i), gr_hypreLower(lb,:), dirichlet_multiplier, faceAreas(i,:,:,:), solnVec, &
                   blockID)
           end if
           dir = dir + 1
        end do
     end do

     
     call Grid_releaseBlkPtr(blockID, solnVec)
!     call Grid_releaseBlkPtr(blockID, solnData)

     !! For anisotropic thermal conduction along magnetic fields
!     if (gr_hypreAnisoDiffusion) then
!        if (NDIM > 1) then
!           ! Note that facevariables are only defined for multidimensional MHD
!           call Grid_releaseBlkPtr(blockID, facexdata, FACEX)
!           call Grid_releaseBlkPtr(blockID, faceydata, FACEY)
!           if (NDIM > 2) call Grid_releaseBlkPtr(blockID, facezdata, FACEZ)
!        endif
!     endif
          
!print *, 'all good 2', lb
     
  end do !! block

     deallocate (faceAreas)       
     deallocate (BoxVal)

  !!-----------------------------------------------------------------------
  !!         THIS IS A GLOBAL CALL.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixAssemble(gr_hypreMatA, ierr)    
  
  call Timers_stop("gr_hypreCreateMatrix") 
  
  return
  
end subroutine gr_hypreCreateMatrix
