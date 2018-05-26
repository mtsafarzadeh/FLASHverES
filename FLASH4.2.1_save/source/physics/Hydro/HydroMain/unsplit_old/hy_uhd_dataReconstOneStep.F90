!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_dataReconstOneStep
!!
!! NAME
!!
!!  hy_uhd_dataReconstOneStep
!!
!! SYNOPSIS
!!
!!  hy_uhd_dataReconstOneStep(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(:,:),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(IN)    :: dt,
!!                            real(IN)    :: del(MDIM),
!!                            real(IN)    :: ogravX(:,:,:),
!!                            real(IN)    :: ogravY(:,:,:),
!!                            real(IN)    :: ogravZ(:,:,:),
!!                            real(IN)    :: V0   (HY_VARINUMMAX),
!!                            real(IN)    :: Vxp  (HY_VARINUMMAX),
!!                            real(IN)    :: Vxn  (HY_VARINUMMAX),
!!                            real(IN)    :: Vyp  (HY_VARINUMMAX),
!!                            real(IN)    :: Vyn  (HY_VARINUMMAX),
!!                            real(IN)    :: Vzp  (HY_VARINUMMAX),
!!                            real(IN)    :: Vzn  (HY_VARINUMMAX),
!!                            real(IN)    :: Vxpp (HY_VARINUMMAX),
!!                            real(IN)    :: Vxnn (HY_VARINUMMAX),
!!                            real(IN)    :: Vypp (HY_VARINUMMAX),
!!                            real(IN)    :: Vynn (HY_VARINUMMAX),
!!                            real(IN)    :: Vzpp (HY_VARINUMMAX),
!!                            real(IN)    :: Vznn (HY_VARINUMMAX),
!!                            real(IN)    :: Vxnnn(HY_VARINUMMAX),
!!                            real(IN)    :: Vynnn(HY_VARINUMMAX),
!!                            real(IN)    :: Vznnn(HY_VARINUMMAX),
!!                            real(IN)    :: FlatCoeff(:,:,:,:),
!!                            logical(IN) :: TransX_updateOnly,
!!                            logical(IN) :: TransY_updateOnly,
!!                            logical(IN) :: TransZ_updateOnly,
!!                            real(OUT)   :: Wxp(HY_VARINUMMAX),
!!                            real(OUT)   :: Wxn(HY_VARINUMMAX),
!!                            real(OUT)   :: Wyp(HY_VARINUMMAX),
!!                            real(OUT)   :: Wyn(HY_VARINUMMAX),
!!                            real(OUT)   :: Wzp(HY_VARINUMMAX),
!!                            real(OUT)   :: Wzn(HY_VARINUMMAX),
!!                            real(OUT)   :: sig(HY_VARINUMMAX,NDIM),
!!                            real(OUT)   :: lambda(HY_WAVENUM,NDIM),
!!                            real(OUT)   :: leig(HY_WAVENUM,HY_VARINUM,NDIM),
!!                            real(OUT)   :: reig(HY_VARINUM,HY_WAVENUM,NDIM))
!!
!! ARGUMENTS
!!
!!  blockID     - local block ID
!!  blkLimitsGC - block limits including guardcells
!!  ix,iy,iz    - local indices
!!  dt          - timestep
!!  del         - deltas in each {x,y,z} direction
!!  ogravX,ogravY,ogravZ - gravity components at n step in each direction
!!  V0          - array containing primitive variables + gamc + game at cell center
!!  Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Vxpp,Vxnn,Vypp,Vynn,Vzpp,Vznn,Vxnnn,Vynnn,Vznnn - data in neighboring cells
!!  FlatCoeff   - flattening parameters
!!  TransX_updateOnly - a switch for a selective transverse update in x direction
!!  TransY_updateOnly - a switch for a selective transverse update in y direction
!!  TransZ_updateOnly - a switch for a selective transverse update in z direction
!!  Wxp,Wxn,Wyp,Wyn,Wzp,Wzn - left(m) and right(p) states
!!  sig    - a transverse flux term
!!  lambda - eigenvalue
!!  leig   - left eigenvector
!!  reig   - right eigenvector
!!
!!
!! DESCRIPTION
!!
!!  This onestep data reconstruction routine evolves the cell centered 
!!  values by dt/2 time step at each cell interface using characteristic 
!!  tracing method based on 2nd order MUSCL-Hancock, 3rd order PPM,
!!  or 5th order WENO schemes.
!!  For MHD, the cell interface values are reconstructed in multidimensional
!!  fashion with inclusions of MHD multidimensional terms that are
!!  proportional to divB for the unsplit MHD implementation.
!!
!! REFERENCES
!!
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Stone, Gardiner, Teuben, Hawley, Simon, "Athena: A new code for astrophysical MHD"
!!    arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
!!  * Colella and Woodward, 54, 174 (1984), JCP
!!  * Colella, 87, 171-200 (1990), JCP
!!
!!***

Subroutine hy_uhd_dataReconstOneStep(blockID,blkLimitsGC,ix,iy,iz, &
                                     dt,del,ogravX,ogravY,ogravZ,V0,     &
                                     Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn,  &
                                     Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn, &
                                     Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn,&
                                     FlatCoeff, &
                                     TransX_updateOnly,&
                                     TransY_updateOnly,&
                                     TransZ_updateOnly,&
                                     Wxp, Wxn, Wyp, Wyn, Wzp, Wzn, &
                                     sig,lambda,leig,reig )

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: blockID
  integer,intent(IN),dimension(LOW:HIGH,MDIM):: blkLimitsGC
  integer,intent(IN) :: ix,iy,iz
  real,   intent(IN) :: dt
  real,   intent(IN), dimension(MDIM) :: del
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: ogravX,ogravY,ogravZ
  real, dimension(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) :: FlatCoeff
#else
  real, dimension(blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: ogravX,ogravY,ogravZ
  real, dimension(NDIM,blkLimitsGC(HIGH,IAXIS),blkLimitsGC(HIGH,JAXIS),blkLimitsGC(HIGH,KAXIS)), &
       intent(IN) :: FlatCoeff
#endif
  logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly
  real, intent(INOUT),  dimension(HY_VARINUMMAX) :: V0, Vxp,  Vxn,  Vyp,  Vyn,  Vzp,  Vzn, &
                                                    Vxpp, Vxnn, Vypp, Vynn, Vzpp, Vznn,&
                                                    Vxppp,Vxnnn,Vyppp,Vynnn,Vzppp,Vznnn
  real, intent(OUT), dimension(HY_VARINUMMAX)    :: Wxp, Wxn, Wyp, Wyn, Wzp, Wzn
  real, intent(OUT), dimension(HY_VARINUMMAX,NDIM) :: sig
  real, intent(OUT), dimension(HY_WAVENUM,NDIM) :: lambda
  real, intent(OUT), dimension(HY_WAVENUM,HY_VARINUM,NDIM) :: leig
  real, intent(OUT), dimension(HY_VARINUM,HY_WAVENUM,NDIM) :: reig

End subroutine hy_uhd_dataReconstOnestep
