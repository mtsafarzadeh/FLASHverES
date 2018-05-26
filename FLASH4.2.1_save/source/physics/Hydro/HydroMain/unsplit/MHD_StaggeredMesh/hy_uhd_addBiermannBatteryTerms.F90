!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_addBiermannBatteryTerms
!!
!! NAME
!!
!!  hy_uhd_addBiermannBatteryTerms
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addBiermannBatteryTerms(blockID,blkLimitsGC,ix,iy,iz,Flux,magVisc,sweepDir)
!!
!!  hy_uhd_addBiermannBatteryTerms(integer(IN) :: blockID,
!!                            integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM),
!!                            integer(IN) :: ix,
!!                            integer(IN) :: iy,
!!                            integer(IN) :: iz,
!!                            real(IN)    :: Flux,
!!                            integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds Biermann battery terms to total MHD fluxes. This is a stub. We are not supporting this yet. 
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section 
!!                of block with the guard cells 
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing MHD fluxes
!!  sweepDir    - direction of sweep
!!
!!***

!!REORDER(4): U

Subroutine hy_uhd_addBiermannBatteryTerms(blockID,blkLimitsGC,ix,iy,iz,Flux,sweepDir)

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas
  use Eos_interface,  ONLY : Eos_wrapped, Eos_getAbarZbar
  use Hydro_data,     ONLY : hy_avogadro, hy_useBiermann, hy_qele, hy_biermannCoef, hy_speedOfLight
  use hy_uhd_slopeLimiters, ONLY : checkMedian, mc, minmod, vanLeer

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC 
  real, dimension(HY_VARINUM), intent(INOUT) :: Flux
  integer, INTENT(IN) :: sweepDir
  !! ----------------------------------------------------------------------

  return
  
End Subroutine hy_uhd_addBiermannBatteryTerms


Subroutine get_upwind(vel,pe_L,pe_R,pe_Up)

  use hy_uhd_slopeLimiters, ONLY : signum

  implicit none
  real, intent(IN)  :: vel,pe_L,pe_R
  real, intent(OUT) :: pe_Up

  real :: velP,velN

  velP = 0.5*(1.+signum(vel)) !*abs(velo)
  velN = 0.5*(1.-signum(vel)) !*abs(velo)

  pe_Up = velP*pe_L + velN+pe_R

end Subroutine get_upwind
