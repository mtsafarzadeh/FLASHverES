!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up a
!!  Sod-like problem in spherical coordinates to test whether a planar
!!  shock in spherical coordinates stays planar.  This effectively tests the 
!!  fictitious forces in force().  If the forces are setup right, then the planar
!!  shock should stay planar.
!!
!!  Right now, this is setup to do the problem in 2-d spherical coordinates.  
!!  sim_idir = 1 is the x-direction.  sim_idir = 2 is the z-direction.  
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!  
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY: Eos_wrapped
  use Eos_data, ONLY: eos_singleSpeciesA


  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: ii, jj, kk

  
  integer, intent(in) :: blockID
  

  real, allocatable, dimension(:) :: yl, yr

  integer, parameter :: nsub = 5
  integer :: nsubj

  logical :: high_state
      
  real :: x_zone, y_zone, z_zone

  real :: rl_zone, r_zone, rr_zone
  real :: thetal_zone, theta_zone, thetar_zone
  real :: phi_zone


  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  


  real :: xx, yy,  zz, xxL, xxR
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone
  
  logical :: gcell = .true.
  integer :: istat 
  real    :: rtilde  !ES

  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX),stat=istat)
  allocate(xRight(sizeX),stat=istat)
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(yl(sizeY),stat=istat)
  allocate(yr(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  yl(:) = 0.0
  yr(:) = 0.0
  
  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) then
      call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
      call Grid_getCellCoords&
                      (JAXIS, blockId, LEFT_EDGE,gcell, yl, sizeY)
      call Grid_getCellCoords&
                      (JAXIS, blockId, RIGHT_EDGE,gcell, yr, sizeY)
      nsubj = nsub
  else
     blkLimits(LOW,JAXIS) = 1
     blkLimits(HIGH,JAXIS) = 1
     yCoord(1) = PI / 2
     yl(1) = yCoord(1) - 1.0e-10
     yr(1) = yCoord(1) + 1.0e-10
     nsubj = 1
  end if

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)


!-----------------------------------------------------------------------------
! loop over all the zones in the current block and set the temperature,
! density, and thermodynamics variables.
!-----------------------------------------------------------------------------

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     axis(KAXIS) = k
     zz = zCoord(k)
    
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        axis(JAXIS) = j
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS) = i
           
           rtilde = xCenter(i)/sim_rvir

           if(rtilde.le.1) then             
              rhoZone = sim_rho0*exp(-sim_Ac*(1-alog(1+sim_nfwc*rtilde)/(sim_nfwc*rtilde)))
           else
              rhoZone = sim_rho0*exp(-sim_Ac*(1-alog(1+sim_nfwc)/sim_nfwc))
              rhoZone = rhoZone*exp(-(2-2/rtilde))
           endif
           presZone = rhoZone/1.6726219E-24/eos_singleSpeciesA*1.38E-16*sim_tvir
           if(xCenter(i) .le. sim_rblast) &
              presZone = presZone+(2./3.)*rhoZone*sim_eblast/sim_mblast/1.98E33
           velxZone = 0.0
           velyZone = 0.0
           velzZone = 0.0


           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone) !used in Eos
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone) !used in Eos
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, 0.0) !dummy use in Eos_wrapped
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone) !used in Eos_wrapped
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone) !used in Eos_wrapped
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone) !used in Eos_wrapped

           !put in value of default species. Ignored in Gamma Eos.
!           if (NSPECIES > 0) then 
!             call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
!                        axis, 1.0e0-(NSPECIES-1)*sim_smallX)
!             !if there is only 1 species, this loop will not execute
!              do n = SPECIES_BEGIN+1,SPECIES_END
!                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
!                      axis, sim_smallX)
!              enddo
!           end if 


#ifdef UNDEFINED
           call Grid_putPointData(blockId, CENTER, VOLU_VAR, EXTERIOR, axis, 0.0)
#endif
#ifdef TEMP_VAR
!           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, temp_zone) !set in Eos
#endif

#ifdef ENER_VAR
!           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone) !set in Eos
#endif
#ifdef GAME_VAR          
!           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma) !set in Eos
#endif
#ifdef GAMC_VAR
!           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma) !set in Eos
#endif

        enddo
     enddo
  enddo

  call Eos_wrapped(MODE_DENS_PRES, blkLimits, blockID)

! Cleanup
  deallocate(xLeft,stat=istat)
  deallocate(xRight,stat=istat)
  deallocate(xCenter,stat=istat)
  deallocate(yCoord,stat=istat)
  deallocate(yl,stat=istat)
  deallocate(yr,stat=istat)
  deallocate(zCoord,stat=istat)

 
  return
end subroutine Simulation_initBlock










