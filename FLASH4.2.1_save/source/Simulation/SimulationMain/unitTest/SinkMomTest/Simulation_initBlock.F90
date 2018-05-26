

subroutine Simulation_initBlock (blockId)

  use Simulation_data, ONLY : sim_xMin, sim_yMin, sim_zMin, &
       sim_xMax, sim_yMax, sim_zMax, sim_xcenter, sim_ycenter, sim_zcenter, &
       sim_radius, sim_dens, sim_cs, sim_vx, sim_vy, sim_vz
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_putPointData
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype

  use Eos_interface, ONLY : Eos_wrapped
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  
  integer,intent(IN) :: blockId
  
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  integer, dimension(MDIM) :: axis
  real, dimension(:), allocatable :: x, y, z
  logical :: gcell = .true.
  integer :: i, j, k, istat, MyPE
   
  real :: radius, rho, vx, vy, vz, pres
  real :: xdist, ydist, zdist

  logical, parameter :: Debug = .false.

  if (Debug) print *, 'Simulation_initBlock entering.'

  call Driver_getMype(GLOBAL_COMM, MyPE)

  call RuntimeParameters_get("sim_xcenter", sim_xcenter)
  call RuntimeParameters_get("sim_ycenter", sim_ycenter)
  call RuntimeParameters_get("sim_zcenter", sim_zcenter)
  call RuntimeParameters_get("sim_vx", sim_vx)
  call RuntimeParameters_get("sim_vy", sim_vy)
  call RuntimeParameters_get("sim_vz", sim_vz)
  call RuntimeParameters_get("sim_radius", sim_radius)
  call RuntimeParameters_get("sim_dens", sim_dens)  
  call RuntimeParameters_get("sim_cs", sim_cs)

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
  allocate(x(sizeX), stat=istat)
  allocate(y(sizeY), stat=istat)
  allocate(z(sizeZ), stat=istat)
  call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, z, sizeZ)
  call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, y, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, x, sizeX)

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

           xdist = x(i) - sim_xcenter
           ydist = y(j) - sim_ycenter
           zdist = z(k) - sim_zcenter
    
           radius = sqrt( xdist**2 + ydist**2 + zdist**2)
                      
           if (radius .le. sim_radius) then
              rho = sim_dens
              pres = sim_cs*sim_cs*rho
           else
              rho = 1e-2*sim_dens
              pres = sim_cs*sim_cs*sim_dens
           end if

           vx = sim_vx
           vy = sim_vy
           vz = sim_vz
    
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k
           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, pres)

        end do
     end do
  end do
  
  call Eos_wrapped(MODE_DENS_PRES, blklimitsGC, blockId)

  deallocate(x)
  deallocate(y)
  deallocate(z)


  if (Debug) print *, 'Simulation_initBlock done.'

  return

end subroutine Simulation_initBlock
