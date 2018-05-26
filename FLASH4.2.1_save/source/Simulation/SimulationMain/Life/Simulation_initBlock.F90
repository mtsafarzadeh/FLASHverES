!!****if* source/Simulation/SimulationMain/Life/Simulation_initBlock
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
!!  a specified block.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_putPointData,      &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getDeltas
  use Driver_interface, ONLY: Driver_abortFlash
  use RadTrans_interface, ONLY: RadTrans_mgdEFromT

  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  
  integer :: i, j, k
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  real, allocatable :: xcent(:), ycent(:), zcent(:)
  real :: tradActual
  real :: rho, tele, trad, tion
  integer :: species
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData

  integer :: ndx, ndy, ndz, nsctot
  integer :: nsct, nscc, l, m, n
  real :: xs, ys, zs, xsc, ysc, zsc, targFrac, chamFrac, norm
  real :: dxsc, dysc, dzsc
  real :: delta(MDIM)
  real :: r, phi

  real :: cham_vol_frac
  real :: targ_vol_frac

  ndx = sim_ndiv
  ndy = 1
  ndz = 1
  if(NDIM >= 2) ndy = sim_ndiv
  if(NDIM == 3) ndz = sim_ndiv
  nsctot = ndx * ndy * ndz

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  allocate(xcent(blkLimitsGC(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., &
       xcent, blkLimitsGC(HIGH, IAXIS))
  allocate(ycent(blkLimitsGC(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., &
       ycent, blkLimitsGC(HIGH, JAXIS))
  allocate(zcent(blkLimitsGC(HIGH, KAXIS)))
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., &
       zcent, blkLimitsGC(HIGH, KAXIS))

  call Grid_getDeltas(blockId, delta)

  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Loop over cells and set the initial state
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           xs = xcent(i) - delta(IAXIS)/2
           ys = ycent(j) - delta(JAXIS)/2
           zs = zcent(k) - delta(KAXIS)/2

           dxsc = delta(IAXIS) / ndx
           dysc = delta(JAXIS) / ndy
           dzsc = delta(KAXIS) / ndz

           nsct = 0
           do l = 1, ndx
              do m = 1, ndy
                 do n = 1, ndz

                    xsc = xs + dxsc * (2*l - 1) / 2.0
                    ysc = ys + dysc * (2*m - 1) / 2.0
                    zsc = zs + dzsc * (2*n - 1) / 2.0

                    call findSpec(xsc, ysc, zcent(k), species)
                    if(species == TARG_SPEC) nsct = nsct + 1
                 end do
              end do
           end do

           nscc = nsctot - nsct
           targ_vol_frac = real(nsct) / real(nsctot)

           ! ! Set volume fraction using radius:
           ! r = sqrt(xcent(i)**2 + ycent(j)**2)
           ! targ_vol_frac = 0.5 * (-tanh(100.0 * (r - sim_targetRadius)/sim_targetRadius) + 1.0)

           cham_vol_frac = 1.0 - targ_vol_frac

           ! Average Cell Density:                                                               
           rho = sim_rhoTarg * targ_vol_frac + sim_rhoCham * cham_vol_frac

           ! Mass Fractions:                                                                     
           targFrac = (sim_rhoTarg*targ_vol_frac)/rho
           chamFrac = 1.0 - targFrac

           ! Set temperatures:
           tele = sim_teleTarg * targ_vol_frac + sim_teleCham * cham_vol_frac
           tion = sim_tionTarg * targ_vol_frac + sim_tionCham * cham_vol_frac
           trad = sim_tradTarg * targ_vol_frac + sim_tradCham * cham_vol_frac

           ! Set all initial conditions:
           solnData(DENS_VAR,i,j,k) = rho
           solnData(TEMP_VAR,i,j,k) = tele
           solnData(TION_VAR,i,j,k) = tion
           solnData(TELE_VAR,i,j,k) = tele

           targFrac = targFrac + sim_smallX
           chamFrac = chamFrac + sim_smallX
           solnData(TARG_SPEC,i,j,k) = targFrac/(targFrac + chamFrac)
           solnData(CHAM_SPEC,i,j,k) = chamFrac/(targFrac + chamFrac)

           if(nsct == 0) then
              ! This cell only contains chamber material
              solnData(VELX_VAR,i,j,k) = sim_velxCham
              solnData(VELY_VAR,i,j,k) = sim_velyCham
           end if

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, (/i,j,k/), trad, tradActual)
           solnData(TRAD_VAR,i,j,k) = tradActual

        enddo
     enddo
  enddo


  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

contains

  subroutine findSpec(x,y,z, spec)
    implicit none

    real, intent(in) :: x
    real, intent(in) :: y
    real, intent(in) :: z
    
    integer, intent(out) :: spec
    
    spec = CHAM_SPEC

    select case (NDIM)
    case (1)

       select case (sim_geometry)
       case (SPHERICAL)
          if (sim_targetGeom == "sphere") then
             if( x < sim_targetRadius ) &
                  spec = TARG_SPEC
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case default
          call Driver_abortFlash("[Simulation_initBlock] Unsupported geometry")
       end select

    case (2)       

       select case (sim_geometry)
       case (CARTESIAN)

          if(sim_targetGeom == "cylinder") then
             
             r = sqrt(x**2 + (y-sim_targetZOffset)**2)
             phi = atan(y/x)

             if(r <= sim_targetRadius) then
                spec = TARG_SPEC
             end if
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case (CYLINDRICAL)

          
          if(sim_targetGeom == "cylinder") then
             
             if(abs(x-sim_targetOffset) <= sim_targetRadius .and. &
                  abs(y-sim_targetZOffset) <= sim_targetHeight/2 ) then
                spec = TARG_SPEC
             end if

          elseif (sim_targetGeom == "sphere") then
             
             if(sqrt(x**2 + (y-sim_targetZOffset)**2) <= sim_targetRadius) then
                spec = TARG_SPEC
             end if
             
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case default 
          call Driver_abortFlash("[Simulation_initBlock] Unsupported geometry")
       end select

    case (3)

       select case (sim_geometry)
       case (CARTESIAN)
          if (sim_targetGeom == "sphere") then
             if(sqrt(x**2 + y**2 + z**2) <= sim_targetRadius) &
                  spec = TARG_SPEC
          elseif (sim_targetGeom == "cylinder") then
             if(sqrt(x**2 + y**2) < sim_targetRadius .and. &
                  abs(z) <= sim_targetHeight ) &
                  spec = TARG_SPEC
          else
             call Driver_abortFlash("[Simulation_initBlock] Unsupported sim_targetGeom")
          end if

       case default 
          call Driver_abortFlash("[Simulation_initBlock] Unsupported geometry")
       end select

    end select
           
  end subroutine findSpec


end subroutine Simulation_initBlock
