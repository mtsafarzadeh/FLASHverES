!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_computeDt
!!
!! NAME
!!  
!!  Diffuse_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Diffuse_computeDt ( integer(IN) : blockID, 
!!                  real(IN):  xCenter(:), 
!!                  real(IN):  xLeft(:), 
!!                  real(IN):  xRight(:), 
!!                  real(IN): dx(:), 
!!                  real(IN): uxgrid(:),
!!                  real(IN):  yCenter(:), 
!!                  real(IN):  yLeft(:), 
!!                  real(IN):  yRight(:), 
!!                  real(IN): dy(:), 
!!                  real(IN): uygrid(:), 
!!                  real(IN):  zCenter(:), 
!!                  real(IN):  zLeft(:), 
!!                  real(IN):  zRight(:), 
!!                  real(IN): dz(:), 
!!                  real(IN): uzgrid(:), 
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for diffusion source term solver.
!! 
!!  The current implementation may be very conservative, especially with
!!  respect to the viscosity term.  Users may want to change the implementation
!!  to be less conservative, and/or tweak the time step by tweaking the
!!  dt_diff_factor runtime parameter.
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  xCenter         X coordinates at the center of the cell
!!  xLeft           X coordinates at the left edge of the cell
!!  xRight          X coordinates at the right edge of the cell
!!  yCenter         Y coordinates at the center of the cell
!!  yLeft           Y coordinates at the left edge of the cell
!!  yRight          Y coordinates at the right edge of the cell
!!  zCenter         Z coordinates at the center of the cell
!!  zLeft           Z coordinates at the left edge of the cell
!!  zRight          Z coordinates at the right edge of the cell
!!  d*              deltas in each {*=x, y z} directions
!!  u*grid          velocity of grid expansion in {*=x, y z} directions
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***


subroutine Diffuse_computeDt (blockID, &
                              xCenter,xLeft,xRight, dx, uxgrid, &
                              yCenter,yLeft,yRight, dy, uygrid, &
                              zCenter,zLeft,zRight, dz, uzgrid, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )


  use Diffuse_data, ONLY : dt_diff_factor, diffusion_cutoff_density, dirGeom, &
                           useDiffuseSpecies,useDiffuseVisc,useDiffuseTherm,  &
                           useDiffuseMagneticResistivity,useDiffuse,          &
                           diff_useEleCond, &
                           diff_meshMe


  use Driver_interface, ONLY : Driver_abortFlash
  use Conductivity_interface, ONLY : Conductivity
  use MassDiffusivity_interface, ONLY : MassDiffusivity
  use Viscosity_interface, ONLY : Viscosity
  use MagneticResistivity_interface, ONLY : MagneticResistivity 
  use Opacity_interface, ONLY: Opacity
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer, dimension(:,:,:,:) :: solnData

  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: &
       xCenter,xLeft,xRight
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: &
       yCenter,yLeft,yRight
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) ::&
        zCenter,zLeft,zRight
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: &
       dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: &
       dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: &
       dz, uzgrid


  real, PARAMETER :: TINY = 1.e-30

  integer :: i, j, k, n

  integer :: temploc(5)

  real, dimension(MDIM) :: del
  real    ::  dt_temp, dt_ltemp
  ! storage for the 1d mass fractions in a zone
  real :: massfrac(NSPECIES), cond_zone_unusedHere, xtemp, xdens, diff_coeff, &
          visc_zone_unusedHere, viscKinematic, magResist,&
          mass_diff_zone, Raddiff_coeff
#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,&
                 GRID_KLO_GC:GRID_KHI_GC)::max_diffusivity
#else
  real,allocatable::max_diffusivity(:,:,:)
#endif
  real :: deltax,deltay,deltaz
  integer :: isize,jsize,ksize

  if(.not.useDiffuse) return !! Simply return if diffusion is turned off

  dt_temp    = huge(1.0)
  temploc(:) = 0
  
  deltax = dx(1)
  deltay = dy(1)
  deltaz = dz(1)
     
  isize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jsize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  ksize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     
#ifndef FIXEDBLOCKSIZE
  
  allocate(max_diffusivity(isize,jsize,ksize))
#endif
!-----------------------------------------------------------------------------
! start by getting the maximum of the thermal, mass, and viscous diffusivity
! -- the largest value is what will set the timestep for a zone.
!-----------------------------------------------------------------------------


  do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          
           ! load the mass fractions
           do n = 1, NSPECIES
              massfrac(n) = solnData(SPECIES_BEGIN-1+n,i,j,k)
           enddo
           xdens = solnData(DENS_VAR,i,j,k)
           xtemp = solnData(TEMP_VAR,i,j,k)

           diff_coeff    = 0.0 

#ifdef TELE_VAR
           if (diff_useEleCond) then
               call Conductivity(solnData(:,i,j,k), diffCoeff=diff_coeff, component=2)
           endif
#else
           xtemp = solnData(TEMP_VAR,i,j,k)
           call Conductivity(solnData(:,i,j,k),cond_zone_unusedHere, diffCoeff=diff_coeff, component=2)
#endif
           max_diffusivity(i,j,k) = TINY        
              
! diffusion breaks down at low densities, when we free stream.  We use a 
! cutoff density here to turn off the diffusion (if desired).  These zones
! should not affect the timestep, so leave the diffusion coefficient really 
! small in this case
              if (diff_coeff > max_diffusivity(i,j,k)) then
                 if (xdens > diffusion_cutoff_density) then
                    max_diffusivity(i,j,k) = diff_coeff
                 endif
              endif
              
           
           if (useDiffuseVisc) then
              call Viscosity(xtemp, xdens, massfrac, visc_zone_unusedHere, viscKinematic)
              
!!$              if (viscKinematic > max_diffusivity(i,j,k)) then
!!$                 max_diffusivity(i,j,k) = viscKinematic * 2.0
!!$              endif
              max_diffusivity(i,j,k) = max_diffusivity(i,j,k) + viscKinematic * ((8*NDIM-4)/3.0)
              
           endif
           
           
           if (useDiffuseSpecies) then
              
              call MassDiffusivity(xtemp, xdens, massfrac, mass_diff_zone)
              
              if (mass_diff_zone > max_diffusivity(i,j,k)) then
                 max_diffusivity(i,j,k) = mass_diff_zone
              endif
              
           endif


           if (useDiffuseMagneticResistivity) then
              call MagneticResistivity(xtemp,xdens,massfrac,magResist)

              if (magResist > max_diffusivity(i,j,k)) then
                 max_diffusivity(i,j,k) = magResist
              endif

           end if
           
        enddo
     enddo
  enddo


!--------------------------------------------------------------------------
! 1-dimensional
!--------------------------------------------------------------------------
     
  if (NDIM == 1) then
     
     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        
        dt_ltemp = 0.5*deltax**2/max_diffusivity(i,1,1)
        
        if (dt_ltemp < dt_temp) then
           dt_temp    = dt_ltemp
           temploc(1) = i
           temploc(2) = 1
           temploc(3) = 1
           temploc(4) = blockID
           temploc(5) = diff_meshMe
        endif
        
     enddo
     
     
!--------------------------------------------------------------------------
! 2-dimensional
!--------------------------------------------------------------------------

! angular coordinates may be present

  elseif (NDIM == 2) then
     
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        
        if (dirGeom(JAXIS) == XYZ) then ! Cartesian 
           
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              dt_ltemp = 0.5*min&
                   (deltax**2,deltay**2)/max_diffusivity(i,j,1)
              
              if (dt_ltemp < dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = blockID
                 temploc(5) = diff_meshMe
              endif
           enddo
           
        else ! Angular
           
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              deltay = (xCenter(i) * &
                   (yRight(j) - yLeft(j)))
              
              dt_ltemp = 0.5*min(deltax**2,deltay**2)/max_diffusivity(i,j,1)
              
              if (dt_ltemp < dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = blockID
                 temploc(5) = diff_meshMe
              endif
           enddo
        endif
     enddo
        

!--------------------------------------------------------------------------
! 3-dimensional
!--------------------------------------------------------------------------

     ! here we only allow Cartesian and spherical right now
     
  else
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        
        if (dirGeom(IAXIS) == XYZ) then    ! Cartesian
           
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 dt_ltemp = 0.5*min&
                      (deltax**2,deltay**2,deltaz**2)/&
                      max_diffusivity(i,j,k)
                 
                 if (dt_ltemp < dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = blockID
                    temploc(5) = diff_meshMe
                 endif
              enddo
           enddo
        else
           
           ! spherical
           
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

! dy and dz are dependent on the radius
                 deltaz = (xCenter(i) *sin(yCenter(j)) * &
                      (zRight(k) - zLeft(k)))
                 
                 deltay = (xCenter(i) * &
                      (yRight(j) - yLeft(j)))
                 
                 dt_ltemp = 0.5*min&
                      (deltax**2,deltay**2,deltaz**2)/&
                      max_diffusivity(i,j,k)
                 
                 if (dt_ltemp < dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = blockID
                    temploc(5) = diff_meshMe
                 endif
              enddo
           enddo
        endif
     enddo
  endif

!-----------------------------------------------------------------------------
! if the diffusion timestep for this block is the smallest, then store it
!-----------------------------------------------------------------------------

  dt_temp = dt_diff_factor * dt_temp
  
  if (dt_temp < dt_check) then
     dt_check = dt_temp
     dt_minloc = temploc
  endif
  

#ifndef FIXEDBLOCKSIZE
  deallocate(max_diffusivity)
#endif

  if(dt_check <= 0.0) call Driver_abortFlash("[Diffuse]: computed dt is not positive! Aborting!")

  

  return
end subroutine Diffuse_computeDt


