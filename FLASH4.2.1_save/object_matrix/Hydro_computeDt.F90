!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Hydro_computeDt(integer(IN)    :: blockID,
!!                  real(IN)       :: x(GRID_ILO_GC:GRID_IHI_GC),
!!                  real(IN)       :: dx(GRID_ILO_GC:GRID_IHI_GC),
!!                  real(IN)       :: uxgrid(GRID_ILO_GC:GRID_IHI_GC),
!!                  real(IN)       :: y(GRID_JLO_GC:GRID_JHI_GC),
!!                  real(IN)       :: dy(GRID_JLO_GC:GRID_JHI_GC),
!!                  real(IN)       :: uygrid(GRID_JLO_GC:GRID_JHI_GC),
!!                  real(IN)       :: z(GRID_KLO_GC:GRID_KHI_GC),
!!                  real(IN)       :: dz(GRID_KLO_GC:GRID_KHI_GC),
!!                  real(IN)       :: uzgrid(GRID_KLO_GC:GRID_KHI_GC),
!!                  integer(IN)    :: blkLimits(2,MDIM)
!!                  integer(IN)    :: blkLimitsGC(2,MDIM)
!!                  real,pointer   :: U(:,:,:,:),
!!                  real(INOUT)    :: dtCheck,
!!                  integer(INOUT) :: dtMinLoc(5),
!!                  real(INOUT), optional :: extraInfo)
!!  
!!
!! DESCRIPTION
!!
!!  This routine computes the timestep limiter for the Unsplit Hydro solver.
!!  The Courant-Fredrichs-Lewy criterion is used.  The sound
!!  speed is computed and together with the velocities, is used to constrain
!!  the timestep such that no information can propagate more than one zone
!!  per timestep.
!!  Note that this routine only accounts for computing advection time step in hyperbolic
!!  system of equations.
!!
!! ARGUMENTS
!!
!!  blockID       -  local block ID
!!  x, y, z       -  three, directional coordinates
!!  dx,dy,dz      -  distances in each {*=x, y z} directions
!!  uxgrid        -  velocity of grid expansion in x directions
!!  uygrid        -  velocity of grid expansion in y directions
!!  uzgrid        -  velocity of grid expansion in z directions
!!  blkLimits     -  the indices for the interior endpoints of the block
!!  blkLimitsGC   -  the indices for endpoints including the guardcells
!!  U             -  the physical, solution data from grid
!!  dtCheck       -  variable to hold timestep constraint
!!  dtMinLoc(5)   -  array to hold location of cell responsible for minimum dt:
!!                   dtMinLoc(1) = i index
!!                   dtMinLoc(2) = j index
!!                   dtMinLoc(3) = k index
!!                   dtMinLoc(4) = blockID
!!                   dtMinLoc(5) = hy_meshMe
!!  extraInfo     -  Driver_computeDt can provide extra info to the caller
!!                   using this argument.
!!
!!***

!!REORDER(4): U

Subroutine Hydro_computeDt( blockID,       &
                            x, dx, uxgrid, &
                            y, dy, uygrid, &
                            z, dz, uzgrid, &
                            blkLimits, blkLimitsGC, &
                            U,  dtCheck, dtMinLoc,  &
                            extraInfo)


#include "Flash.h"
#include "constants.h"

  use Hydro_data
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

  !! Arguments type declaration ------------------------------------------
  integer, intent(IN) :: blockID 
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC

#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif

  real, pointer         :: U(:,:,:,:)
  real,   intent(INOUT) :: dtCheck
  integer,intent(INOUT) :: dtMinLoc(5)
  real, OPTIONAL,intent(INOUT) :: extraInfo
  !! ----------------------------------------------------------------------
  integer :: i, j, k, temploc(5)
  real    :: sndspd2, delxinv, delyinv, delzinv, dt_temp, dt_ltemp
  real    :: cfx2,cfy2,cfz2,bbx2,bby2,bbz2,b2

  !! Case 1: we exit this routine if not needed.
  if ((.not. hy_useHydro) .or. (.not. hy_updateHydroFluxes)) return

  !! Case 2: we simply pass the already computed dt information to Driver_computeDt.
  !!         In this case, hydro dt gets computed in either 
  !!          (i) hy_uhd_energyFix (hy_hydroComputeDtOption=0), or 
  !!         (ii) hy_uhd_getFaceFlux (hy_hydroComputeDtOption=1)
  if ((hy_hydroComputeDtOption .ne. -1) .and. &
      (.not. hy_hydroComputeDtFirstCall)) then
     if ( hy_dtmin < dtCheck ) then
        dtCheck  = hy_dtmin
        dtMinLoc = hy_dtminloc
     endif

     if (present(extraInfo)) then
        extraInfo = 0.
        if (hy_useVaryingCFL) then
           extraInfo = hy_cfl_original
           if (hy_cfl <= extraInfo) then
              extraInfo = hy_cfl
           endif
        else
           extrainfo = 0.
        endif
     endif

  else
     !! Case 3: hy_hydroComputeDtOption = -1
     !!         We perform the global loop to compute hydro dt in the old way.
     !!         This implementation provides the full geometry supports for computing dt.
     !! NOTE: we always call this global compute dt for the very first step.
     if (hy_hydroComputeDtFirstCall) hy_hydroComputeDtFirstCall = .false.
     dt_temp    = 0.
     temploc(:) = 0

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        !! Conversion for unitSystem if needed ---------------------------------------
        U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
        U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
        U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
        U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)/hy_vref
#if defined(FLASH_USM_MHD)
        U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)/hy_bref
#endif
        !! ---------------------------------------------------------------------------
     endif


     ! First, set default grid delta values in Cartesian;
     ! Otherwise, they will be reset for other geometries later.
     ! Case1: (hy_geometry == CARTESIAN) then
     delyinv = 1.
     delzinv = 1.

     delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
     if (NDIM > 1) &
     delyinv = 1.0/dy(blkLimits(LOW,JAXIS))
     if (NDIM > 2) &
     delzinv = 1.0/dz(blkLimits(LOW,KAXIS))

     if (hy_geometry == POLAR) & !Polar in 3D (that's a no no)
          call Driver_abortFlash("[Hydro_computeDt] ERROR: Polar geometry not supported in 3D")

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef BDRY_VAR /* Do not compute time step dt when in the solid boundary cells */
              if (U(BDRY_VAR,i,j,k) < 0.0) then
#endif
                 sndspd2 = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 cfx2    = sndspd2
                 cfy2    = sndspd2
                 cfz2    = sndspd2

#if defined(FLASH_USM_MHD) /*compute additional magneto-acoustic speeds for MHD */
                 bbx2 = U(MAGX_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 bby2 = U(MAGY_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 bbz2 = U(MAGZ_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 b2   = bbx2 + bby2 + bbz2
                 sndspd2= U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)

                 cfx2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bbx2))
                 cfy2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bby2))
                 cfz2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bbz2))
#endif

                 ! For other geometry supports
                 if (hy_geometry == CYLINDRICAL) then
#if NDIM > 2
                    delzinv = 1.0/(x(i)*dz(k))           ! z is phi
#endif
                 elseif (hy_geometry == SPHERICAL) then
#if NDIM > 1
                    delyinv = 1.0/(x(i)*dy(j))           ! y is theta
                    delzinv = 1.0/(x(i)*sin(y(j))*dz(k)) ! z is phi
#endif
                 endif


                 dt_ltemp = (abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(cfx2))*delxinv
                 if (NDIM > 1) dt_ltemp = max(dt_ltemp,(abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(cfy2))*delyinv)
                 if (NDIM > 2) dt_ltemp = max(dt_ltemp,(abs(U(VELZ_VAR,i,j,k)-uzgrid(j))+sqrt(cfz2))*delzinv)

                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = blockID
                    temploc(5) = hy_meshMe
                 endif
#ifdef BDRY_VAR
              endif
#endif
           enddo
        enddo
     enddo



     dt_temp = hy_cfl / dt_temp
     if (dt_temp < dtCheck) then
        dtCheck = dt_temp
        dtMinLoc = temploc
     endif

 
     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        !! Conversion for unitSystem if needed ---------------------------------------
        U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)*hy_dref
        U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)*hy_eref
        U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)*hy_pref
        U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)*hy_vref
#if defined(FLASH_USM_MHD)
        U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)*hy_bref
#endif
        !! ---------------------------------------------------------------------------
     endif

  endif ! end if of if (hy_hydroComputeDtOption .ne. -1) then



  !! For the purpose of having screen output for varying CFL
  if (present(extraInfo)) then
     extraInfo = 0.
     if (hy_useVaryingCFL) then
        extraInfo = hy_cfl_original
        if (hy_cfl <= extraInfo) then
           extraInfo = hy_cfl
        endif
     else
        extrainfo = 0.
     endif
  endif

  if(dtCheck <= 0.0) then
     print*,dtCheck
     call Driver_abortFlash("[Hydro]: Computed dt is not positive! Aborting!")
  endif
  return

End Subroutine Hydro_computeDt


