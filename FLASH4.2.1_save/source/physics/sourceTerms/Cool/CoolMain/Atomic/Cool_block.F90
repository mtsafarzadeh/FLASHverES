!!****f* source/source_terms/Cool/CoolMain/Cool_block
!!
!! NAME
!!  
!!  Cool_block
!!
!!
!! SYNOPSIS
!! 
!!  call Cool_block(blockID)
!!  
!! DESCRIPTION
!!
!!  Apply the radloss source term operator to a block of zones.
!!  The radiative losses rate is used to update the
!!  internal energy in the zone. The radiative losses from an 
!!  optically thin plasma are from Sarazin (1986) Rev.Mod.Phys.
!!  and from Raymond (1976).
!!  (see also Peres et al. 1982, ApJ 252, 791)
!!
!!  After we call radloss, call the eos to update the pressure 
!!  and temperature based on the radiative losses.
!!
!!
!!***

 subroutine Cool_block (blockID)
      use Cool_data
!      use Heat_data, ONLY: rbubble,taububble,zbubble,sphericalBomb
      use Driver_interface, ONLY: Driver_getDt,Driver_getSimTime
      use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getBlkPtr, &
        Grid_releaseBlkPtr
      use Eos_interface, ONLY : Eos_wrapped

      implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

      integer, intent(IN) :: blockID

      real :: radia,sdot,sdot_total,cool_out

      real, pointer :: solnData(:,:,:,:)

!     These are counters
      integer :: i, j, k

!!    A flag that gets returned when you have memory problems
      integer :: istat

!     This is the time step
      real :: dt
      real :: dtprime
      real :: dtmax

!!    These are for measuing the limits of the block
      integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

!!    These are the limits of the block
      integer  :: lowX,lowY,lowZ,highX,highY,highZ

!!    A bunch of local variables
      real :: tmp, rho, den_H, den_e, ei, eiold, ek,pres
      real :: center
      real :: metal 

      real :: indexfloat
      integer :: index

!! Here we allocate arrays with the different coordinates
      real, allocatable, dimension(:)         :: xCoord, yCoord, zCoord!    A bunch of local variables
      real :: time
      real :: enhance

      logical :: getGuardCells

! Mark zones that suffer cooling
      logical :: rad_zone

      external radloss
!
!==============================================================================
!
! get the current timestep
      call Driver_getDt(dt)
      call Driver_getSimTime(time)

      call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
      call Grid_getBlkPtr(blockID,solnData)


! get the array limits
     lowX  = blkLimits(LOW,IAXIS)
     lowY  = blkLimits(LOW,JAXIS)
     lowZ  = blkLimits(LOW,KAXIS)
     highX = blkLimits(HIGH,IAXIS)
     highY = blkLimits(HIGH,JAXIS)
     highZ = blkLimits(HIGH,KAXIS)

      rad_zone = .FALSE.

! sweep over all the zones

     allocate(xCoord(highX),stat=istat)
     allocate(yCoord(highY),stat=istat)
     allocate(zCoord(highZ),stat=istat)

     getGuardCells = .true.
     call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,highX)
     call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,highY)
     call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,highZ)

! sweep over all the zones
      do k = lowZ,highZ
         do j = lowY,highY
            do i = lowX,highX
               tmp = solnData(TEMP_VAR,i,j,k)
               rho = solnData(DENS_VAR,i,j,k)
               pres = solnData(PRES_VAR,i,j,k)


               den_H = rho*massfracH/amu
               den_e = rho*(massfracH+(2.*(1.-massfracH)/4.))/amu 

               sdot = 0.0e0
               sdot_total = 0.0e0

! if the temperature is between the limits and
! the density is between the limits

!	print *, tmp, den_H*amu,tradmin, tradmax,dradmin,dradmax

               if ( (tmp >= tradmin .AND. tmp <= tradmax) .AND.  &         
                  (rho >= dradmin) .AND. (rho <= dradmax)) then
!                  .AND. (ycoord(j) >= 4.*3.08d21 )) then 
                  rad_zone = .TRUE.

! radiative losses from an optically thin plasma
                  radia = 0.
                  if(tmp>=tradmin) then
                    metal = 1.0 !0.3 !50.*(solnData(METL_MSCALAR,i,j,k)+solnData(MTSN_MSCALAR,i,j,k))
                    indexfloat = 353.*(log(tmp)/log(10.)-2.)/7.
                    index = indexfloat
                    if(index>=352) then 
                      index = 352
                    endif
                    radia = hhe_cooling(index)+metal_cooling(index)*metal

                  endif

                  sdot = -(den_H*den_e*radia)/rho
              
                  solnData(RADI_VAR,i,j,k) = sdot

!                  print *, radia, tmp, den_H,index

! change in internal energy due to radiative losses
                  ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                              solnData(VELY_VAR,i,j,k)**2 +  &
                              solnData(VELZ_VAR,i,j,k)**2)

                  ei = solnData(EINT_VAR,i ,j,k)

! ES This is probably wrong and so commented out for the moment

! Set nocool to zero to turn of cooling, but keep the variable around for plotting
!                  eiold = ei
!                  ei = ei + dt*sdot*nocool
!                  if(ei .le. (eiold*0.9)) then
!                         if(eiold .ge. 0) then
!			 ei = eiold*0.9
!                   endif
!
!                  endif



! New stuff: This will loop over time to make sure we don't cool too fast and get wrong answers

		  dtprime = dt
		  dtmax = 0.1*ei/(-sdot*nocool+1.E-30)

		 do while(dtprime>0.) 
		     if(dtprime < dtmax) then
		     	ei = ei + dtprime*sdot*nocool
		        dtprime = 0.
                        sdot_total = sdot_total + dtprime*sdot*nocool
!		     endif
		     else
			ei = 0.9*ei
			tmp = 0.9*tmp
                        sdot_total = sdot_total - ei*0.1
			dtprime = dtprime-dtmax
			den_H = rho*massfracH/amu
               		den_e = rho*(massfracH+(2.*(1.-massfracH)/4.))/amu
			
			radia = 0.
                        if(tmp>=tradmin) then
                           indexfloat = 353.*(log(tmp)/log(10.)-2.)/7.
                           index = indexfloat
                           if(index>=352) then 
                              index = 352
                           endif
                           radia = hhe_cooling(index)+metal_cooling(index)*metal
!                           radia = radia+molecular(tmp)
                         endif
                      sdot = -(den_H*den_e*radia)/rho

                      dtmax = 0.1*ei/(-sdot*nocool+1E-30)
		  endif
		 enddo



! update the global thermodynamic quantities due to the radiative losses
                  solnData(ENER_VAR,i,j,k) = ei + ek
                  solnData(EINT_VAR,i,j,k) = ei

! store the radiative losses rate -- we store it along with the
! solution variables since it is useful to be able to refine on enuc.



               endif

            enddo
         enddo
      enddo

! if we called the RADLOSS on any zones in this block, then crank the
! eos out on the entire block
      if (rad_zone) then 

        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)

      end if

     call Grid_releaseBlkPtr(blockID,solnData)
     return
end subroutine Cool_block
