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
      use Driver_interface, ONLY: Driver_getDt
      use Driver_data, ONLY: dr_nstep
      use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getBlkPtr, &
        Grid_releaseBlkPtr
      use Eos_interface, ONLY : Eos_wrapped

      implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

      integer, intent(IN) :: blockID

      real :: radia,sdot

      real, pointer :: solnData(:,:,:,:)

!     These are counters
      integer :: i, j, k

!     This is the time step
      real :: dt

!!    These are for measuing the limits of the block
      integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

!!    These are the limits of the block
      integer  :: lowX,lowY,lowZ,highX,highY,highZ

!!    A bunch of local variables
      real :: tmp, rho, den_H, den_e, ei, eiold, ek, ini_t, ini_d

! Mark zones that suffer cooling
      logical :: rad_zone

      external radloss
!
!==============================================================================
!
! get the current timestep
      call Driver_getDt(dt)

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

      do k = lowZ,highZ
         do j = lowY,highY
            do i = lowX,highX
!
               tmp = solnData(TEMP_VAR,i,j,k)
               rho = solnData(DENS_VAR,i,j,k)

! derive the Hydrogen and the electron number density

               den_H = rho*massfracH/amu
               den_e = rho*(massfracH+(2.*(1.-massfracH)/4.))/amu 

               sdot = 0.0e0

! if the temperature is between the limits and
! the density is between the limits

!                  print *, radia, tmp, den_H

               if (dr_nstep.eq.1) solnData(INIT_VAR,i,j,k) = tmp
               if (dr_nstep.eq.1) solnData(INID_VAR,i,j,k) = rho

	       ini_t = solnData(INIT_VAR,i,j,k) 
               ini_d = solnData(INID_VAR,i,j,k)

               if ( (tmp >= tradmin) .AND. (tmp <= tradmax) .AND.         &
     &        (den_H*amu >= dradmin) .AND. (den_H*amu <= dradmax) .AND. &
     &        (tmp.gt.ini_t).and.(rho.lt.ini_d) .and.(dr_nstep.gt.1) ) then

                  rad_zone = .TRUE.

! radiative losses from an optically thin plasma
                  call radloss(tmp,radia)

                  sdot = -(den_H*den_e*radia)/rho


! change in internal energy due to radiative losses
                  ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                              solnData(VELY_VAR,i,j,k)**2 +  &
                              solnData(VELZ_VAR,i,j,k)**2)

                  ei = solnData(EINT_VAR,i ,j,k)
! Set nocool to zero to turn of cooling, but keep the variable around for plotting

!                  if (dr_nstep.gt.1) ei = ei + dt*sdot*nocool - dt*solnData(RADI_VAR,i,j,k)*rho


                   ei = ei + dt*sdot*nocool

! update the global thermodynamic quantities due to the radiative losses
                  solnData(ENER_VAR,i,j,k) = ei + ek
                  solnData(EINT_VAR,i,j,k) = ei
                  solnData(RADI_VAR,i,j,k) = sdot*nocool

! store the initial radiative losses rate -- need it for the equilibrium heating rate
!                  solnData(RADI_VAR,i,j,k) = sdot


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
