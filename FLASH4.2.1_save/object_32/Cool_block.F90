!
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
!     These are the cooling rates due to the two dust mechanisms
      real :: radia_rec,sdot_rec,radia_rad,sdot_rad
!     This the the Compton HEATING rate
      real :: radia_compton,sdot_compton
!     This is the size of the dust in microns
      real :: a
!     This is xstar as defined in Montier and Giard 2004
      real :: xstar
!     This is ionization rate divided by the recobination rate
      real :: Xi
!     This is ionization fraction a estimated from only phoionization
      real :: xe

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
      real :: tmp, rho, den_B, den_e, den_H, den_tot,dustfrac
      real :: tmp6, deltagam, deltadust, ei, eiold, ek,pres
      real :: center

      real :: indexfloat,indexrhofloat,remain,remainrho
      integer :: index,  indexrho

!! Here we allocate arrays with the different coordinates
      real, allocatable, dimension(:)         :: xCoord, yCoord, zCoord!    A bunch of local variables
      real :: time
      real :: enhance

      logical :: getGuardCells

!     factors used for attenuating the background due to self-sheilding 
      real    :: delta, gammafactor

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
               dustfrac = solnData(DUST_MSCALAR,i,j,k)
               rho = solnData(DENS_VAR,i,j,k)
               pres = solnData(PRES_VAR,i,j,k)

               den_H   = rho*(   massfracH  +0.*(1.-massfracH)/4.)/amu 
               den_e   = rho*(   massfracH  +2.*(1.-massfracH)/4.)/amu 
               den_B   = rho*(   massfracH  +4.*(1.-massfracH)/4.)/amu 
               den_tot = rho*(2.*massfracH  +3.*(1.-massfracH)/4.)/amu 

               sdot = 0.0e0
               sdot_total = 0.0e0

! if the temperature is between the limits and
! the density is between the limits

               if ( (tmp >= tradmin .AND. tmp <= tradmax) .AND.  &         
                  (rho >= dradmin) .AND. (rho <= dradmax)) then
                  rad_zone = .TRUE.

! radiative losses from an optically thin plasma
                  radia = 0.
                  if(tmp>=tradmin) then
                    indexfloat = 160.*(log(tmp)/log(10.)-1.)/8.+1
                    index = indexfloat
                    remain = indexfloat-index
                    if(index>=161) then 
                      index  = 161
                      remain = 1.
                    endif
                    indexrhofloat = log10(den_B/60.)*20.+61.

!  Here's where we put in the attenuation and self-sheiding
!  estimat ionized fractio 
                    if((tmp.gt.1E5)) then 
                      xe = 1.
                    else  
                      Xi = Gammapi/(6.2E-10*tmp**(-0.845)*den_H)
                      xe = 0.5*(sqrt(Xi*Xi+4.*Xi)-Xi)             
                      if(Xi .ge. 100) xe = 1.-1./Xi
                    endif
                    solnData(IONI_VAR,i,j,k) = xe
                    gammafactor = 1.
                    if((tmp.lt.1E5).and.(selfshield)) then
                      deltagam = den_H*(1.-xe)*6.3E-18*shieldlength
                      gammafactor = 0.98*(1+deltagam**1.64)**(-2.28) &
                                   +0.02*(1.+deltagam)**(-0.84)
                    endif
!   THe table is dependent on gamma/rho so we just divide by gammafactor
                    indexrhofloat = log10(den_B/60./gammafactor)*20.+61.
                    indexrho = indexrhofloat
                    remainrho = indexrhofloat-indexrho
                    if(indexrho>=240) then 
                      indexrho  = 240
                      remainrho = 1.
                    endif
                    if(indexrho<=0) then 
                      indexrho = 0.
                      remainrho = 0.
                    endif
                    radia =  total_cooling(indexrho,  index)  *(1.-remainrho)*(1.-remain)&
                            +total_cooling(indexrho,  index+1)*(1.-remainrho)*(remain)   &     
                            +total_cooling(indexrho+1,index)  *(remainrho)   *(1.-remain)&
                            +total_cooling(indexrho+1,index+1)*(remainrho)   *(remain)     
                    radia = radia &
                            -total_heating(indexrho,  index)  *(1.-remainrho)*(1.-remain)&
                            -total_heating(indexrho,  index+1)*(1.-remainrho)*(remain)   &     
                            -total_heating(indexrho+1,index)  *(remainrho)   *(1.-remain)&
                            -total_heating(indexrho+1,index+1)*(remainrho)   *(remain)     
!                    if(radia .le. 0) radia = 0.
                    radia = radia*den_B*den_B

!   This one is dust recombination cooling
                    delta = 0.74*tmp**(-0.068)
                    radia_rec = 4.65E-30*tmp**0.94*(G0*sqrt(tmp)/den_tot)**delta
                    radia_rec = radia_rec*dustfrac*den_tot*den_tot

!   This is on the dust collsional heating
                    a=A0*dustfrac**0.3333
! A0 is in units of microns so we convert it to cm here
                    a=a*1E-4
                    xstar = a**0.6667*(1.26E11/tmp)
                    if (xstar .ge. 4.5) then
                      radia_rad =  5.38E-10*a*a*tmp**1.5
                      else
                      if(xstar .ge. 1.5) then
                        radia_rad = 1.47E-3*a**2.41*tmp**0.88   
                      else
                        radia_rad = 6.48E6*a**3
                      endif
                    endif 
                    if(dustfrac .lt. 1E-6) radia_rad = 0.
! Solar is 0.01 by mass we are assuming solar and taking the
! dustinit ne 1 to be representative in a decrease in grain size and not in
! number density
                    radia_rad = radia_rad*den_e*(0.01*rho/mdust0)

! finally we add the compton heating/cooling
                    radia_compton = den_e*xe*6.6525E-25*Comptonflux*Comptontemp/511. 
                    radia_compton = radia_compton*(1.-tmp/ComptontempK)
                  endif
                  sdot_rec = -(radia_rec)/rho
                  solnData(RADR_VAR,i,j,k) = sdot_rec 
                  sdot_rad = -(radia_rad)/rho
                  solnData(RADD_VAR,i,j,k) = sdot_rad
! note there is no minus sign here
                  sdot_compton = (radia_compton)/rho
                  solnData(RADC_VAR,i,j,k) = sdot_compton
                  sdot     = -(radia)/rho
                  solnData(RADI_VAR,i,j,k) = sdot
                  sdot     = radia_compton/rho-(radia+radia_rec+radia_rad)/rho
 
! DONT LET COOLING ACTUALLY HEAT
                  if(sdot .ge. sdot_compton) sdot = sdot_compton 
 
                  

! change in internal energy due to radiative losses
                  ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                              solnData(VELY_VAR,i,j,k)**2 +  &
                              solnData(VELZ_VAR,i,j,k)**2)

                  ei = solnData(EINT_VAR,i ,j,k)

! New stuff: This will loop over time to make sure we don't cool too fast and get wrong answers

                 dtprime = dt
                 dtmax = 0.1*ei/abs(sdot*nocool+1.E-30)

                 do while(dtprime>0.) 
                     if(dtprime < dtmax) then
                       ei = ei + dtprime*sdot*nocool
                       dtprime = 0.
                       sdot_total = sdot_total + dtprime*sdot*nocool
                     else
                        ei = ei*(1.+0.1*sdot/abs(sdot))
                        tmp =tmp*(1.+0.1*sdot/abs(sdot))
                        sdot_total = sdot_total + ei*0.1*sdot/abs(sdot)
                        dtprime = dtprime-dtmax
                        radia = 0.
                        if(tmp>=tradmin) then
                         indexfloat = 160.*(log(tmp)/log(10.)-1.)/8.+1
                         index = indexfloat
                         remain = indexfloat-index
                         if(index>=161) then
                           index  = 161
                           remain = 1.
                         endif
                         indexrhofloat = log10(den_B/60.)*20.+61.
!  Here's where we put in the attenuation and self-sheiding
!  estimate ionized fractio 
                         if((tmp.gt.1E5)) then
                           xe = 1.
                         else
                           Xi = Gammapi/(6.2E-10*tmp**(-0.845)*den_H)
                           xe = 0.5*(sqrt(Xi*Xi+4.*Xi)-Xi)
                           if(Xi .ge. 100) xe = 1.-1./Xi
                         endif
                         gammafactor = 1.
                         if((tmp.lt.1E5).and.(selfshield)) then
                           deltagam = den_H*(1.-xe)*6.3E-18*shieldlength
                           gammafactor = 0.98*(1+deltagam**1.64)**(-2.28) &
                                       +0.02*(1.+deltagam)**(-0.84)
!                           print *,'Yeah',rho,tmp,gammafactor
                         endif
!   THe table is dependent on gamma/rho so we just divide by gammafactor
                        indexrhofloat = log10(den_B/60./gammafactor)*20.+61.

                         indexrho = indexrhofloat
                         remainrho = indexrhofloat-indexrho
                         if(indexrho>=240) then
                           indexrho  = 240
                           remainrho = 1.
                         endif
                         if(indexrho<=0) then
                           indexrho = 0.
                           remainrho = 0.
                         endif
                         radia =  total_cooling(indexrho,  index) *(1.-remainrho)*(1.-remain)&
                                 +total_cooling(indexrho, index+1)*(1.-remainrho)*(remain)   &
                                 +total_cooling(indexrho+1,index)  *(remainrho) *(1.-remain)&
                                 +total_cooling(indexrho+1,index+1)*(remainrho) *(remain)
                         radia = radia  &
                                 -total_heating(indexrho,  index) *(1.-remainrho)*(1.-remain)&
                                 -total_heating(indexrho, index+1)*(1.-remainrho)*(remain)   &
                                 -total_heating(indexrho+1,index)  *(remainrho) *(1.-remain)&
                                 -total_heating(indexrho+1,index+1)*(remainrho) *(remain)
!                         if(radia .le. 0) radia = 0.
                         radia = radia*den_B*den_B

!   This one is dust recombination cooling
                           delta = 0.74*tmp**(-0.068)
                           radia_rec = 4.65E-30*tmp**0.94*(G0*sqrt(tmp)/den_tot)**delta
                           radia_rec = radia_rec*dustfrac*den_tot*den_tot

!   This is on the dust collsional heating
                           a=A0*dustfrac**0.3333
                           a=a*1E-4
                           xstar = a**0.6667*(1.26E11/tmp)
                           if (xstar .ge. 4.5) then
                             radia_rad =  5.38E-10*a*a*tmp**1.5
                             else
                             if(xstar .ge. 1.5) then
                               radia_rad = 1.47E-3*a**2.41*tmp**0.88
                             else
                               radia_rad = 6.48E6*a**3
                             endif
                           endif
                           if(dustfrac .lt. 1E-6) radia_rad = 0.
                           radia_rad = radia_rad*den_e*(0.01*rho/mdust0)

! finally we add the compton heating/cooling
                           radia_compton = xe*den_e*6.6525E-25*Comptonflux*Comptontemp/511.
                           radia_compton = radia_compton*(1.-tmp/ComptontempK)

                           radia = radia+radia_rec+radia_rad-radia_compton
                         endif
                         sdot = -(radia)/rho
                         sdot_compton = radia_compton/rho
                      if(sdot .ge. sdot_compton) sdot = sdot_compton
                      dtmax = 0.1*ei/abs(sdot*nocool+1E-30)
                 endif
                 enddo

! update the global thermodynamic quantities due to the radiative losses
                 solnData(ENER_VAR,i,j,k) = ei + ek
                 solnData(EINT_VAR,i,j,k) = ei

! here we add the dust destruction term
                  if(dustfrac.gt.0) then
                    tmp6      = tmp/1E6
                    deltadust = -dt*den_tot*(Adust/A0)*(tmp6)**(-.25)*exp(-Bdust/sqrt(tmp6)) 
                    if(-deltadust.le.dustfrac*(1./3.)*0.001) then
                      dustfrac = dustfrac + 3.*dustfrac**(2./3)*deltadust
                   else
                     dustfrac = (dustfrac**(1./3.)+deltadust)**3
                   endif
                   if(dustfrac.le.0) dustfrac = 0
                   solnData(DUST_MSCALAR,i,j,k) = dustfrac 
                   solndata(DELD_VAR,i,j,k) = deltadust/dt*3.*dustfrac**(2./3)
                 endif
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
