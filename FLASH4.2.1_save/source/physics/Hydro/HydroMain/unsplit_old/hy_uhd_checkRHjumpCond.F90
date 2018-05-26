!!****if* source/physics/Hydro/HydroMain/unsplit_old/hy_uhd_checkRHjumpCond
!!
!! NAME
!!
!!  hy_uhd_checkRHjumpCond
!!
!! SYNOPSIS
!!
!!  hy_uhd_checkRHjumpCond(real(IN)      :: V0
!!                         real(INOUT)   :: Wr(HY_VARINUMMAX),
!!                         real(INOUT)   :: Wl(HY_VARINUMMAX))
!!
!! ARGUMENTS
!!
!!  V0     - array containing primitive variables + gamc + game at cell centers
!!  Wr,Wl  - reconstructed right(r) and left(l) Riemann states
!!
!!
!! DESCRIPTION
!!
!!  behind shock |  ahead of shock
!!  (region 2)   |   (region 1)
!!               | 
!!     u2        |     u1 (= shock speed, that is the flow of gas comes to the shock with the shock speed)
!!   <-----      |  <------   
!!     rho2      |    rho1
!!     p2        |     p1
!!               |
!! (shocked gas) |  (undisturbed gas)
!!               |
!!             shock
!!
!!        shock frame (in this shock frame, shock speed = 0)
!!
!! Across the hydro shock, the Rankine-Hugoniot jump relations should satisfy:
!!
!!  (1) rho2/rho1 = [(gamc+1)*M1^2]   / [2+(gamc-1)*M1^2]
!!  (2) u2/u1     = [2+(gamc-1)*M1^2] / [(gamc+1)*M1^2]
!!  (3) p2/p1     = [2*gamc*M1^2 - (gamc-1)] / [gamc+1]
!!
!! Consequences:
!!  (1) M1>=1 (the shock speed (u1) exceeds the sound speed ahead of the shock
!!  (2) u2 <= soundSpeed in region 2 (subsonic behind the shock; supersonic ahead of it)
!!  (3) p2 >= p1 and rho2 >= rho1 (the shock is compressive)
!!  (4) u2 <= u1 and temp2 >= temp1  (the shock slows down the gas and heats it up)
!!  (5) 1<= rho2/rho1 < [gamc+1]/[gamc-1] (the maximum density ratio is [gamc+1]/[gamc-1];
!!                                         the pressure ratio increases with M1^2)
!! In the presence of magnetic fields, the relationship is lightly modified: letting X= rho2/rho1,
!!  (6) u2/u1 = 1/X
!!  (7) the shock is compressive with X>= 1
!!  (8) the effect of magnetic fields is to recude X below its hydro limit
!!  (9) the shock speed (u1) must exceeds the fast magneosonic speed sqrt(C1^2+V_alfven1^2), 
!!      where C1 is sound speed in region1
!!  (10) 1<B2/B1 < [gamc+1]/[gamc-1]
!!
!! REFERENCES
!!
!!  * Kirk, Melrose, Priest, Plasma Astrophysics
!!  * Toro
!!
!!***


Subroutine hy_uhd_checkRHjumpCond(dir,idx,idy,idz,V0,Vxr,Vxl,Vyr,Vyl,Vzr,Vzl,Wr,Wl,SWr,SWl)


  implicit none

#include "Flash.h"
#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer, intent(IN)  :: dir
  real, intent(IN)     :: idx,idy,idz
  real, dimension(HY_VARINUMMAX), intent(IN)    :: V0,Vxr,Vxl,Vyr,Vyl,Vzr,Vzl
  real, dimension(HY_VARINUMMAX), intent(INOUT) :: Wr,Wl
  logical, intent(OUT) :: SWr,SWl
  !!------------------------------------------------------------------------

  real :: Mach2, gammaRatio, densRatio,gamp,gamm
  real :: tinyD,bigD,tinyB,bigB,tinyV,bigV,tinyP,bigP
  integer :: hyVelEnd
  integer :: HY_VELN, HY_MAGN
  real :: divV

  SWl = .false.
  SWr = .false.


  gamp = V0(HY_GAMC)+1.
  gamm = V0(HY_GAMC)-1.
  gammaRatio = gamp/gamm

  if (NDIM == 1) then
     hyVelEnd=HY_VELX
  elseif (NDIM == 2) then
     hyVelEnd=HY_VELY
  elseif (NDIM == 3) then
     hyVelEnd=HY_VELZ
  endif


  if (dir==DIR_X) then
     HY_VELN   = HY_VELX
#ifdef FLASH_USM_MHD
     HY_MAGN   = HY_MAGX
#endif
  elseif (dir==DIR_Y) then
     HY_VELN   = HY_VELY
#ifdef FLASH_USM_MHD
     HY_MAGN   = HY_MAGY
#endif
  elseif (dir==DIR_Z) then
     HY_VELN   = HY_VELZ
#ifdef FLASH_USM_MHD
     HY_MAGN   = HY_MAGZ
#endif
  endif


  ! Mach2 = Mach^2
  Mach2 = dot_product(V0(HY_VELX:hyVelEnd),V0(HY_VELX:hyVelEnd))&
       /(V0(HY_GAMC)*V0(HY_PRES)/V0(HY_DENS))

  if (Mach2 > 1.e-4) then

     ! Mach2 = (2*gamc*Mach^2 - (gamc-1))(gamc-1)
     Mach2 = (2.*V0(HY_GAMC)*Mach2 - gamm)/gamp

     ! (1a) check if the left state satisfies the RH condition
     !      assuming the left state is behind the shock, that is,
     !      the left state is the shocked gas.
     if (Wl(HY_DENS) >    gammaRatio*V0(HY_DENS)) SWl = .true.
     !NOTE: velocity check won't be useful if gravity is included as
     !      gravity components are added to velocity fields in
     !      reconstruction.
     !if (Wl(HY_VELN) > 1./gammaRatio*V0(HY_VELN)) SWl = .true.
!      if (Wl(HY_PRES) >         Mach2*V0(HY_PRES)) SWl = .true.
! #ifdef FLASH_USM_MHD
!      if (Wl(HY_MAGN) >    gammaRatio*V0(HY_MAGN)) SWl = .true.
! #endif

     ! (1b) check if the left state satisfies the RH condition
     !      assuming the left state is ahead of the shock, that is,
     !      the left state is the pre-shocked gas.
     if (   gammaRatio*Wl(HY_DENS) < V0(HY_DENS)) SWl = .true.
     !if (1./gammaRatio*Wl(HY_VELN) < V0(HY_VELN)) SWl = .true.
!      if (        Mach2*Wl(HY_PRES) < V0(HY_PRES)) SWl = .true.
! #ifdef FLASH_USM_MHD
!      if (   gammaRatio*Wl(HY_MAGN) < V0(HY_MAGN)) SWl = .true.
! #endif

     ! (2a) check if the right state satisfies the RH condition
     !      assuming the right state is behind the shock, that is,
     !      the right state is the shocked gas.
     if (Wr(HY_DENS) >    gammaRatio*V0(HY_DENS)) SWr = .true.
     !if (Wr(HY_VELN) > 1./gammaRatio*V0(HY_VELN)) SWr = .true.
!      if (Wr(HY_PRES) >         Mach2*V0(HY_PRES)) SWr = .true.
! #ifdef FLASH_USM_MHD
!      if (Wr(HY_MAGN) >    gammaRatio*V0(HY_MAGN)) SWr = .true.
! #endif

     ! (2b) check if the right state satisfies the RH condition
     !      assuming the right state is ahead of the shock, that is,
     !      the right state is the pre-shocked gas.
     if (   gammaRatio*Wr(HY_DENS) < V0(HY_DENS)) SWr = .true.
     !if (1./gammaRatio*Wr(HY_VELN) < V0(HY_VELN)) SWr = .true.
!      if (        Mach2*Wr(HY_PRES) < V0(HY_PRES)) SWr = .true.
! #ifdef FLASH_USM_MHD
!      if (   gammaRatio*Wr(HY_MAGN) < V0(HY_MAGN)) SWr = .true.
! #endif

  endif

  ! =========================================================================
  ! (3) Now set first-order if Riemann states are over/under shooting 
  !     at the zone under consideration.
  if ((.not. SWl) .and. (.not. SWr)) then
     return
  else
     Wr(HY_DENS:HY_END_VARS) = V0(HY_DENS:HY_END_VARS)
     Wl(HY_DENS:HY_END_VARS) = V0(HY_DENS:HY_END_VARS)
  endif

End Subroutine hy_uhd_checkRHjumpCond
