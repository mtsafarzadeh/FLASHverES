!!****if* source/physics/Hydro/HydroMain/unsplit_old/Hydro_Unsplit/hy_uhd_energyFix
!!
!! NAME
!!
!!  hy_uhd_energyFix
!!
!! SYNOPSIS
!!
!!  hy_uhd_energyFix( integer (IN) :: blockID,
!!                    integer (IN) :: blkLimits(2,MDIM),
!!                    real(IN)     :: dt,
!!                    real(IN)     :: del(MDIM),
!!                    integer(IN)  :: eosMode)
!!
!! DESCRIPTION
!!
!!  This routine corrects energy by using the internal energy evolution
!!  to avoid any negativity states of pressure in the Eos routines.
!!  This routine also takes care of abundances.
!!
!! ARGUMENTS
!!
!!  blockID   - a local block ID
!!  blkLimits - an array that holds the lower and upper indices of the section
!!              of block without the guard cells
!!  dt        - time step
!!  del       - grid deltas in each direction
!!  eosMode   - a mode used in a call to Eos
!!
!!***

!!REORDER(4):U

Subroutine hy_uhd_energyFix(blockID,blkLimits,dt,del,eosMode)

  use Hydro_data,     ONLY : hy_eswitch, hy_irenorm, hy_smallE

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_renormAbundance, Grid_limitAbundance


  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  real, intent(IN) :: dt
  real, dimension(MDIM), intent(IN) :: del
  integer, intent(IN) :: eosMode
  !! -----------------------------------------------------

  integer :: i,j,k
  real    :: ekin,eint
  real, pointer, dimension(:,:,:,:) :: U
  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone

  call Grid_getBlkPtr(blockID,U,CENTER)


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) < 0.0) then
#endif

#ifndef FLASH_UHD_3T
              ! In case 3T is used then the energy updates are already done in hy_uhd_unsplitUpdate
              ! and are not needed here.
              ekin = .5*dot_product(U(VELX_VAR:VELZ_VAR,i,j,k),U(VELX_VAR:VELZ_VAR,i,j,k))&
                   *U(DENS_VAR,i,j,k)

              eint = U(ENER_VAR,i,j,k)-ekin

              if (eint > hy_eswitch*ekin) then
                 U(EINT_VAR,i,j,k) = max(hy_smallE,eint/U(DENS_VAR,i,j,k))
              else
                 U(EINT_VAR,i,j,k) = max(hy_smallE,U(EINT_VAR,i,j,k)/U(DENS_VAR,i,j,k))
              endif
              !! Store specific gas energy ener = ekin + eint
              U(ENER_VAR,i,j,k) = U(EINT_VAR,i,j,k) + ekin/U(DENS_VAR,i,j,k)
#endif

#ifdef DIVV_VAR
#if NDIM == 1
              U(DIVV_VAR,i,j,k) = (U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k))/del(DIR_X)
#elif NDIM == 2
              U(DIVV_VAR,i,j,k) = (U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k))/del(DIR_X)&
                                 +(U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k))/del(DIR_Y)
#elif NDIM == 3
              U(DIVV_VAR,i,j,k) = (U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k))/del(DIR_X)&
                                 +(U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k))/del(DIR_Y)&
                                 +(U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1))/del(DIR_Z)
#endif
#endif

#ifdef BDRY_VAR
           endif
#endif
        enddo
     enddo
  enddo


  if (eosMode==MODE_DENS_PRES) then
     U(PRES_VAR,:,:,:) = U(EINT_VAR,:,:,:)*U(DENS_VAR,:,:,:)*(U(GAME_VAR,:,:,:)-1.)
  endif

  ! Renormalize or limit abundances
  if (hy_irenorm == 1) then
     call Grid_renormAbundance(blockID,blkLimits,U)
  else
     call Grid_limitAbundance(blkLimits,U)
  endif
  

  call Grid_releaseBlkPtr(blockID,U,CENTER)

  return
End Subroutine hy_uhd_energyFix
