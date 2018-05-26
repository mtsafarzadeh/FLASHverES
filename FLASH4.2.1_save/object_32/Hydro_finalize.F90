!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_finalize
!!
!! NAME
!!
!!  Hydro_finalize
!!
!!
!! SYNOPSIS
!!
!!  Hydro_finalize()
!!  
!!
!! DESCRIPTION
!! 
!!  Finalizes the Hydro_unit.  Most likely will deallocate any memory allocated in 
!!  Hydro_init
!!
!!***

subroutine Hydro_finalize()

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
  use Hydro_data, only : hy_SpcR,hy_SpcL,hy_SpcSig
#endif
#endif

  implicit none

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
  if (allocated(hy_SpcR))   deallocate(hy_SpcR)
  if (allocated(hy_SpcL))   deallocate(hy_SpcL)
  if (allocated(hy_SpcSig)) deallocate(hy_SpcSig)
#endif
#endif

end subroutine Hydro_finalize

