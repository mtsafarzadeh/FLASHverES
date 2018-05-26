!!****f* source/physics/Diffuse/Diffuse_fluxLimiter
!!
!!  NAME 
!!
!!  Diffuse_fluxLimiter
!!
!!  SYNOPSIS
!!
!!  call Diffuse_fluxLimiter(integer, intent(in) :: idcoeff, 
!!                           integer, intent(in) :: ifunc,
!!                           integer, intent(in) :: ifl,
!!                           integer, intent(in) :: mode,
!!                           integer, intent(IN) :: blkcnt,
!!                           integer, intent(IN) :: blklst(blkcnt))
!!
!!  DESCRIPTION 
!!      This routine modifes the diffusion coefficient (idcoeff) 
!!      and applies limiting.
!!
!! ARGUMENTS
!!
!!   idcoeff : Coefficient on which limiter is applied.
!!   ifunc   : 
!!   ifl     : 
!!   mode    : Flux limiter mode.
!!   blkcnt  : The number of blocks in the list
!!   blklst  : The list of blocks on which the solution must be updated.
!!
!! SIDE EFFECTS
!!
!! Modifies the diffusion coefficient variable, indicted by the
!! argument idcoef
!!
!!***
subroutine Diffuse_fluxLimiter(idcoef, ifunc, ifl, mode, blkcnt, blklst)
  implicit none
  
  ! Arguments:
  integer, intent(in) :: idcoef
  integer, intent(in) :: ifunc
  integer, intent(in) :: ifl
  integer, intent(in) :: mode
  integer, intent(IN) :: blkcnt
  integer, intent(IN) :: blklst(blkcnt)

end subroutine Diffuse_fluxLimiter
