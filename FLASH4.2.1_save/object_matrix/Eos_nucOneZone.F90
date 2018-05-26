!!****f* source/physics/Eos/Eos_nucOneZone
!!
!! NAME
!!
!!  Eos_nucOneZone
!!
!! SYNOPSIS
!!
!!  call Eos_nucOneZone(real, intent(IN)  :: xdens,
!!                      real, intent(INOUT)  :: xtemp,
!!                      real, intent(IN)  :: xye,
!!                      real, intent(INOUT)  :: xener,
!!                      real, intent(INOUT)  :: xpres,
!!                      real, intent(INOUT)  :: xentr,
!!                      real, intent(OUT)  :: xmunu,
!!                      integer, intent(IN)  :: mode)
!!
!! DESCRIPTION
!!
!!   Stub only.
!!
!! ARGUMENTS
!!
!!   xdens : 
!!
!!   xtemp : 
!!
!!   xye : 
!!
!!   xener : 
!!
!!   xpres : 
!!
!!   xentr : 
!!
!!   xmunu : 
!!
!!   mode : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xMuNu,mode)

  implicit none 

  real, intent(IN) :: xDens, xYe
  real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
  real, intent(OUT) :: xMuNu
  integer, intent(IN) :: mode

  ! STUB

  return
end subroutine Eos_nucOneZone
