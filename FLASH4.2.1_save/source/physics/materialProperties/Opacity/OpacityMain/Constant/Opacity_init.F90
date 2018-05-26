!!****if* source/physics/materialProperties/Opacity/OpacityMain/Constant/Opacity_init
!!
!! NAME
!!
!!  Opacity_init
!!
!!
!! SYNOPSIS
!!
!!  call Opacity_init()
!!
!! Description
!!
!! Initialiazed data for the Constant opacity model using run time
!! parameters.
!!
!! ARGUMENTS
!!
!!***
subroutine Opacity_init()
  use Opacity_data, ONLY : op_emitConst, op_transConst, &
       op_absorbConst
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  call RuntimeParameters_get("op_emitConst", op_emitConst)
  call RuntimeParameters_get("op_transConst", op_transConst)
  call RuntimeParameters_get("op_absorbConst", op_absorbConst)

end subroutine Opacity_init
