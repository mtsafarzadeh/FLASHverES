!!****h* source/physics/materialProperties/Opacity/Opacity_interface
!!
!! NAME
!!
!!  Opacity_interface
!!
!! SYNOPSIS
!!
!!  use Opacity_interface
!!
!! DESCRIPTION
!!
!!  This is the interface file for the Opacity Unit that defines its public interfaces.
!!
!!***

module Opacity_interface

  interface
     subroutine Opacity_init ()
     end subroutine Opacity_init
  end interface

  interface
     subroutine Opacity_finalize ()
     end subroutine Opacity_finalize
  end interface

  interface
     subroutine Opacity (soln, ngrp, opacityAbsorption, opacityEmission, opacityTransport)
       integer, intent(in)  :: ngrp
       real,    intent(out) :: opacityAbsorption
       real,    intent(out) :: opacityEmission
       real,    intent(out) :: opacityTransport
       real,    intent(in), dimension (:) :: soln
     end subroutine Opacity
  end interface

  interface
     subroutine Opacity_unitTest (fileUnit,perfect)
       integer, intent (in)    :: fileUnit
       logical, intent (inout) :: perfect
     end subroutine Opacity_unitTest
  end interface

end module Opacity_interface
