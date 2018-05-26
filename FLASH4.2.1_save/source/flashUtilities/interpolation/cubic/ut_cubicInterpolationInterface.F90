!!****ih* source/flashUtilities/interpolation/cubic/ut_cubicInterpolationInterface
!!
!! NAME
!!
!!  ut_cubicInterpolationInterface
!!
!! SYNOPSIS
!!
!!   use ut_cubicInterpolationInterface
!!
!! DESCRIPTION
!!
!!  This is the module needed for defining the cubic interpolation interface.
!!
!!***

Module ut_cubicInterpolationInterface

  interface
     subroutine ut_monoCubicCoefficients (numberOfLines, a)
       integer, intent (in)    :: numberOfLines
       real,    intent (inout) :: a (1:4,*)
     end subroutine ut_monoCubicCoefficients
  end interface

  interface
     subroutine ut_biCubicCoefficients (numberOfSquares, a)
       integer, intent (in)    :: numberOfSquares
       real,    intent (inout) :: a (1:16,*)
     end subroutine ut_biCubicCoefficients
  end interface

  interface
     subroutine ut_triCubicCoefficients (numberOfCubes, a)
       integer, intent (in)    :: numberOfCubes
       real,    intent (inout) :: a (1:64,*)
     end subroutine ut_triCubicCoefficients
  end interface


end Module ut_cubicInterpolationInterface
