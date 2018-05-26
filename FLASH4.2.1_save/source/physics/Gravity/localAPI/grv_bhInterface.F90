!!****ih* source/physics/Gravity/localAPI/grv_bhInterface
!!
!! This is an interface module for the Gravity unit that defines some
!! private interfaces specific to the BHTree implementation.
!!***

Module grv_bhInterface

  interface grv_elintF
    real function grv_elintF(theta, k)
      implicit none
      real, intent(in) :: theta, k
    end function grv_elintF
  end interface

  interface grv_intCoef1P2I
    real function grv_intCoef1P2I(k,x,y)
      implicit none
      integer, intent(in) :: k
      real, intent(in) :: x,y
    end function grv_intCoef1P2I
  end interface

  interface grv_derfc
    DOUBLE PRECISION FUNCTION DERFC(X)
      DOUBLE PRECISION X
    END FUNCTION DERFC
  end interface

  interface grv_derfcx
    DOUBLE PRECISION FUNCTION DERFCX(X)
      DOUBLE PRECISION X
    END FUNCTION DERFCX
  end interface

  interface grv_derf
    DOUBLE PRECISION FUNCTION DERF(X)
      DOUBLE PRECISION X
    END FUNCTION DERF
  end interface

  interface grv_besj0
    DOUBLE PRECISION FUNCTION BESJ0(X)
      DOUBLE PRECISION X
    END FUNCTION BESJ0
  end interface


end Module grv_bhInterface

