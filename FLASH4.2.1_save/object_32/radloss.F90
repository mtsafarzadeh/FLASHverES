!
!       The cooling function is from Sarazin (1986) Rev.Mod.Phys.
!       and from Raymond 1976
!       
!       ROUTINE FROM PERES ET AL. 1982, APJ 252, 791

!!****f* source/source_terms/cool/radloss/radloss
!!
!! NAME
!!  
!!  radloss
!!
!!
!! SYNOPSIS
!! 
!!  call radloss(T, radia)
!!  call radloss(real, real)
!! 
!! 
!! DESCRIPTION
!!
!!
!!  .
!!
!!
!! ARGUMENTS
!!
!!  T: 
!!
!!  radia: 
!!
!!***

        subroutine radloss(T,radia)
!
        implicit none

        real*4 ts(9)
        real, intent(IN) :: T
        real, intent(OUT) :: radia

        real exp
        integer iir, is, ial

        data ts/4.44e3,8.e3,2.00e4,3.98e4,7.94e4,2.51e5,5.62e5,2.00e6,  &
     &              1.00e7/
        data iir /9/
        data exp/-0.66666666/

!
!
!
        if (T.ge.4.e7) then
            radia = 2.692e-27*sqrt(T)
        else if (T.ge.1.e7) then
            radia = 6.2e-19/T**0.6
        else
        if (T.gt.4.44e3) goto 40
            radia = 5.e-28
            return
40          continue
!                               !binary search
        is  = 8
        ial = 8
35      ial = ial/2
        if (is.le.iir) goto 50
            is = is-ial
            goto 35
50          continue
        if (T.gt.ts(is).or.ial.eq.0) goto 60
            is = is-ial
            goto 35
60      if (T.le.ts(is).or.ial.eq.0) goto 70
            is = is+ial
            goto 35
70      if (T.lt.ts(is).and.ial.eq.0) is = is-1
!
        goto (79,80,81,82,83,84,85,86) is
        stop 'error in radloss'
!
79      radia = (1.0606e-6*T)**11.732
        return
!
80      radia = (1.3972e-8*T)**6.1496
        return
!
81      radia = 1.41e-22
        return
!
82      radia = 1.e-31*T**2
        return
!
83      radia = 6.31e-22
        return
!
84      radia = 3.98e-11/T**2
        return
!
85      radia = 1.15e-22
        return
!
86      radia = 1.86e-18*T**exp
        return
!
        endif
!
!
        return
        end
