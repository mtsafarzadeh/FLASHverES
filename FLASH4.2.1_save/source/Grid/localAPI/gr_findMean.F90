!!****if* source/Grid/localAPI/gr_findMean
!!
!! NAME
!!  gr_findMean
!! 
!! SYNOPSIS
!!  gr_findMean(integer, intent(in)  :: iSrc,
!!                     integer, intent(in)  :: iType,
!!                     logical, intent(in)  :: bGuardcell
!!                        real, intent(out) :: mean)
!!
!! DESCRIPTION
!!  Calculates the mean of a function
!!
!!
!! ARGUMENTS
!!  iSrc -- the index (e.g. DENS_VAR) into the unk routine to calculate over
!!  iType -- the type of mean.  Valid values will be  (feel free to implement more)
!!     1 = L1 norm = arithmetic average
!!  bGuardcell -- logical indicating whether guard cells should be included in the calculation
!!  mean -- the requested mean
!!
!!
!!***

subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
  
  implicit none

  integer, intent(in) :: iSrc, iType
  logical, intent(in) :: bGuardcell
  real, intent(out) :: mean


end subroutine gr_findMean

