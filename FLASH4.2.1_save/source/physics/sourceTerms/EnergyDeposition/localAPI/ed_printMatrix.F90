!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_printMatrix
!!
!! NAME
!!
!!  ed_printMatrix
!!
!! SYNOPSIS
!!
!!  call ed_printMatrix (integer,           intent (in) :: fileUnit,
!!                       character (len=*), intent (in) :: title,
!!                       integer,           intent (in) :: rowMin,
!!                       integer,           intent (in) :: rowMax,
!!                       integer,           intent (in) :: colMin,
!!                       integer,           intent (in) :: colMax,
!!                       real,              intent (in) :: matrix (:,:))
!!
!! DESCRIPTION
!!
!!  General matrix printing routine. The matrix is printed, together with the descriptive
!!  title, to the file associated with the unit number 'fileUnit'. The dimensions of the
!!  matrix need not to be equal to the actual matrix section printed, however the size
!!  of the matrix will be checked with the index printing range.
!!
!! ARGUMENTS
!!
!!  fileUnit : the unit number for the printout file
!!  title    : the printout title
!!  rowMin   : the minimum row index of the matrix section to be printed
!!  rowMax   : the maximum row index of the matrix section to be printed
!!  colMin   : the minimum column index of the matrix section to be printed
!!  colMax   : the maximum column index of the matrix section to be printed
!!  matrix   : the matrix of assumed shape
!!
!! NOTES
!!
!!  The columns of the matrix are printed such that only 10 columns are printed
!!  simultaneously on each line.
!!
!!***

subroutine ed_printMatrix (fileUnit, title, rowMin, rowMax, colMin, colMax, matrix)

  use Driver_interface,       ONLY : Driver_abortFlash

  implicit none

  integer,           intent (in) :: fileUnit
  character (len=*), intent (in) :: title
  integer,           intent (in) :: rowMin, rowMax
  integer,           intent (in) :: colMin, colMax
  real,              intent (in) :: matrix (:,:)

  return
end subroutine ed_printMatrix
