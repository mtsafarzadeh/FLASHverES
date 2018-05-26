!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printMatrix
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
!!  of the matrix will be checked against the requested index printing range.
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

  integer :: block
  integer :: col, colBeg, colEnd
  integer :: maxMatrixRow, maxMatrixCol
  integer :: nColBlocks, nColTotal
  integer :: row

  integer, parameter :: nColPerBlock = 10
!
!
!     ...Print out the title.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') title
  write (fileUnit,'(/)')
!
!
!     ...Check, if printing index ranges are ok.
!
!
  maxMatrixRow = size (matrix,1)
  maxMatrixCol = size (matrix,2)

  if (rowMax > maxMatrixRow .or. colMax > maxMatrixCol) then
      call Driver_abortFlash ("ed_printMatrix: Matrix printing range too large ")
  end if

  if (rowMin < 1 .or. colMin < 1) then
      call Driver_abortFlash ("ed_printMatrix: Matrix printing range < 1 ")
  end if

  if (rowMin > rowMax .or. colMin > colMax) then
      call Driver_abortFlash ("ed_printMatrix: No matrix printing range!")
  end if
!
!
!     ...Everything ok. Start the printing.
!
!
  nColTotal  = colMax - colMin + 1
  nColBlocks = nColTotal / nColPerBlock
!
!
!     ...Print all full column blocks.
!
!
  colEnd = colMin - 1

  do block = 1,nColBlocks

     colBeg = colEnd + 1
     colEnd = colEnd + nColPerBlock

     write (fileUnit,'(/,1x,10i14)') (col , col = colBeg , colEnd)
     write (fileUnit,'(/)')

     do row = rowMin,rowMax
        write (fileUnit,'(i6,10es14.6)') row, (matrix (row,col), col = colBeg , colEnd)
     end do

  end do
!
!
!     ...Print remaining partial column block (if any).
!
!
  colBeg = colEnd + 1
  
  if (colBeg > colMax) return

  write (fileUnit,'(/,1x,10i14)') (col , col = colBeg , colMax)
  write (fileUnit,'(/)')

  do row = rowMin,rowMax
     write (fileUnit,'(i6,10es14.6)') row, (matrix (row,col), col = colBeg , colMax)
  end do
!
!
!     ...Ready!
!
!
  return
end subroutine ed_printMatrix
