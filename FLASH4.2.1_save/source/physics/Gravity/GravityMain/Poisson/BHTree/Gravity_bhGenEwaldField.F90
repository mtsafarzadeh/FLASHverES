!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhGenEwaldField
!!
!! NAME
!!
!!  Gravity_bhGenEwaldField
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhGenEwaldField()
!!
!! DESCRIPTION
!!
!!  Generates the Ewald field. It is done in parallel, the whole field is
!!  distributed among all cpus. The field is eventually saved in file 
!!  grv_bhEwaldFName. In case grv_bhEwaldAlwaysGenerate is FALSE and
!!  the file with the Ewald field exists, the Ewald field is read from it.
!!  Ewald field is stored on a nested grid (with the highest density at small
!!  distances), one level of the EF grid is calculated by
!!  grv_bhGenEwaldFieldLevel which is called by this subroutine.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***




subroutine Gravity_bhGenEwaldField()
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Gravity_data, ONLY : grv_bhEwaldFieldNx, grv_bhEwaldFieldNy, grv_bhEwaldFieldNz, &
    grv_bhEwaldFName, grv_bhEwaldAlwaysGenerate, &
    grv_bhLx, grv_bhLy, grv_bhLz, grv_bhEwaldLMax, &
    grv_bhDxI, grv_bhDyI, grv_bhDzI, grv_bhMinEFSize, &
    grv_meshMe, grv_meshComm, grv_bhTreeEwald, &
    grv_bhEwaldSeriesN, grv_bhEwaldNRef, grv_bhDirectionQuad
  use Gravity_interface, ONLY : Gravity_bhGenEwaldFieldLevel
  use Grid_interface, ONLY : Grid_getMinCellSizes
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer :: i, j, k, l, ierr, nref
  integer :: ewald_periodicity
  real :: xmax, xmin, ymax, ymin, zmax, zmin, dum
  real :: x, y, z, ewald_xmax, ewald_ymax, ewald_zmax
  real, dimension(MDIM) :: mcs
  real, allocatable :: field_prep(:,:,:)
  logical :: file_exists
  logical :: grv_bhLinearInterpolOnly
  character(len=MAX_STRING_LENGTH) :: strBuff, grav_boundary_type_x &
  & , grav_boundary_type_y, grav_boundary_type_z

  call RuntimeParameters_get("grav_boundary_type_x", grav_boundary_type_x)
  call RuntimeParameters_get("grav_boundary_type_y", grav_boundary_type_y)
  call RuntimeParameters_get("grav_boundary_type_z", grav_boundary_type_z)

  ! return if periodic boundaries are not used in any direction
  if (    (grav_boundary_type_x .ne. "periodic") &
  & .and. (grav_boundary_type_y .ne. "periodic") &
  & .and. (grav_boundary_type_z .ne. "periodic")) return

  ! read runtime parameters
  call RuntimeParameters_get("grv_bhEwaldFieldNx", grv_bhEwaldFieldNx)
  call RuntimeParameters_get("grv_bhEwaldFieldNy", grv_bhEwaldFieldNy)
  call RuntimeParameters_get("grv_bhEwaldFieldNz", grv_bhEwaldFieldNz)
  call RuntimeParameters_get("grv_bhEwaldSeriesN", grv_bhEwaldSeriesN)
  call RuntimeParameters_get("grv_bhEwaldNRef", grv_bhEwaldNRef)
  call RuntimeParameters_get("grv_bhEwaldFName", grv_bhEwaldFName)
  call RuntimeParameters_get("grv_bhEwaldAlwaysGenerate", grv_bhEwaldAlwaysGenerate)
  call RuntimeParameters_get("xmin", xmin)
  call RuntimeParameters_get("xmax", xmax)
  call RuntimeParameters_get("ymin", ymin)
  call RuntimeParameters_get("ymax", ymax)
  call RuntimeParameters_get("zmin", zmin)
  call RuntimeParameters_get("zmax", zmax)
  call RuntimeParameters_get("grv_bhLinearInterpolOnly", grv_bhLinearInterpolOnly)

  ! dimensions of the computational domain
  grv_bhLx = xmax - xmin
  grv_bhLy = ymax - ymin
  grv_bhLz = zmax - zmin
  grv_bhEwaldLMax = max(grv_bhLx,grv_bhLy,grv_bhLz)

  ! dimensions of one grid-cell of the largest Ewald field,
  ! used in grv_bhEwald.F90
  grv_bhDxI = (grv_bhEwaldFieldNx-1) / grv_bhEwaldLMax
  grv_bhDyI = (grv_bhEwaldFieldNy-1) / grv_bhEwaldLMax
  grv_bhDzI = (grv_bhEwaldFieldNz-1) / grv_bhEwaldLMax

  ! number of refinement levels of the Ewald field
  if (grv_bhEwaldNRef <= 0) then
    call Grid_getMinCellSizes(mcs)
    grv_bhEwaldNRef = 1
    do
      !if (grv_meshMe == MASTER_PE)print *, "NRef: ", grv_bhEwaldNRef, mcs, &
      !& "|", 1/grv_bhDxI, 1/grv_bhDyI, 1/grv_bhDzI
      grv_bhMinEFSize = grv_bhEwaldLMax/(2**(grv_bhEwaldNRef-1))
      if (    (grv_bhMinEFSize < 0.5*mcs(IAXIS)) &
      & .and. (grv_bhMinEFSize < 0.5*mcs(JAXIS)) &
      & .and. (grv_bhMinEFSize < 0.5*mcs(KAXIS))) exit
      grv_bhEwaldNRef = grv_bhEwaldNRef + 1
    enddo
  else
    grv_bhMinEFSize = grv_bhEwaldLMax/(2**(grv_bhEwaldNRef-1))
  endif

  write (strBuff, '("Ewald refinement level determined: ",i3)') grv_bhEwaldNRef 
  call Logfile_stamp( strBuff, "[BHTree]")

  ! periodicity and orientation of the problem
  ewald_periodicity = 0
  if (grav_boundary_type_x == "periodic") ewald_periodicity = ewald_periodicity + 1
  if (grav_boundary_type_y == "periodic") ewald_periodicity = ewald_periodicity + 2
  if (grav_boundary_type_z == "periodic") ewald_periodicity = ewald_periodicity + 4

  ! choose between linear and quadratic (quadratic only in one direction) 
  ! interpolation method (for routine Gravity_bhEwald)
  if (grv_bhLinearInterpolOnly) then
    grv_bhDirectionQuad = 0
  else
    if (ewald_periodicity == 6) then
      grv_bhDirectionQuad = 1
    else if (ewald_periodicity == 5) then
      grv_bhDirectionQuad = 2
    else if (ewald_periodicity == 3) then
      grv_bhDirectionQuad = 3
    else
      grv_bhDirectionQuad = 0
    endif
  endif

  ! allocate array for the Ewald field
  allocate(grv_bhTreeEwald(0:grv_bhEwaldNRef-1,  -1:grv_bhEwaldFieldNx, &
  &                        -1:grv_bhEwaldFieldNy, -1:grv_bhEwaldFieldNz), stat=ierr)
  if (ierr /= 0) call Driver_abortFlash("could not allocate grv_bhTreeEwald in Gravity_bhGenEwaldField.F90")


  ! check the existence of the ewald_field file with correct grv_bhEwaldNRef
  if (grv_meshMe == MASTER_PE) then
    inquire(file=grv_bhEwaldFName, exist=file_exists)
    if (file_exists .and. (.not. grv_bhEwaldAlwaysGenerate)) then
      open(unit=53, file=grv_bhEwaldFName, status='old')
      read(53,*) nref
      if (nref /= grv_bhEwaldNRef) then
        file_exists = .false.
        close(unit=53)
      endif
    else
      file_exists = .false.
    endif
  endif
  call MPI_Bcast(file_exists, 1, FLASH_LOGICAL, MASTER_PE, grv_meshComm, ierr)

  if (file_exists) then
    ! ewald_field file with correct grv_bhEwaldNRef is present and there is 
    ! no user request to regenerate it => read the ewald field from the file
    write (strBuff, '("Reading Ewald field from file: ",a20)') grv_bhEwaldFName
    call Logfile_stamp( strBuff, "[BHTree]")
    if (grv_meshMe == MASTER_PE) then
      do l = 0, grv_bhEwaldNRef-1
        do k = -1,grv_bhEwaldFieldNz
          do j = -1,grv_bhEwaldFieldNy
            do i = -1,grv_bhEwaldFieldNx
              read(53,*) x, y, z, grv_bhTreeEwald(l,i,j,k), dum
            enddo
            read(53,*)
          enddo
          read(53,*)
        enddo
        read(53,*)
      enddo
      close(unit=53)
    endif
    call MPI_Bcast(grv_bhTreeEwald, grv_bhEwaldFieldNx*grv_bhEwaldFieldNy*grv_bhEwaldFieldNz &
    &  , FLASH_REAL, MASTER_PE, grv_meshComm, ierr)
    write (strBuff, '("Ewald correction field read and communicated.")')
    call Logfile_stamp( strBuff, "[BHTree]")

  else

    ! Ewald field will be generated
    write (strBuff, '("Generating Ewald correction field")')
    call Logfile_stamp( strBuff, "[BHTree]")
    allocate(field_prep(-1:grv_bhEwaldFieldNx, -1:grv_bhEwaldFieldNy, &
    & -1:grv_bhEwaldFieldNz), stat=ierr)
    if (ierr /= 0) call &
    &  Driver_abortFlash("could not allocate grv_bhTreeEwald in Gravity_bhGenEwaldField.F90")
    do l = 0, grv_bhEwaldNRef-1
      ewald_xmax = grv_bhEwaldLMax/(2**l)
      ewald_ymax = grv_bhEwaldLMax/(2**l)
      ewald_zmax = grv_bhEwaldLMax/(2**l)

      ! call subroutine which generate ewald field with appropriate boundary conditions
      call Gravity_bhGenEwaldFieldLevel(grv_bhEwaldFieldNx, grv_bhEwaldFieldNy, &
           & grv_bhEwaldFieldNz, ewald_xmax, ewald_ymax, ewald_zmax, &
           & ewald_periodicity, field_prep)

      grv_bhTreeEwald(l,:,:,:) = field_prep(:,:,:)

    enddo
    deallocate(field_prep)


    ! write the ewald field into the file
    if (grv_meshMe == MASTER_PE) then
 
      open(unit=53, file=grv_bhEwaldFName, status='replace')
      write(53,*) grv_bhEwaldNRef
 
      do l = 0, grv_bhEwaldNRef-1
        ewald_xmax = grv_bhEwaldLMax/(2**l)
        ewald_ymax = grv_bhEwaldLMax/(2**l)
        ewald_zmax = grv_bhEwaldLMax/(2**l)
        do k = -1,grv_bhEwaldFieldNz
          z = (k * ewald_zmax) / (grv_bhEwaldFieldNz-1)
          do j = -1,grv_bhEwaldFieldNy
            y = (j *  ewald_ymax) / (grv_bhEwaldFieldNy-1)
            do i = -1,grv_bhEwaldFieldNx
              x = (i * ewald_xmax) / (grv_bhEwaldFieldNx-1)
              write(53,'(5(e17.10, 2x))') x, y, z, grv_bhTreeEwald(l,i,j,k) &
              , 1.0/(sqrt(x*x+y*y+z*z)+1d-99)
            enddo
            write(53,*)
          enddo
          write(53,*)
        enddo
        write(53,*)
      enddo
      close(unit=53)
      write (strBuff, '("Ewald correction field written into file:", a20)') grv_bhEwaldFName
      call Logfile_stamp( strBuff, "[BHTree]")
    endif ! MASTER_PE

  endif ! ewald field generated

  return
end subroutine Gravity_bhGenEwaldField

