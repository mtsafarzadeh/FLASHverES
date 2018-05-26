!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!  
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  This modules stores the data for the Gravity unit
!!
!!***

module Gravity_data

#include "constants.h"
#include "Flash.h"

  ! Common Gravity module variables
  character(len=MAX_STRING_LENGTH), save :: grav_boundary_type !string boundary condition
  integer, save :: grav_boundary(6)  !integer boundary condition

  integer, save :: grav_geometry  !mesh geometry
  integer, save :: grv_meshMe, grv_meshNumProcs, grv_meshComm
  integer, save :: grv_commSize=1

  logical, save :: useGravity, updateGravity
  logical, save :: grav_temporal_extrp !extrapolate or otherwise rescale
  logical, save :: grav_unjunkPden

  real,    save :: grav_poisfact

  ! Variables specific to the BHTree
  character(len=MAX_STRING_LENGTH), save :: grv_bhMAC

  ! configuration parameters
  logical, save :: grv_bhPhysMACTW, grv_bhPhysMACComm
  logical, save :: grv_bhUseRelAccErr
  logical, save :: grv_bhUseEwald

  ! integer vars (mainly degree of the multipole expansion for APE MAC)
  integer, save :: grv_bhMPDegree, grv_bhMACNum
  integer, save :: grv_bhMPDegree_p2, grv_bhMPDegree_half, grv_bhMPDegree_halfp1

  ! real vars (mainly error in acceleration)
  real,    save :: grv_bhMeanBlockAccErrInv
  real,    save :: grv_bhAccErr, grv_bhAccErrInv

  ! integer constants (mainly indeces to obtain given quantities from arrays)
  integer, save :: grv_bhNSIZE = 0, grv_bhBNSIZE = 0
  integer, save :: grv_bhIM, grv_bhIX, grv_bhIY, grv_bhIZ
  integer, save :: grv_bhNODE5, grv_bhIDMIN, grv_bhIB2, grv_bhIB3
  integer, save :: grv_bhA1DIR, grv_bhA2DIR, grv_bhA3DIR
  integer, save :: grv_defaultGpotVar = GPOT_VAR

  ! real constants
  real,    save :: grv_bhPiGhalf
  real,    save :: grv_bhNewton
  real,    save :: grv_bhSContrGeom ! self contribution geometry

  ! variables related to periodic boundaries and the Ewald field
  character(len=MAX_STRING_LENGTH), save :: grv_bhEwaldFName
  logical, save :: grv_bhEwaldAlwaysGenerate
  integer, save :: grv_bhEwaldFieldNx, grv_bhEwaldFieldNy, grv_bhEwaldFieldNz
  integer, save :: grv_bhEwaldNRef, grv_bhEwaldSeriesN
  integer, save :: grv_bhDirectionQuad
  real,    save :: grv_bhEwaldLMax, grv_bhDxI, grv_bhDyI, grv_bhDzI
  real,    save :: grv_bhMinEFSize, grv_bhLx,  grv_bhLy,  grv_bhLz
  real, save, allocatable :: grv_bhTreeEwald(:,:,:,:)

  ! constants representing individual MACs
  integer, parameter :: grv_bhMAC_APE = 1
  integer, parameter :: grv_bhMAC_MPE = 2
  integer, parameter :: grv_bhMAC_SS  = 3

  ! meaning of the fifth field in tree nodes
  integer, parameter :: grv_bhN5_NONE = 0
  integer, parameter :: grv_bhN5_DMIN = 1 ! minimum distance
  integer, parameter :: grv_bhN5_B2   = 2 ! second order moment

  ! B3 moment of a uniform unit cube with rho = 1
  ! obtained numerically, because I didn't find the analytical integral :(
  !real, parameter :: GRAV_B3_CONST = 0.137674463674

end module Gravity_data
