!!****if* source/Simulation/SimulationMain/Cool_Test/old_Driver_computeDt
!!
!! NAME
!!
!!  Driver_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDt(integer(IN) :: myPE,  
!!                  integer(IN) :: numProcs, 
!!                  integer(IN) :: nbegin,
!!                  integer(IN) :: nstep, 
!!                  real(IN)    :: simTime,
!!                  real(IN)    :: dtOld,
!!                  real(OUT)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Determine the stability-limited time step.
!!  This timestep is determined using information from the included
!!  physics modules - many different timestep limiters are polled.
!!
!!  The global driver might use a different (hopefully smaller) time
!!  step, to match a file write time (tplot or trstr) or if the
!!  simulation end time has been reached; such possibilities are
!!  not considered here.
!!
!! ARGUMENTS
!!
!!  myPE - current processor number
!!  numProcs - number of processors in the run 
!!  nbegin - first step of the simulation (nbegin is only used
!!              to determine if a label header should be written to
!!              the screen)
!!  nstep - current step of the simulation
!!  simTime - current simulation time of the run
!!  dtOld - the dt from the timestep that we just finished 
!!         (it's old because we be using dtOld to calculate 
!!          and return the dt for the next timestep (dtNew)
!!  dtNew - returned value of the dt calculated for the next timestep
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_myPE or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!! The calls to units currently not included in the code are commented out.
!!
!!
!!
!!*** 
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_computeDt(myPE, numProcs,  &
                    nbegin, nstep, &
                    simTime, dtOld, dtNew)

  use Driver_data, ONLY : dr_dtMin,dr_dtMax, dr_tstepChangeFactor, &
       dr_redshift, dr_useRedshift
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getCellCoords, Grid_getDeltas, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getSingleCellCoords
  use Hydro_interface, ONLY : Hydro_computeDt
  use Stir_interface, ONLY: Stir_computeDt
  use Cosmology_interface, ONLY: Cosmology_computeDt
  use Cool_interface, ONLY : Cool_computeDt
  use Heat_interface, ONLY : Heat_computeDt
  use Diffuse_interface, ONLY : Diffuse_computeDt 
  use Burn_interface, ONLY : Burn_computeDt
  use PrimordialChemistry_interface, ONLY : PrimordialChemistry_computeDt
  implicit none

#include "constants.h"
#include "Flash.h"
  include "Flash_mpi.h"

  integer, intent(IN) :: myPE, numProcs
  integer, intent(IN) :: nbegin, nstep
  real,    intent(IN) :: simTime    !! current simulation time
  real, intent(IN) :: dtOld      !! last time step we used
  real, intent(OUT):: dtNew      !! the new timestep we get. to be returned.
 

  ! Local variables and functions
  integer :: i, error, blockID, numLeafBlocks, iout

  !! This is arbitrarily fixed. Users that add more units that influence the time
  !! should change this.

  integer, parameter :: nUnits = 11
  real, PARAMETER :: MAX_TSTEP = huge(1.0)

  
  real    :: dtModule(2,nUnits), dtLocal(2,nUnits)
  integer :: dtMinLoc(5), lminloc(5,nUnits), ngmin, pgmin
  integer :: status(MPI_Status_Size)

  logical :: gcell = .true.
  real, DIMENSION(MDIM) :: coords

  integer, dimension(MAXBLOCKS) :: blockList

  real, dimension(nUnits) :: tstepOutput
  character (len=20), DIMENSION(nUnits) :: &
                         limiterName, limiterNameOutput

  !!prepatory data structures for passing coords to timestep routines
  real, dimension(MDIM) :: del
  integer, dimension(MDIM) :: index

#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_ILO_GC:GRID_IHI_GC) :: xLeft,xRight,xCenter
  real,dimension(GRID_JLO_GC:GRID_JHI_GC) :: yLeft,yRight,yCenter
  real,dimension(GRID_KLO_GC:GRID_KHI_GC) :: zLeft,zRight,zCenter
  real, dimension(GRID_ILO_GC:GRID_IHI_GC) ::  dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC) ::  dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC) ::  dz, uzgrid
#else
  real, allocatable,dimension(:)::&
       dx,uxgrid,dy,uygrid,dz,uzgrid
  real, allocatable,dimension(:)::xLeft,xRight,xCenter
  real, allocatable,dimension(:)::yLeft,yRight,yCenter
  real, allocatable,dimension(:)::zLeft,zRight,zCenter

#endif

  !arrays which hold the starting and ending indicies of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer :: solnData(:,:,:,:)
  integer :: isize,jsize,ksize

  logical :: printTStepLoc
  integer :: itempLimit = 0
  integer, parameter :: HYDRO=1,BURN=2,GRAV=3,HEAT=4,COOL=5,TEMP=6,&
                        PART=7,DIFF=8,COSMO=9,STIR=10,CHEM=11


  data limiterName(HYDRO) /'dt_hydro'/
  data limiterName(HEAT) /'dt_Heat'/
  data limiterName(PART) /'dt_Part '/
  data limiterName(BURN) /'dt_Burn '/
  data limiterName(COOL) /'dt_Cool '/
  data limiterName(TEMP) /'dt_Temp '/
  data limiterName(DIFF) /'dt_Diff '/
  data limiterName(COSMO) /'dt_Cosm'/
  data limiterName(STIR) /'dt_Stir'/
  data limiterName(CHEM) /'dt_Chem'/



  !            Find the local minimum timestep among the included physics
  !            modules for locally stored blocks.
  
  !            Initialize all timestep variables.




  printTStepLoc=.false.
  
  dtMinLoc(:) = 0
  lminloc(:,:) = 0

  do i = 1, nUnits
     dtLocal(1,i) = MAX_TSTEP
     dtLocal(2,i) = real(myPE)
  enddo

  !            Loop over all local leaf-node blocks
  
  call Grid_getListOfBlocks(LEAF,blockList, numLeafBlocks)
  
  do i = 1, numLeafBlocks
     !!Get the coordinate information for all the
     call Grid_getBlkIndexLimits(blockList(i),blkLimits,blkLimitsGC)
     isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
#ifndef FIXEDBLOCKSIZE
     allocate(xLeft(isize))
     allocate(xRight(isize))
     allocate(xCenter(isize))
     allocate(dx(isize))
     allocate(uxgrid(isize))
     allocate(yLeft(jsize))
     allocate(yRight(jsize))
     allocate(yCenter(jsize))
     allocate(dy(jsize))
     allocate(uygrid(jsize))
     allocate(zLeft(ksize))
     allocate(zRight(ksize))
     allocate(zCenter(ksize))
     allocate(dz(ksize))
     allocate(uzgrid(ksize))
#endif
#ifdef DEBUG_DRIVER
     print*,'before calling get coordinates',isize,gcell
#endif
     call Grid_getCellCoords(IAXIS,blockList(i),CENTER,gcell,xCenter,isize)
     call Grid_getCellCoords(IAXIS,blockList(i),LEFT_EDGE,gcell,xLeft,isize)
     call Grid_getCellCoords(IAXIS,blockList(i),RIGHT_EDGE,gcell,xRight,isize)

#ifdef DEBUG_DRIVER
     print*,'before calling get coordinates',jsize,gcell
#endif
     call Grid_getCellCoords(JAXIS,blockList(i),CENTER,gcell,yCenter,jsize)
     call Grid_getCellCoords(JAXIS,blockList(i),LEFT_EDGE,gcell,yLeft,jsize)
     call Grid_getCellCoords(JAXIS,blockList(i),RIGHT_EDGE,gcell,yRight,jsize)



#ifdef DEBUG_DRIVER
     print*,'before calling get coordinates',ksize,gcell
#endif
     call Grid_getCellCoords(KAXIS,blockList(i),CENTER,gcell,zCenter,ksize)
     call Grid_getCellCoords(KAXIS,blockList(i),LEFT_EDGE,gcell,zLeft,ksize)
     call Grid_getCellCoords(KAXIS,blockList(i),RIGHT_EDGE,gcell,zRight,ksize)


     call Grid_getDeltas(blockList(i), del)
     dx(:) = del(1)
     dy(:) = del(2)
     dz(:) = del(3)

     uxgrid(:) = 0
     uygrid(:) = 0
     uzgrid(:) = 0

     call Grid_getBlkPtr(blockList(i),solnData)

     ! hydro
#ifdef DEBUG_DRIVER
     print*,'going to call Hydro timestep'
#endif


     call Hydro_computeDt ( blockList(i), myPE, &
                           xCenter, dx, uxgrid, &
                           yCenter, dy, uygrid, &
                           zCenter, dz, uzgrid, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,HYDRO), lminloc(:,HYDRO) )


     call Stir_computeDt ( blockList(i), myPE, &
                           blkLimits,blkLimitsGC,  &
                           solnData,               &
                          dtLocal(1,STIR), lminloc(:,STIR) )
   

#ifdef DEBUG_DRIVER
     print*,'returned from hydro timestep'
#endif

     call Burn_computeDt ( blockList(i), myPE, &
                           blkLimits,blkLimitsGC,  &
                           solnData,               &
                           dtLocal(1,BURN), lminloc(:,BURN) )
     
!!$     call Gravity_computeDt ( blockList(i), myPE, &
!!$                           xCenter,xLeft,xRight, dx, uxgrid, &
!!$                           yCenter,yLeft,yRight, dy, uygrid, &
!!$                           zCenter,zLeft,zRight, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,GRAV), lminloc(:,GRAV) )
     
     
     call Heat_computeDt ( blockList(i), myPE, &
                           xCenter, dx, uxgrid, &
                           yCenter, dy, uygrid, &
                           zCenter, dz, uzgrid, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,HEAT), lminloc(:,HEAT) )

     
     call Cool_computeDt ( blockList(i), myPE, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,COOL), lminloc(:,COOL) )

     call PrimordialChemistry_computeDt ( blockList(i), myPE, &
                                blkLimits,blkLimitsGC, &
                                solnData,    &
                                dtLocal(1,CHEM), lminloc(:,CHEM))


     call Particles_computeDt &
          ( blockList(i),myPE, dtLocal(1,PART), lminloc(:,PART))
                           

     call Diffuse_computeDt ( blockList(i), myPE, &
                           xCenter,xLeft,xRight, dx, uxgrid, &
                           yCenter,yLeft,yRight, dy, uygrid, &
                           zCenter,zLeft,zRight, dz, uzgrid, &
                           blkLimits,blkLimitsGC,  &
                           solnData,      &
                           dtLocal(1,DIFF), lminloc(:,DIFF) )

!! WARNING -- should there be a Stir_computeDT here? There isn't one in the unit.

     
!!$     call Cosmo_timestep ( blockList(i), myPE, &
!!$                           xCenter,xLeft,xRight, dx, uxgrid, &
!!$                           yCenter,yLeft,yRight, dy, uygrid, &
!!$                           zCenter,zLeft,zRight, dz, uzgrid, &
!!$                           blkLimits,blkLimitsGC,  &
!!$                           solnData,      &
!!$                           dtLocal(1,COSMO), lminloc(:,COSMO) )
      call Cosmology_computeDt(dtLocal(1,COSMO))
     
#ifndef FIXEDBLOCKSIZE
     deallocate(xCenter)
     deallocate(xLeft)
     deallocate(xRight)
     deallocate(dx)
     deallocate(uxgrid)
     deallocate(yCenter)
     deallocate(yLeft)
     deallocate(yRight)
     deallocate(dy)
     deallocate(uygrid)
     deallocate(zCenter)
     deallocate(zLeft)
     deallocate(zRight)
     deallocate(dz)
     deallocate(uzgrid)
#endif
#ifdef DEBUG_DRIVER
     print*,'release blockpointer'
#endif

     call Grid_releaseBlkPtr(blockList(i),solnData)
  enddo



  ! we disabled temperature timestep limiter for now.  
  ! The old temperature was not updated with the refinement, 
  ! so dT/T was precomputed and the number of blocks may not be
  ! the same as there are now.
!!$  if (itempLimit == 1) then
!!$     do blockID = 1,  MAXBLOCKS
!!$       call Driver_computeDtTemp(myPE, dtOld, dtLocal(1,6), &
!!$             lminloc(1,6), blockID)
!!$     enddo
!!$  endif
  


  !            Find the minimum timestep across all processors and all
  !            modules.

  call MPI_AllReduce (dtLocal(1,1), dtModule(1,1), nUnits, & 
       MPI_2Double_Precision, MPI_MinLoc, MPI_Comm_World, error)

  dtNew    = huge(dtNew)             ! dt will hold the minimum timestep
  ngmin = 1                 ! ngmin will hold the winning module #
  pgmin = MASTER_PE          ! pgmin will hold the winning PE #
  
  do i = 1, nUnits
     if (dtModule(1,i) < dtNew) then
        dtNew    = dtModule(1,i)
        pgmin = dtModule(2,i)
        ngmin = i
     endif
  enddo
      
  ! have the processor that is determining the timestep limit broadcast the
  ! proc number, block number, and i,j,k of the zone that set the timestep
  ! to all processors

  dtMinLoc(:) = lminloc(:,ngmin)

  call MPI_Bcast(dtMinLoc(1), 5, MPI_INTEGER, pgmin, MPI_COMM_WORLD, error)

  ! limit the timestep to increase by at most a factor of dr_tstepChangeFactor

  dtNew = min( dtNew, dtOld*dr_tstepChangeFactor )
  
  if (nstep .GE. 50) then   !! This is where Cellular starts to fail
!     print *, 'nstep = ',nstep
  endif

  ! Use dr_dtmin and dr_dtmax to limit the timestep.  If this makes the code
  ! unstable, it's not our fault.
  dtNew = min( max( dtNew, dr_dtMin ), dr_dtMax )
  



  if (printTStepLoc) then
         
     ! convert the dtMinLoc array into a physical position (x,y,z) where the
     ! timestep is being set.  dtMinLoc(5) is the processor number, dtMinLoc(4)
     ! is the blockID on that proc.
     coords(:) = 0.0

     if (myPE == dtMinLoc(5)) then

        index(:)=dtMinLoc(1:MDIM)
        call Grid_getSingleCellCoords(index,dtMinLoc(4),CENTER, EXTERIOR,coords)

        ! send this to the master processor
        if (myPE /= MASTER_PE) then

           call MPI_Send (coords(1), 3, MPI_INTEGER, MASTER_PE, & 
                0, MPI_COMM_WORLD, error)
           
        endif
        
     elseif (myPE == MASTER_PE) then
        
        call MPI_Recv (coords(1), 3, MPI_INTEGER, dtMinLoc(5), 0, & 
             MPI_COMM_WORLD, status, error)            
        
     endif
     
  endif
  
  ! Print out the time and next timestep.
  
  ! only print out the timestep from the limiters that are active
  iout = 0
  do i = 1, nUnits
     if (dtModule(1,i) /= MAX_TSTEP) then
        iout = iout + 1
        tstepOutput(iout) = dtModule(1,i)
        limiterNameOutput(iout) = limiterName(i)
     endif
  enddo
  
  
  if (myPE == MASTER_PE) then
     
     if (printTStepLoc) then
        
        if (nstep == nbegin) then

           if (.not. dr_useRedshift) then
              write (*,803) 'n', 't', 'dt', 'x', 'y', 'z', &
                   (limiterNameOutput(i),i=1,iout)
           else
              write (*,804) 'n', 't', 'z', 'dt', 'x', 'y', 'z', &
                   (limiterNameOutput(i),i=1,iout)
           endif
              
        endif
        
        if (.not. dr_useRedshift) then
           write(*,801) nstep, simTime, dtNew, coords(1), coords(2), &
                coords(3), (tstepOutput(i),i=1,iout)
        else
           write(*,802) nstep, simTime, dr_redshift, dtNew, coords(1), &
                coords(2), coords(3), (tstepOutput(i),i=1,iout)
        endif

     else
        
        if (nstep .eq. nbegin) then
           
           if (.not. dr_useRedshift) then
              write (*,903) 'n', 't', 'dt', (limiterNameOutput(i),i=1,iout)
           else
              write (*,904) 'n', 't', 'z', 'dt', (limiterNameOutput(i),i=1,iout)
           endif
           
        endif
        
        if (.not. dr_useRedshift) then
           write(*,901) nstep, simTime, dtNew, (tstepOutput(i),i=1,iout)
        else
           write(*,902) nstep, simTime, dr_redshift, dtNew, &
                (tstepOutput(i),i=1,iout)
        endif

     endif

  end if
      
801 format (1X, I7, 1x, ES10.4, 1x, ES10.4, 2x, '(', ES9.3, ', ', &
            ES 9.3, ', ', ES9.3, ')', ' | ', 11(1X, :, ES9.3))
802 format (1X, I7, 1x, ES10.4, 1x, F8.3, 1x, ES10.4, 2x, '(', ES9.3, ', ', &
            ES 9.3, ', ', ES9.3, ')', ' | ', 11(1X, :, ES9.3))
803 format (1X, A7, 1x, A10, 1x, A10, 2x, '(', A9, ', ', A9, ', ', A9, ')', &
            ' | ', 11(1X, :, A9))
804 format (1X, A7, 1x, A10, 1x, A7, 1x, A10, 2x, '(', A9, ', ', A9, ', ', &
         A9, ')', ' | ', 11(1X, :, A9))

901 format (1X, I7, 1X, ES10.4, 1x, ES10.4, ' | ', 11(1X, :, ES11.5))
902 format (1X, I7, 1X, ES10.4, 1x, F8.3, 1x, ES10.4, ' | ', 11(1X, :, ES11.5))
903 format (1X, A7, 1x, A10   , 1x, A10,    ' | ', 11(1X, :, A11))
904 format (1X, A7, 1x, A10   , 1x, A7, 1x, A10,    ' | ', 11(1X, :, A11))

  return
end subroutine Driver_computeDt
