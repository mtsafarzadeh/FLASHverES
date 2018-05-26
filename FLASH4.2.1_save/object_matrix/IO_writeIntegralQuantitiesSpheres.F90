!!****if* source/Simulation/SimulationMain/StirTurbHydro/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities() 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4): solnData, scr

subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalMe
  use Grid_interface, ONLY : Grid_computeUserVars, Grid_getCellCoords, &
    Grid_getListOfBlocks, Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getSingleCellVol, Grid_getDeltas, Grid_getBlkData
  use Runtimeparameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Simulation_data, ONLY: tone,ttwo,tfour,teight,npdfstart, &
                             n_one,ntwo,nfour,neight,&
                             sim_rhoambient,sim_writematrix,&
                             steps_one,steps_two,steps_four,steps_eight, &
                             sim_init_sten, sim_rms_mach_target

  use Driver_data, ONLY: dr_nstep
  use Stir_data,   ONLY: st_energy,st_OUvar,st_decay

  implicit none

!!#include "mpif.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  

  real, intent(in) :: simTime
  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer :: size(MDIM)

  integer, parameter ::  nGlobalSum = 25   ! Number of globally-summed quantities
  real :: gsum(nGlobalSum)                 ! Global summed quantities
  real :: lsum(nGlobalSum)                 ! Local summed quantities

  integer :: i, j, k                       ! Counters over cells
  real :: dvol, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real, DIMENSION(:,:,:,:), POINTER :: scr

  integer :: point(MDIM)

  real, save :: sim_pi, sim_grv_const, sim_gamma
  logical :: gcell = .true.

!! 
!! SS : variables whose time series we output 
!!
  real :: lmin_dens, lmax_dens, gmin_dens, gmax_dens 
  real :: rms_velh_dw, rms_velv_dw, rms_vel_dw, rms_mach, rms_force
  real :: mean_pres, mean_temp
!! ES  :: this is for modying the kicks to maintain a mach number
  real :: sten

!! ES for computing the mach number at 2 4 and 8 grid scales
  real :: rms_mach2,rms_mach4,rms_mach8
  real :: ddd,dvx,dvy,dvz                  ! accumulate density, rho*vx, rho*vy, rho*vz

!! ES indicies for past masses
  integer :: ieul                          ! index over eularian past density variables
  integer :: imas                          ! index over lagrangian past density variables           
!! ES these are the pdfarry, local and global versions of the pdf and transition
  real :: lpdfone(201),      lpdftwo(201)
  real :: lpdffour(201),lpdfeight(201)
  real :: gpdfone(201),      gpdftwo(201)
  real :: gpdffour(201),gpdfeight(201)
  real :: pdfonetot,    pdfonetotS2,  pdfonetotS4,  pdfonetotS8
  real :: pdftwotot,    pdftwototS2,  pdftwototS4,  pdftwototS8
  real :: pdffourtot,   pdffourtotS2, pdffourtotS4, pdffourtotS8
  real :: pdfeighttot,  pdfeighttotS2,pdfeighttotS4,pdfeighttotS8

!! these are the arrays for the mean of delta s 
  real :: lsone(201),  lstwo(201)
  real :: lsfour(201), lseight(201)
  real :: gsone(201),  gstwo(201)
  real :: gsfour(201), gseight(201)
  real :: lsMone(201), lsMtwo(201)
  real :: lsMfour(201),lsMeight(201)
  real :: gsMone(201), gsMtwo(201)
  real :: gsMfour(201),gsMeight(201)

!! ES these are the arrays for deltas^2
  real :: ls2one(201),  ls2two(201)
  real :: ls2four(201), ls2eight(201)
  real :: gs2one(201),  gs2two(201)
  real :: gs2four(201), gs2eight(201)
  real :: ls2Mone(201), ls2Mtwo(201)
  real :: ls2Mfour(201),ls2Meight(201)
  real :: gs2Mone(201), gs2Mtwo(201)
  real :: gs2Mfour(201),gs2Meight(201)


!! ES these are smoothed versions of the pdfs
  real :: lpdfoneS(603), lpdftwoS(603)
  real :: lpdffourS(603),lpdfeightS(603)
  real :: gpdfoneS(603), gpdftwoS(603)
  real :: gpdffourS(603),gpdfeightS(603)

!! ES these are the arrays for smoothed version of 
!!   the mean of delta s 
  real :: lsoneS(603),  lstwoS(603)
  real :: lsfourS(603), lseightS(603)
  real :: gsoneS(603),  gstwoS(603)
  real :: gsfourS(603), gseightS(603)
  real :: lsMoneS(603), lsMtwoS(603)
  real :: lsMfourS(603),lsMeightS(603)
  real :: gsMoneS(603), gsMtwoS(603)
  real :: gsMfourS(603),gsMeightS(603)

!! ES these are the arrays for smoothed version of
!! the arrays for deltas^2
  real :: ls2oneS(603),  ls2twoS(603)
  real :: ls2fourS(603), ls2eightS(603)
  real :: gs2oneS(603),  gs2twoS(603)
  real :: gs2fourS(603), gs2eightS(603)
  real :: ls2MoneS(603), ls2MtwoS(603)
  real :: ls2MfourS(603),ls2MeightS(603)
  real :: gs2MoneS(603), gs2MtwoS(603)
  real :: gs2MfourS(603),gs2MeightS(603)
 

! ES these are arrays of the transition matrix
! 401 (in the delta s direction ) *
! 201 (in the s direction) =  80601
  real :: lTMone(80601),      lTMtwo(80601)
  real :: lTMfour(80601),    lTMeight(80601)
  real :: gTMone(80601),      gTMtwo(80601)
  real :: gTMfour(80601),    gTMeight(80601)
  real :: lMTMone(80601),     lMTMtwo(80601)
  real :: lMTMfour(80601),   lMTMeight(80601)
  real :: gMTMone(80601),     gMTMtwo(80601)
  real :: gMTMfour(80601),   gMTMeight(80601)
  real :: TMonetot, TMtwotot, TMfourtot, TMeighttot
  real :: MTMonetot,  MTMonetotS2,  MTMonetotS4,  MTMonetotS8
  real :: MTMtwotot,   MTMtwototS2,   MTMtwototS4,   MTMtwototS8
  real :: MTMfourtot,MTMfourtotS2,MTMfourtotS4,MTMfourtotS8
  real :: MTMeighttot,MTMeighttotS2,MTMeighttotS4,MTMeighttotS8

! ES these are the transition matrixes for the smoothed variables
! 3* 80601 =  241803
  real :: lMTMoneS(241803),  lMTMtwoS(241803)
  real :: lMTMfourS(241803),lMTMeightS(241803)
  real :: gMTMoneS(241803),  gMTMtwoS(241803)
  real :: gMTMfourS(241803),gMTMeightS(241803)
  real :: MTMonetotS,MTMtwototS,MTMfourtotS,MTMeighttotS
  real :: pdfonetotS,pdftwototS,pdffourtotS,pdfeighttotS

  real    :: sold,  soldm,  snew
  integer :: indnew, inddif, inddifm
  integer :: csmooth                 ! a counter over smoothing sizes
  integer :: istep                   ! radius of sphere in cells 
  integer :: istsq,isq,jsq,ksq       ! squared indexes
  integer :: offset1,offset2         ! offsets for smoothed pdf and transition matrix
  integer :: ip,jp,kp                ! these count over cells in the smoothing volume
  real    :: dsmoothnew,dsmoothold,dsmootholdm ! smoothed densities
  real    :: totzones                ! total number of zones

! ES  this saves the smoothed values of the densities 
! I just make this one 80*80*80, if we use smaller blocks then it is bigger
! than needed
  character (len=4) ::  fnumStr
  character (len=MAX_STRING_LENGTH) :: TMFileName

! ES  This is the multiplier for the transition matrix
  real    ::  TM_multiplier

  TM_multiplier = 20


  call RuntimeParameters_get('gamma', sim_gamma)
  call PhysicalConstants_get( 'Newton', sim_grv_const)               !! Gravitational constant G
  call PhysicalConstants_get( 'Pi', sim_pi)                          !! Value of pi

  ! Make sure the vorticity is up-to-date
  call Grid_computeUserVars()

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  gmin_dens = 1.0
  gmax_dens = 1.0
  lmin_dens = 1.0
  lmax_dens = 1.0

  call Grid_getListOfBlocks(LEAF, blockList, count)
  

!! SS : start looping over the blocks 

  do lb = 1, count

     ! get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)


     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
             
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 

              lmin_dens = min(lmin_dens, solnData(DENS_VAR,i,j,k))
              lmax_dens = max(lmin_dens, solnData(DENS_VAR,i,j,k))

#endif           
              ! Total volume 
              lsum(2) = lsum(2) + dvol 

              ! Pressure 
              lsum(3) = lsum(3) + solnData(PRES_VAR,i,j,k)*dvol 


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! x-momentum
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              ! y-momentum
              lsum(5) = lsum(5) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              ! z-momentum
              lsum(6) = lsum(6) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol

              !! horizontal velocity dispersion squared - density weighted
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k)*(solnData(VELX_VAR,i,j,k)**2 + solnData(VELY_VAR,i,j,k)**2)*dvol

              !! vertical velocity dispersion squared - density weighted 

               lsum(8) = lsum(8) + solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)**2*dvol 
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(9) = lsum(9) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol 
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy - density weighted 
              lsum(10) = lsum(10) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol     

              ! density weighted velocity squared 

              lsum(11) = lsum(11) + solnData(DENS_VAR,i,j,k)*(solnData(VELX_VAR,i,j,k)**2 + & 
                                                              solnData(VELY_VAR,i,j,k)**2+solnData(VELZ_VAR,i,j,k)**2)*dvol      

     
              ! Sound speed squared 

              lsum(12) = lsum(12) + sim_gamma*solnData(PRES_VAR,i,j,k)*dvol
#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(13) = lsum(13) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

               ! Quantities related to the artificial forcing

               ! Forcing squared 

               lsum(14) = lsum(14) + (solnData(ACCX_VAR,i,j,k)**2 + solnData(ACCY_VAR,i,j,k)**2 + solnData(ACCZ_VAR,i,j,k)**2)*dvol

               !  rho*accx
               lsum(15) = lsum(15) + solnData(DENS_VAR,i,j,k)* &
                      &                         solnData(ACCX_VAR,i,j,k)*dvol

              ! rho*accy
               lsum(16) = lsum(16) + solnData(DENS_VAR,i,j,k)* &
                        &                       solnData(ACCY_VAR,i,j,k)*dvol

              ! rho*accz
               lsum(17) = lsum(17) + solnData(DENS_VAR,i,j,k)* &
                         &                      solnData(ACCZ_VAR,i,j,k)*dvol

              ! rho * accx *velx
               lsum(18) = lsum(18) + 2.0*solnData(DENS_VAR,i,j,k)* &
                          solnData(VELX_VAR,i,j,k)*solnData(ACCX_VAR,i,j,k)*dvol

              ! rho * accy * vely
               lsum(19) = lsum(19) + 2.0*solnData(DENS_VAR,i,j,k)* &
                          solnData(VELY_VAR,i,j,k)*solnData(ACCY_VAR,i,j,k)*dvol 

              ! rho * accz * velz 
               lsum(20) = lsum(20) + 2.0*solnData(DENS_VAR,i,j,k)* &
                          solnData(VELZ_VAR,i,j,k)*solnData(ACCZ_VAR,i,j,k)*dvol 
 
              ! Temperature 
               lsum(21) = lsum(21) + solnData(TEMP_VAR,i,j,k)*dvol           


           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)
  enddo     !! SS : end of the loop over all the blocks 
  do lb = 1, count

   ! get a pointer to the current block of data
   call Grid_getBlkPtr(blockList(lb), solnData)

  !! SS : Calculate the sum of vorticity over all cells in the different blocks
     lsum(22) = 0.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k
              lsum(22) = lsum(22) + (solnData(MVRT_SCRATCH_CENTER_VAR,i,j,k) * dvol)/2
           end do
        end do
     end do

 !   ES compute the rms velocity at the size of 2 cells,4, and 8 cells
     do csmooth = 1,3 
       istep = 2**(csmooth)
       lsum(22+csmooth) = 0.
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS),istep
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS),istep
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS),istep
                ddd  = 0.
                dvx = 0.
                dvy = 0.
                dvz = 0.
                do kp=k,k+istep-1
                   do jp=j,j+istep-1
                      do ip=i,i+istep-1
                        ddd = ddd+solnData(DENS_VAR,ip,jp,kp) 
                        dvx = dvx+solnData(DENS_VAR,ip,jp,kp)*solnData(VELX_VAR,ip,jp,kp) 
                        dvy = dvy+solnData(DENS_VAR,ip,jp,kp)*solnData(VELY_VAR,ip,jp,kp) 
                        dvz = dvz+solnData(DENS_VAR,ip,jp,kp)*solnData(VELZ_VAR,ip,jp,kp) 
                     end do
                   end do
                end do
                lsum(22+csmooth) = lsum(22+csmooth) + (dvx*dvx+dvy*dvy+dvz*dvz)*dvol/ddd
             end do
            end do
       end do
     end do

     call Grid_releaseBlkPtr(blockList(lb), solnData)
  enddo     !! SS : end of the loop over all the blocks 


  ! we have to do an MPI call to gather it together
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, &
       &                MASTER_PE, MPI_Comm_World, error)


  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)
  
!! SS : Reduce calls for the min and max density
  call MPI_Reduce (lmin_dens, gmin_dens, 1, MPI_Double_Precision, MPI_Min, &
       &                MASTER_PE, MPI_Comm_World, error)

  call MPI_Reduce (lmax_dens, gmax_dens, 1, MPI_Double_Precision, MPI_Max, &
       &                MASTER_PE, MPI_Comm_World, error)


  if (io_globalMe  == MASTER_PE) then

                 mean_pres     = gsum(3)/gsum(2)
                 mean_temp     = gsum(21)/gsum(2)
                 rms_velh_dw   = sqrt(gsum(7)/gsum(1))     
                 rms_velv_dw   = sqrt(gsum(8)/gsum(1))     
                 rms_vel_dw    = sqrt(gsum(11)/gsum(1))
                 rms_mach      = sqrt(gsum(11)/gsum(12))
                 rms_mach2     = sqrt(gsum(23)/gsum(12))
                 rms_mach4     = sqrt(gsum(24)/gsum(12))
                 rms_mach8     = sqrt(gsum(25)/gsum(12))
                 rms_force     = sqrt(gsum(14)/gsum(2))

     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time(1)                    ', &
                'mass(2)                      ', &
                'x-momentum(3)                ', &
                'y-momentum(4)                ', & 
                'z-momentum(5)                ', &
                'E_total(6)                   ', &
                'E_kinetic(7)                 ', &
                'E_internal(8)                ', &
                'rms_vel_horz(9)              ', &
                'rms_vel_vert(10)              ', &
                'rms_vel_dw(11)                ', &
                'rms_mach(12)                  ', &
                'rms_mach2(13)                 ', &
                'rms_mach4(14)                 ', &
                'rms_mach8(15)                 ', &
                'rms_force(16)                 ', &
                'mean_pres(17)                 ', &
                'mean_temp(18)                 ', &
                'Min. dens(19)                 ', &
                'Max. dens(20)                 ', &
                'Mag_vorticity(21)             '
           
10         format (2x,50(a25, :, 1X))
           
        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12)     & ! Write the global sums to the file.
              simtime,     &                 ! Time 
              gsum(1),     &                 ! Mass 
              gsum(4),     &                 ! x-momentum 
              gsum(5),     &                 ! y-momentum 
              gsum(6),     &                 ! z-momentum 
              gsum(9),     &                 ! E_total
              gsum(10),    &                 ! E_kin 
              gsum(13),    &                 ! E_int
              rms_velh_dw, &                 ! rms value of horz. vel. 
              rms_velv_dw, &                 ! rms value of vert. vel. 
              rms_vel_dw,  &                 ! rms velocity - 3d 
              rms_mach,    &                 ! rms Mach number 
              rms_mach2,   &                 ! rms Mach number smoothed 2 cells
              rms_mach4,   &                 ! rms Mach number smoothed 4 cells
              rms_mach8,   &                 ! rms Mach number smoothed 8 cells
              rms_force,   &                 ! rms value of the forcing 
              mean_pres,   &                 ! Mean pressure 
              mean_temp,   &                 ! Mean temperature
              gmin_dens,   &                 ! Min. density
              gmax_dens,   &                 ! Max. density
              gsum(21)                       ! Mag. of vorticity

! if we are past the initial startup phase 
! then we modify the strength of the kicks to maintain
! the mach number we are seeking             
              sten = sim_init_sten
              if(dr_nstep.ge.npdfstart)  then
                sten = sim_init_sten*(sim_rms_mach_target/rms_mach)**10
              endif
              st_energy = sten
              st_OUvar = sqrt(st_energy/st_decay) 
              print *,'mach',sim_rms_mach_target,rms_mach
              print *,'sten',st_energy,st_OUvar,sim_init_sten
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!ES
!!
!!  Begin part 2 of the routine writing out the data on the transition
!!  matrixes, shifts, and pdfs (repeats 4 times the same code) 
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call MPI_Barrier (MPI_Comm_World, error)
  call MPI_Bcast(st_energy ,1,FLASH_REAL,MASTER_PE,MPI_Comm_World,error)
  call MPI_Bcast(st_OUvar  ,1,FLASH_REAL,MASTER_PE,MPI_Comm_World,error)
  call MPI_Barrier(MPI_Comm_World,error) 
    

!! ES in this first case the delay is all one time step
!! if it has been 1 cycle since the last time we wrote the pdf file
!! then we write the pdf file again 
  if((dr_nstep .ge. steps_one*n_one).and.(dr_nstep .ge. npdfstart)) then

    ! initialize data
    lpdfone  = 0.
    lsone    = 0.
    ls2one   = 0. 
    lsMone   = 0. 
    ls2Mone  = 0.
    lTMone   = 0.
    lMTMone  = 0.

    ! initialize smoothed versions of the data
    lpdfoneS = 0.
    lsoneS   = 0.
    ls2oneS  = 0.
    lsMoneS  = 0.
    ls2MoneS = 0.
    lMTMoneS = 0.

    !loop over the block and accumulate the transition matrix
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DONE_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDON_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*TM_multiplier +201
             inddifm = (snew-soldm)*TM_multiplier+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdfone(indnew) = lpdfone(indnew)+1.
                 lsone(indnew)   = lsone(indnew)+(snew-sold)
                 ls2one(indnew)  = ls2one(indnew)+(snew-sold)*(snew-sold)
                 lsMone(indnew)  = lsMone(indnew)+(snew-soldm)
                 ls2Mone(indnew) = ls2Mone(indnew)+(snew-soldm)*(snew-soldm)
             endif
             if(((inddif.ge.1).and.(indnew.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lTMone((inddif-1)*201+indnew)   = lTMone((inddif-1)*201+indnew)+1.
             endif
             if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lMTMone((inddifm-1)*201+indnew)  = lMTMone((inddifm-1)*201+indnew)+1.
             endif
             solnData(DONE_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDON_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do

!   now do the case where you smooth over spheres of various sizes
     do csmooth = 1,3
       istep = 1+csmooth
       istsq = istep*istep
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)

       if (csmooth.eq.1) then
         ieul = TWON_VAR
         imas = MTWO_MSCALAR
       endif
       if (csmooth.eq.2) then
         ieul = THON_VAR
         imas = MTHO_MSCALAR
       endif
       if (csmooth.eq.3) then
         ieul = FRON_VAR
         imas = MFRO_MSCALAR
       endif

       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!           compute the new and old smoothed densities 
              dsmoothnew  = 0.
              dsmoothold  = 0.
              dsmootholdm = 0.
              totzones    = 0.
              do kp=k-istep,k+istep
                ksq = (kp-k)*(kp-k)
                do jp=j-istep,j+istep
                  jsq = (jp-j)*(jp-j)
                  do ip=i-istep,i+istep
                    isq = (ip-i)*(ip-i)
                    if(ksq+jsq+isq.le.istsq) then
                      dsmoothnew  = dsmoothnew +solnData(DENS_VAR,ip,jp,kp)
                      totzones    = totzones + 1
                    endif
                  end do
                end do
              end do
              dsmoothold  = solnData(ieul,i,j,k)
              dsmootholdm = solnData(imas,i,j,k)
              dsmoothnew  = dsmoothnew/totzones
              solnData(ieul,i,j,k)     = dsmoothnew
              solnData(imas,i,j,k)     = dsmoothnew

! add to the pdf and the transition matrix
              snew  = alog(dsmoothnew/ sim_rhoAmbient)
              sold  = alog(dsmoothold/ sim_rhoAmbient)
              soldm = alog(dsmootholdm/sim_rhoAmbient)
              indnew  = snew *20+101
              inddif  = (snew-sold)*TM_multiplier+201
              inddifm = (snew-soldm)*TM_multiplier+201
              if((indnew.ge.1).and.(indnew.le.201)) then
                  lpdfoneS(offset1+indnew)         = lpdfoneS(offset1+indnew)+1.
                  lsoneS(offset1+indnew)           = lsoneS  (offset1+indnew)+(snew-sold)
                  ls2oneS(offset1+indnew)          = ls2oneS (offset1+indnew)+(snew-sold)*(snew-sold)
                  lsMoneS(offset1+indnew)          = lsMoneS (offset1+indnew)+(snew-soldm)
                  ls2MoneS(offset1+indnew)         = ls2MoneS(offset1+indnew)+(snew-soldm)*(snew-soldm)
              endif
              if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                  lMTMoneS(offset2+(inddifm-1)*201+indnew)  = lMTMoneS(offset2+(inddifm-1)*201+indnew)+1.
              endif
             end do
          end do
        end do
      end do 
      call Grid_releaseBlkPtr(blockList(lb), solnData)
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdfone,  gpdfone, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdfoneS, gpdfoneS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsone,   gsone, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMone,  gsMone, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsoneS, gsoneS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMoneS, gsMoneS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2one,  gs2one, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mone, gs2Mone, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2oneS, gs2oneS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2MoneS, gs2MoneS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMone,  gTMone,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMone, gMTMone, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMoneS, gMTMoneS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)


    ! nomalize the pdf and the transision matrix
    if (io_globalMe  == MASTER_PE) then
      pdfonetot    = 0.
      pdfonetotS2  = 0.
      pdfonetotS4  = 0.
      pdfonetotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMonetot     = 0.
        MTMonetot    = 0.
        MTMonetotS2  = 0.
        MTMonetotS4  = 0.
        MTMonetotS8  = 0.
        do inddif=1,401
         TMonetot  = TMonetot +gTMone((inddif-1)*201+indnew)
         MTMonetot = MTMonetot+gMTMone((inddif-1)*201+indnew)
         MTMonetotS2 = MTMonetotS2+gMTMoneS((inddif-1)*201+indnew)
         MTMonetotS4 = MTMonetotS4+gMTMoneS(80601+(inddif-1)*201+indnew)
         MTMonetotS8 = MTMonetotS8+gMTMoneS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMonetot    .eq. 0.) TMonetot  = 1.
        if(MTMonetot   .eq. 0.) MTMonetot = 1.
        if(MTMonetotS2 .eq. 0.) MTMonetotS2 = 1.
        if(MTMonetotS4 .eq. 0.) MTMonetotS4 = 1.
        if(MTMonetotS8 .eq. 0.) MTMonetotS8 = 1.
        do inddif=1,401
         gTMone((inddif-1)*201+indnew)= & 
           gTMone((inddif-1)*201+indnew)/TMonetot
         gMTMone((inddif-1)*201+indnew)= & 
           gMTMone((inddif-1)*201+indnew)/MTMonetot
         gMTMoneS((inddif-1)*201+indnew)= & 
           gMTMoneS((inddif-1)*201+indnew)/MTMonetotS2
         gMTMoneS(80601+(inddif-1)*201+indnew)= & 
           gMTMoneS(80601+(inddif-1)*201+indnew)/MTMonetotS4
         gMTMoneS(161202+(inddif-1)*201+indnew)= & 
           gMTMoneS(161202+(inddif-1)*201+indnew)/MTMonetotS8
        enddo
       pdfonetot = pdfonetot+gpdfone(indnew)
       pdfonetotS2 = pdfonetotS2+gpdfoneS(indnew)
       pdfonetotS4 = pdfonetotS4+gpdfoneS(201+indnew)
       pdfonetotS8 = pdfonetotS8+gpdfoneS(402+indnew)
      enddo 

      ! write the file
      if((sim_writematrix).or.(mod(n_one,16).eq.1)) then
        write (fnumStr, '(i4.4)') (n_one-npdfstart/steps_one)
        TMFileName = 'TM1_' // fnumStr // '.dat'
        open (funit, file=trim(TMFileName))
        write(funit,21) dr_nstep 
        write(funit,20) simTime
        write(funit,20) simTime-tone
        write(funit,20) rms_mach
        write(funit,20) rms_mach2
        write(funit,20) rms_mach4
        write(funit,20) rms_mach8
!       write the pdfs
        do indnew=1,201
         write(funit,20) gpdfone(indnew)/pdfonetot
        enddo
        do indnew=1,201
         write(funit,20) gpdfoneS(indnew)/pdfonetotS2
        enddo
        do indnew=202,402
         write(funit,20) gpdfoneS(indnew)/pdfonetotS4
        enddo
        do indnew=403,603
         write(funit,20) gpdfoneS(indnew)/pdfonetotS8
        enddo
!       write the transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gMTMone((inddif-1)*201+indnew)
          enddo
        enddo 
!       write the smoothed transition matrix
        do inddif=1,1203
          do indnew=1,201
            write(funit,20) gMTMoneS((inddif-1)*201+indnew)
          enddo
        enddo
!       write the eulerian transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gTMone((inddif-1)*201+indnew)
          enddo
        enddo
        close(funit)
      endif !! Close the file that contains the transition matrix, if we are
            !! writing it. 

!     write the file with s_s2 in it
      write (fnumStr, '(i4.4)') (n_one-npdfstart/steps_one)
      TMFileName = 's_s2_1_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-tone
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdfone(indnew)/pdfonetot
       gpdfone(indnew) = gpdfone(indnew) + 1E-10
      enddo
      do indnew=1,201
       write(funit,20) gpdfoneS(indnew)/pdfonetotS2
       gpdfoneS(indnew) = gpdfoneS(indnew) + 1E-10
      enddo
      do indnew=202,402
       write(funit,20) gpdfoneS(indnew)/pdfonetotS4
       gpdfoneS(indnew) = gpdfoneS(indnew) + 1E-10
      enddo
      do indnew=403,603
       write(funit,20) gpdfoneS(indnew)/pdfonetotS8
       gpdfoneS(indnew) = gpdfoneS(indnew) + 1E-10
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gsone(indnew)/pdfonetot
      enddo
      do indnew=1,201
       write(funit,20) gsoneS(indnew)/pdfonetotS2
      enddo
      do indnew=202,402
       write(funit,20) gsoneS(indnew)/pdfonetotS4
      enddo
      do indnew=403,603
       write(funit,20) gsoneS(indnew)/pdfonetotS8
      enddo
      do indnew=1,201
       write(funit,20) gsMone(indnew)/pdfonetot
      enddo
      do indnew=1,201
       write(funit,20) gsMoneS(indnew)/pdfonetotS2
      enddo
      do indnew=202,402
       write(funit,20) gsMoneS(indnew)/pdfonetotS4
      enddo
      do indnew=403,603
       write(funit,20) gsMoneS(indnew)/pdfonetotS8
      enddo
!     write the s2 values
      do indnew=1,201
       write(funit,20) (gs2one(indnew)-gsone(indnew)**2/gpdfone(indnew))/pdfonetot
      enddo
      do indnew=1,201
       write(funit,20) (gs2oneS(indnew)-gsoneS(indnew)**2/gpdfoneS(indnew))/pdfonetotS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2oneS(indnew)-gsoneS(indnew)**2/gpdfoneS(indnew))/pdfonetotS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2oneS(indnew)-gsoneS(indnew)**2/gpdfoneS(indnew))/pdfonetotS8
      enddo
      do indnew=1,201
       write(funit,20) (gs2Mone(indnew)-gsMone(indnew)**2/gpdfone(indnew))/pdfonetot
      enddo
      do indnew=1,201
       write(funit,20) (gs2MoneS(indnew)-gsMoneS(indnew)**2/gpdfoneS(indnew))/pdfonetotS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2MoneS(indnew)-gsMoneS(indnew)**2/gpdfoneS(indnew))/pdfonetotS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2MoneS(indnew)-gsMoneS(indnew)**2/gpdfoneS(indnew))/pdfonetotS8
      enddo
      close(funit)


    endif
    n_one = n_one+1
    tone = simTime
  endif
 
20  format (' ',E15.7)
21  format (' ',I5.2)
  !! ES if it has been 2 cycles since the last time we wrote the pdf file
  !! then we write the pdf file again 
  if((dr_nstep .ge. steps_two*ntwo).and.(dr_nstep .ge. npdfstart)) then

    ! initialize data
    lpdftwo  = 0.
    lstwo    = 0.
    ls2two   = 0.
    lsMtwo   = 0.
    ls2Mtwo  = 0.
    lTMtwo   = 0.
    lMTMtwo  = 0.

    ! initialize smoothed versions of the data
    lpdftwoS = 0.
    lstwoS   = 0.
    ls2twoS  = 0.
    lsMtwoS  = 0.
    ls2MtwoS = 0.
    lMTMtwoS = 0.

    !loop over the block and accumulate the transition matrix
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DTWO_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDTW_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*TM_multiplier +201
             inddifm = (snew-soldm)*TM_multiplier+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdftwo(indnew) = lpdftwo(indnew)+1.
                 lstwo(indnew)   = lstwo(indnew)+(snew-sold)
                 ls2two(indnew)  = ls2two(indnew)+(snew-sold)*(snew-sold)
                 lsMtwo(indnew)  = lsMtwo(indnew)+(snew-soldm)
                 ls2Mtwo(indnew) = ls2Mtwo(indnew)+(snew-soldm)*(snew-soldm)
             endif
             if(((inddif.ge.1).and.(inddif.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lTMtwo((inddif-1)*201+indnew)   = lTMtwo((inddif-1)*201+indnew)+1.
             endif
             if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lMTMtwo((inddifm-1)*201+indnew)  = lMTMtwo((inddifm-1)*201+indnew)+1.
             endif
             solnData(DTWO_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDTW_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do

!   now do the case where you smooth over spheres of various sizes
     do csmooth = 1,3
       istep = 1+csmooth
       istsq = istep*istep
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)

       if (csmooth.eq.1) then
         ieul = TWTW_VAR
         imas = MTWT_MSCALAR
       endif
       if (csmooth.eq.2) then
         ieul = THTW_VAR
         imas = MTHT_MSCALAR
       endif
       if (csmooth.eq.3) then
         ieul = FRTW_VAR
         imas = MFRT_MSCALAR
       endif

       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!           compute the new and old averaged densities 
              dsmoothnew  = 0.
              dsmoothold  = 0.
              dsmootholdm = 0.
              totzones    = 0.
              do kp=k-istep,k+istep
                ksq = (kp-k)*(kp-k)
                do jp=j-istep,j+istep
                  jsq = (jp-j)*(jp-j)
                  do ip=i-istep,i+istep
                    isq = (ip-i)*(ip-i)
                    if(ksq+jsq+isq.le.istsq) then
                      dsmoothnew  = dsmoothnew +solnData(DENS_VAR,ip,jp,kp)
                      totzones    = totzones + 1
                    endif
                  end do
                end do
              end do
              dsmoothold  = solnData(ieul,i,j,k)
              dsmootholdm = solnData(imas,i,j,k)
              dsmoothnew  = dsmoothnew/totzones
              solnData(ieul,i,j,k)     = dsmoothnew
              solnData(imas,i,j,k)     = dsmoothnew

! add to the pdf and the transition matrix
              snew  = alog(dsmoothnew/ sim_rhoAmbient)
              sold  = alog(dsmoothold/ sim_rhoAmbient)
              soldm = alog(dsmootholdm/sim_rhoAmbient)
              indnew  = snew *20+101
              inddif  = (snew-sold)*TM_multiplier+201
              inddifm = (snew-soldm)*TM_multiplier+201
              if((indnew.ge.1).and.(indnew.le.201)) then
                  lpdftwoS(offset1+indnew)         = lpdftwoS(offset1+indnew)+1.
                  lstwoS(offset1+indnew)           = lstwoS(offset1+indnew)+(snew-sold)
                  ls2twoS(offset1+indnew)          = ls2twoS(offset1+indnew)+(snew-sold)*(snew-sold)
                  lsMtwoS(offset1+indnew)          = lsMtwoS(offset1+indnew)+(snew-soldm)
                  ls2MtwoS(offset1+indnew)         = ls2MtwoS(offset1+indnew)+(snew-soldm)*(snew-soldm)
              endif
              if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                  lMTMtwoS(offset2+(inddifm-1)*201+indnew)  = lMTMtwoS(offset2+(inddifm-1)*201+indnew)+1.
              endif
            end do
          end do
       end do
      end do 
      call Grid_releaseBlkPtr(blockList(lb), solnData)
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdftwo,  gpdftwo, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdftwoS, gpdftwoS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lstwo,   gstwo, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMtwo,  gsMtwo, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lstwoS, gstwoS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMtwoS, gsMtwoS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2two,  gs2two, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mtwo, gs2Mtwo, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2twoS, gs2twoS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2MtwoS, gs2MtwoS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMtwo,  gTMtwo,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMtwo, gMTMtwo, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMtwoS, gMTMtwoS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)



    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdftwotot    = 0.
      pdftwototS2  = 0.
      pdftwototS4  = 0.
      pdftwototS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMtwotot     = 0.
        MTMtwotot    = 0.
        MTMtwototS2  = 0.
        MTMtwototS4  = 0.
        MTMtwototS8  = 0.
        do inddif =1,401
         TMtwotot  = TMtwotot +gTMtwo((inddif-1)*201+indnew)
         MTMtwotot = MTMtwotot+gMTMtwo((inddif-1)*201+indnew)
         MTMtwototS2 = MTMtwototS2+gMTMtwoS((inddif-1)*201+indnew)
         MTMtwototS4 = MTMtwototS4+gMTMtwoS(80601+(inddif-1)*201+indnew)
         MTMtwototS8 = MTMtwototS8+gMTMtwoS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMtwotot    .eq. 0.) TMtwotot  = 1.
        if(MTMtwotot   .eq. 0.) MTMtwotot = 1.
        if(MTMtwototS2 .eq. 0.) MTMtwototS2 = 1.
        if(MTMtwototS4 .eq. 0.) MTMtwototS4 = 1.
        if(MTMtwototS8 .eq. 0.) MTMtwototS8 = 1.
        do inddif=1,401
         gTMtwo((inddif-1)*201+indnew)= & 
           gTMtwo((inddif-1)*201+indnew)/TMtwotot
         gMTMtwo((inddif-1)*201+indnew)= & 
           gMTMtwo((inddif-1)*201+indnew)/MTMtwotot
         gMTMtwoS((inddif-1)*201+indnew)= & 
           gMTMtwoS((inddif-1)*201+indnew)/MTMtwototS2
         gMTMtwoS(80601+(inddif-1)*201+indnew)= & 
           gMTMtwoS(80601+(inddif-1)*201+indnew)/MTMtwototS4
         gMTMtwoS(161202+(inddif-1)*201+indnew)= & 
           gMTMtwoS(161202+(inddif-1)*201+indnew)/MTMtwototS8
        enddo
       pdftwotot = pdftwotot+gpdftwo(indnew)
       pdftwototS2 = pdftwototS2+gpdftwoS(indnew)
       pdftwototS4 = pdftwototS4+gpdftwoS(201+indnew)
       pdftwototS8 = pdftwototS8+gpdftwoS(402+indnew)
      enddo 

      ! write the file
      if((sim_writematrix).or.(mod(ntwo,8).eq.1)) then
        write (fnumStr, '(i4.4)') (ntwo-npdfstart/steps_two)
        TMFileName = 'TM2_' // fnumStr // '.dat'
        open (funit, file=trim(TMFileName))
        write(funit,21) dr_nstep 
        write(funit,20) simTime
        write(funit,20) simTime-ttwo
        write(funit,20) rms_mach
        write(funit,20) rms_mach2
        write(funit,20) rms_mach4
        write(funit,20) rms_mach8
     !     write the pdfs
        do indnew=1,201
          write(funit,20) gpdftwo(indnew)/pdftwotot
        enddo
        do indnew=1,201
          write(funit,20) gpdftwoS(indnew)/pdftwototS2
        enddo
        do indnew=202,402
          write(funit,20) gpdftwoS(indnew)/pdftwototS4
        enddo
        do indnew=403,603
          write(funit,20) gpdftwoS(indnew)/pdftwototS8
        enddo
!       write the transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gMTMtwo((inddif-1)*201+indnew)
          enddo
        enddo 
!       write the smoothed transition matrix
        do inddif=1,1203
          do indnew=1,201
            write(funit,20) gMTMtwoS((inddif-1)*201+indnew)
          enddo
        enddo
!       write the eulerian transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gTMtwo((inddif-1)*201+indnew)
          enddo
        enddo
        close(funit)
      endif

!     write the file with s and s2 in it
      write (fnumStr, '(i4.4)') (ntwo-npdfstart/steps_two)
      TMFileName = 's_s2_2_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-ttwo
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
        write(funit,20) gpdftwo(indnew)/pdftwotot
        gpdftwo(indnew) = gpdftwo(indnew) + 1E-10
      enddo
      do indnew=1,201
        write(funit,20) gpdftwoS(indnew)/pdftwototS2
        gpdftwoS(indnew) = gpdftwoS(indnew) + 1E-10
      enddo
      do indnew=202,402
        write(funit,20) gpdftwoS(indnew)/pdftwototS4
        gpdftwoS(indnew) = gpdftwoS(indnew) + 1E-10
      enddo
      do indnew=403,603
        write(funit,20) gpdftwoS(indnew)/pdftwototS8
        gpdftwoS(indnew) = gpdftwoS(indnew) + 1E-10
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gstwo(indnew)/pdftwotot
      enddo
      do indnew=1,201
       write(funit,20) gstwoS(indnew)/pdftwototS2
      enddo
      do indnew=202,402
       write(funit,20) gstwoS(indnew)/pdftwototS4
      enddo
      do indnew=403,603
       write(funit,20) gstwoS(indnew)/pdftwototS8
      enddo
      do indnew=1,201
       write(funit,20) gsMtwo(indnew)/pdftwotot
      enddo
      do indnew=1,201
       write(funit,20) gsMtwoS(indnew)/pdftwototS2
      enddo
      do indnew=202,402
       write(funit,20) gsMtwoS(indnew)/pdftwototS4
      enddo
      do indnew=403,603
       write(funit,20) gsMtwoS(indnew)/pdftwototS8
      enddo
      do indnew=1,201
       write(funit,20) (gs2two(indnew)-gstwo(indnew)**2/gpdftwo(indnew))/pdftwotot
      enddo
      do indnew=1,201
       write(funit,20) (gs2twoS(indnew)-gstwoS(indnew)**2/gpdftwoS(indnew))/pdftwototS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2twoS(indnew)-gstwoS(indnew)**2/gpdftwoS(indnew))/pdftwototS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2twoS(indnew)-gstwoS(indnew)**2/gpdftwoS(indnew))/pdftwototS8
      enddo
      do indnew=1,201
       write(funit,20) (gs2Mtwo(indnew)-gsMtwo(indnew)**2/gpdftwo(indnew))/pdftwotot
      enddo
      do indnew=1,201
       write(funit,20) (gs2MtwoS(indnew)-gsMtwoS(indnew)**2/gpdftwoS(indnew))/pdftwototS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2MtwoS(indnew)-gsMtwoS(indnew)**2/gpdftwoS(indnew))/pdftwototS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2MtwoS(indnew)-gsMtwoS(indnew)**2/gpdftwoS(indnew))/pdftwototS8
      enddo
      close(funit)

    endif
    ntwo = ntwo+1
    ttwo = simTime
  endif

  !! ES if it has been 4 then cycles since the last time we wrote the pdf file
  !! then we write the pdf file again 
  if((dr_nstep .ge. steps_four*nfour).and.(dr_nstep .ge. npdfstart)) then

    ! initialize data
    lpdffour  = 0.
    lsfour    = 0.
    ls2four   = 0.
    lsMfour   = 0.
    ls2Mfour  = 0.
    lTMfour   = 0.
    lMTMfour  = 0.

    ! initialize smoothed versions of the data
    lpdffourS = 0.
    lsfourS   = 0.
    ls2fourS  = 0.
    lsMfourS  = 0.
    ls2MfourS = 0.
    lMTMfourS = 0.

    !loop over the block and accumulate the transition matrix
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DFOR_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDFR_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*TM_multiplier +201
             inddifm = (snew-soldm)*TM_multiplier+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdffour(indnew) = lpdffour(indnew)+1.
                 lsfour(indnew)   = lsfour(indnew)+(snew-sold)
                 ls2four(indnew)  = ls2four(indnew)+(snew-sold)*(snew-sold)
                 lsMfour(indnew)  = lsMfour(indnew)+(snew-soldm)
                 ls2Mfour(indnew) = ls2Mfour(indnew)+(snew-soldm)*(snew-soldm)
             endif
             if(((inddif.ge.1).and.(inddif.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lTMfour((inddif-1)*201+indnew)   = lTMfour((inddif-1)*201+indnew)+1.
             endif
             if(((inddifm.ge.1).and.(inddifm.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lMTMfour((inddifm-1)*201+indnew)  = lMTMfour((inddifm-1)*201+indnew)+1.
             endif
             solnData(DFOR_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDFR_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do

!   now do the case where you smooth over spheres of various sizes
     do csmooth = 1,3
       istep = 1+csmooth
       istsq = istep*istep
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)

       if (csmooth.eq.1) then
         ieul = TWFR_VAR
         imas = MTWF_MSCALAR
       endif
       if (csmooth.eq.2) then
         ieul = THFR_VAR
         imas = MTHF_MSCALAR
       endif
       if (csmooth.eq.3) then
         ieul = FRFR_VAR
         imas = MFRF_MSCALAR
       endif

       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!           compute the new and old averaged densities
              dsmoothnew  = 0.
              dsmoothold  = 0.
              dsmootholdm = 0.
              totzones    = 0.
              do kp=k-istep,k+istep
                ksq = (kp-k)*(kp-k)
                do jp=j-istep,j+istep
                  jsq = (jp-j)*(jp-j)
                  do ip=i-istep,i+istep
                    isq = (ip-i)*(ip-i)
                    if(ksq+jsq+isq.le.istsq) then
                      dsmoothnew  = dsmoothnew +solnData(DENS_VAR,ip,jp,kp)
                      totzones    = totzones + 1
                    endif
                  end do
                end do
              end do
              dsmoothold  = solnData(ieul,i,j,k)
              dsmootholdm = solnData(imas,i,j,k)
              dsmoothnew  = dsmoothnew/totzones
              solnData(ieul,i,j,k)     = dsmoothnew
              solnData(imas,i,j,k)     = dsmoothnew

! add to the pdf and the transition matrix
              snew  = alog(dsmoothnew/ sim_rhoAmbient)
              sold  = alog(dsmoothold/ sim_rhoAmbient)
              soldm = alog(dsmootholdm/sim_rhoAmbient)
              indnew  = snew *20+101
              inddif  = (snew-sold)*TM_multiplier+201
              inddifm = (snew-soldm)*TM_multiplier+201
              if((indnew.ge.1).and.(indnew.le.201)) then
                  lpdffourS(offset1+indnew)       = lpdffourS(offset1+indnew)+1.
                  lsfourS(offset1+indnew)         = lsfourS(offset1+indnew)+(snew-sold)
                  ls2fourS(offset1+indnew)        = ls2fourS(offset1+indnew)+(snew-sold)*(snew-sold)
                  lsMfourS(offset1+indnew)        = lsMfourS(offset1+indnew)+(snew-soldm)
                  ls2MfourS(offset1+indnew)       = ls2MfourS(offset1+indnew)+(snew-soldm)*(snew-soldm)
              endif
              if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                  lMTMfourS(offset2+(inddifm-1)*201+indnew)  = lMTMfourS(offset2+(inddifm-1)*201+indnew)+1.
              endif
           end do
          end do
        end do
      end do 
      call Grid_releaseBlkPtr(blockList(lb), solnData)
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdffour,  gpdffour, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdffourS, gpdffourS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsfour,   gsfour, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMfour,  gsMfour, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsfourS, gsfourS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMfourS, gsMfourS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2four,  gs2four, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mfour, gs2Mfour, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2fourS, gs2fourS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2MfourS, gs2MfourS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMfour,  gTMfour,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMfour, gMTMfour, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMfourS, gMTMfourS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)


    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdffourtot    = 0.
      pdffourtotS2  = 0.
      pdffourtotS4  = 0.
      pdffourtotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMfourtot     = 0.
        MTMfourtot    = 0.
        MTMfourtotS2  = 0.
        MTMfourtotS4  = 0.
        MTMfourtotS8  = 0.
        do inddif=1,401
         TMfourtot  = TMfourtot +gTMfour((inddif-1)*201+indnew)
         MTMfourtot = MTMfourtot+gMTMfour((inddif-1)*201+indnew)
         MTMfourtotS2 = MTMfourtotS2+gMTMfourS((inddif-1)*201+indnew)
         MTMfourtotS4 = MTMfourtotS4+gMTMfourS(80601+(inddif-1)*201+indnew)
         MTMfourtotS8 = MTMfourtotS8+gMTMfourS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMfourtot    .eq. 0.) TMfourtot  = 1.
        if(MTMfourtot   .eq. 0.) MTMfourtot = 1.
        if(MTMfourtotS2 .eq. 0.) MTMfourtotS2 = 1.
        if(MTMfourtotS4 .eq. 0.) MTMfourtotS4 = 1.
        if(MTMfourtotS8 .eq. 0.) MTMfourtotS8 = 1.
        do inddif=1,401
         gTMfour((inddif-1)*201+indnew)= & 
           gTMfour((inddif-1)*201+indnew)/TMfourtot
         gMTMfour((inddif-1)*201+indnew)= & 
           gMTMfour((inddif-1)*201+indnew)/MTMfourtot
         gMTMfourS((inddif-1)*201+indnew)= & 
           gMTMfourS((inddif-1)*201+indnew)/MTMfourtotS2
         gMTMfourS(80601+(inddif-1)*201+indnew)= & 
           gMTMfourS(80601+(inddif-1)*201+indnew)/MTMfourtotS4
         gMTMfourS(161202+(inddif-1)*201+indnew)= & 
           gMTMfourS(161202+(inddif-1)*201+indnew)/MTMfourtotS8
        enddo
       pdffourtot = pdffourtot+gpdffour(indnew)
       pdffourtotS2 = pdffourtotS2+gpdffourS(indnew)
       pdffourtotS4 = pdffourtotS4+gpdffourS(201+indnew)
       pdffourtotS8 = pdffourtotS8+gpdffourS(402+indnew)
      enddo 

      ! write the file
      if((sim_writematrix).or.(mod(nfour,4).eq.1)) then
        print *,'writing the file',io_globalMe 
        write (fnumStr, '(i4.4)') (nfour-npdfstart/steps_four)
        TMFileName = 'TM4_' // fnumStr // '.dat'
        open (funit, file=trim(TMFileName))
        write(funit,21) dr_nstep 
        write(funit,20) simTime
        write(funit,20) simTime-tfour
        write(funit,20) rms_mach
        write(funit,20) rms_mach2
        write(funit,20) rms_mach4
        write(funit,20) rms_mach8
     !     write the pdfs
        do indnew=1,201
         write(funit,20) gpdffour(indnew)/pdffourtot
        enddo
        do indnew=1,201
         write(funit,20) gpdffourS(indnew)/pdffourtotS2
        enddo
        do indnew=202,402
         write(funit,20) gpdffourS(indnew)/pdffourtotS4
        enddo
        do indnew=403,603
         write(funit,20) gpdffourS(indnew)/pdffourtotS8
        enddo
!       write the transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gMTMfour((inddif-1)*201+indnew)
          enddo
        enddo 
!       write the smoothed transition matrix
        do inddif=1,1203
          do indnew=1,201
            write(funit,20) gMTMfourS((inddif-1)*201+indnew)
          enddo
        enddo
!       write the eulerian transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gTMfour((inddif-1)*201+indnew)
          enddo
        enddo
        close(funit)
      endif

!     write the file with s and s2 in it
      write (fnumStr, '(i4.4)') (nfour-npdfstart/steps_four)
      TMFileName = 's_s2_4_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-tfour
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
        write(funit,20) gpdffour(indnew)/pdffourtot
        gpdffour(indnew) = gpdffour(indnew) + 1E-10
      enddo
      do indnew=1,201
        write(funit,20) gpdffourS(indnew)/pdffourtotS2
        gpdffourS(indnew) = gpdffourS(indnew) + 1E-10
      enddo
      do indnew=202,402
        write(funit,20) gpdffourS(indnew)/pdffourtotS4
        gpdffourS(indnew) = gpdffourS(indnew) + 1E-10
      enddo
      do indnew=403,603
        write(funit,20) gpdffourS(indnew)/pdffourtotS8
        gpdffourS(indnew) = gpdffourS(indnew) + 1E-10
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gsfour(indnew)/pdffourtot
      enddo
      do indnew=1,201
       write(funit,20) gsfourS(indnew)/pdffourtotS2
      enddo
      do indnew=202,402
       write(funit,20) gsfourS(indnew)/pdffourtotS4
      enddo
      do indnew=403,603
       write(funit,20) gsfourS(indnew)/pdffourtotS8
      enddo
      do indnew=1,201
       write(funit,20) gsMfour(indnew)/pdffourtot
      enddo
      do indnew=1,201
       write(funit,20) gsMfourS(indnew)/pdffourtotS2
      enddo
      do indnew=202,402
       write(funit,20) gsMfourS(indnew)/pdffourtotS4
      enddo
      do indnew=403,603
       write(funit,20) gsMfourS(indnew)/pdffourtotS8
      enddo
!     write the s2 values
      do indnew=1,201
       write(funit,20) (gs2four(indnew)-gsfour(indnew)**2/gpdffour(indnew))/pdffourtot
      enddo
      do indnew=1,201
       write(funit,20) (gs2fourS(indnew)-gsfourS(indnew)**2/gpdffourS(indnew))/pdffourtotS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2fourS(indnew)-gsfourS(indnew)**2/gpdffourS(indnew))/pdffourtotS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2fourS(indnew)-gsfourS(indnew)**2/gpdffourS(indnew))/pdffourtotS8
      enddo
      do indnew=1,201
       write(funit,20) (gs2Mfour(indnew)-gsMfour(indnew)**2/gpdffour(indnew))/pdffourtot
      enddo
      do indnew=1,201
       write(funit,20) (gs2MfourS(indnew)-gsMfourS(indnew)**2/gpdffourS(indnew))/pdffourtotS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2MfourS(indnew)-gsMfourS(indnew)**2/gpdffourS(indnew))/pdffourtotS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2MfourS(indnew)-gsMfourS(indnew)**2/gpdffourS(indnew))/pdffourtotS8
      enddo
      close(funit)

    endif
    nfour = nfour+1
    tfour = simTime
  endif


  !! ES if it has been 8 then cycles since the last time we wrote the pdf file
  !! then we write the pdf file again 
  if((dr_nstep .ge. steps_eight*neight).and.(dr_nstep .ge. npdfstart)) then

    ! initialize data
    lpdfeight  = 0.
    lseight    = 0.
    ls2eight   = 0.
    lsMeight   = 0.
    ls2Meight  = 0.
    lTMeight   = 0.
    lMTMeight  = 0.

    ! initialize smoothed versions of the data
    lpdfeightS = 0.
    lseightS   = 0.
    ls2eightS  = 0.
    lsMeightS  = 0.
    ls2MeightS = 0.
    lMTMeightS = 0.

    !loop over the block and accumulate the transition matrix
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DEIT_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDEI_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*TM_multiplier +201
             inddifm = (snew-soldm)*TM_multiplier+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdfeight(indnew) = lpdfeight(indnew)+1.
                 lseight(indnew)   = lseight(indnew)+(snew-sold)
                 ls2eight(indnew)  = ls2eight(indnew)+(snew-sold)*(snew-sold)
                 lsMeight(indnew)  = lsMeight(indnew)+(snew-soldm)
                 ls2Meight(indnew) = ls2Meight(indnew)+(snew-soldm)*(snew-soldm)
             endif
             if(((inddif.ge.1).and.(indnew.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lTMeight((inddif-1)*201+indnew)   = lTMeight((inddif-1)*201+indnew)+1.
             endif
             if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lMTMeight((inddifm-1)*201+indnew)  = lMTMeight((inddifm-1)*201+indnew)+1.
             endif
             solnData(DEIT_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDEI_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do

!   now do the case where you smooth over spheres of various sizes
     do csmooth = 1,3
       istep = 1+csmooth
       istsq = istep*istep
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)

       if (csmooth.eq.1) then
         ieul = TWEI_VAR
         imas = MTWE_MSCALAR
       endif
       if (csmooth.eq.2) then
         ieul = THEI_VAR
         imas = MTHE_MSCALAR
       endif
       if (csmooth.eq.3) then
         ieul = FREI_VAR
         imas = MFRE_MSCALAR
       endif

       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!           compute the new and old averaged densities 
              dsmoothnew  = 0.
              dsmoothold  = 0.
              dsmootholdm = 0.
              totzones    = 0.
              do kp=k-istep,k+istep
                ksq = (kp-k)*(kp-k)
                do jp=j-istep,j+istep
                  jsq = (jp-j)*(jp-j)
                  do ip=i-istep,i+istep
                    isq = (ip-i)*(ip-i)
                    if(ksq+jsq+isq.le.istsq) then
                      dsmoothnew  = dsmoothnew +solnData(DENS_VAR,ip,jp,kp)
                      totzones    = totzones + 1
                    endif
                  end do
                end do
              end do
              dsmoothold  = solnData(ieul,i,j,k)
              dsmootholdm = solnData(imas,i,j,k)
              dsmoothnew  = dsmoothnew/totzones
              solnData(ieul,i,j,k)     = dsmoothnew
              solnData(imas,i,j,k)     = dsmoothnew

! add to the pdf and the transition matrix
              snew  = alog(dsmoothnew/ sim_rhoAmbient)
              sold  = alog(dsmoothold/ sim_rhoAmbient)
              soldm = alog(dsmootholdm/sim_rhoAmbient)
              indnew  = snew *20+101
              inddif  = (snew-sold)*TM_multiplier+201
              inddifm = (snew-soldm)*TM_multiplier+201
              if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdfeightS(offset1+indnew)        = lpdfeightS(offset1+indnew)+1.
                  lseightS(offset1+indnew)         = lseightS(offset1+indnew)+(snew-sold)
                  ls2eightS(offset1+indnew)        = ls2eightS(offset1+indnew)+(snew-sold)*(snew-sold)
                  lsMeightS(offset1+indnew)        = lsMeightS(offset1+indnew)+(snew-soldm)
                  ls2MeightS(offset1+indnew)       = ls2MeightS(offset1+indnew)+(snew-soldm)*(snew-soldm)
               endif
               if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                  lMTMeightS(offset2+(inddifm-1)*201+indnew)  = lMTMeightS(offset2+(inddifm-1)*201+indnew)+1.
              endif
            end do
          end do
        end do
      end do 
      call Grid_releaseBlkPtr(blockList(lb), solnData)
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdfeight,  gpdfeight, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdfeightS, gpdfeightS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lseight,   gseight, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMeight,  gsMeight, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lseightS, gseightS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMeightS, gsMeightS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2eight,  gs2eight, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Meight, gs2Meight, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2eightS, gs2eightS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2MeightS, gs2MeightS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMeight,  gTMeight,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMeight, gMTMeight, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMeightS, gMTMeightS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)


    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdfeighttot    = 0.
      pdfeighttotS2  = 0.
      pdfeighttotS4  = 0.
      pdfeighttotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMeighttot     = 0.
        MTMeighttot    = 0.
        MTMeighttotS2  = 0.
        MTMeighttotS4  = 0.
        MTMeighttotS8  = 0.
        do inddif=1,401
         TMeighttot  = TMeighttot +gTMeight((inddif-1)*201+indnew)
         MTMeighttot = MTMeighttot+gMTMeight((inddif-1)*201+indnew)
         MTMeighttotS2 = MTMeighttotS2+gMTMeightS((inddif-1)*201+indnew)
         MTMeighttotS4 = MTMeighttotS4+gMTMeightS(80601+(inddif-1)*201+indnew)
         MTMeighttotS8 = MTMeighttotS8+gMTMeightS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMeighttot    .eq. 0.) TMeighttot  = 1.
        if(MTMeighttot   .eq. 0.) MTMeighttot = 1.
        if(MTMeighttotS2 .eq. 0.) MTMeighttotS2 = 1.
        if(MTMeighttotS4 .eq. 0.) MTMeighttotS4 = 1.
        if(MTMeighttotS8 .eq. 0.) MTMeighttotS8 = 1.
        do inddif=1,401
         gTMeight((inddif-1)*201+indnew)= & 
           gTMeight((inddif-1)*201+indnew)/TMeighttot
         gMTMeight((inddif-1)*201+indnew)= & 
           gMTMeight((inddif-1)*201+indnew)/MTMeighttot
         gMTMeightS((inddif-1)*201+indnew)= & 
           gMTMeightS((inddif-1)*201+indnew)/MTMeighttotS2
         gMTMeightS(80601+(inddif-1)*201+indnew)= & 
           gMTMeightS(80601+(inddif-1)*201+indnew)/MTMeighttotS4
         gMTMeightS(161202+(inddif-1)*201+indnew)= & 
           gMTMeightS(161202+(inddif-1)*201+indnew)/MTMeighttotS8
        enddo
       pdfeighttot = pdfeighttot+gpdfeight(indnew)
       pdfeighttotS2 = pdfeighttotS2+gpdfeightS(indnew)
       pdfeighttotS4 = pdfeighttotS4+gpdfeightS(201+indnew)
       pdfeighttotS8 = pdfeighttotS8+gpdfeightS(402+indnew)
      enddo 

      ! write the matrix file
      ! always write the one once every two cycles
      if((sim_writematrix).or.(mod(neight,2).eq.1)) then
        write (fnumStr, '(i4.4)') (neight-npdfstart/steps_eight)
        TMFileName = 'TM8_' // fnumStr // '.dat'
        open (funit, file=trim(TMFileName))
        write(funit,21) dr_nstep 
        write(funit,20) simTime
        write(funit,20) simTime-teight
        write(funit,20) rms_mach
        write(funit,20) rms_mach2
        write(funit,20) rms_mach4
        write(funit,20) rms_mach8
     !     write the pdfs
        do indnew=1,201
         write(funit,20) gpdfeight(indnew)/pdfeighttot
        enddo
        do indnew=1,201
         write(funit,20) gpdfeightS(indnew)/pdfeighttotS2
        enddo
        do indnew=202,402
         write(funit,20) gpdfeightS(indnew)/pdfeighttotS4
        enddo
        do indnew=403,603
         write(funit,20) gpdfeightS(indnew)/pdfeighttotS8
        enddo
!       write the transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gMTMeight((inddif-1)*201+indnew)
          enddo
        enddo 
!       write the smoothed transition matrix
        do inddif=1,1203
          do indnew=1,201
            write(funit,20) gMTMeightS((inddif-1)*201+indnew)
          enddo
        enddo
!       write the eulerian transition matrix
        do inddif=1,401
          do indnew=1,201
            write(funit,20) gTMeight((inddif-1)*201+indnew)
          enddo
        enddo

        close(funit)
      endif

!     write the file with s and s2 in it
      write (fnumStr, '(i4.4)') (neight-npdfstart/steps_eight)
      TMFileName = 's_s2_8_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-teight
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdfeight(indnew)/pdfeighttot
       gpdfeight(indnew) = gpdfeight(indnew) + 1E-10
      enddo
      do indnew=1,201
        write(funit,20) gpdfeightS(indnew)/pdfeighttotS2
        gpdfeightS(indnew) = gpdfeightS(indnew) + 1E-10
      enddo
      do indnew=202,402
        write(funit,20) gpdfeightS(indnew)/pdfeighttotS4
        gpdfeightS(indnew) = gpdfeightS(indnew) + 1E-10
      enddo
      do indnew=403,603
        write(funit,20) gpdfeightS(indnew)/pdfeighttotS8
        gpdfeightS(indnew) = gpdfeightS(indnew) + 1E-10
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gseight(indnew)/pdfeighttot
      enddo
      do indnew=1,201
       write(funit,20) gseightS(indnew)/pdfeighttotS2
      enddo
      do indnew=202,402
       write(funit,20) gseightS(indnew)/pdfeighttotS4
      enddo
      do indnew=403,603
       write(funit,20) gseightS(indnew)/pdfeighttotS8
      enddo
      do indnew=1,201
       write(funit,20) gsMeight(indnew)/pdfeighttot
      enddo
      do indnew=1,201
       write(funit,20) gsMeightS(indnew)/pdfeighttotS2
      enddo
      do indnew=202,402
       write(funit,20) gsMeightS(indnew)/pdfeighttotS4
      enddo
      do indnew=403,603
       write(funit,20) gsMeightS(indnew)/pdfeighttotS8
      enddo
     do indnew=1,201
       write(funit,20) (gs2eight(indnew)-gseight(indnew)**2/gpdfeight(indnew))/pdfeighttot
      enddo
      do indnew=1,201
       write(funit,20) (gs2eightS(indnew)-gseightS(indnew)**2/gpdfeightS(indnew))/pdfeighttotS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2eightS(indnew)-gseightS(indnew)**2/gpdfeightS(indnew))/pdfeighttotS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2eightS(indnew)-gseightS(indnew)**2/gpdfeightS(indnew))/pdfeighttotS8
      enddo
      do indnew=1,201
       write(funit,20) (gs2Meight(indnew)-gsMeight(indnew)**2/gpdfeight(indnew))/pdfeighttot
      enddo
      do indnew=1,201
       write(funit,20) (gs2MeightS(indnew)-gsMeightS(indnew)**2/gpdfeightS(indnew))/pdfeighttotS2
      enddo
      do indnew=202,402
       write(funit,20) (gs2MeightS(indnew)-gsMeightS(indnew)**2/gpdfeightS(indnew))/pdfeighttotS4
      enddo
      do indnew=403,603
       write(funit,20) (gs2MeightS(indnew)-gsMeightS(indnew)**2/gpdfeightS(indnew))/pdfeighttotS8
      enddo
      close(funit)

    endif
    neight = neight+1
    teight = simTime
  endif

end subroutine IO_writeIntegralQuantities
