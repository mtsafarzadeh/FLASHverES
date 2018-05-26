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
  use Simulation_data, ONLY: tfive,tten,ttwenty,tfourty,npdfstart, &
                             nfive,nten,ntwenty,nfourty,sim_rhoambient
  use Driver_data, ONLY: dr_nstep

  implicit none

#include "mpif.h"
#include "constants.h"
#include "Flash.h"
  

  real, intent(in) :: simTime
  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer :: size(MDIM)

  integer, parameter ::  nGlobalSum = 25         ! Number of globally-summed quantities
  real :: gsum(nGlobalSum)                        ! Global summed quantities
  real :: lsum(nGlobalSum)                        ! Local summed quantities

  integer :: i, j, k
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
!! ES for computing the mach number at 2 4 and 8 grid scales
  real :: rms_mach2,rms_mach4,rms_mach8
  real :: ddd,dvx,dvy,dvz   ! accumulate density, rho*vx, rho*vy, rho*vz

!! ES these are the pdfarry, local and global versions of the pdf and transition
  real :: lpdffive(201),      lpdften(201)
  real :: lpdftwenty(201),lpdffourty(201)
  real :: gpdffive(201),      gpdften(201)
  real :: gpdftwenty(201),gpdffourty(201)
  real :: pdffivetot,  pdffivetotS2,  pdffivetotS4,  pdffivetotS8
  real :: pdftentot,   pdftentotS2,   pdftentotS4,   pdftentotS8
  real :: pdftwentytot,pdftwentytotS2,pdftwentytotS4,pdftwentytotS8
  real :: pdffourtytot,pdffourtytotS2,pdffourtytotS4,pdffourtytotS8

!! these are the arrays for the mean of delta s 
  real :: lsfive(201),      lsten(201)
  real :: lstwenty(201),lsfourty(201)
  real :: gsfive(201),      gsten(201)
  real :: gstwenty(201),gsfourty(201)
  real :: lsMfive(201),      lsMten(201)
  real :: lsMtwenty(201),lsMfourty(201)
  real :: gsMfive(201),      gsMten(201)
  real :: gsMtwenty(201),gsMfourty(201)


!! these are the arrays for deltas^2
  real :: ls2five(201),      ls2ten(201)
  real :: ls2twenty(201),ls2fourty(201)
  real :: gs2five(201),      gs2ten(201)
  real :: gs2twenty(201),gs2fourty(201)
  real :: ls2Mfive(201),      ls2Mten(201)
  real :: ls2Mtwenty(201),ls2Mfourty(201)
  real :: gs2Mfive(201),      gs2Mten(201)
  real :: gs2Mtwenty(201),gs2Mfourty(201)


!! ES these are smoothed versions of the pdfs
  real :: lpdffiveS(603),      lpdftenS(603)
  real :: lpdftwentyS(603),lpdffourtyS(603)
  real :: gpdffiveS(603),      gpdftenS(603)
  real :: gpdftwentyS(603),gpdffourtyS(603)
 

!! these are the arrays for the mean of delta s 
  real :: lsfiveS(601),      lstenS(601)
  real :: lstwentyS(601),lsfourtyS(601)
  real :: gsfiveS(601),      gstenS(601)
  real :: gstwentyS(601),gsfourtyS(601)


!! these are the arrays for deltas^2
  real :: ls2fiveS(601),      ls2tenS(601)
  real :: ls2twentyS(601),ls2fourtyS(601)
  real :: gs2fiveS(601),      gs2tenS(601)
  real :: gs2twentyS(601),gs2fourtyS(601)
 

! 401*201 =  80601
  real :: lTMfive(80601),      lTMten(80601)
  real :: lTMtwenty(80601),    lTMfourty(80601)
  real :: gTMfive(80601),      gTMten(80601)
  real :: gTMtwenty(80601),    gTMfourty(80601)
  real :: lMTMfive(80601),     lMTMten(80601)
  real :: lMTMtwenty(80601),   lMTMfourty(80601)
  real :: gMTMfive(80601),     gMTMten(80601)
  real :: gMTMtwenty(80601),   gMTMfourty(80601)
  real :: TMfivetot, TMtentot, TMtwentytot, TMfourtytot
  real :: MTMfivetot,  MTMfivetotS2,  MTMfivetotS4,  MTMfivetotS8
  real :: MTMtentot,   MTMtentotS2,   MTMtentotS4,   MTMtentotS8
  real :: MTMtwentytot,MTMtwentytotS2,MTMtwentytotS4,MTMtwentytotS8
  real :: MTMfourtytot,MTMfourtytotS2,MTMfourtytotS4,MTMfourtytotS8

! box average over two, four, and eight  cells
! 3* 80601 =  241803
  real :: lMTMfiveS(241803),  lMTMtenS(241803)
  real :: lMTMtwentyS(241803),lMTMfourtyS(241803)
  real :: gMTMfiveS(241803),  gMTMtenS(241803)
  real :: gMTMtwentyS(241803),gMTMfourtyS(241803)
  real :: MTMfivetotS,MTMtentotS,MTMtwentytotS,MTMfourtytotS
  real :: pdffivetotS,pdftentotS,pdftwentytotS,pdffourtytotS

  real    :: sold,  soldm,  snew
  integer :: indnew, inddif, inddifm
  integer :: csmooth ! a counter over smoothing sizes
  integer :: istep   ! 2, 4, or 8
  integer :: offset1,offset2 ! offsets for smoothed pdf and transitiona matrix
  integer :: ip,jp,kp   ! these count over cells in the smoothing volume
  real    :: dsmoothnew,dsmoothold ! smoothed densities
  real    :: totzones

  character (len=4) ::  fnumStr
  character (len=MAX_STRING_LENGTH) :: TMFileName

  
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

 !   compute the rms velocity at the size of 2 cells
     istep = 2
     lsum(23) = 0.
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
              lsum(23) = lsum(23) + (dvx*dvx+dvy*dvy+dvz*dvz)*dvol/ddd
           end do
        end do
     end do


 !   compute the rms velocity at the size of 4 cells
     istep = 4
     lsum(24) = 0.
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
              lsum(24) = lsum(24) + (dvx*dvx+dvy*dvy+dvz*dvz)*dvol/ddd
           end do
        end do
     end do

 !   compute the rms velocity at the size of 8 cells
     istep = 8
     lsum(25) = 0.
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
              lsum(25) = lsum(25) + (dvx*dvx+dvy*dvy+dvz*dvz)*dvol/ddd
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(lb), solnData)
  enddo     !! SS : end of the loop over all the blocks 


  ! we have to do an MPI call to gather it together
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, &
       &                MASTER_PE, MPI_Comm_World, error)

  ! then we write a file called pdf_n

  
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

12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
 
  !! ES if it has been 2 then cycles since the last time we wrote the pdf file 
  !! then we write the pdf file again 
  if((dr_nstep .ge. 1*nfive).and.(dr_nstep .ge. npdfstart)) then

    lpdffive  = 0.
    lsfive    = 0.
    ls2five   = 0. 
    lsMfive   = 0. 
    ls2Mfive  = 0.
    lTMfive   = 0.
    lMTMfive  = 0.

    !loop over the block and accumulate the transition matrix 
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DFIV_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDFV_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*1000 +201
             inddifm = (snew-soldm)*1000+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdffive(indnew) = lpdffive(indnew)+1.
             endif
             if(((inddif.ge.1).and.(indnew.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lTMfive((inddif-1)*201+indnew)   = lTMfive((inddif-1)*201+indnew)+1.
                 lsfive(indnew)   = lsfive(indnew)+(snew-sold)
                 ls2five(indnew)  = ls2five(indnew)+(snew-sold)*(snew-sold)
             endif
             if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lMTMfive((inddifm-1)*201+indnew)  = lMTMfive((inddifm-1)*201+indnew)+1.
                 lsMfive(indnew)  = lsMfive(indnew)+(snew-soldm)
                 ls2Mfive(indnew) = ls2Mfive(indnew)+(snew-soldm)*(snew-soldm)
             endif
             solnData(DFIV_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDFV_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(lb), scr,SCRATCH_CTR)

!   now do the case where you box smooth over 2 points
     do csmooth = 1,3
       istep = 2**csmooth
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)
       totzones   = istep*istep*istep 
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS),istep
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS),istep
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS),istep
!           compute the new and old averaged densities 
              dsmoothnew = 0.
              dsmoothold = 0.
              if (csmooth.eq.1) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(TDFV_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
              endif
              if (csmooth.eq.2) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(FDFV_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
               endif
               if (csmooth.eq.3) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp)
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(EDFV_MSCALAR,ip,jp,kp)
                   end do
                 end do
               end do 
               endif

               dsmoothold = dsmoothold/dsmoothnew
               dsmoothnew = dsmoothnew/totzones

! add to the pdf and the transition matrix
               snew  = alog(dsmoothnew/sim_rhoAmbient)
               soldm = alog(dsmoothold/sim_rhoAmbient)
               indnew  = snew *20+101
               inddifm = (snew-soldm)*1000+201
               if((indnew.ge.1).and.(indnew.le.401)) then
                   lpdffiveS(offset1+indnew)         = lpdffiveS(offset1+indnew)+1.
               endif
               if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                   lMTMfiveS(offset2+(inddifm-1)*201+indnew)  = lMTMfiveS(offset2+(inddifm-1)*201+indnew)+1.
                   lsfiveS(offset1+indnew)           = lsfiveS(offset1+indnew)+(snew-soldm)
                   ls2fiveS(offset1+indnew)          = ls2fiveS(offset1+indnew)+(snew-soldm)*(snew-soldm)
               endif
! write the old one 
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     if (csmooth .eq. 1 ) solnData(TDFV_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 2 ) solnData(FDFV_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 3 ) solnData(EDFV_MSCALAR,ip,jp,kp) = dsmoothnew
                   end do
                 end do  
                end do 
             end do
          end do
        end do
      end do 
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdffive,  gpdffive, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdffiveS, gpdffiveS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsfive,   gsfive, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMfive,  gsMfive, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsfiveS, gsfiveS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2five,  gs2five, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mfive, gs2Mfive, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2fiveS, gs2fiveS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMfive,  gTMfive,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMfive, gMTMfive, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMfiveS, gMTMfiveS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)


    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdffivetot    = 0.
      pdffivetotS2  = 0.
      pdffivetotS4  = 0.
      pdffivetotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
       ! idnew held fixed and inddif varied
        TMfivetot     = 0.
        MTMfivetot    = 0.
        MTMfivetotS2  = 0.
        MTMfivetotS4  = 0.
        MTMfivetotS8  = 0.
        do inddif=1,401
         TMfivetot  = TMfivetot +gTMfive((inddif-1)*201+indnew)
         MTMfivetot = MTMfivetot+gMTMfive((inddif-1)*201+indnew)
         MTMfivetotS2 = MTMfivetotS2+gMTMfiveS((inddif-1)*201+indnew)
         MTMfivetotS4 = MTMfivetotS4+gMTMfiveS(80601+(inddif-1)*201+indnew)
         MTMfivetotS8 = MTMfivetotS8+gMTMfiveS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMfivetot    .eq. 0.) TMfivetot  = 1.
        if(MTMfivetot   .eq. 0.) MTMfivetot = 1.
        if(MTMfivetotS2 .eq. 0.) MTMfivetotS2 = 1.
        if(MTMfivetotS4 .eq. 0.) MTMfivetotS4 = 1.
        if(MTMfivetotS8 .eq. 0.) MTMfivetotS8 = 1.
        do inddif=1,401
         gTMfive((inddif-1)*201+indnew)= & 
           gTMfive((inddif-1)*201+indnew)/TMfivetot
         gMTMfive((inddif-1)*201+indnew)= & 
           gMTMfive((inddif-1)*201+indnew)/MTMfivetot
         gMTMfiveS((inddif-1)*201+indnew)= & 
           gMTMfiveS((inddif-1)*201+indnew)/MTMfivetotS2
         gMTMfiveS(80601+(inddif-1)*201+indnew)= & 
           gMTMfiveS(80601+(inddif-1)*201+indnew)/MTMfivetotS4
         gMTMfiveS(161202+(inddif-1)*201+indnew)= & 
           gMTMfiveS(161202+(inddif-1)*201+indnew)/MTMfivetotS8
        enddo
       pdffivetot = pdffivetot+gpdffive(indnew)
       pdffivetotS2 = pdffivetotS2+gpdffiveS(indnew)
       pdffivetotS4 = pdffivetotS4+gpdffiveS(201+indnew)
       pdffivetotS8 = pdffivetotS8+gpdffiveS(402+indnew)
      enddo 
      ! write the file
      write (fnumStr, '(i4.4)') nfive
      TMFileName = 'TM1_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep 
      write(funit,20) simTime
      write(funit,20) simTime-tfive
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdffive(indnew)/pdffivetot
      enddo
      do indnew=1,201
       write(funit,20) gpdffiveS(indnew)/pdffivetotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdffiveS(indnew)/pdffivetotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdffiveS(indnew)/pdffivetotS8
      enddo
!     write the transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gMTMfive((inddif-1)*201+indnew)
        enddo
      enddo 
!     write the smoothed transition matrix
      do inddif=1,1203
        do indnew=1,201
          write(funit,20) gMTMfiveS((inddif-1)*201+indnew)
        enddo
      enddo
!     write the eulerian transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gTMfive((inddif-1)*201+indnew)
        enddo
      enddo
      close(funit)

!     write the file with s_s2 in it
      write (fnumStr, '(i4.4)') nfive
      TMFileName = 's_s2_1_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-tfive
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdffive(indnew)/pdffivetot
      enddo
      do indnew=1,201
       write(funit,20) gpdffiveS(indnew)/pdffivetotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdffiveS(indnew)/pdffivetotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdffiveS(indnew)/pdffivetotS8
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gsfive(indnew)/pdffivetot
      enddo
      do indnew=1,201
       write(funit,20) gsMfive(indnew)/pdffivetot
      enddo
      do indnew=1,201
       write(funit,20) gsfiveS(indnew)/pdffivetotS2
      enddo
      do indnew=202,402
       write(funit,20) gsfiveS(indnew)/pdffivetotS4
      enddo
      do indnew=403,603
       write(funit,20) gsfiveS(indnew)/pdffivetotS8
      enddo
!     write the s2 values
      do indnew=1,201
       write(funit,20) gs2five(indnew)/pdffivetot
      enddo
      do indnew=1,201
       write(funit,20) gs2Mfive(indnew)/pdffivetot
      enddo
      do indnew=1,201
       write(funit,20) gs2fiveS(indnew)/pdffivetotS2
      enddo
      do indnew=202,402
       write(funit,20) gs2fiveS(indnew)/pdffivetotS4
      enddo
      do indnew=403,603
       write(funit,20) gs2fiveS(indnew)/pdffivetotS8
      enddo
      close(funit)


    endif
    nfive = nfive+1
    tfive = simTime
  endif
 
20  format (' ',E15.7)
21  format (' ',I5.2)
  !! ES if it has been 10 then cycles since the last time we wrote the pdf file 
  !! then we write the pdf file again 
  if((dr_nstep .ge. 2*nten).and.(dr_nstep .ge. npdfstart)) then

    lpdften  = 0.
    lsten    = 0.
    ls2ten   = 0.  
    lsMten   = 0.  
    ls2Mten  = 0.
    lTMten   = 0.
    lMTMten  = 0.

    !loop over the block and accumulate the transition matrix 
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DTEN_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDTN_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*1000 +201
             inddifm = (snew-soldm)*1000+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdften(indnew) = lpdften(indnew)+1.
             endif
             if(((inddif.ge.1).and.(inddif.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lTMten((inddif-1)*201+indnew)   = lTMten((inddif-1)*201+indnew)+1.
                 lsten(indnew)   = lsten(indnew)+(snew-sold)
                 ls2ten(indnew)  = ls2ten(indnew)+(snew-sold)*(snew-sold)
             endif
             if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lMTMten((inddifm-1)*201+indnew)  = lMTMten((inddifm-1)*201+indnew)+1.
                 lsMten(indnew)  = lsMten(indnew)+(snew-soldm)
                 ls2Mten(indnew) = ls2Mten(indnew)+(snew-soldm)*(snew-soldm)
             endif
             solnData(DTEN_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDTN_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(lb), scr,SCRATCH_CTR)

!   now do the case where you box smooth over 2 points
     do csmooth = 1,3
       istep = 2**csmooth
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)
       totzones   = istep*istep*istep 
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS),istep
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS),istep
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS),istep
!           compute the new and old averaged densities 
              dsmoothnew = 0.
              dsmoothold = 0.
              if (csmooth.eq.1) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(TDTN_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
              endif
              if (csmooth.eq.2) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(FDTN_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
               endif
               if (csmooth.eq.3) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp)
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(EDTN_MSCALAR,ip,jp,kp)
                   end do
                 end do
               end do 
               endif

               dsmoothold = dsmoothold/dsmoothnew
               dsmoothnew = dsmoothnew/totzones

! add to the pdf and the transition matrix
               snew  = alog(dsmoothnew/sim_rhoAmbient)
               soldm = alog(dsmoothold/sim_rhoAmbient)
               indnew  = snew *20+101
               inddifm = (snew-soldm)*1000+201
               if((indnew.ge.1).and.(indnew.le.201)) then
                   lpdftenS(offset1+indnew)         = lpdftenS(offset1+indnew)+1.
               endif
               if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                   lMTMtenS(offset2+(inddifm-1)*201+indnew)  = lMTMtenS(offset2+(inddifm-1)*201+indnew)+1.
                   lstenS(offset1+indnew)           = lstenS(offset1+indnew)+(snew-soldm)
                   ls2tenS(offset1+indnew)          = ls2tenS(offset1+indnew)+(snew-soldm)*(snew-soldm)
               endif

! write the old one 
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     if (csmooth .eq. 1 ) solnData(TDTN_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 2 ) solnData(FDTN_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 3 ) solnData(EDTN_MSCALAR,ip,jp,kp) = dsmoothnew
                   end do
                 end do  
                end do 
             end do
          end do
        end do
      end do 
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdften,  gpdften, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdftenS, gpdftenS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsten,   gsten, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMten,  gsMten, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lstenS, gstenS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2ten,  gs2ten, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mten, gs2Mten, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2tenS, gs2tenS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMten,  gTMten,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMten, gMTMten, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMtenS, gMTMtenS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)



    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdftentot    = 0.
      pdftentotS2  = 0.
      pdftentotS4  = 0.
      pdftentotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMtentot     = 0.
        MTMtentot    = 0.
        MTMtentotS2  = 0.
        MTMtentotS4  = 0.
        MTMtentotS8  = 0.
        do inddif =1,401
         TMtentot  = TMtentot +gTMten((inddif-1)*201+indnew)
         MTMtentot = MTMtentot+gMTMten((inddif-1)*201+indnew)
         MTMtentotS2 = MTMtentotS2+gMTMtenS((inddif-1)*201+indnew)
         MTMtentotS4 = MTMtentotS4+gMTMtenS(80601+(inddif-1)*201+indnew)
         MTMtentotS8 = MTMtentotS8+gMTMtenS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMtentot    .eq. 0.) TMtentot  = 1.
        if(MTMtentot   .eq. 0.) MTMtentot = 1.
        if(MTMtentotS2 .eq. 0.) MTMtentotS2 = 1.
        if(MTMtentotS4 .eq. 0.) MTMtentotS4 = 1.
        if(MTMtentotS8 .eq. 0.) MTMtentotS8 = 1.
        do inddif=1,401
         gTMten((inddif-1)*201+indnew)= & 
           gTMten((inddif-1)*201+indnew)/TMtentot
         gMTMten((inddif-1)*201+indnew)= & 
           gMTMten((inddif-1)*201+indnew)/MTMtentot
         gMTMtenS((inddif-1)*201+indnew)= & 
           gMTMtenS((inddif-1)*201+indnew)/MTMtentotS2
         gMTMtenS(80601+(inddif-1)*201+indnew)= & 
           gMTMtenS(80601+(inddif-1)*201+indnew)/MTMtentotS4
         gMTMtenS(161202+(inddif-1)*201+indnew)= & 
           gMTMtenS(161202+(inddif-1)*201+indnew)/MTMtentotS8
        enddo
       pdftentot = pdftentot+gpdften(indnew)
       pdftentotS2 = pdftentotS2+gpdftenS(indnew)
       pdftentotS4 = pdftentotS4+gpdftenS(201+indnew)
       pdftentotS8 = pdftentotS8+gpdftenS(402+indnew)
      enddo 
      print *,'asdf',pdftentotS2,pdftentotS4,pdftentotS8

      ! write the file
      print *,'writing the file',io_globalMe 
      write (fnumStr, '(i4.4)') nten
      TMFileName = 'TM2_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep 
      write(funit,20) simTime
      write(funit,20) simTime-tten
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdften(indnew)/pdftentot
      enddo
      do indnew=1,201
       write(funit,20) gpdftenS(indnew)/pdftentotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdftenS(indnew)/pdftentotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdftenS(indnew)/pdftentotS8
      enddo
!     write the transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gMTMten((inddif-1)*201+indnew)
        enddo
      enddo 
!     write the smoothed transition matrix
      do inddif=1,1203
        do indnew=1,201
          write(funit,20) gMTMtenS((inddif-1)*201+indnew)
        enddo
      enddo
!     write the eulerian transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gTMten((inddif-1)*201+indnew)
        enddo
      enddo
      close(funit)

!     write the file with s and s2 in it
      write (fnumStr, '(i4.4)') nten
      TMFileName = 's_s2_2_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-tten
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdften(indnew)/pdftentot
      enddo
      do indnew=1,201
       write(funit,20) gpdftenS(indnew)/pdftentotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdftenS(indnew)/pdftentotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdftenS(indnew)/pdftentotS8
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gsten(indnew)/pdftentot
      enddo
      do indnew=1,201
       write(funit,20) gsMten(indnew)/pdftentot
      enddo
      do indnew=1,201
       write(funit,20) gstenS(indnew)/pdftentotS2
      enddo
      do indnew=202,402
       write(funit,20) gstenS(indnew)/pdftentotS4
      enddo
      do indnew=403,603
       write(funit,20) gstenS(indnew)/pdftentotS8
      enddo
!     write the s2 values
      do indnew=1,201
       write(funit,20) gs2ten(indnew)/pdftentot
      enddo
      do indnew=1,201
       write(funit,20) gs2Mten(indnew)/pdftentot
      enddo
      do indnew=1,201
       write(funit,20) gs2tenS(indnew)/pdftentotS2
      enddo
      do indnew=202,402
       write(funit,20) gs2tenS(indnew)/pdftentotS4
      enddo
      do indnew=403,603
       write(funit,20) gs2tenS(indnew)/pdftentotS8
      enddo
      close(funit)

    endif
    nten = nten+1
    tten = simTime
  endif
  !! ES if it has been 20 then cycles since the last time we wrote the pdf file 
  !! then we write the pdf file again 
  if((dr_nstep .ge. 4*ntwenty).and.(dr_nstep .ge. npdfstart)) then

    lpdftwenty  = 0.
    lstwenty    = 0.
    ls2twenty   = 0.
    lsMtwenty   = 0.
    ls2Mtwenty  = 0.
    lTMtwenty   = 0.
    lMTMtwenty  = 0.

    !loop over the block and accumulate the transition matrix 
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DTWT_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDTW_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*1000 +201
             inddifm = (snew-soldm)*1000+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdftwenty(indnew) = lpdftwenty(indnew)+1.
             endif
             if(((inddif.ge.1).and.(inddif.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lstwenty(indnew)   = lstwenty(indnew)+(snew-sold)
                 ls2twenty(indnew)  = ls2twenty(indnew)+(snew-sold)*(snew-sold)
                 lTMtwenty((inddif-1)*201+indnew)   = lTMtwenty((inddif-1)*201+indnew)+1.
             endif
             if(((inddifm.ge.1).and.(inddifm.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lsMtwenty(indnew)  = lsMtwenty(indnew)+(snew-soldm)
                 ls2Mtwenty(indnew) = ls2Mtwenty(indnew)+(snew-soldm)*(snew-soldm)
                 lMTMtwenty((inddifm-1)*201+indnew)  = lMTMtwenty((inddifm-1)*201+indnew)+1.
             endif
             solnData(DTWT_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDTW_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(lb), scr,SCRATCH_CTR)

!   now do the case where you box smooth over 2 points
     do csmooth = 1,3
       istep = 2**csmooth
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)
       totzones   = istep*istep*istep 
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS),istep
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS),istep
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS),istep
!           compute the new and old averaged densities 
              dsmoothnew = 0.
              dsmoothold = 0.
              if (csmooth.eq.1) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(TDTW_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
              endif
              if (csmooth.eq.2) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(FDTW_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
               endif
               if (csmooth.eq.3) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp)
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(EDTW_MSCALAR,ip,jp,kp)
                   end do
                 end do
               end do 
               endif

               dsmoothold = dsmoothold/dsmoothnew
               dsmoothnew = dsmoothnew/totzones

! add to the pdf and the transition matrix
               snew  = alog(dsmoothnew/sim_rhoAmbient)
               soldm = alog(dsmoothold/sim_rhoAmbient)
               indnew  = snew *20+101
               inddifm = (snew-soldm)*1000+201
               if((indnew.ge.1).and.(indnew.le.201)) then
                   lpdftwentyS(offset1+indnew)                  = lpdftwentyS(offset1+indnew)+1.
               endif
               if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                   lstwentyS(offset1+indnew)         = lstwentyS(offset1+indnew)+(snew-soldm)
                   ls2twentyS(offset1+indnew)        = ls2twentyS(offset1+indnew)+(snew-soldm)*(snew-soldm)
                   lMTMtwentyS(offset2+(inddifm-1)*201+indnew)  = lMTMtwentyS(offset2+(inddifm-1)*201+indnew)+1.
               endif
! write the old one 
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     if (csmooth .eq. 1 ) solnData(TDTW_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 2 ) solnData(FDTW_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 3 ) solnData(EDTW_MSCALAR,ip,jp,kp) = dsmoothnew
                   end do
                 end do  
                end do 
             end do
          end do
        end do
      end do 
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdftwenty,  gpdftwenty, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdftwentyS, gpdftwentyS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lstwenty,   gstwenty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMtwenty,  gsMtwenty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lstwentyS, gstwentyS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2twenty,  gs2twenty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mtwenty, gs2Mtwenty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2twentyS, gs2twentyS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMtwenty,  gTMtwenty,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMtwenty, gMTMtwenty, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMtwentyS, gMTMtwentyS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)


    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdftwentytot    = 0.
      pdftwentytotS2  = 0.
      pdftwentytotS4  = 0.
      pdftwentytotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMtwentytot     = 0.
        MTMtwentytot    = 0.
        MTMtwentytotS2  = 0.
        MTMtwentytotS4  = 0.
        MTMtwentytotS8  = 0.
        do inddif=1,401
         TMtwentytot  = TMtwentytot +gTMtwenty((inddif-1)*201+indnew)
         MTMtwentytot = MTMtwentytot+gMTMtwenty((inddif-1)*201+indnew)
         MTMtwentytotS2 = MTMtwentytotS2+gMTMtwentyS((inddif-1)*201+indnew)
         MTMtwentytotS4 = MTMtwentytotS4+gMTMtwentyS(80601+(inddif-1)*201+indnew)
         MTMtwentytotS8 = MTMtwentytotS8+gMTMtwentyS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMtwentytot    .eq. 0.) TMtwentytot  = 1.
        if(MTMtwentytot   .eq. 0.) MTMtwentytot = 1.
        if(MTMtwentytotS2 .eq. 0.) MTMtwentytotS2 = 1.
        if(MTMtwentytotS4 .eq. 0.) MTMtwentytotS4 = 1.
        if(MTMtwentytotS8 .eq. 0.) MTMtwentytotS8 = 1.
        do inddif=1,401
         gTMtwenty((inddif-1)*201+indnew)= & 
           gTMtwenty((inddif-1)*201+indnew)/TMtwentytot
         gMTMtwenty((inddif-1)*201+indnew)= & 
           gMTMtwenty((inddif-1)*201+indnew)/MTMtwentytot
         gMTMtwentyS((inddif-1)*201+indnew)= & 
           gMTMtwentyS((inddif-1)*201+indnew)/MTMtwentytotS2
         gMTMtwentyS(80601+(inddif-1)*201+indnew)= & 
           gMTMtwentyS(80601+(inddif-1)*201+indnew)/MTMtwentytotS4
         gMTMtwentyS(161202+(inddif-1)*201+indnew)= & 
           gMTMtwentyS(161202+(inddif-1)*201+indnew)/MTMtwentytotS8
        enddo
       pdftwentytot = pdftwentytot+gpdftwenty(indnew)
       pdftwentytotS2 = pdftwentytotS2+gpdftwentyS(indnew)
       pdftwentytotS4 = pdftwentytotS4+gpdftwentyS(201+indnew)
       pdftwentytotS8 = pdftwentytotS8+gpdftwentyS(402+indnew)
      enddo 
      print *,'asdf',pdftwentytotS2,pdftwentytotS4,pdftwentytotS8

      ! write the file
      print *,'writing the file',io_globalMe 
      write (fnumStr, '(i4.4)') ntwenty
      TMFileName = 'TM4_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep 
      write(funit,20) simTime
      write(funit,20) simTime-ttwenty
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdftwenty(indnew)/pdftwentytot
      enddo
      do indnew=1,201
       write(funit,20) gpdftwentyS(indnew)/pdftwentytotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdftwentyS(indnew)/pdftwentytotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdftwentyS(indnew)/pdftwentytotS8
      enddo
!     write the transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gMTMtwenty((inddif-1)*201+indnew)
        enddo
      enddo 
!     write the smoothed transition matrix
      do inddif=1,1203
        do indnew=1,201
          write(funit,20) gMTMtwentyS((inddif-1)*201+indnew)
        enddo
      enddo
!     write the eulerian transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gTMtwenty((inddif-1)*201+indnew)
        enddo
      enddo
      close(funit)

!     write the file with s and s2 in it
      write (fnumStr, '(i4.4)') ntwenty
      TMFileName = 's_s2_4_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-ttwenty
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdftwenty(indnew)/pdftwentytot
      enddo
      do indnew=1,201
       write(funit,20) gpdftwentyS(indnew)/pdftwentytotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdftwentyS(indnew)/pdftwentytotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdftwentyS(indnew)/pdftwentytotS8
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gstwenty(indnew)/pdftwentytot
      enddo
      do indnew=1,201
       write(funit,20) gsMtwenty(indnew)/pdftwentytot
      enddo
      do indnew=1,201
       write(funit,20) gstwentyS(indnew)/pdftwentytotS2
      enddo
      do indnew=202,402
       write(funit,20) gstwentyS(indnew)/pdftwentytotS4
      enddo
      do indnew=403,603
       write(funit,20) gstwentyS(indnew)/pdftwentytotS8
      enddo
!     write the s2 values
      do indnew=1,201
       write(funit,20) gs2twenty(indnew)/pdftwentytot
      enddo
      do indnew=1,201
       write(funit,20) gs2Mtwenty(indnew)/pdftwentytot
      enddo
      do indnew=1,201
       write(funit,20) gs2twentyS(indnew)/pdftwentytotS2
      enddo
      do indnew=202,402
       write(funit,20) gs2twentyS(indnew)/pdftwentytotS4
      enddo
      do indnew=403,603
       write(funit,20) gs2twentyS(indnew)/pdftwentytotS8
      enddo
      close(funit)

    endif
    ntwenty = ntwenty+1
    ttwenty = simTime
  endif
  !! ES if it has been 40 then cycles since the last time we wrote the pdf file 
  !! then we write the pdf file again 
  if((dr_nstep .ge. 8*nfourty).and.(dr_nstep .ge. npdfstart)) then

    lpdffourty  = 0.
    lsfourty    = 0. 
    ls2fourty   = 0. 
    lsMfourty   = 0. 
    ls2Mfourty  = 0. 
    lTMfourty   = 0.
    lMTMfourty  = 0.

    !loop over the block and accumulate the transition matrix 
    do lb = 1, count
     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             snew = alog(solnData(DENS_VAR,i,j,k)/sim_rhoAmbient)
             sold = alog(solnData(DFTY_VAR,i,j,k)/sim_rhoAmbient)
             soldm = alog(solnData(MDFY_MSCALAR,i,j,k)/sim_rhoAmbient)
             indnew  = snew *20+101
             inddif  = (snew-sold)*1000 +201
             inddifm = (snew-soldm)*1000+201
             if((indnew.ge.1).and.(indnew.le.201)) then
                 lpdffourty(indnew) = lpdffourty(indnew)+1.
             endif
             if(((inddif.ge.1).and.(indnew.ge.1)).and.((inddif.le.401).and.(indnew.le.201))) then
                 lsfourty(indnew)   = lsfourty(indnew)+(snew-sold)
                 ls2fourty(indnew)  = ls2fourty(indnew)+(snew-sold)*(snew-sold)
                 lTMfourty((inddif-1)*201+indnew)   = lTMfourty((inddif-1)*201+indnew)+1.
             endif
             if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                 lsMfourty(indnew)  = lsMfourty(indnew)+(snew-soldm)
                 ls2Mfourty(indnew) = ls2Mfourty(indnew)+(snew-soldm)*(snew-soldm)
                 lMTMfourty((inddifm-1)*201+indnew)  = lMTMfourty((inddifm-1)*201+indnew)+1.
             endif
             solnData(DFTY_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)
             solnData(MDFY_MSCALAR,i,j,k) = solnData(DENS_VAR,i,j,k)
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockList(lb), scr,SCRATCH_CTR)

!   now do the case where you box smooth over 2 points
     do csmooth = 1,3
       istep = 2**csmooth
       offset1 = 201*(csmooth-1)
       offset2 = 80601*(csmooth-1)
       totzones   = istep*istep*istep 
       do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS),istep
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS),istep
             do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS),istep
!           compute the new and old averaged densities 
              dsmoothnew = 0.
              dsmoothold = 0.
              if (csmooth.eq.1) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(TDFY_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
              endif
              if (csmooth.eq.2) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp) 
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(FDFY_MSCALAR,ip,jp,kp) 
                   end do 
                 end do 
               end do 
               endif
               if (csmooth.eq.3) then
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     dsmoothnew = dsmoothnew+solnData(DENS_VAR,ip,jp,kp)
                     dsmoothold = dsmoothold+solnData(DENS_VAR,ip,jp,kp)*solnData(EDFY_MSCALAR,ip,jp,kp)
                   end do
                 end do
               end do 
               endif

               dsmoothold = dsmoothold/dsmoothnew
               dsmoothnew = dsmoothnew/totzones

! add to the pdf and the transition matrix
               snew  = alog(dsmoothnew/sim_rhoAmbient)
               soldm = alog(dsmoothold/sim_rhoAmbient)
               indnew  = snew *20+101
               inddifm = (snew-soldm)*1000+201
               if((indnew.ge.1).and.(indnew.le.201)) then
                   lpdffourtyS(offset1+indnew)                  = lpdffourtyS(offset1+indnew)+1.
               endif
               if(((inddifm.ge.1).and.(indnew.ge.1)).and.((inddifm.le.401).and.(indnew.le.201))) then
                   lMTMfourtyS(offset2+(inddifm-1)*201+indnew)  = lMTMfourtyS(offset2+(inddifm-1)*201+indnew)+1.
                   lsfourtyS(offset1+indnew)         = lsfourtyS(offset1+indnew)+(snew-soldm)
                   ls2fourtyS(offset1+indnew)        = ls2fourtyS(offset1+indnew)+(snew-soldm)*(snew-soldm)
               endif
! write the old one 
               do kp=k,k+istep-1
                 do jp=j,j+istep-1
                   do ip=i,i+istep-1
                     if (csmooth .eq. 1 ) solnData(TDFY_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 2 ) solnData(FDFY_MSCALAR,ip,jp,kp) = dsmoothnew
                     if (csmooth .eq. 3 ) solnData(EDFY_MSCALAR,ip,jp,kp) = dsmoothnew
                   end do
                 end do  
                end do 
             end do
          end do
        end do
      end do 
    end do  !! end of loop over all the blocks

    ! combine the local arrays in to a global array 
    call MPI_Barrier (MPI_Comm_World, error)
    call MPI_Reduce (lpdffourty,  gpdffourty, 201, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lpdffourtyS, gpdffourtyS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsfourty,   gsfourty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsMfourty,  gsMfourty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lsfourtyS, gsfourtyS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2fourty,  gs2fourty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2Mfourty, gs2Mfourty, 201, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (ls2fourtyS, gs2fourtyS, 603, MPI_Double_Precision, MPI_Sum, &
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lTMfourty,  gTMfourty,  80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMfourty, gMTMfourty, 80601, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)
    call MPI_Reduce (lMTMfourtyS, gMTMfourtyS, 241803, MPI_Double_Precision, MPI_Sum, &  
                MASTER_PE, MPI_Comm_World, error)


    ! nomalize the pdf and print it to a file
    if (io_globalMe  == MASTER_PE) then
      pdffourtytot    = 0.
      pdffourtytotS2  = 0.
      pdffourtytotS4  = 0.
      pdffourtytotS8  = 0.
      do indnew=1,201
       ! each row of the transion matrix is normalized
        TMfourtytot     = 0.
        MTMfourtytot    = 0.
        MTMfourtytotS2  = 0.
        MTMfourtytotS4  = 0.
        MTMfourtytotS8  = 0.
        do inddif=1,401
         TMfourtytot  = TMfourtytot +gTMfourty((inddif-1)*201+indnew)
         MTMfourtytot = MTMfourtytot+gMTMfourty((inddif-1)*201+indnew)
         MTMfourtytotS2 = MTMfourtytotS2+gMTMfourtyS((inddif-1)*201+indnew)
         MTMfourtytotS4 = MTMfourtytotS4+gMTMfourtyS(80601+(inddif-1)*201+indnew)
         MTMfourtytotS8 = MTMfourtytotS8+gMTMfourtyS(161202+(inddif-1)*201+indnew)
        enddo
        if(TMfourtytot    .eq. 0.) TMfourtytot  = 1.
        if(MTMfourtytot   .eq. 0.) MTMfourtytot = 1.
        if(MTMfourtytotS2 .eq. 0.) MTMfourtytotS2 = 1.
        if(MTMfourtytotS4 .eq. 0.) MTMfourtytotS4 = 1.
        if(MTMfourtytotS8 .eq. 0.) MTMfourtytotS8 = 1.
        do inddif=1,401
         gTMfourty((inddif-1)*201+indnew)= & 
           gTMfourty((inddif-1)*201+indnew)/TMfourtytot
         gMTMfourty((inddif-1)*201+indnew)= & 
           gMTMfourty((inddif-1)*201+indnew)/MTMfourtytot
         gMTMfourtyS((inddif-1)*201+indnew)= & 
           gMTMfourtyS((inddif-1)*201+indnew)/MTMfourtytotS2
         gMTMfourtyS(80601+(inddif-1)*201+indnew)= & 
           gMTMfourtyS(80601+(inddif-1)*201+indnew)/MTMfourtytotS4
         gMTMfourtyS(161202+(inddif-1)*201+indnew)= & 
           gMTMfourtyS(161202+(inddif-1)*201+indnew)/MTMfourtytotS8
        enddo
       pdffourtytot = pdffourtytot+gpdffourty(indnew)
       pdffourtytotS2 = pdffourtytotS2+gpdffourtyS(indnew)
       pdffourtytotS4 = pdffourtytotS4+gpdffourtyS(201+indnew)
       pdffourtytotS8 = pdffourtytotS8+gpdffourtyS(402+indnew)
      enddo 
      print *,'asdf',pdffourtytotS2,pdffourtytotS4,pdffourtytotS8

      ! write the file
      print *,'writing the file',io_globalMe 
      write (fnumStr, '(i4.4)') nfourty
      TMFileName = 'TM8_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep 
      write(funit,20) simTime
      write(funit,20) simTime-tfourty
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdffourty(indnew)/pdffourtytot
      enddo
      do indnew=1,201
       write(funit,20) gpdffourtyS(indnew)/pdffourtytotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdffourtyS(indnew)/pdffourtytotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdffourtyS(indnew)/pdffourtytotS8
      enddo
!     write the transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gMTMfourty((inddif-1)*201+indnew)
        enddo
      enddo 
!     write the smoothed transition matrix
      do inddif=1,1203
        do indnew=1,201
          write(funit,20) gMTMfourtyS((inddif-1)*201+indnew)
        enddo
      enddo
!     write the eulerian transition matrix
      do inddif=1,401
        do indnew=1,201
          write(funit,20) gTMfourty((inddif-1)*201+indnew)
        enddo
      enddo

      close(funit)

!     write the file with s and s2 in it
      write (fnumStr, '(i4.4)') nfourty
      TMFileName = 's_s2_8_' // fnumStr // '.dat'
      open (funit, file=trim(TMFileName))
      write(funit,21) dr_nstep
      write(funit,20) simTime
      write(funit,20) simTime-tfourty
      write(funit,20) rms_mach
      write(funit,20) rms_mach2
      write(funit,20) rms_mach4
      write(funit,20) rms_mach8
!     write the pdfs
      do indnew=1,201
       write(funit,20) gpdffourty(indnew)/pdffourtytot
      enddo
      do indnew=1,201
       write(funit,20) gpdffourtyS(indnew)/pdffourtytotS2
      enddo
      do indnew=202,402
       write(funit,20) gpdffourtyS(indnew)/pdffourtytotS4
      enddo
      do indnew=403,603
       write(funit,20) gpdffourtyS(indnew)/pdffourtytotS8
      enddo
!     write the s values
      do indnew=1,201
       write(funit,20) gsfourty(indnew)/pdffourtytot
      enddo
      do indnew=1,201
       write(funit,20) gsMfourty(indnew)/pdffourtytot
      enddo
      do indnew=1,201
       write(funit,20) gsfourtyS(indnew)/pdffourtytotS2
      enddo
      do indnew=202,402
       write(funit,20) gsfourtyS(indnew)/pdffourtytotS4
      enddo
      do indnew=403,603
       write(funit,20) gsfourtyS(indnew)/pdffourtytotS8
      enddo
!     write the s2 values
      do indnew=1,201
       write(funit,20) gs2fourty(indnew)/pdffourtytot
      enddo
      do indnew=1,201
       write(funit,20) gs2Mfourty(indnew)/pdffourtytot
      enddo
      do indnew=1,201
       write(funit,20) gs2fourtyS(indnew)/pdffourtytotS2
      enddo
      do indnew=202,402
       write(funit,20) gs2fourtyS(indnew)/pdffourtytotS4
      enddo
      do indnew=403,603
       write(funit,20) gs2fourtyS(indnew)/pdffourtytotS8
      enddo
      close(funit)

    endif
    nfourty = nfourty+1
    tfourty = simTime
  endif

end subroutine IO_writeIntegralQuantities
