!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_DataReconstructNormalDir_WENO5
!!
!! NAME
!!
!!  hy_uhd_DataReconstructNormalDir_WENO5
!!
!! SYNOPSIS
!!
!! hy_uhd_DataReconstructNormalDir_WENO5(integer(IN) :: dir,
!!                                    real(IN)    :: dt,
!!                                    real(IN)    :: delta,
!!                                    pointer(IN) :: Data1D(:,:),
!!                                    pointer(IN) :: DataGrav1D(:),
!!                                    real(IN)    :: FlatCoeff,
!!                                    logical(IN) :: TransUpdateOnly,
!!                                    real(OUT)   :: lambda0(:),
!!                                    real(OUT)   :: leig0(:,:),
!!                                    real(OUT)   :: reig0(:,:),
!!                                    real(OUT)   :: Wp(:),
!!                                    real(OUT)   :: Wm(:),
!!                                    real(OUT)   :: sig(:),
!!                                    real(IN),optional :: dnBnormal,
!!                                    real(IN),optional :: aBnormal(:),
!!                                    real(OUT),optional :: Sr(:),
!!                                    real(OUT),optional :: Sl(:),
!!                                    real(OUT),optional :: SpcSig(:),
!!
!! ARGUMENTS
!!
!!  dir         - normal direction along which the reconstuction is performed
!!  dt          - timestep
!!  delta       - deltas in each {x,y,z} direction
!!  Data1D      - pointer array holding neighboring stencil hydro/MHD data for reconstruction
!!  DataGrav1D  - pointer array holding neighboring stencil gravity data for reconstruction
!!  FlatCoeff   - flattening parameters, primarily for PPM
!!  TransUpdateOnly - a switch for a selective transverse flux update in the normal direction
!!  lambda      - eigenvalue
!!  leig        - left eigenvector
!!  reig        - right eigenvector
!!  Wp,Wm       - left(minus) and right(plus) Riemann states
!!  sig         - transverse flux term
!!  dnBnormal   - MHD multidimensional term for constrained-transport (i.e., USM) MHD scheme
!!  aBnormal    - MHD multidimensional term for constrained-transport (i.e., USM) MHD scheme
!!  Sr          - right Riemann state of species and mass scalars 
!!  Sl          - left  Riemann state of species and mass scalars 
!!  SpcSig      - transverse flux term for spieces and mass scalars
!!
!!
!! DESCRIPTION
!!
!!  This subroutine provides a fifth-order spatially accurate WENO5 for
!!  reconstruction.
!!
!! REFERENCES
!!
!!  * G. Jiang and C.-W. Shu, Efficient implementation of weighted ENO schemes, JCP, 126:202-228, 1996.
!!  * Colella and Woodward, 54, 174 (1984), JCP
!!  * Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics, Springer, 1997
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!  * Stone, Gardiner, Teuben, Hawley, Simon, "Athena: A new code for astrophysical MHD"
!!    arXiv:0804.0402v1 [astro-ph] 2 Apr 2008
!!  * Colella, 87, 171-200 (1990), JCP
!!
!!***

Subroutine hy_uhd_DataReconstructNormalDir_WENO5&
     (dir,dt,delta,Data1D,DataGrav1D,&
      FlatCoeff,TransUpdateOnly, &
      lambda0,leig0,reig0,&
      Wp,Wm,sig,&
      dnBnormal,aBnormal,&     !These are optional
      Sr,Sl,SpcSig)            !These are optional

  use Hydro_data,           ONLY : hy_charLimiting,   &
                                   hy_eswitch,        &
                                   hy_RiemannSolver,  &
                                   hy_tiny,           &
                                   hy_entropy,        &
                                   hy_flattening,     &
                                   hy_transOrder,     &
                                   hy_ContactSteepening,&
                                   hy_upwindTVD,      &
                                   hy_3Torder,        &
                                   hy_useAuxEintEqn

  use hy_uhd_interface,     ONLY : hy_uhd_TVDslope,       &
                                   hy_uhd_TVDslopeUpwind, & 
                                   hy_uhd_upwindTransverseFlux,&
                                   hy_uhd_eigenParameters, &
                                   hy_uhd_eigenValue,      &
                                   hy_uhd_eigenVector

  use hy_uhd_slopeLimiters, ONLY : checkMedian, minmod, mc, vanLeer
  use Timers_interface,     ONLY : Timers_start, Timers_stop

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !!-----Arguments---------------------------------------------------------
  integer,intent(IN) :: dir
  real,   intent(IN) :: dt,delta
  real, pointer, dimension(:,:) :: Data1D
  real, pointer, dimension(:)   :: DataGrav1D
  real,    intent(IN) :: FlatCoeff
  logical, intent(IN) :: TransUpdateOnly
  real, dimension(HY_WAVENUM),intent(OUT) :: lambda0
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
  real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
  real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
  real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
  !optional arguments for MHD-----
  real,intent(IN), optional :: dnBnormal
  real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal
  !-------------------------------
  real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl,SpcSig
  !!------------------------------------------------------------------------

  integer :: n, nVar, iBeg, hyEndVar,hyEndPrimVar
  real    :: dtn,hdtn,qdtn,factor
  real    :: constA,constB,constC,constD,lambdaMax,lambdaMin
  real    :: temp1, temp2, temp3, eta_steep, del2rhoR, del2rhoL, Flattening
  real, dimension(HY_VARINUMMAX) :: vec,sigL,sigR
  real, dimension(HY_VARINUMMAX) :: vecL,vecR,delW,W6,delbar0,delbarP,delbarN
  real, dimension(HY_VARINUMMAX),target :: Vmmm,Vmm,Vm,Vc,Vp,Vpp,Vppp
  real, PARAMETER :: eta1=20.E0, eta2=0.05E0,epsln=0.01E0,K0=0.1E0


  ! EIG SYSTEM
  logical :: cons=.false.
  real    :: cf,uN,cs,ca,as,af,hyp
  real, pointer, dimension(:) :: vm1Ptr,vc0Ptr,vp1Ptr
  real, dimension(MDIM) :: beta
  real, dimension(HY_WAVENUM)   :: lambdaP,lambdaN,lambdaPP,lambdaNN
  real, dimension(HY_VARINUM,HY_WAVENUM) :: leigP,leigN,reigP,reigN


#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
  ! Species and mass scalars
  integer :: isph, ispu !isph for hydro scope index, ispu for unk index
  real, dimension(HY_NSPEC), target :: Sc,Sp,Sm,Spp,Smm
  real, dimension(HY_NSPEC) :: delbarSp0,delbarSpP,delbarSpN,Sp6
#endif
#endif


  ! WENO parameters
  integer,parameter :: wenoExp=1 !(for the expontent of weno5) 
  integer :: NVAR_DO
  real, parameter :: tiny=1.e-36
  real :: sumAlpha,errorCheck,errorCheck1
  real, dimension(3,3) :: coeff1m,coeff1p
  real, dimension(3)   :: coeff2m,coeff2p,W5p,W5m,betaWeno,Alpha5p,Alpha5m,omega,omegaBar,V5p,V5m
  real, dimension(5)   :: Intw
  real, dimension(6)   :: Intp,Intm
  real, dimension(HY_VARINUMMAX) :: Wc0,Wp0,Wpp0,Wm0,Wmm0,Wppp0,Wmmm0,&  ! characteristic variables
                                    vecL_temp, vecR_temp, &
                                    vecL_temp1,vecR_temp1,&
                                    vecL_temp2,vecR_temp2


  ! Debugging mode compilation complains when these are not initialized.
  if (present(Sr)) Sr = 0.
  if (present(Sl)) Sl = 0.
  if (present(SpcSig)) SpcSig = 0.


  ! Set index range depending on hydro or MHD
  ! default for for hydro
  hyEndVar  = HY_ENER
#ifdef FLASH_USM_MHD /* for USM-MHD */
  hyEndVar  = HY_MAGZ
#endif
  hyEndPrimVar = hyEndVar

  if (hy_useAuxEintEqn) hyEndVar = HY_EINT
#ifdef FLASH_UHD_3T
  hyEndVar = HY_END_VARS
#endif
#ifdef GRAVITY
  hyEndVar = HY_GRAV
#endif

  iBeg = 0
  If (TransUpdateOnly) Then
     iBeg = iBeg-1
  Endif
  !! ----------------------------------------------------------------------
  !! 1D array under consideration
  !! ----------------------------------------------------------------------

  !! [A] Usual hydro variables
  IF (.not. TransUpdateOnly) THEN
#if NGUARD > 4
  iBeg = 1

  !Always initialize with zero before storing arrays!
  Vmmm=0; Vppp=0.
  Vmmm(HY_DENS)         = Data1D(DENS_VAR,           iBeg)
  Vmmm(HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,  iBeg)
  Vmmm(HY_PRES)         = Data1D(PRES_VAR,           iBeg)
  Vmmm(HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,  iBeg)
  Vmmm(HY_EINT)         = Data1D(EINT_VAR,           iBeg)

  Vppp(HY_DENS)         = Data1D(DENS_VAR,         6+iBeg)
  Vppp(HY_VELX:HY_VELZ) = Data1D(VELX_VAR:VELZ_VAR,6+iBeg)
  Vppp(HY_PRES)         = Data1D(PRES_VAR,         6+iBeg)
  Vppp(HY_GAMC:HY_GAME) = Data1D(GAMC_VAR:GAME_VAR,6+iBeg)
  Vppp(HY_EINT)         = Data1D(EINT_VAR,         6+iBeg)
#endif

  Vmm=0.;Vpp=0.
  Vmm(HY_DENS)          = Data1D(DENS_VAR,         1+iBeg)
  Vmm(HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,1+iBeg)
  Vmm(HY_PRES)          = Data1D(PRES_VAR,         1+iBeg)
  Vmm(HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,1+iBeg)
  Vmm(HY_EINT)          = Data1D(EINT_VAR,         1+iBeg)

  Vpp(HY_DENS)          = Data1D(DENS_VAR,         5+iBeg)
  Vpp(HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,5+iBeg)
  Vpp(HY_PRES)          = Data1D(PRES_VAR,         5+iBeg)
  Vpp(HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,5+iBeg)
  Vpp(HY_EINT)          = Data1D(EINT_VAR,         5+iBeg)
  ENDIF

  Vm=0.;Vc=0.;Vp=0.
  Vm (HY_DENS)          = Data1D(DENS_VAR,         2+iBeg)
  Vm (HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,2+iBeg)
  Vm (HY_PRES)          = Data1D(PRES_VAR,         2+iBeg)
  Vm (HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,2+iBeg)
  Vm (HY_EINT)          = Data1D(EINT_VAR,         2+iBeg)

  Vc (HY_DENS)          = Data1D(DENS_VAR,         3+iBeg)
  Vc (HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,3+iBeg)
  Vc (HY_PRES)          = Data1D(PRES_VAR,         3+iBeg)
  Vc (HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,3+iBeg)
  Vc (HY_EINT)          = Data1D(EINT_VAR,         3+iBeg)

  Vp (HY_DENS)          = Data1D(DENS_VAR,         4+iBeg)
  Vp (HY_VELX:HY_VELZ)  = Data1D(VELX_VAR:VELZ_VAR,4+iBeg)
  Vp (HY_PRES)          = Data1D(PRES_VAR,         4+iBeg)
  Vp (HY_GAMC:HY_GAME)  = Data1D(GAMC_VAR:GAME_VAR,4+iBeg)
  Vp (HY_EINT)          = Data1D(EINT_VAR,         4+iBeg)



  !! [B] magnetic fields for MHD
#ifdef FLASH_USM_MHD
  IF (.not. TransUpdateOnly) THEN
  Vmm (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,1+iBeg)
  Vpp (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,5+iBeg)
#if NGUARD > 4
  Vmmm(HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,  iBeg)
  Vppp(HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,6+iBeg)
#endif
  ENDIF
  Vm  (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,2+iBeg)
  Vc  (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,3+iBeg)
  Vp  (HY_MAGX:HY_MAGZ) = Data1D(MAGX_VAR:MAGZ_VAR,4+iBeg)
#endif

  !! [C] 3T variables
#ifdef FLASH_UHD_3T
IF (.not. TransUpdateOnly) THEN
  Vmm(HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,1+iBeg),Data1D(EION_VAR,1+iBeg),Data1D(ERAD_VAR,1+iBeg)/)
  Vpp(HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,5+iBeg),Data1D(EION_VAR,5+iBeg),Data1D(ERAD_VAR,5+iBeg)/)
#if NGUARD > 4
  Vmmm(HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,  iBeg),Data1D(EION_VAR,  iBeg),Data1D(ERAD_VAR,  iBeg)/)
  Vppp(HY_EELE:HY_ERAD) = (/Data1D(EELE_VAR,6+iBeg),Data1D(EION_VAR,6+iBeg),Data1D(ERAD_VAR,6+iBeg)/)
#endif
Endif
  Vm (HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,2+iBeg),Data1D(EION_VAR,2+iBeg),Data1D(ERAD_VAR,2+iBeg)/)
  Vc (HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,3+iBeg),Data1D(EION_VAR,3+iBeg),Data1D(ERAD_VAR,3+iBeg)/)
  Vp (HY_EELE:HY_ERAD)  = (/Data1D(EELE_VAR,4+iBeg),Data1D(EION_VAR,4+iBeg),Data1D(ERAD_VAR,4+iBeg)/)
#endif

  !! [D] Gravity component
#ifdef GRAVITY
  IF (.not. TransUpdateOnly) THEN
  Vmm(HY_GRAV)  = DataGrav1D(1+iBeg)
  Vpp(HY_GRAV)  = DataGrav1D(5+iBeg)
#if NGUARD > 4
  Vmmm(HY_GRAV) = DataGrav1D(  iBeg)
  Vppp(HY_GRAV) = DataGrav1D(6+iBeg)
#endif
  Endif
  Vm (HY_GRAV)  = DataGrav1D(2+iBeg)
  Vc (HY_GRAV)  = DataGrav1D(3+iBeg)
  Vp (HY_GRAV)  = DataGrav1D(4+iBeg)
#endif

  !! [E] Species
#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
  do ispu = SPECIES_BEGIN, MASS_SCALARS_END !SPECIES_END
     isph= ispu-NPROP_VARS
     IF (.not. TransUpdateOnly) THEN
     Smm(isph)  = Data1D(ispu,1+iBeg)
     Spp(isph)  = Data1D(ispu,5+iBeg)
     Endif
     Sm (isph)  = Data1D(ispu,2+iBeg)
     Sc (isph)  = Data1D(ispu,3+iBeg)
     Sp (isph)  = Data1D(ispu,4+iBeg)
  enddo
#endif
#endif


  !!**********************!!
  !! BEGIN RECONSTRUCTION !!
  !!**********************!!

  !! half dt & delta
  dtn=dt/delta
  hdtn=0.5*dtn
  qdtn=0.25*dtn

  !! initialize arrays with zero
  sigL=0.;    sigR=0.
  vecL=0.;    vecR=0.
  Wp=0.;      Wm=0.


  call hy_uhd_eigenParameters&
       (Vc(HY_DENS:HY_GAME),dir,uN,cf&
#if defined(FLASH_USM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )
  call hy_uhd_eigenValue&
       (lambda0,uN,cf&
#if defined(FLASH_USM_MHD)
       ,C_alfn=ca,C_slow=cs&
#endif
       )
  call hy_uhd_eigenVector&
       (leig0,reig0,Vc(HY_DENS:HY_GAME),dir,.false.,cf&
#if defined(FLASH_USM_MHD)
       ,C_alfn=ca,C_slow=cs,A_f=af,A_s=as,B_beta=beta&
#endif
       )

  !! -------------------------------------------------------------------------------------!
  !! [1] First we compute transverse flux if we are interested in running ----------------!
  !!     multi-dimensional cases ---------------------------------------------------------!
  !! -------------------------------------------------------------------------------------!
  if (NDIM > 1) then

     vm1Ptr => Vm (:)
     vc0Ptr => Vc (:)
     vp1Ptr => Vp (:)

     call hy_uhd_upwindTransverseFlux&
          (dir,hy_transOrder,vm1Ptr,vc0Ptr,vp1Ptr,lambda0,leig0,reig0,HY_END_VARS,sig)

#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0
     vm1Ptr => Sm (:)
     vc0Ptr => Sc (:)
     vp1Ptr => Sp (:)

     call hy_uhd_upwindTransverseFlux&
          (dir,hy_transOrder,vm1Ptr,vc0Ptr,vp1Ptr,lambda0,leig0,reig0,HY_NSPEC,SpcSig,speciesScalar=.true.)
#endif
#endif
  endif ! NDIM > 1



  !! -------------------------------------------------------------------------------------!
  !! [2] Apply TVD slope limiter for normal gradients ------------------------------------!
  !! -------------------------------------------------------------------------------------!
  IF (.not. TransUpdateOnly) THEN
     call hy_uhd_TVDslope(dir,Vm, Vc, Vp, lambda0,leig0,delbar0)

     !! First initialize flattening coefficients
     if (hy_flattening) then
        Flattening = FlatCoeff
     else
        Flattening = 0.
     endif


     !! -------------------------------------------------------------------------------!
     !! [3] Begin polynomial interpolation for WENO -----------------------------------!
     !! -------------------------------------------------------------------------------!
     !! (a) WENO5 interpolation at the left and right cell interfaces
     if (hy_charLimiting) then
        ! Initialize characteristic variables W*
        Wc0  (hyEndPrimVar+1:) = Vc  (hyEndPrimVar+1:) 
        Wp0  (hyEndPrimVar+1:) = Vp  (hyEndPrimVar+1:) 
        Wm0  (hyEndPrimVar+1:) = Vm  (hyEndPrimVar+1:) 
        Wpp0 (hyEndPrimVar+1:) = Vpp (hyEndPrimVar+1:) 
        Wmm0 (hyEndPrimVar+1:) = Vmm (hyEndPrimVar+1:) 
        Wppp0(hyEndPrimVar+1:) = Vppp(hyEndPrimVar+1:) 
        Wmmm0(hyEndPrimVar+1:) = Vmmm(hyEndPrimVar+1:)         
        NVAR_DO = HY_WAVENUM
     else
        Wc0  = Vc
        Wp0  = Vp
        Wm0  = Vm
        Wpp0 = Vpp
        Wmm0 = Vmm
        Wppp0= Vppp
        Wmmm0= Vmmm
        NVAR_DO = HY_VARINUM3
     endif

     !! Set WENO5 coefficients once and for all

     coeff1p(1,1:3) = (/ 2., -7., 11./) !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
     coeff1p(2,1:3) = (/-1.,  5.,  2./) !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
     coeff1p(3,1:3) = (/ 2.,  5., -1./) !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
     coeff1p        = coeff1p/6.
     coeff2p(1:3)   = (/0.1, 0.6, 0.3/) !=(gamma1,gamma2,gamma3)

     coeff1m(1,1:3) = (/-1.,  5.,  2./) !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
     coeff1m(2,1:3) = (/ 2.,  5., -1./) !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
     coeff1m(3,1:3) = (/ 11.,-7.,  2./) !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
     coeff1m        = coeff1m/6.
     coeff2m(1:3)   = (/0.3, 0.6, 0.1/) !=(gamma1,gamma2,gamma3)

     do n=1,hyEndVar 
 	if (n .le. NVAR_DO) then      
           if (hy_charLimiting) then
           ! Note: It is important to use eigenvectors of the cell Vc at which the reconstruction of
           !       the two left and right states, Wm and Wp, are considered 
              Wmmm0(n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vmmm(HY_DENS:hyEndPrimVar)) 
              Wmm0 (n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vmm (HY_DENS:hyEndPrimVar)) 
              Wm0  (n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vm  (HY_DENS:hyEndPrimVar)) 
              Wc0  (n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vc  (HY_DENS:hyEndPrimVar)) 
              Wp0  (n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vp  (HY_DENS:hyEndPrimVar)) 
              Wpp0 (n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vpp (HY_DENS:hyEndPrimVar)) 
              Wppp0(n) = dot_product(leig0(HY_DENS:hyEndPrimVar,n), Vppp(HY_DENS:hyEndPrimVar)) 
            endif        
        endif
        
        ! Interpolation stencil for weno
        Intw(1:5) = (/Wmm0(n),Wm0(n),Wc0(n),Wp0(n),Wpp0(n)/)


        !! Calculate interface values at i+1/2
        W5p(1) = dot_product(coeff1p(1,1:3),Intw(1:3))
        W5p(2) = dot_product(coeff1p(2,1:3),Intw(2:4))
        W5p(3) = dot_product(coeff1p(3,1:3),Intw(3:5))

        !! Calculate interface values at i-1/2
        W5m(1) = dot_product(coeff1m(1,1:3),Intw(1:3))
        W5m(2) = dot_product(coeff1m(2,1:3),Intw(2:4))
        W5m(3) = dot_product(coeff1m(3,1:3),Intw(3:5))

        !! Calculate smoothness indicators at i+1/2
        betaWeno(1) = 13./12.*(Intw(1)-2.*Intw(2)+Intw(3))**2 + 0.25*(   Intw(1)-4.*Intw(2)+3.*Intw(3))**2
        betaWeno(2) = 13./12.*(Intw(2)-2.*Intw(3)+Intw(4))**2 + 0.25*(   Intw(2)              -Intw(4))**2
        betaWeno(3) = 13./12.*(Intw(3)-2.*Intw(4)+Intw(5))**2 + 0.25*(3.*Intw(3)-4.*Intw(4)   +Intw(5))**2

        !! Calculate weights at i+1/2
        Alpha5p(1) = coeff2p(1)/(tiny+betaWeno(1))**wenoExp 
        Alpha5p(2) = coeff2p(2)/(tiny+betaWeno(2))**wenoExp 
        Alpha5p(3) = coeff2p(3)/(tiny+betaWeno(3))**wenoExp 
 
        !! Normalize weights at i+1/2
        sumAlpha = Alpha5p(1)+Alpha5p(2)+Alpha5p(3)
        omega(1) = Alpha5p(1)/sumAlpha
        omega(2) = Alpha5p(2)/sumAlpha
        omega(3) = Alpha5p(3)/sumAlpha

        !! Compute interface value at i+1/2
        vecR(n) = dot_product(omega(1:3), W5p(1:3))

        !! Calculate weights at i-1/2
        Alpha5m(1) = coeff2m(1)/(tiny+betaWeno(1))**wenoExp 
        Alpha5m(2) = coeff2m(2)/(tiny+betaWeno(2))**wenoExp 
        Alpha5m(3) = coeff2m(3)/(tiny+betaWeno(3))**wenoExp 

        !! Normalize weights at i-1/2
        sumAlpha = Alpha5m(1)+Alpha5m(2)+Alpha5m(3)
        omega(1) = Alpha5m(1)/sumAlpha
        omega(2) = Alpha5m(2)/sumAlpha
        omega(3) = Alpha5m(3)/sumAlpha

        !! Compute interface value at i+1/2
        vecL(n) = dot_product(omega(1:3), W5m(1:3))


     enddo


     if (hy_charLimiting) then
        ! characteristic variables at i+1/2 ---------------------------
        vecR(HY_DENS:hyEndPrimVar) = &
             reig0(HY_DENS:hyEndPrimVar,HY_FASTLEFT)*vecR(HY_FASTLEFT)+&
#ifdef FLASH_USM_MHD
             reig0(HY_DENS:hyEndPrimVar,HY_ALFNLEFT)*vecR(HY_ALFNLEFT)+&
#endif
             reig0(HY_DENS:hyEndPrimVar,HY_SLOWLEFT)*vecR(HY_SLOWLEFT)+&
             reig0(HY_DENS:hyEndPrimVar,HY_ENTROPY )*vecR(HY_ENTROPY )+&
             reig0(HY_DENS:hyEndPrimVar,HY_SLOWRGHT)*vecR(HY_SLOWRGHT)+&
#ifdef FLASH_USM_MHD
             reig0(HY_DENS:hyEndPrimVar,HY_ALFNRGHT)*vecR(HY_ALFNRGHT)+&
#endif
             reig0(HY_DENS:hyEndPrimVar,HY_FASTRGHT)*vecR(HY_FASTRGHT)

        ! characteristic variables at i-1/2 ---------------------------
        vecL(HY_DENS:hyEndPrimVar) = &
             reig0(HY_DENS:hyEndPrimVar,HY_FASTLEFT)*vecL(HY_FASTLEFT)+&
#ifdef FLASH_USM_MHD
             reig0(HY_DENS:hyEndPrimVar,HY_ALFNLEFT)*vecL(HY_ALFNLEFT)+&
#endif
             reig0(HY_DENS:hyEndPrimVar,HY_SLOWLEFT)*vecL(HY_SLOWLEFT)+&
             reig0(HY_DENS:hyEndPrimVar,HY_ENTROPY )*vecL(HY_ENTROPY )+&
             reig0(HY_DENS:hyEndPrimVar,HY_SLOWRGHT)*vecL(HY_SLOWRGHT)+&
#ifdef FLASH_USM_MHD
             reig0(HY_DENS:hyEndPrimVar,HY_ALFNRGHT)*vecL(HY_ALFNRGHT)+&
#endif
             reig0(HY_DENS:hyEndPrimVar,HY_FASTRGHT)*vecL(HY_FASTRGHT)
     endif



     !! (b) PPM interpolate species and mass scalars
#ifndef FLASH_UHD_NEED_SCRATCHVARS
#if (NSPECIES+NMASS_SCALARS) > 0

     do isph= 1, HY_NSPEC

         ! Interpolation stencil for weno 
         Intw(1:5) = (/Smm(isph),Sm(isph),Sc(isph),Sp(isph),Spp(isph)/) 
 	 
  
         !! Calculate interface values at i+1/2 
         W5p(1) = dot_product(coeff1p(1,1:3),Intw(1:3)) 
         W5p(2) = dot_product(coeff1p(2,1:3),Intw(2:4)) 
         W5p(3) = dot_product(coeff1p(3,1:3),Intw(3:5)) 
 
         !! Calculate interface values at i-1/2 
         W5m(1) = dot_product(coeff1m(1,1:3),Intw(1:3)) 
         W5m(2) = dot_product(coeff1m(2,1:3),Intw(2:4)) 
         W5m(3) = dot_product(coeff1m(3,1:3),Intw(3:5)) 
  
         !! Calculate smoothness indicators at i+1/2 
         betaWeno(1) = 13./12.*(Intw(1)-2.*Intw(2)+Intw(3))**2 + 0.25*(   Intw(1)-4.*Intw(2)+3.*Intw(3))**2 
         betaWeno(2) = 13./12.*(Intw(2)-2.*Intw(3)+Intw(4))**2 + 0.25*(   Intw(2)              -Intw(4))**2 
         betaWeno(3) = 13./12.*(Intw(3)-2.*Intw(4)+Intw(5))**2 + 0.25*(3.*Intw(3)-4.*Intw(4)   +Intw(5))**2 
	 
         !! Calculate weights at i+1/2 
         Alpha5p(1) = coeff2p(1)/(tiny+betaWeno(1))**wenoExp 
         Alpha5p(2) = coeff2p(2)/(tiny+betaWeno(2))**wenoExp 
         Alpha5p(3) = coeff2p(3)/(tiny+betaWeno(3))**wenoExp 
 
        !! Normalize weights at i+1/2 
         sumAlpha = Alpha5p(1)+Alpha5p(2)+Alpha5p(3) 
         omega(1) = Alpha5p(1)/sumAlpha 
         omega(2) = Alpha5p(2)/sumAlpha 
         omega(3) = Alpha5p(3)/sumAlpha 
 
         !! Compute interface value at i+1/2 
         Sr(isph) = dot_product(omega(1:3), W5p(1:3)) 
 	 
         !! Calculate weights at i-1/2 
         Alpha5m(1) = coeff2m(1)/(tiny+betaWeno(1))**wenoExp 
         Alpha5m(2) = coeff2m(2)/(tiny+betaWeno(2))**wenoExp 
         Alpha5m(3) = coeff2m(3)/(tiny+betaWeno(3))**wenoExp 
 	 
         !! Normalize weights at i-1/2 
         sumAlpha = Alpha5m(1)+Alpha5m(2)+Alpha5m(3) 
         omega(1) = Alpha5m(1)/sumAlpha 
         omega(2) = Alpha5m(2)/sumAlpha 
         omega(3) = Alpha5m(3)/sumAlpha 
 
  
        !! Compute interface value at i+1/2 
        Sl(isph) = dot_product(omega(1:3), W5m(1:3)) 
	 

        if (6.*(Sr(isph)-Sl(isph))*(Sc(isph)-0.5*(Sr(isph)+Sl(isph))) &
             > (Sr(isph)-Sl(isph))**2) then
           Sl(isph) = 3.*Sc(isph) - 2.*Sr(isph)
        endif
        if (6.*(Sr(isph)-Sl(isph))*(Sc(isph)-0.5*(Sr(isph)+Sl(isph))) &
             <-(Sr(isph)-Sl(isph))**2) then
           Sr(isph) = 3.*Sc(isph) - 2.*Sl(isph)
        endif
        delbarSp0(isph) = Sr(isph)-Sl(isph)
        Sp6(isph) = 6.*(Sc(isph)-0.5*(Sr(isph)+Sl(isph)))

        Sl(isph) = Sl(isph)-min(lambda0(HY_ENTROPY),0.)*hdtn*&
       (delbarSp0(isph)+(1.+min(lambda0(HY_ENTROPY),0.)*4./3.*hdtn)*Sp6(isph))

        Sr(isph) = Sr(isph)-max(lambda0(HY_ENTROPY),0.)*hdtn*&
       (delbarSp0(isph)-(1.-max(lambda0(HY_ENTROPY),0.)*4./3.*hdtn)*Sp6(isph))
     enddo

     if (hy_flattening) then 
        Sl(:) = Flattening*Sc(:) + (1.0-Flattening)*Sl(:) 
        Sr(:) = Flattening*Sc(:) + (1.0-Flattening)*Sr(:) 
     endif      

#endif /*  (NSPECIES+NMASS_SCALARS) > 0 */
#endif /* ifndef FLASH_UHD_NEED_SCRATCHVARS */
     !! End of polynomial interpolation for PPM


     !! -------------------------------------------------------------------------------------!
     !! [4] Contact steepening for PPM  -----------------------------------------------------!
     !! -------------------------------------------------------------------------------------!
!!     if (hy_ContactSteepening) then
!!        temp1 = Vp(HY_DENS) - Vm(HY_DENS)
!!        if (abs(temp1) > hy_tiny) then
!!           ! Eqn 1.17 : Second derivatives
!!           del2rhoR = (Vpp(HY_DENS)-2.*Vp(HY_DENS)+ Vc(HY_DENS))/(6.*delta*delta)
!!           del2rhoL = ( Vc(HY_DENS)-2.*Vm(HY_DENS)+Vmm(HY_DENS))/(6.*delta*delta)

           ! Third derivative
!!           eta_steep = (del2rhoL-del2rhoR)*delta**2/temp1
!!           if (del2rhoR*del2rhoL >= 0.) then
!!              eta_steep = 0.
!!           endif
!!           if (epsln*min(Vp(HY_DENS),Vm(HY_DENS))-abs(Vp(HY_DENS) - Vm(HY_DENS)) >= 0.) then
!!              eta_steep = 0.
!!           endif

           ! Eqn 1.16
!!           eta_steep = max(0., min(1., eta1*(eta_steep - eta2)))

           ! Eqn 3.2
!!           temp2 = abs(Vp(HY_PRES)-Vm(HY_PRES))/min(Vp(HY_PRES),Vm(HY_PRES))
!!           temp3 = abs(Vp(HY_DENS)-Vm(HY_DENS))/min(Vp(HY_DENS),Vm(HY_DENS))

!!           if (Vc(HY_GAME)*K0*temp3-temp2 < 0.0) then
!!              eta_steep = 0.
!!           endif

           ! Eqn 1.15
!!           vecL(HY_DENS) = vecL(HY_DENS)*(1.-eta_steep) + (Vm(HY_DENS)+0.5*delbarN(HY_DENS))*eta_steep
!!           vecR(HY_DENS) = vecR(HY_DENS)*(1.-eta_steep) + (Vp(HY_DENS)-0.5*delbarP(HY_DENS))*eta_steep
!!        endif
!!    endif
     !! End of Contact steepening for PPM

     !! -------------------------------------------------------------------------------------!
     !! [5] Flattening for PPM  -------------------------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     if (hy_flattening) then
        vecL(:) = Flattening*Vc(:) + (1.0-Flattening)*vecL(:)
        vecR(:) = Flattening*Vc(:) + (1.0-Flattening)*vecR(:)
     endif

     !! -------------------------------------------------------------------------------------!
     !! [6] Monotonicity check for PPM  -----------------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     ! Ensure that the interpolated values lie between the cell-centered values
     ! Limit according to Colella-Woodward Eqn 1.10
     do n=HY_DENS,hyEndVar
        if ( 6.*(vecR(n)-vecL(n))*(Vc(n)-0.5*(vecL(n)+vecR(n))) > (vecR(n) - vecL(n))**2  ) then
           vecL(n) = 3.*Vc(n) - 2.*vecR(n)
        endif
        if ( 6.*(vecR(n)-vecL(n))*(Vc(n)-0.5*(vecL(n)+vecR(n))) < -(vecR(n) - vecL(n))**2  ) then
           vecR(n) = 3.*Vc(n) - 2.*vecL(n)
        endif
     enddo
     !! End of Contact steepeing, Flattening, and Monotonicity constraint for PPM


     !! -------------------------------------------------------------------------------------!
     !! [7] Take initial guesses for the left and right states-------------------------------!
     !! -------------------------------------------------------------------------------------!
     !! PPM coefficients for parabolic interpolations
     delW(HY_DENS:hyEndVar) = vecR(HY_DENS:hyEndVar)-vecL(HY_DENS:hyEndVar)
     W6(HY_DENS:hyEndVar)   = 6.*(Vc(HY_DENS:hyEndVar)&
                             -0.5*(vecR(HY_DENS:hyEndVar)+vecL(HY_DENS:hyEndVar)))

     !! [7-a] Right states
     !! Primary variables first
     lambdaMax =max(lambda0(HY_FASTRGHT),0.)
     Wp(HY_DENS:hyEndPrimVar) = vecR(HY_DENS:hyEndPrimVar) - lambdaMax*hdtn &
          *(delW(HY_DENS:hyEndPrimVar) - (1.0 - lambdaMax*hdtn*4./3.)*W6(HY_DENS:hyEndPrimVar))

     !! Secondary variables, gamc, game, eint, 3T vars, grav
     lambdaMax =max(lambda0(HY_ENTROPY),0.)
     Wp(hyEndPrimVar+1:hyEndVar) = vecR(hyEndPrimVar+1:hyEndVar) - lambdaMax*hdtn &
          *(delW(hyEndPrimVar+1:hyEndVar) - (1.0 - lambdaMax*hdtn*4./3.)*W6(hyEndPrimVar+1:hyEndVar))

     !! [7-b] Left states
     !! Primary variables first
     lambdaMin = -min(lambda0(HY_FASTLEFT),0.)
     Wm(HY_DENS:hyEndPrimVar) = vecL(HY_DENS:hyEndPrimVar) + lambdaMin*hdtn &
          *(delW(HY_DENS:hyEndPrimVar) + (1.0 - lambdaMin*hdtn*4./3.)*W6(HY_DENS:hyEndPrimVar))

     !! Secondary variables, gamc, game, eint, 3T vars, grav
     lambdaMin = -min(lambda0(HY_ENTROPY),0.)
     Wm(hyEndPrimVar+1:hyEndVar) = vecL(hyEndPrimVar+1:hyEndVar) + lambdaMin*hdtn &
          *(delW(hyEndPrimVar+1:hyEndVar) + (1.0 - lambdaMin*hdtn*4./3.)*W6(hyEndPrimVar+1:hyEndVar))

     !! [7-c] Apply constraints
     !! Force constant state if simple gamma laws
#ifdef FLASH_EOS_GAMMA
     Wp(HY_GAMC:HY_GAME)=Vc(HY_GAMC:HY_GAME)
     Wm(HY_GAMC:HY_GAME)=Vc(HY_GAMC:HY_GAME)
#endif

     !! [7-d] Gravity component shoud only be spatially reconstructed, 
     !!       without being characteristically traced in time along with fluid's wave
#ifdef GRAVITY
     Wp(HY_GRAV) = vecR(HY_GRAV)
     Wm(HY_GRAV) = vecL(HY_GRAV)
#endif
     !! End of initial guesses
     !! End of high-order polynomial interpolations for PPM interface values


     !! -------------------------------------------------------------------------------------!
     !! [8] Advance the above interpolated interface values by 1/2 time step using ----------!
     !!     characteristic tracing method     -----------------------------------------------!
     !! -------------------------------------------------------------------------------------!
     do n=1,HY_WAVENUM
        
        constA = dot_product(leig0(HY_DENS:hyEndPrimVar,n), delW(HY_DENS:hyEndPrimVar))
        constB = dot_product(leig0(HY_DENS:hyEndPrimVar,n),  -W6(HY_DENS:hyEndPrimVar))

        if (hy_RiemannSolver == ROE) then
           if (lambda0(n) < 0.) then
              ! PPM step 10
              !! Left states:
              vecL(HY_DENS:hyEndPrimVar) =  &
                .5*(-1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.+2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

              sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

           elseif (lambda0(n) > 0.) then
              !! Right states:
              vecR(HY_DENS:hyEndPrimVar) =  &
                .5*( 1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.-2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

              sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)

           endif
        else
           !! Left and right states for HLL* type solvers:
           !! Apply monotone slope limiting for normal flux
           !! PPM step 10

           !! Left states:
           vecL(HY_DENS:hyEndPrimVar) =  &
                .5*(-1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.+2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

           sigL(HY_DENS:hyEndPrimVar) = sigL(HY_DENS:hyEndPrimVar) + vecL(HY_DENS:hyEndPrimVar)

           !! Right states:
           vecR(HY_DENS:hyEndPrimVar) =  &
                .5*( 1.-   dtn*lambda0(n)                          )*reig0(HY_DENS:hyEndPrimVar,n)*constA &
              +.25*( 1.-2.*dtn*lambda0(n)+4./3.*(dtn*lambda0(n))**2)*reig0(HY_DENS:hyEndPrimVar,n)*constB

           sigR(HY_DENS:hyEndPrimVar) = sigR(HY_DENS:hyEndPrimVar) + vecR(HY_DENS:hyEndPrimVar)

        endif
     enddo ! do n=1,HY_WAVENUM

 
     !! -------------------------------------------------------------------------------------!
     !! [9] Consider lower order schemes for 3T variables if requested ----------------------!
     !! -------------------------------------------------------------------------------------!
#ifdef FLASH_UHD_3T
     IF (hy_3Torder .ne. 3) THEN
        if (hy_3Torder == 1) then
           Wp(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)
           Wm(HY_EINT:HY_ERAD) = Vc(HY_EINT:HY_ERAD)

        elseif (hy_3Torder == 2) then


           Wp(HY_EINT:HY_ERAD)=Vc(HY_EINT:HY_ERAD)&
                +0.5*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)

           Wm(HY_EINT:HY_ERAD)=Vc(HY_EINT:HY_ERAD)&
                -0.5*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)

           do nVar=HY_EINT,HY_ERAD
              Wm(nVar) = max(min(Vc(nVar),Wm(nVar)),min(max(Vc(nVar),Wm(nVar)),Wm(nVar)))
              Wp(nVar) = max(min(Vc(nVar),Wp(nVar)),min(max(Vc(nVar),Wp(nVar)),Wp(nVar)))
           enddo

           Wp(HY_EINT:HY_ERAD)=Wp(HY_EINT:HY_ERAD)&
                -max(lambda0(HY_ENTROPY),0.)*hdtn*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)

           Wm(HY_EINT:HY_ERAD)=Wm(HY_EINT:HY_ERAD)&
                -min(lambda0(HY_ENTROPY),0.)*hdtn*delbar0(HY_EINT:HY_ERAD)*(1.-Flattening)
        endif
     ENDIF
#endif


     !! -------------------------------------------------------------------------------------!
     !! [10] Finalize the Riemann states in 1D normal direction -----------------------------!
     !! -------------------------------------------------------------------------------------!
     !! Riemann states in normal direction
     !! Note that gamc, game, grav, & 3T variables are treated separately in the above
     Wm(HY_DENS:hyEndPrimVar) = Vc(HY_DENS:hyEndPrimVar)+W6(HY_DENS:hyEndPrimVar)/12.+sigL(HY_DENS:hyEndPrimVar)
     Wp(HY_DENS:hyEndPrimVar) = Vc(HY_DENS:hyEndPrimVar)+W6(HY_DENS:hyEndPrimVar)/12.+sigR(HY_DENS:hyEndPrimVar)

#if defined(FLASH_USM_MHD)
     Wm(HY_DENS:hyEndPrimVar) = Wm(HY_DENS:hyEndPrimVar)-hdtn*aBnormal(HY_DENS:hyEndPrimVar)*dnBnormal
     Wp(HY_DENS:hyEndPrimVar) = Wp(HY_DENS:hyEndPrimVar)-hdtn*aBnormal(HY_DENS:hyEndPrimVar)*dnBnormal
#endif

  ENDIF ! end of IF (.not. TransUpdateOnly) THEN

End Subroutine hy_uhd_DataReconstructNormalDir_WENO5
