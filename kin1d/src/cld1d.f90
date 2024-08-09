subroutine columnmodel

!----------------------------------------------------------------------------!
!  1-D kinematic cloud model driver for testing microphysics scheme.
!  Reads in temperature and humidity from a sounding, prescribes an updraft
!  profile, and begins a time integration for a simple cloud simulation.
!
!  The following procedures are performed at each time step:
!
!  - update vertical velocity
!  - compute column-integrated mass (and number) for the beginning of step
!  - update model variables due to advection, compressibility, and divergence
!  - force global concervation of mass (and number)
!  - compute total integrated mass before call to microphysics
!  - call microphysics
!  - compute total integrated mass after call to microphysics
!  - add low-level moisture (to prevent depletion)
!  - write output files
!
!--------------------------------------------------------------------------!
! Variable names:
! -----   ----
! cld1d    p3
! -----   ----
!  QI     qitot  - total (deposition + rime) ice mass mixing ratio
!  QG     qrim   - rime ice mass mixing ratio
!  QL     qiliq_in  - liquid on ice mass mixing ratio (optional)
!  NI     nitot  - ice number mixing ratio
!  BG     birim  - rime volume mixing ratio
!  ZI     zitot  - 6th moment mixing ratio
!
!--------------------------------------------------------------------------!
!  Author:         Jason Milbrandt
!  Last modified:  2023-05-28
!--------------------------------------------------------------------------!

      use subs_cld1d
      use microphy_p3

      implicit none

      integer, parameter :: n_iceCat     =  2
      logical, parameter :: liqFrac      = .true.
      logical, parameter :: trplMomIce   = .true.

      logical, parameter :: scpf_on      = .false.  ! switch for cloud fraction parameterization (SCPF)
      real,    parameter :: scpf_pfrac   = 1.        ! precipitation fraction factor (SCPF)
      real,    parameter :: scpf_resfact = 1.        ! model resolution factor (SCPF)
      logical, parameter :: prog_nc_ssat = .true.
     !logical, parameter :: nk_BOTTOM    = .true.   !.T. --> nk at bottom
      logical, parameter :: debug_on     = .false.   ! switch for run-time check-values in p3_main
      real, parameter    :: clbfact_dep = 1.0       !calibration factor for deposition
      real, parameter    :: clbfact_sub = 1.0       !calibration factor for sublimation

      integer      :: airmass,SCHEME,nk,outfreq,ttotmin,iice
      real         :: AMPA,AMPB,Htop0,tscale1,tscale2,dt
      logical      :: microON,EVOLVING,TDMIN
      character*10 :: sndcase

      parameter (sndcase = 'ALBERTA')
      parameter (microON = .true. )     ! call microphysics SCHEME
      parameter (TDMIN   = .true. )     ! prevent low-level moisture depletion
      parameter (EVOLVING= .true. )     ! switch for evolving updraft
      parameter (AMPA    = 2.     )     ! initial central updraft speed [m s-1] (evolving only)
      parameter (AMPB    = 5.     )     ! maximum central updraft speed [m s-1]
      parameter (Htop0   = 5000.  )     ! initial height of cloud top   [m]
      parameter (tscale1 = 5400.  )     ! period for evolving AMPB      [s]
      parameter (tscale2 = 5400.  )     ! period for evolving Hcld      [s]
      parameter (nk      =  41    )     ! number of vertical levels
!     parameter (nk      =  62    )     ! number of vertical levels
!     parameter (nk      =  86    )     ! number of vertical levels
      parameter (outfreq =  1     )     ! output every 'OUTFREQ' minutes
      parameter (dt      = 10.    )     ! time step                     [s]
      parameter (ttotmin = 90     )     ! total integration time	[min]

      character(len=16), parameter :: model = 'KIN1D'
!     character(len=16), parameter :: model = 'WRF'  !for level tests
      logical, parameter           :: abort_on_err = .false.
      logical, parameter           :: dowr = .true.

      character(len=1024), parameter :: LT_path  = './lookup_tables'
!     character(len=1024), parameter :: LT_path = '/users/milbrand/mp_p3/lookupTables/tables'  ! override default


!---------------------------------------------------------------------------------!
!#include "consphy.cdk"  (necessary parameters only)
      real, parameter :: TRPL     =.27316e+3          !K; triple point of water
      real, parameter :: EPS1     =.62194800221014    ! ; RGASD/RGASV
      real, parameter :: EPS2     =.3780199778986     !; 1 - EPS1
      real, parameter :: PI       =.314159265359e+1   ! PI constant = ACOS(-1)
      real, parameter :: GRAV     =.980616e+1         ! M s-2; gravitational acceleration

!#include "dintern.cdk"  (necessary variables only)
      real   :: TTT,PRS,QQQ
      real*8 :: FOEW,FOQST

!------------------------------------------------------------------------------!
!#include "fintern.cdk"
!   DEFINITION DES FONCTIONS THERMODYNAMIQUES DE BASE
!   POUR LES CONSTANTES, UTILISER LE COMMON /CONSPHY/
!     NOTE: TOUTES LES FONCTIONS TRAVAILLENT AVEC LES UNITES S.I.
!     FONCTION DE TENSION DE VAPEUR SATURANTE (TETENS) - EW OU EI SELON TT
      FOEW(TTT) = 610.78D0*DEXP( DMIN1(DSIGN(17.269D0,                     &
       DBLE(TTT)-DBLE(TRPL)),DSIGN                                         &
       (21.875D0,DBLE(TTT)-DBLE(TRPL)))*DABS(DBLE(TTT)-DBLE(TRPL))/        &
       (DBLE(TTT)-35.86D0+DMAX1(0.D0,DSIGN                                 &
       (28.2D0,DBLE(TRPL)-DBLE(TTT)))))

!     FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT)
      FOQST(TTT,PRS) = DBLE(EPS1)/(DMAX1(1.D0,DBLE(PRS)/FOEW(TTT))-        &
       DBLE(EPS2))
!------------------------------------------------------------------------------!


      integer :: i,j,k,k1,k2,nkcld,n,step,lv,nlvs,tmin
      integer :: its,ite,kts,kte

      real PMv0,PMc0,PMr0,PMi0,PMg0,PMs0,PMh0,PMv1,PMc1,PMr1,PMi1,PMg1,PMs1,PMh1,     &
           PMv3,PMc3,PMr3,PMi3,PMg3,PMs3,PMh3,rhoQmax,depth,gam1,gam3,lamda,No,cnt,   &
           dm,alpha,cxh,cmx,Nfact1,Dmx(nk),ref(nk),Zplt,H,H0,Hcld,ZePlt,HcldTop,      &
           wmaxH,c1,c2,esat,alfa,alfa0,walfa,Kdiff,dz,dzsq,INTthr,dum,HcldBase,       &
           LAMr,ZZ,AMPL,AMPL0,Dc,Dr,Di,Ds,Dg,Dh,cmr,cmi,cms,cmg,cmh,thrd,sig,Ltot,    &
           M1,M2,Aii,BASE,RATIOc,RATIOv,RATIOr,RATIOi,RATIOg,RATIOs,RATIOh,tsec,      &
           tminr,af,bf,cm,Cx,Nox,rhox,PR

      real, parameter    :: eps    = 0.622
      real, parameter    :: g_cld  = 9.81
      real, parameter    :: Rd_cld = 287.0
      real, parameter    :: T0     = 273.15
      real, parameter    :: cp_cld = 1005.

      integer, parameter :: ni     = 1
      integer, parameter :: lnmax  = 100
      integer, parameter :: nt     = nint(ttotmin*60/dt)

      real, dimension(nk)    :: z,p,tt,td,tt2,w1,DIV,Tdry,zcld,w1cld,rho,alfa2
      real, dimension(ni,nk) :: tt0,tt1,Qv0,Qv1,w,SIGMA,Qsat,womega,th2d0,th2d1,p2d,dz2d
      real, dimension(lnmax) :: Pin,Zin,TTin,TDin

    ! Prognostic hydrometeor variables:
      !liquid:
      real, dimension(ni,nk) :: Qc0,Qc1,Nc0,Nc1,Qr0,Qr1,Nr0,Nr1,ssat1  !ssat0
      !ice:
      real, dimension(ni,nk,n_iceCat) :: Qi0,Qi1,Qg0,Qg1,Ni0,Zi0,Ni1,Bg0,Bg1,Zi1,Ql1,Ql0

    ! Source-Sink term arrays:
      integer, parameter               :: n_diag_2d = 20
      integer, parameter               :: n_diag_3d = 20
      real, dimension(ni,n_diag_2d)    :: diag_2d     !user-defined 2D diagnostic arrays (for output)
      real, dimension(ni,nk,n_diag_3d) :: diag_3d     !user-defined 3D diagnostic arrays (for output)
      real, dimension(ni,nk,n_iceCat)  :: diag_reffi,diag_vmi,diag_di,diag_rhoi,diag_dhmax,diag_vni,diag_vzi
      real, dimension(ni,nk)           :: diag_reffc

    ! Precipitation rates:
!       real, dimension(ni)      :: prt_liq,prt_sol,prt_drzl,prt_rain,prt_crys,prt_snow, &
!                                   prt_grpl,prt_pell,prt_hail,prt_sndp

      real, dimension(ni)      :: prt_liq       ! precip rate, total liquid     m s-1
      real, dimension(ni)      :: prt_sol       ! precip rate, total solid      m s-1
      real, dimension(ni)      :: prt_drzl      ! precip rate, drizzle          m s-1
      real, dimension(ni)      :: prt_rain      ! precip rate, rain             m s-1
      real, dimension(ni)      :: prt_crys      ! precip rate, ice cystals      m s-1
      real, dimension(ni)      :: prt_snow      ! precip rate, snow             m s-1
      real, dimension(ni)      :: prt_grpl      ! precip rate, graupel          m s-1
      real, dimension(ni)      :: prt_pell      ! precip rate, ice pellets      m s-1
      real, dimension(ni)      :: prt_hail      ! precip rate, hail             m s-1
      real, dimension(ni)      :: prt_sndp      ! precip rate, unmelted snow    m s-1
      real, dimension(ni)      :: prt_wsnow     ! precip rate, very wet snow    m s-1

    ! Diagnostics, etc.
      real, dimension(ni,nk)   :: diag_ZET,Qvinit,GZ,scpf_cldfrac
      real, dimension(nk)      :: COMP
      real, dimension(ni)      :: p0,p0m,LR,SR,diag_ZEC
      real                     :: delz,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real, dimension(ni,nk,6) :: qi_type

      integer                  :: nk_read,kk,stat
      logical, parameter       :: log_predictNc = .true.
      integer, dimension(200)  :: kskiplevs

      real                     :: time1,time2,time3,time4,time5 ! for timing tests

!---------------------------------------------------------------------------------------!

      print*
      print*, '** Remember to use compiler debug options for testing new code (modfify Makefile appropriately) **'
      print*

      diag_dhmax = -1.

      open (unit=30, file='out_p3.dat')
! for bit-pattern
!      open (unit=30, file='out_p3_bit.dat')

      its = 1
      ite = 1
      kts = 1
      kte = nk

      ssat1 = 0.  !for P3_v2.3.3, prognostic supersaturation is not available

      PR = 0. ! accumulated precipitation
      time3 = 0. ! for total timing test

!------- INITIALIZE w, tt, td, p AND Q[x]  PROFILES: -------

      H= 12000.
!     open(unit=12, file='./soundings/snd_input.Alta_3.data')   !grnd=  0m
!     open(unit=12, file='./Tsfc2.txt')
      open(unit=12, file='./soundings/snd_input.KOUN_00z1june2008.data')

!  Input uninterpolated sounding data:
      read(12,*) nlvs
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      do lv=1,nlvs
         !read(12,*) dum,Pin(lv),Zin(lv),TTin(lv),TDin(lv)
         read(12,*) Pin(lv),Zin(lv),TTin(lv),TDin(lv)   !new
         TTin(lv) = TTin(lv) + T0		        !convert from C to Kelvins
	 TDin(lv) = TDin(lv) + T0
         !if(sndcase=='ALBERTA') TDin(lv) = TDin(lv) + T0 !        " "
         !if(sndcase=='HARP')    TDin(lv) = TDin(lv)      !'TD' is actually Qv
         Pin(lv)  = Pin(lv)*100                          !convert from mb to Pa
      enddo
      close(12)

!-- Set up vertical level (evenly-spaced z-levels)
!       RE-CODE if necessary

!!-- Read in levels from file:
!!     open(unit=20, file='./levels/levs_86.dat')
!!     open(unit=20, file='./levels/levs_62.dat')
      open(unit=20, file='./levels/levs_41.dat')

      read(20,*) nk_read
      if (nk /= nk_read) then
         print*, '*** Abort in CLD1D ***'
         print*, 'Mismatch in nk specified for arrays and nk from table of levels: ',nk,nk_read
         stop
      endif
      read(20,*)
      read(20,*)
      do kk = 1,nk_read
         read(20,*) k,dum1,z(k),dum1,dum2
      enddo
      close(20)

!-- Compute dz2d(1,k) for upwind sedimentation.
!    note: this overrides the "DELTA_Z" in the file, which is not appropriate for P3
      if (trim(model)=='WRF') then
         do k = 1,nk-1
            dz2d(1,k) = z(k+1) - z(k)
         enddo
         dz2d(1,nk) = dz2d(1,nk-1)  !copy second highest DZ to highest
      else
         do k = 2,nk
            dz2d(1,k) = z(k-1) - z(k)
         enddo
         dz2d(1,1) = dz2d(1,2)  !copy second highest DZ to highest
      endif

!  Set up vertical level (evenly-spaced z-levels)
!      z(nk) = Zin(1)
!      delz = (H-Zin(1))/(nk-1)
!      do k=1,nk
!         z(k) = Zin(1) + delz*(nk-k)
!         dz2d(1,k) = delz
!      enddo
!      dz = z(1)-z(2)
!      dzsq = dz**2.
!      H0= z(nk)
!      GZ(1,:)= z*GRAV

!  Interpolate p,tt,td from sounding data and initialize Q[x]:
      do k=1,nk
         call vertint2b( p(k),z(k),Pin, Zin,nk,nlvs)
         call vertint2b(tt(k),z(k),TTin,Zin,nk,nlvs)
         call vertint2b(td(k),z(k),TDin,Zin,nk,nlvs)

         rho(k)= p(k)/(Rd_cld*tt(k))
!        rho(k)= 1.
         p2d(1,k) = p(k)

         Qsat(1,k)= FOQST(tt(k),p(k))
         if(sndcase=='ALBERTA') Qv0(1,k) = FOQST(td(k),p(k))  !Alberta (Td read in)
         Qvinit(1,k) = Qv0(1,k)

         Qc0(1,k) = 0.
         Nc0(1,k) = 0.
         Qr0(1,k) = 0.
         Nr0(1,k) = 0.
!        ssat0(1,k) = 0.   !supersaturation mixing ratio (qv-qv_sat)

         Qi0(1,k,:) = 0.   !qitot (total ice mass mixing ratio)
         Qg0(1,k,:) = 0.   !qirim (riming mass mixing ratio)
         Ql0(1,k,:) = 0.   !qiliq (liquid on ice mass mixing ratio)
         Ni0(1,k,:) = 0.   !nitot (number mixing ratio)
         Bg0(1,k,:) = 0.   !birim (rime volume mixing ratio)
         Zi0(1,k,:) = 0.   !zitot (reflectivity mixing ratio)

         tt0(1,k) = tt(k)
         tt1(1,k) = tt(k)

         Dmx(k) = 0.
      enddo

      p0(1)   = 100000.    !surface air pressure [Pa]
      p0m(1)  = p0(1)
      do k = 1,nk
        SIGMA(1,k) = p(k)/p(nk)
      enddo

!  Define initial w-profiles (and profiles related to w):
      H0 = 0.
      Hcld = H
      nkcld= nk
      if (.not.EVOLVING) then
        call NewWprof2(w,w1,AMPB,H,DIV,COMP,rho,dt,z,nk)
        do k=1,nk
          womega(:,k)= -w(:,k)*rho(k)*GRAV
        enddo
      endif
!**  z, w, p, tt, td and Qx0 are now in arrays from 1 to nk  **


!------ Sedimentation test:
!    !prescribe initial hydrometeor values at selected levels

!      k = 2
!      Qc0(1,k) = 0.2e-2;  Nc0(1,k) = 0.1e+8


!    Qr0(1,10) = 0.2e-2;  Nr0(1,10) = 0.1e+6
!    Qr0(1,11) = 0.2e-2;  Nr0(1,11) = 0.1e+5
!    Qr0(1,12) = 0.2e-2;  Nr0(1,12) = 0.1e+6

!       k=2
!       Qr0(1,k) = 0.1e-2;  Nr0(1,k) = 0.1e+5
!       k=3
!       Qr0(1,k) = 0.2e-2;  Nr0(1,k) = 0.1e+5
!       k=4
!       Qr0(1,k) = 0.2e-2;  Nr0(1,k) = 0.1e+5
!       k=5
!       Qr0(1,k) = 0.1e-2;  Nr0(1,k) = 0.1e+5


!       k=3
!       Qr0(1,k) = 0.1e-4;  Nr0(1,k) = 0.1e+6


!       k = 2
!       Qi0(1,k,1) = 5.e-3
!       Ni0(1,k,1) = 1.e+3
!       Qg0(1,k,1) = Qi0(1,k,1)
!       Bg0(1,k,1) = Qg0(1,k,1)/900.
! !     Zi0(1,k,1) = 1.e-12
!       Zi0(1,k,1) = 1. !use a large initial value; zi-limit will adjust it such that mu_i = 0
! !     Zi0(1,k,1) = 0. !use a tiny (or zero) initial value; zi-limit will adjust it such that mu_i = 8 (maximium)

! full column:
!       Qi0(1,15:nk-20,1) = 5.e-3
!       Ni0(1,15:nk-20,1) = 1.e+3
!       Qg0(1,15:nk-20,1) = Qi0(1,15:nk-20,1)
!       Bg0(1,15:nk-20,1) = Qg0(1,15:nk-20,1)/900.
! !     Zi0(1,15:nk-20,1) = 1.e-12
!       Zi0(1,15:nk-20,1) = 1. !use a large initial value; zi-limit will adjust it such that mu_i = 0
! !     Zi0(1,15:nk-20,1) = 0. !use a tiny (or zero) initial value; zi-limit will adjust it such that mu_i = 8 (maximium)

!         k = 3
!         Qi0(1,k,1) = 5.e-3
!
!         k = 12
!         Qi0(1,k,1) = 5.e-3

!       k = 3
!       Qi0(1,k,1) = 5.e-3
!       Ni0(1,k,1) = 1.e-6
!       Qg0(1,k,1) = Qi0(1,k,1)
!       Bg0(1,k,1) = Qg0(1,k,1)/900.
!======

      do k = 1,nk
         Qsat(1,k)= 0.
         Qv1(1,k)= Qv0(1,k)
         Qc1(1,k)= Qc0(1,k)
         Nc1(1,k)= Nc0(1,k)
         Qr1(1,k)= Qr0(1,k)
         Nr1(1,k)= Nr0(1,k)
!        ssat1(1,k) = ssat0(1,k)

         Qi1(1,k,:)= Qi0(1,k,:)
         Qg1(1,k,:)= Qg0(1,k,:)
         Ql1(1,k,:)= Ql0(1,k,:)
         Ni1(1,k,:)= Ni0(1,k,:)
         Bg1(1,k,:)= Bg0(1,k,:)
         Zi1(1,k,:)= Zi0(1,k,:)
      enddo

!-------------------------------------------------------------------------!
!  END OF INITIAL SET-UP

! !     !-- write header for output file
! !       write(30,'(1f4.0)') dt
! !       write(30,'(1i4)')   nt
! !       write(30,'(1i4)')   nk
! !       write(30,'(1i4)')   n_iceCat
! !     !==

!     call P3_INIT('./lookup_tables/',n_iceCat)                  !v2.8.4
!     call P3_INIT('./lookup_tables/',n_iceCat,stat)             !v2.9.1, v2.10.1
!     call P3_INIT('./lookup_tables/',n_iceCat,1.,1.)            !v3.0.0
!     call P3_INIT('./lookup_tables/',n_iceCat,stat)             !v2.9.1, v2.10.1, v3.1.0, v3.1.1
!     call P3_INIT('./lookup_tables/',n_iceCat,trplMomIce,stat)  !v4.0.0
!     call P3_INIT('./lookup_tables/',override_path=my_LT_path,n_iceCat,trplMomIce,stat)  !v4.0.0_b38

      call P3_INIT(LT_path,n_iceCat,trplMomIce,liqFrac,model,stat,abort_on_err,dowr)      !v4.0.0_b42


      do k=1,nk
         write(100,'(1f10.0, 1f8.1, 8e16.6, 1f10.2, 1f10.1, 1f10.3)') z(k),diag_ZET(1,k),          &
                    Qc1(1,k), Nc1(1,k), Qr1(1,k), Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1),  &
                    Ql1(1,k,1), w(1,k)  , tt1(1,k), scpf_cldfrac(1,k)

!        write(200,'(1i6,1f10.0,1f8.1,6e16.6,2f10.2,2e16.6)') k,z(k),diag_ZET(1,k),Qr1(1,k), Nr1(1,k)
!        write(300,'(1f10.0,1f8.1,20e16.6)') z(k),diag_ZET(1,k),Qi1(1,k,1),Zi1(1,k,1)
       enddo


!=========================================================================!
!  MAIN TIME LOOP:
      print*, 'Starting main loop...'

      DO step = 1,nt

        tsec = step*dt         !integration time [s]    (real)
        tminr= tsec/60.        !integration time [min]  (real)
        tmin = nint(tminr)     !integration time [min]  (integer)

! !         write(30,'(1f8.0,3f10.3)') tsec,prt_liq(1),prt_sol(1),maxval(w)

! ******* CHO ***********
! Initialize de la neige au top (qi0)
!	do k = 2,2
!	   Qi1(1,k,:) = 0.000265
!	   Qg1(1,k,:) = 0.00017
!	   Bg1(1,k,:) = Qg1(1,k,:)/500.
!	   Ni1(1,k,:) = 5000
!
!	   Qi0(1,k,:) = Qi1(1,k,:)
!	   Qg0(1,k,:) = Qg1(1,k,:)
!	   Bg0(1,k,:) = Bg1(1,k,:)
!	   Ni0(1,k,:) = Ni1(1,k,:)
!	enddo
! ******* CHO ***********

!  tsec-dependent profiles: (evolving updraft only)
        if (EVOLVING) then
          wmaxH= AMPA + (AMPB-AMPA)*0.5*(cos((tsec/tscale1)*2.*pi+pi)+1.)
          if (tsec .gt. 3600.) wmaxH= 0.

          Hcld=  Htop0+ (H-Htop0)*0.5*(cos((tsec/tscale2)*2.*pi+pi)+1.)
          k=1  !Find nkcld:
          do while(Hcld<z(k))
            nkcld= k
            k=k+1
          enddo
          nkcld= min(nk,nk+1-nkcld)
!         call NewWprof2(w,w1,wmaxH,H,DIV,COMP,rho,dt,z,nk)     !t-dep W only
          call NewWprof2(w,w1,wmaxH,Hcld,DIV,COMP,rho,dt,z,nk)  !t-dep W & top
          do k=1,nk
            womega(:,k)= -w(:,k)*rho(k)*GRAV
          enddo
        endif

!       goto 300  !skip advection (testing only)

        call advec(tt1,w,H,H0,Hcld,nk,nkcld,dt) !T-adv
        call advec(Qv1,w,H,H0,Hcld,nk,nkcld,dt)
        call advec(Qc1,w,H,H0,Hcld,nk,nkcld,dt)
        call advec(Qr1,w,H,H0,Hcld,nk,nkcld,dt)
        call advec(Nc1,w,H,H0,Hcld,nk,nkcld,dt)
        call advec(Nr1,w,H,H0,Hcld,nk,nkcld,dt)
!       call advec(ssat1,w,H,H0,Hcld,nk,nkcld,dt)

        do iice = 1,n_iceCat
           call advec(Qi1(:,:,iice),w,H,H0,Hcld,nk,nkcld,dt)
           call advec(Qg1(:,:,iice),w,H,H0,Hcld,nk,nkcld,dt)
           call advec(Ql1(:,:,iice),w,H,H0,Hcld,nk,nkcld,dt)
           call advec(Ni1(:,:,iice),w,H,H0,Hcld,nk,nkcld,dt)
           call advec(Bg1(:,:,iice),w,H,H0,Hcld,nk,nkcld,dt)
           call advec(Zi1(:,:,iice),w,H,H0,Hcld,nk,nkcld,dt)
        enddo


        do k=1,nk
          Qv1(1,k)= amax1(0., Qv1(1,k)-Qv0(1,k)*(DIV(k)+COMP(k))*dt)
          Qc1(1,k)= amax1(0., Qc1(1,k)-Qc0(1,k)*(DIV(k)+COMP(k))*dt)
          Qr1(1,k)= amax1(0., Qr1(1,k)-Qr0(1,k)*(DIV(k)+COMP(k))*dt)
          Nc1(1,k)= amax1(0., Nc1(1,k)-Nc0(1,k)*(DIV(k)+COMP(k))*dt)
          Nr1(1,k)= amax1(0., Nr1(1,k)-Nr0(1,k)*(DIV(k)+COMP(k))*dt)
!         ssat1(1,k)= amax1(0., ssat1(1,k)-ssat0(1,k)*(DIV(k)+COMP(k))*dt)
!         ssat1(1,k)= ssat1(1,k)-ssat0(1,k)*(DIV(k)+COMP(k))*dt


          do iice = 1,n_iceCat
             Qi1(1,k,iice)= amax1(0., Qi1(1,k,iice)-Qi0(1,k,iice)*(DIV(k)+COMP(k))*dt)
             Qg1(1,k,iice)= amax1(0., Qg1(1,k,iice)-Qg0(1,k,iice)*(DIV(k)+COMP(k))*dt)
             Ql1(1,k,iice)= amax1(0., Ql1(1,k,iice)-Ql0(1,k,iice)*(DIV(k)+COMP(k))*dt)
             Bg1(1,k,iice)= amax1(0., Bg1(1,k,iice)-Bg0(1,k,iice)*(DIV(k)+COMP(k))*dt)
             Ni1(1,k,iice)= amax1(0., Ni1(1,k,iice)-Ni0(1,k,iice)*(DIV(k)+COMP(k))*dt)
             Zi1(1,k,iice)= amax1(0., Zi1(1,k,iice)-Zi0(1,k,iice)*(DIV(k)+COMP(k))*dt)
          enddo

!       Adiabatic (cooling due to ascent):
          tt1(1,k)= tt1(1,k) + dt*(-g_cld/cp_cld*w(1,k))
        enddo

 300    continue

!------------------------------------------------------------------------!
       IF (microON) THEN

          th2d0(1,:) = tt0(1,:)*(1.e+5/p(:))**0.286
          th2d1(1,:) = tt1(1,:)*(1.e+5/p(:))**0.286

!v5.2.0

          call cpu_time(time1)

          if (.not. trplMomIce) Zi1(:,:,:) = 0.
          if (.not. liqFrac)    Ql1(:,:,:) = 0.

          CALL P3_MAIN(Qc1,Nc1,Qr1,Nr1,th2d0,th2d1,Qv0,Qv1,dt,Qi1,Qg1,Ql1,Ni1,Bg1,Zi1,         &
                             ssat1,w,p2d,dz2d,step,prt_liq,prt_sol,its,ite,kts,kte,n_iceCat,   &
                             diag_ZET,diag_reffc,diag_reffi,diag_vmi,diag_di,diag_rhoi,        &
                             n_diag_2d,diag_2d,n_diag_3d,diag_3d,log_predictNc,                &
                             trim(model),clbfact_dep,clbfact_sub,debug_on,                     &
                             scpf_on,scpf_pfrac,scpf_resfact,scpf_cldfrac,trplMomIce,liqFrac,  &
                             prt_drzl = prt_drzl,  &
                             prt_rain = prt_rain,  &
                             prt_crys = prt_crys,  &
                             prt_snow = prt_snow,  &
                             prt_grpl = prt_grpl,  &
                             prt_pell = prt_pell,  &
                             prt_hail = prt_hail,  &
                             prt_sndp = prt_sndp,  &
                             prt_wsnow = prt_wsnow, &
                             qi_type  = qi_type,    &
                             diag_dhmax = diag_dhmax)

          call cpu_time(time2)

          tt1(1,:) = th2d1(1,:)*(p(:)*1.e-5)**0.286

         !convert precipitation rates to units mm h-1 (from m s-1)
         prt_liq  = prt_liq *3.6e+6  !total liquid
         prt_sol  = prt_sol *3.6e+6  !total solid
         prt_drzl = prt_drzl*3.6e+6
         prt_rain = prt_rain*3.6e+6
         prt_crys = prt_crys*3.6e+6
         prt_snow = prt_snow*3.6e+6
         prt_grpl = prt_grpl*3.6e+6
         prt_pell = prt_pell*3.6e+6
         prt_hail = prt_hail*3.6e+6
         prt_sndp = prt_sndp*3.6e+6
         prt_wsnow = prt_wsnow*3.6e+6



!*-----------------------------------------------------------------------*!

      ENDIF   !(microON)

!  Rearrange arrays:
         DO k = 1,nk
            Qv0(1,k) = Qv1(1,k)
            tt0(1,k) = tt1(1,k)
            Qc0(1,k) = Qc1(1,k)
            Qr0(1,k) = Qr1(1,k)
            Nc0(1,k) = Nc1(1,k)
            Nr0(1,k) = Nr1(1,k)
            Qi0(1,k,:) = Qi1(1,k,:)
            Qg0(1,k,:) = Qg1(1,k,:)
            Ql0(1,k,:) = Ql1(1,k,:)
            Ni0(1,k,:) = Ni1(1,k,:)
            Bg0(1,k,:) = Bg1(1,k,:)
            Zi0(1,k,:) = Zi1(1,k,:)
!           ssat0(1,k) = ssat1(1,k)
         ENDDO
         p0m(1) = p0(1)

!--  Prevent low-level moisture depletion:     ----------------!
!    (Impose Td >= Td_initial after  xx min. [..< xx)] )
!        if (TDMIN .and. step<70) then
         if (TDMIN) then
           do k = 1,nk    ! or, e.g., k=nk-10,nk
             if (z(k)<1000.) then
!              Qv0(1,k)= amax1(Qv1(1,k),Qvinit(1,k))
               Qv0(1,k)= amax1(Qv1(1,k),Qvinit(1,k)*0.4) !test
               Qv1(1,k)= Qv0(1,k)
             endif
           enddo
         endif
!----------------------------------------------------------------!

         PR = PR + (prt_liq(1) + prt_sol(1))*dt*1.e-3  !accumulated precipitation, mm
         time3 =time3+(time2-time1)

!  Output to files:

!    Precipitation rates at lowest level (surface):
         write(41,*) tminr, prt_liq(1), prt_sol(1), PR,time3

! if (tminr.eq.1) then
!     print*, Zi1(1,2,1),Qi1(1,2,1),Ni1(1,2,1)
! endif

! if (tminr.eq.60) then
! 	print*, Qr1(1,50)*100/(Qi1(1,50,1)+Qr1(1,50)), Qi1(1,50,1)*100/(Qi1(1,50,1)+Qr1(1,50))
! endif


!	do k=2,50
!	if (tminr.eq.30) then
!		write(26,*) max(0.,Qg1(1,2,1)/(Qi1(1,2,1))), k, &
!				max(0.,Qg1(1,k,1)/(Qi1(1,k,1))), &
!				tt1(1,k)-T0,diag_vmi(1,k,1),diag_rhoi(1,k,1),diag_di(1,k,1)
!	endif
!	enddo

!	do k=2,50
!	if (tminr.eq.45) then
!		write(36,*) max(0.,Qg1(1,2,1)/(Qi1(1,2,1))), k,  &
!				max(0.,Qg1(1,k,1)/(Qi1(1,k,1))), &
!				tt1(1,k)-T0,diag_vmi(1,k,1),diag_rhoi(1,k,1),diag_di(1,k,1)
!	endif
!	enddo


!	do k=2,50
!	if (tminr.eq.60) then
!		write(28,*) max(0.,Qg1(1,2,1)/(Qi1(1,2,1))), k,  &
!				max(0.,Qg1(1,k,1)/(Qi1(1,k,1))), &
!				tt1(1,k)-T0,diag_vmi(1,k,1),diag_rhoi(1,k,1),diag_di(1,k,1)
!	endif
!	enddo

         if (mod(tminr,float(1)) < 1.e-5) print*, 'time (min): ',tminr,(prt_liq(1) + prt_sol(1)),PR, time3
        ! if (mod(tminr,float(5)) < 1.e-5) print*, 'time (min): ',tminr,prt_sol(1),PR

        ! if (mod(tminr,float(1)) < 1.e-5) print*, Qg1(1,50,1)*100/(Qi1(1,50,1)+Qr1(1,50)), &
	!			Qr1(1,50)*100/(Qi1(1,50,1)+Qr1(1,50)), 			&
	!			(Qi1(1,50,1)-Qg1(1,50,1))*100/(Qi1(1,50,1)+Qr1(1,50))


         IF (mod(tminr,float(outfreq)) < 1.e-5) THEN  !output only every OUTFREQ min   !for TESING


	!	write(27,*), tminr, Qg1(1,2,1)/(Qi1(1,2,1)), 	&
	!			Qr1(1,50)*100/(Qi1(1,50,1)+Qr1(1,50)), 		&
	!			Qi1(1,50,1)*100/(Qi1(1,50,1)+Qr1(1,50)),	&
	!			(Qg1(1,50,1))*100/(Qi1(1,50,1)+Qr1(1,50))

	   do k=1,nk

!    Hydrometeors mass content profiles:

! for bit-matching
!              write(20,'(1i5,2f7.0,11f15.5)') k,tminr,z(k),rho(k),p2d(1,k),	        &
!                        tt1(1,k)-T0, &
!			Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,	   		   	&
!			Qi1(1,k,1)*10**8,Qv1(1,k),diag_vmi(1,k,1),                      &
!			diag_rhoi(1,k,1),Qg1(1,k,1)/(Qi1(1,k,1))


!              write(20,'(1i5,1f7.0,16f15.5)') k,z(k),rho(k),p2d(1,k),	        &
!                        tt1(1,k)-T0, &
!			Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,Ql1(1,k,1)*10**8,	   	&
!			Qi1(1,k,1)*10**8,Qv1(1,k),diag_vmi(1,k,1),                      &
!			diag_rhoi(1,k,1),Qg1(1,k,1)/(Qi1(1,k,1)-Ql1(1,k,1)),            &
!	                Ql1(1,k,1)/Qi1(1,k,1),Qc1(1,k)

!	      write(24,'(1i5,1f7.0,5f10.5)') k,z(k),  		    	&
!		     Qc1(1,k), Qr1(1,k), Qi1(1,k,1),		    	&
!		     Qg1(1,k,1),Ql1(1,k,1)

!diag_di(1,k,1),diag_rhoi(1,k,1),  &
!                     Zi1(1,k,1)*10**8,diag_mui(1,k,1),diag_lami(1,k,1),       &
!                     diag_ZET(1,k),Qg1(1,k,1)



        !if (mod(tminr,float(1)) < 1.e-5) print*, 'time (min): ',tminr,(prt_liq(1) + prt_sol(1)),PR
        !if (mod(tminr,float(5)) < 1.e-5) print*, 'time (min): ',tminr,(prt_liq(1) + prt_sol(1)),PR

 !       IF (mod(tminr,float(outfreq)) < 1.e-5) THEN  !output only every OUTFREQ min   !for TESING
 !       IF (.true.) then   !for output to 'plot_1d.pro'

!	   do k=1,nk

!              write(300+tmin,'(i4,2f8.0,2e16.6,10f10.2)') k,diag_ss3d(1,k,10),diag_ss3d(1,k,11),Qr1(1,k),Nr1(1,k), &
!                 diag_ss3d(1,k,1),diag_ss3d(1,k,2),diag_ss3d(1,k,3)

!              write(300+tmin,'(i4,2f8.0,2e16.6,10f10.2)') k,diag_ss3d(1,k,10),diag_ss3d(1,k,11),Qi1(1,k,1),Ni1(1,k,1), &
!                 diag_ss3d(1,k,1),diag_ss3d(1,k,2),diag_ss3d(1,k,3)



!    Hydrometeors mass content profiles:
!              write(100+tmin,'(1f10.0, 1f8.1, 7e16.6, 1f10.2, 1f10.1, 1f10.3)') z(k),diag_ZET(1,k),             &
!                     Qc1(1,k), Nc1(1,k), Qr1(1,k), Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1),  &
!                     w(1,k)  , tt1(1,k), scpf_cldfrac(1,k)

! from P3_main:
! diag_3d(i,k,1) = mu_i
! diag_3d(i,:,2) = diag_vmi(i,:,iice)
! diag_3d(i,:,3) = diag_di(i,:,iice)     diag_di
! diag_3d(i,:,4) = diag_rhoi(i,:,iice)
! diag_3d(i,k,5) = V_zit(k)
! diag_3d(i,k,6) = V_qit(k)
! diag_3d(i,k,7) = V_nit(k)
! diag_3d(i,k,11) = f1pr22   ! lambda_i
! diag_3d(i,k,12) = f1pr23   ! mu_i

!               write(200+tmin,'(1f10.0, 1f8.1, 7e16.6, 1f10.2, 7e16.6)') z(k),diag_ZET(1,k),             &
!                      Qc1(1,k), Qr1(1,k), Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1), Zi1(1,k,1), diag_3d(1,k,1), &
!                      diag_3d(1,k,2),diag_3d(1,k,3),diag_3d(1,k,4) !,diag_3d(1,k,5),diag_3d(1,k,6),diag_3d(1,k,7)

!               write(200+tmin,'(1f10.0, 1f8.1, 7e16.6, 1f10.2, 2f10.3, 1f10.1, 17e16.6)') z(k),diag_ZET(1,k), &
!                      Qc1(1,k), Qr1(1,k), Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1), Zi1(1,k,1),         &
!                      1000.*diag_3d(1,k,3),diag_3d(1,k,2),diag_3d(1,k,4)

!              write(200+tmin,'(1f10.0, 1f8.1, 20e16.6)') z(k),diag_ZET(1,k), &
!                     Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1), Zi1(1,k,1), &
!                     diag_3d(1,k,1),diag_3d(1,k,2),diag_3d(1,k,3),diag_3d(1,k,4)

!                    Qi1(1,k,2), Qg1(1,k,2), Bg1(1,k,2), Ni1(1,k,2), Zi1(1,k,2)

!!compute n0:
!  dum = Ni1(1,k,1)*diag_3d(1,k,11)**(diag_3d(1,k,12)+1.)/gamma(diag_3d(1,k,12)+1.)

!              write(300+tmin,'(1f10.0, 1f8.1, 5e16.6, 2f10.2, 1e16.6, 3f10.2)') z(k),diag_ZET(1,k), &
!                     Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1),    &
!                     dum,                      &  ! n0_i
!                     diag_3d(1,k,1),           &  ! mu_i (3mom only)
!                     diag_mui(1,k,1),          &  ! mu_i, table
!                     diag_lami(1,k,1),         &  ! lambda_i
!                     diag_di(1,k,1) *1000. ,   &  ! di
!                     diag_dhmax(1,k,1)*1000. , &  ! Dh_max
!                     Qg1(1,k,1)/Qi1(1,k,1) ! , &  ! Frime



!               write(100+tmin,'(1f10.0,1f8.1,6e16.6,2f10.2)') z(k),diag_ZET(1,k),          &
!                       (1,k), Qr1(1,k), Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1),  &
!                      w(1,k)  , tt1(1,k) !, Qv1(1,k)

!               write(100+tmin,'(1f10.0,1f8.1,6e16.6,2f10.2,2e16.6)') z(k),diag_ZET(1,k),     &
!                      Qi1(1,k,1), Qg1(1,k,1), Bg1(1,k,1), Ni1(1,k,1)
!
!               write(200+tmin,'(1i6,1f10.0,1f8.1,6e16.6,2f10.2,2e16.6)') k,z(k),diag_ZET(1,k),Qr1(1,k), Nr1(1,k)

!               write(300+tmin,'(1f10.0,1f8.1,6e16.6,2f10.2,2e16.6)') z(k),diag_ZET(1,k),Qc1(1,k), Nc1(1,k)


!             if (prog_nc_ssat) write(200+tmin,'(1f10.0,3e16.6)') z(k),Qc1(1,k), Nc1(1,k), ssat1(1,k)

!               write(200+tmin,'(1f10.0,7e16.6)') z(k),diag_3d(1,k,11),diag_3d(1,k,12)
!
!               write(500+tmin,'(1f10.0,f9.3,2f8.1,9e16.6)') z(k),diag_3d(1,k,05),diag_3d(1,k,06)*1.e+6,diag_3d(1,k,07),diag_3d(1,k,13), &
!                                       diag_3d(1,k,14),diag_3d(1,k,15),diag_3d(1,k,16),diag_3d(1,k,17),diag_3d(1,k,18)

              !write(200+tmin,'(1f10.0,2e16.6)') z(k),Qr1(1,k), Nr1(1,k)

!                 write(200+tmin,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1, 1e18.4,1e13.3,1f9.3,3f9.1)')  &
! !               write(2000+step,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1, 1e18.4,1e13.3,1f9.3,3f9.1)') &
!                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
!                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
!                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2)

! from P3_main:
! diag_di(i,:,iice)
! diag_vmi(i,:,iice)
! diag_rhoi(i,:,iice)
! diag_3d(i,k,1) = mu_i
! diag_3d(i,k,5) = V_zit(k)
! diag_3d(i,k,6) = V_qit(k)
! diag_3d(i,k,7) = V_nit(k)

! for bit-pattern
!! 1-category:
!            if (n_iceCat.eq.1) then
!              dum1 = 0.
!              if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)  !rime fraction
!              dum3 = 0.
!              if (Nr1(1,k).gt.0) dum3 = (6.*Qr1(1,k)/(1000.*Nr1(1,k)*pi))**(1./3.)
!              write(30,'(1f10.0, 2f8.1, 4e13.4, 1f8.2, 4e13.4, 2f8.2, 1e13.4, 8e13.4)')   &
!                     z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),   &
!                      Qi1(1,k,1),Ni1(1,k,1),  &
!                      dum1,              &  ! F_rime
!                      diag_vmi(1,k,1),   &
!                      diag_rhoi(1,k,1),  &
!                      diag_di(1,k,1),    &  ! Di
!                      diag_dhmax(1,k,1), &  ! Dh_max
!                      diag_3d(1,k,3) ,   &  ! mu_i     (calcualted)
!                      diag_mui(1,k,1),   &  ! mu_i     (table)
!                      diag_lami(1,k,1),  &  ! lambda_i (table)
!                      tt1(1,k)-T0,       &  ! temperature (C)
!                      Nc1(1,k),Nr1(1,k), &
!                      Ni1(1,k,1),Bg1(1,k,1)*10**8, &
!                      Zi1(1,k,1),dum3,   &
!                      Qv1(1,k)
!             endif

              dum1 = 0.
              dum2 = 0.
              dum3 = 0.
              if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/(Qi1(1,k,1)-Ql1(1,k,1))  !rime fraction
              if (Qi1(1,k,1).gt.0.) dum2 = Ql1(1,k,1)/Qi1(1,k,1)  !liquid fraction
              if (Nr1(1,k).gt.0.) dum3 = (6.*Qr1(1,k)/(pi*1000.*Nr1(1,k)))**(1./3.)  !Drm
            !for nCat=1: (misc)
              write(30,'(26e13.4)')      &
                      z(k),              &  !col 0 (for python plot script)
                      w(1,k),            &  !col 1
                      diag_ZET(1,k),     &  !col 2
                      Qc1(1,k),          &  !col 3
                      Qr1(1,k),          &  !col 4
!                     Qi1(1,k,1),        &  !col 5
                      sum(Qi1(1,k,:)),   &  !col 5
                      0.,   &  !col 6
!                     diag_mui(1,k,1),   &  !col 6
                      diag_rhoi(1,k,1),  &  !col 7
                      diag_di(1,k,1),    &  !col 8
                      diag_dhmax(1,k,1), &  !col 9 Dh_max
                      prt_liq(1),        &  !col 10
                      prt_sol(1),        &   !col 11
                      Nc1(1,k),          &  !col 12
                      Nr1(1,k),          &  !col 13
                      sum(Ni1(1,k,:)),   &  !col 14
                      dum1,              &       !col 15
                      dum2,     &  !col 16
                      dum3,      & !col 17
                      tt1(1,k)-T0  ! col 18

! real outputs
!! 1-category:
!            if (n_iceCat.eq.1) then
!              dum1 = 0.
!              if ((Qi1(1,k,1)-Ql1(1,k,1)).gt.1.e-14) dum1 = Qg1(1,k,1)/(Qi1(1,k,1)-Ql1(1,k,1))  !rime fraction
!              dum2 = 0.
!              if (Qi1(1,k,1).gt.1.e-14) dum2 = Ql1(1,k,1)/Qi1(1,k,1)  !liquid fraction
!              dum3 = 0.
!              if (Nr1(1,k).gt.1.e-16) dum3 = (6.*Qr1(1,k)/(1000.*Nr1(1,k)*pi))**(1./3.)
!              write(30,'(1f10.0, 44e13.4)')   &
!                     z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),   &
!                     Qi1(1,k,1),Ni1(1,k,1),  &
!                      dum1,              &  ! F_rime
!                      diag_vmi(1,k,1),   &
!                      diag_rhoi(1,k,1),  &
!                      diag_di(1,k,1),    &  ! Di
!                      diag_dhmax(1,k,1), &  ! Dh_max
!                      Qg1(1,k,1) ,   &  !
!                      diag_mui(1,k,1),   &  ! mu_i     (table)
!                      diag_lami(1,k,1),  &  ! lambda_i (table)
!                      tt1(1,k)-T0,       &  ! temperature (C)
!                      Nc1(1,k),Nr1(1,k), &
!                      Ni1(1,k,1),(Bg1(1,k,1))*10**8, &
!                      Zi1(1,k,1),dum3,   &
!                      Qv1(1,k),Ql1(1,k,1),dum2, &
!                      diag_3d(1,k,1), &
!                      diag_3d(1,k,2),diag_3d(1,k,3), &
!                      diag_3d(1,k,4),diag_3d(1,k,5), diag_3d(1,k,6), &
!                      diag_3d(1,k,7),diag_3d(1,k,8), &
!                      diag_3d(1,k,9),diag_3d(1,k,10), &
!                      diag_3d(1,k,11),diag_3d(1,k,12), &
!                      diag_3d(1,k,13),diag_3d(1,k,14), &
!                      diag_3d(1,k,15),diag_3d(1,k,16), &
!                      diag_3d(1,k,17),diag_3d(1,k,18), &
!                      diag_3d(1,k,19),diag_3d(1,k,20)
!!             endif

! real outputs
!! 1-category:
!            if (n_iceCat.eq.1) then
!              dum1 = 0.
!              if ((Qi1(1,k,1)-Ql1(1,k,1)).gt.0.) dum1 = Qg1(1,k,1)/(Qi1(1,k,1)-Ql1(1,k,1))  !rime fraction
!              dum2 = 0.
!              if (Qi1(1,k,1).gt.0.) dum2 = Ql1(1,k,1)/Qi1(1,k,1)  !liquid fraction
!              dum3 = 0.
!              if (Nr1(1,k).gt.0) dum3 = (6.*Qr1(1,k)/(1000.*Nr1(1,k)*pi))**(1./3.)
!              write(30,'(1f10.0, 2f8.1, 4e13.4, 1f8.2, 4e13.4, 2f8.2, 1e13.4, 47e13.4)')   &
!                     z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),   &
!                      Qi1(1,k,1),Ni1(1,k,1),  &
!                      dum1,              &  ! F_rime
!                      diag_vmi(1,k,1),   &
!                      diag_rhoi(1,k,1),  &
!                      diag_di(1,k,1),    &  ! Di
!                      diag_dhmax(1,k,1), &  ! Dh_max
!                      diag_3d(1,k,3) ,   &  ! mu_i     (calcualted)
!                      diag_mui(1,k,1),   &  ! mu_i     (table)
!                      diag_lami(1,k,1),  &  ! lambda_i (table)
!                      tt1(1,k)-T0,       &  ! temperature (C)
!                      Nc1(1,k),Nr1(1,k), &
!                      Ni1(1,k,1),Bg1(1,k,1)*10**8, &
!                      Zi1(1,k,1),dum3,   &
!                      Qv1(1,k),diag_3d(1,k,4), &
!                      diag_3d(1,k,5),diag_3d(1,k,6), &
!                      diag_3d(1,k,7),diag_3d(1,k,8), &
!                      diag_3d(1,k,9),diag_3d(1,k,10), &
!                      diag_3d(1,k,11),diag_3d(1,k,12), &
!                      diag_3d(1,k,13),diag_3d(1,k,14), &
!                      diag_3d(1,k,15),diag_3d(1,k,16), &
!                      diag_3d(1,k,17),diag_3d(1,k,18), &
!                      diag_3d(1,k,19),diag_3d(1,k,20), &
!                      diag_3d(1,k,21),diag_3d(1,k,22), &
!                      diag_3d(1,k,23),diag_3d(1,k,24), &
!                      diag_3d(1,k,25),diag_3d(i,k,26), &
!                      diag_3d(1,k,27),diag_3d(1,k,28), &
!                      diag_3d(1,k,29),diag_3d(1,k,30), &
!                      diag_3d(1,k,31),diag_3d(i,k,32), &
!                      diag_3d(1,k,33),diag_3d(i,k,34), &
!                      Ql1(1,k,1),dum2,diag_3d(1,k,1),  &
!                      diag_3d(1,k,2),diag_3d(i,k,3),   &
!                      diag_3d(1,k,4),diag_3d(i,k,5)
!             endif


! ! ! ! 2-category:
! ! !             if (n_iceCat.eq.2) then
! ! !               dum1 = 0.;  dum2 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2)
! ! !              endif
! ! !
! ! !
! ! ! ! 3-category:
! ! !             if (n_iceCat.eq.3) then
! ! !               dum1 = 0.;  dum2 = 0.;  dum3 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               if (Qi1(1,k,3).gt.0.) dum3 = Qg1(1,k,3)/Qi1(1,k,3)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2),        &
! ! !                      Qi1(1,k,3),Ni1(1,k,3),dum3,diag_vmi(1,k,3),diag_di(1,k,3)*1.e+6,diag_rhoi(1,k,3)
! ! !              endif
! ! !
! ! ! ! 4-category:
! ! !             if (n_iceCat.eq.4) then
! ! !               dum1 = 0.;  dum2 = 0.;  dum3 = 0.;  dum4 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               if (Qi1(1,k,3).gt.0.) dum3 = Qg1(1,k,3)/Qi1(1,k,3)
! ! !               if (Qi1(1,k,4).gt.0.) dum4 = Qg1(1,k,4)/Qi1(1,k,4)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2),        &
! ! !                      Qi1(1,k,3),Ni1(1,k,3),dum3,diag_vmi(1,k,3),diag_di(1,k,3)*1.e+6,diag_rhoi(1,k,3),        &
! ! !                      Qi1(1,k,4),Ni1(1,k,4),dum3,diag_vmi(1,k,4),diag_di(1,k,4)*1.e+6,diag_rhoi(1,k,4)
! ! !              endif
! ! !
! ! ! ! 5-category:
! ! !             if (n_iceCat.eq.5) then
! ! !               dum1 = 0.;  dum2 = 0.;  dum3 = 0.;  dum4 = 0.;   dum5 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               if (Qi1(1,k,3).gt.0.) dum3 = Qg1(1,k,3)/Qi1(1,k,3)
! ! !               if (Qi1(1,k,4).gt.0.) dum4 = Qg1(1,k,4)/Qi1(1,k,4)
! ! !               if (Qi1(1,k,5).gt.0.) dum5 = Qg1(1,k,5)/Qi1(1,k,5)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2),        &
! ! !                      Qi1(1,k,3),Ni1(1,k,3),dum3,diag_vmi(1,k,3),diag_di(1,k,3)*1.e+6,diag_rhoi(1,k,3),        &
! ! !                      Qi1(1,k,4),Ni1(1,k,4),dum3,diag_vmi(1,k,4),diag_di(1,k,4)*1.e+6,diag_rhoi(1,k,4),        &
! ! !                      Qi1(1,k,5),Ni1(1,k,5),dum3,diag_vmi(1,k,5),diag_di(1,k,5)*1.e+6,diag_rhoi(1,k,5)
! ! !              endif
! ! !
! ! ! ! 6-category:
! ! !             if (n_iceCat.eq.6) then
! ! !               dum1 = 0.;  dum2 = 0.;  dum3 = 0.;  dum4 = 0.;   dum5 = 0.;   dum6 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               if (Qi1(1,k,3).gt.0.) dum3 = Qg1(1,k,3)/Qi1(1,k,3)
! ! !               if (Qi1(1,k,4).gt.0.) dum4 = Qg1(1,k,4)/Qi1(1,k,4)
! ! !               if (Qi1(1,k,5).gt.0.) dum5 = Qg1(1,k,5)/Qi1(1,k,5)
! ! !               if (Qi1(1,k,6).gt.0.) dum6 = Qg1(1,k,6)/Qi1(1,k,6)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2),        &
! ! !                      Qi1(1,k,3),Ni1(1,k,3),dum3,diag_vmi(1,k,3),diag_di(1,k,3)*1.e+6,diag_rhoi(1,k,3),        &
! ! !                      Qi1(1,k,4),Ni1(1,k,4),dum3,diag_vmi(1,k,4),diag_di(1,k,4)*1.e+6,diag_rhoi(1,k,4),        &
! ! !                      Qi1(1,k,5),Ni1(1,k,5),dum3,diag_vmi(1,k,5),diag_di(1,k,5)*1.e+6,diag_rhoi(1,k,5),        &
! ! !                      Qi1(1,k,6),Ni1(1,k,6),dum3,diag_vmi(1,k,6),diag_di(1,k,6)*1.e+6,diag_rhoi(1,k,6)
! ! !              endif
! ! !
! ! ! ! 7-category:
! ! !             if (n_iceCat.eq.7) then
! ! !               dum1 = 0.;  dum2 = 0.;  dum3 = 0.;  dum4 = 0.;   dum5 = 0.;   dum6 = 0.;   dum7 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               if (Qi1(1,k,3).gt.0.) dum3 = Qg1(1,k,3)/Qi1(1,k,3)
! ! !               if (Qi1(1,k,4).gt.0.) dum4 = Qg1(1,k,4)/Qi1(1,k,4)
! ! !               if (Qi1(1,k,5).gt.0.) dum5 = Qg1(1,k,5)/Qi1(1,k,5)
! ! !               if (Qi1(1,k,6).gt.0.) dum6 = Qg1(1,k,6)/Qi1(1,k,6)
! ! !               if (Qi1(1,k,7).gt.0.) dum6 = Qg1(1,k,7)/Qi1(1,k,7)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2),        &
! ! !                      Qi1(1,k,3),Ni1(1,k,3),dum3,diag_vmi(1,k,3),diag_di(1,k,3)*1.e+6,diag_rhoi(1,k,3),        &
! ! !                      Qi1(1,k,4),Ni1(1,k,4),dum3,diag_vmi(1,k,4),diag_di(1,k,4)*1.e+6,diag_rhoi(1,k,4),        &
! ! !                      Qi1(1,k,5),Ni1(1,k,5),dum3,diag_vmi(1,k,5),diag_di(1,k,5)*1.e+6,diag_rhoi(1,k,5),        &
! ! !                      Qi1(1,k,6),Ni1(1,k,6),dum3,diag_vmi(1,k,6),diag_di(1,k,6)*1.e+6,diag_rhoi(1,k,6),        &
! ! !                      Qi1(1,k,7),Ni1(1,k,7),dum3,diag_vmi(1,k,7),diag_di(1,k,7)*1.e+6,diag_rhoi(1,k,7)
! ! !              endif
! ! !
! ! ! ! 8-category:
! ! !             if (n_iceCat.eq.8) then
! ! !               dum1 = 0.;  dum2 = 0.;  dum3 = 0.;  dum4 = 0.;   dum5 = 0.;   dum6 = 0.;   dum7 = 0.;   dum8 = 0.
! ! !               if (Qi1(1,k,1).gt.0.) dum1 = Qg1(1,k,1)/Qi1(1,k,1)
! ! !               if (Qi1(1,k,2).gt.0.) dum2 = Qg1(1,k,2)/Qi1(1,k,2)
! ! !               if (Qi1(1,k,3).gt.0.) dum3 = Qg1(1,k,3)/Qi1(1,k,3)
! ! !               if (Qi1(1,k,4).gt.0.) dum4 = Qg1(1,k,4)/Qi1(1,k,4)
! ! !               if (Qi1(1,k,5).gt.0.) dum5 = Qg1(1,k,5)/Qi1(1,k,5)
! ! !               if (Qi1(1,k,6).gt.0.) dum6 = Qg1(1,k,6)/Qi1(1,k,6)
! ! !               if (Qi1(1,k,7).gt.0.) dum6 = Qg1(1,k,7)/Qi1(1,k,7)
! ! !               if (Qi1(1,k,8).gt.0.) dum6 = Qg1(1,k,8)/Qi1(1,k,8)
! ! !               write(30,'(1f10.0,2f8.1,2e13.3, 1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1,    &
! ! !                                               1e18.4,1e13.3,1f9.3,3f9.1)')  &
! ! !                      z(k),w(1,k),diag_ZET(1,k),Qc1(1,k),Qr1(1,k),                                              &
! ! !                      Qi1(1,k,1),Ni1(1,k,1),dum1,diag_vmi(1,k,1),diag_di(1,k,1)*1.e+6,diag_rhoi(1,k,1),        &
! ! !                      Qi1(1,k,2),Ni1(1,k,2),dum2,diag_vmi(1,k,2),diag_di(1,k,2)*1.e+6,diag_rhoi(1,k,2),        &
! ! !                      Qi1(1,k,3),Ni1(1,k,3),dum3,diag_vmi(1,k,3),diag_di(1,k,3)*1.e+6,diag_rhoi(1,k,3),        &
! ! !                      Qi1(1,k,4),Ni1(1,k,4),dum3,diag_vmi(1,k,4),diag_di(1,k,4)*1.e+6,diag_rhoi(1,k,4),        &
! ! !                      Qi1(1,k,5),Ni1(1,k,5),dum3,diag_vmi(1,k,5),diag_di(1,k,5)*1.e+6,diag_rhoi(1,k,5),        &
! ! !                      Qi1(1,k,6),Ni1(1,k,6),dum3,diag_vmi(1,k,6),diag_di(1,k,6)*1.e+6,diag_rhoi(1,k,6),        &
! ! !                      Qi1(1,k,7),Ni1(1,k,7),dum3,diag_vmi(1,k,7),diag_di(1,k,7)*1.e+6,diag_rhoi(1,k,7),        &
! ! !                      Qi1(1,k,8),Ni1(1,k,8),dum3,diag_vmi(1,k,8),diag_di(1,k,8)*1.e+6,diag_rhoi(1,k,8)
! ! !              endif


!*---------------------------------------------------------*!
            enddo  ! k-loop

         ENDIF    ! outfreq

      ENDDO  !main time loop

      close(30)

!----------------------------------------------------------!
!            ---  End of time integreation  ---            !
!----------------------------------------------------------!
      print*
      print*, 'DONE'

end subroutine columnmodel
!=============================================================================!
