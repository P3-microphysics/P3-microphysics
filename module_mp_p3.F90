MODULE MP_P3

!__________________________________________________________________________________________
! This module contains the Predicted Particle Property (P3) bulk microphysics scheme.      !
!                                                                                          !
! This code was originally written by H. Morrison,  MMM Division, NCAR (Dec 2012).         !
! Modification were made by J. Milbrandt, RPN, Environment Canada (July 2014).             !
!                                                                                          !
! For details see: Morrison and Milbrandt (2015) [J. Atmos. Sci., 72, 287-311]             !
!                  Milbrandt and Morrison (2016) [J. Atmos. Sci., 73, 975-995]             !
!                                                                                          !
! For questions or bug reports, please contact:                                            !
!    Hugh Morrison   (morrison@ucar.edu)        or                                         !
!    Jason Milbrandt (jason.milbrandt@canda.ca)                                            !
!                                                                                          !
! --------------------------                                                               !
! Version:        v2.0                                                                     !
! Last modified:  2015-07-05                                                               !
!__________________________________________________________________________________________!

 implicit none

 real, parameter :: pi = 3.14159265

 public  :: mp_p3_wrapper_wrf,mp_p3_wrapper_gem,p3_main
 public  :: polysvp1
 private :: gamma,derf,pi,isnan

! integer switch for warm rain autoconversion/accretion schemes
 integer, private :: iparam

! ice microphysics lookup table array dimensions
 integer, parameter :: isize     = 20
 integer, parameter :: jsize     = 20
 integer, parameter :: iisize    = 12
 integer, parameter :: jjsize    = 12
 integer, parameter :: densize   =  5
 integer, parameter :: rimsize   =  4
 integer, parameter :: rcollsize = 30

! number of ice microphysical quantities used from lookup table
 integer, parameter :: tabsize   = 12

! number of ice-rain collection microphysical quantities used from lookup table
 integer, parameter :: colltabsize = 2

! number of ice-ice collection microphysical quantities used from lookup table
 integer, parameter :: collitabsize = 2

 private :: isize,jsize,densize,rimsize,rcollsize,tabsize,colltabsize,collitabsize,      &
            iisize,jjsize

!ice lookup table values
 real, private, dimension(densize,rimsize,isize,jsize,tabsize) :: itab

!ice lookup table values for ice-rain collision/collection
 double precision, private, dimension(densize,rimsize,isize,jsize,rcollsize,colltabsize)    :: itabcoll
! below is separated into itabcolli1 and itabcolli2, because there is a max of
! 7 dimensions for arrays on some fortran compilers
 double precision, private, dimension(iisize,rimsize,densize,iisize,jjsize,rimsize,densize) :: itabcolli1
 double precision, private, dimension(iisize,rimsize,densize,iisize,jjsize,rimsize,densize) :: itabcolli2

! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
! warm rain autoconversion/accretion option only (iparam = 1)
 real, private, dimension(16) :: dnu

! physical constants
 real, private :: rhosur,rhosui,ar,br,f1r,f2r,ecr,qsmall,rhow,kr,kc,bimm,aimm,rin,mi0,   &
                  eci,eri,eii,bcn,cpw,e0,cons1,cons2,cons3,cons4,cons5,cons6,cons7,      &
                  nccnst,inv_rhow

! physical constants
 real, private :: cp,g,rd,rv,ep_2,inv_cp

! aerosol and droplet activation parameters
 real, private :: mw,osm,vi,epsm,rhoa,map,ma,rr,bact,inv_rm1,sig1,nanew1,f11,f21,       &
                  inv_rm2,sig2,nanew2,f12,f22,inv_bact

! lookup table values for rain shape parameter mur
 real, private, dimension(150) :: mur_table

! lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
 real, private, dimension(300,10) :: vn_table,vm_table,revap_table

 contains


!==========================================================================================!
 SUBROUTINE p3_init(lookup_file_1,lookup_file_2)

!------------------------------------------------------------------------------------------!
! This subroutine initializes all physical constants and parameters needed by the P3       !
! scheme.  The subroutine 'P3_INIT' must be called at the first model time step from there !
! wrapper subroutine, prior to first call to the main scheme subroutine, 'P3_MAIN'.        !
!------------------------------------------------------------------------------------------!

 implicit none

! Passed arguments:
 character*(*), intent(in) :: lookup_file_1    !lookup table for main processes
 character*(*), intent(in) :: lookup_file_2    !lookup table for ice category interactions

! Local variables/parameters:
 integer         :: i,j,k,ii,jj,kk,jjj,jjjj,jjj2,jjjj2
 real            :: lamr,mur,lamold,dum,initlamr,dm,dum1,dum2,dum3,dum4,dum5,dd,amg,vt,  &
                    dia,vn,vm
 real, parameter :: thrd     = 1./3.
 real, parameter :: sxth     = 1./6.
 real, parameter :: piov6rho = pi/6.*997.

!------------------------------------------------------------------------------------------!

 print*
 print*, 'P3_INIT (READING/CREATING P3 LOOK-UP TABLES):'

! switch for warm-rain parameterization
! = 1 seifert and beheng 2001
! = 2 beheng 1994
! = 3 khairoutdinov and kogan 2000
 iparam = 3

! droplet concentration (m-3)
 nccnst = 250.e6

! parameters for Seifert and Beheng (2001) autoconversion/accretion
 kc     = 9.44e9
 kr     = 5.78e3

! physical constants
 cp     = 1005.
 inv_cp = 1./cp
 g      = 9.816
 rd     = 287.15
 rv     = 461.51
 ep_2   = 0.622
 rhosur = 100000./(rd*273.15)
 rhosui = 60000./(rd*253.15)
 ar     = 841.99667
 br     = 0.8
 f1r    = 0.78
 f2r    = 0.32
 ecr    = 1.
 qsmall = 1.e-14
 rhow   = 997.
 cpw    = 4218.

 inv_rhow = 1.e-3  !inverse of (max.) density of liquid water

! Bigg (1953)
!bimm   = 100.
!aimm   = 0.66
! Barklie and Gokhale (1959)
 bimm   = 2.
 aimm   = 0.65
 rin    = 0.1e-6
 mi0    = 4./3.*pi*900.*1.e-18

 eci    = 0.5
 eri    = 1.
 eii    = 0.1
 bcn    = 2.

! saturation pressure at T = 0 C
 e0     = polysvp1(273.15,0)
 cons1  = pi/6.*rhow
 cons2  = 4./3.*pi*rhow
 cons3  = 1./(cons2*(25.e-6)**3)
 cons4  = 1./((550.e-6)**3*pi*rhow)
 cons5  = pi/6.*bimm
 cons6  = pi*pi/36.*rhow*bimm
 cons7  = 4./3.*pi*rhow*(1.e-6)**3

! aerosol/droplet activation parameters
 mw     = 0.018
 osm    = 1.
 vi     = 3.
 epsm   = 0.9
 rhoa   = 1777.
 map    = 0.132
 ma     = 0.0284
 rr     = 8.3187
 bact   = vi*osm*epsm*mw*rhoa/(map*rhow)
 inv_bact = (map*rhow)/(vi*osm*epsm*mw*rhoa)

! mode 1
 inv_rm1= 2.e+8
 sig1   = 2.0
 nanew1 = 300.e6
 f11    = 0.5*exp(2.5*(log(sig1))**2)
 f21    = 1.+0.25*log(sig1)

! mode 2
 inv_rm2= 1./1.3e-6
 sig2   = 2.5
 nanew2 = 0.
 f12    = 0.5*exp(2.5*(log(sig2))**2)
 f22    = 1.+0.25*log(sig2)

! parameters for droplet mass spectral shape, used by Seifert and Beheng (2001)
! warm rain scheme only (iparam = 1)
 dnu(1)  = 0.
 dnu(2)  = -0.557
 dnu(3)  = -0.430
 dnu(4)  = -0.307
 dnu(5)  = -0.186
 dnu(6)  = -0.067
 dnu(7)  = 0.050
 dnu(8)  = 0.167
 dnu(9)  = 0.282
 dnu(10) = 0.397
 dnu(11) = 0.512
 dnu(12) = 0.626
 dnu(13) = 0.739
 dnu(14) = 0.853
 dnu(15) = 0.966
 dnu(16) = 0.966

!------------------------------------------------------------------------------------------!
! read in ice microphysics table

 open(unit=10,file=lookup_file_1, status='old')

 do jj = 1,densize
    do ii = 1,rimsize
       do i = 1,isize
          do k = 1,jsize
             read(10,*) dum,dum,dum,dum,itab(jj,ii,i,k,1),itab(jj,ii,i,k,2),             &
                  itab(jj,ii,i,k,3),itab(jj,ii,i,k,4),itab(jj,ii,i,k,5),                 &
                  itab(jj,ii,i,k,6),itab(jj,ii,i,k,7),itab(jj,ii,i,k,8),dum,             &
                  itab(jj,ii,i,k,9),itab(jj,ii,i,k,10),itab(jj,ii,i,k,11),               &
                  itab(jj,ii,i,k,12)
           enddo
        enddo
! read in table for ice-rain collection
       do i = 1,isize
          do k = 1,jsize
             do j = 1,rcollsize
                read(10,*) dum,dum,dum,dum,dum,itabcoll(jj,ii,i,k,j,1),                  &
                 itabcoll(jj,ii,i,k,j,2),dum
                 itabcoll(jj,ii,i,k,j,1) = dlog10(itabcoll(jj,ii,i,k,j,1))
                 itabcoll(jj,ii,i,k,j,2) = dlog10(itabcoll(jj,ii,i,k,j,2))
             enddo
          enddo
       enddo
    enddo
 enddo

 close(unit=10)

! read in ice-ice collision lookup table

 open(unit=10,file=lookup_file_2,status='old')

 do i = 1,iisize
    do jjj = 1,rimsize
       do jjjj = 1,densize
          do ii = 1,iisize
             do jj = 1,jjsize
                do jjj2 = 1,rimsize
                   do jjjj2 = 1,densize
                      read(10,*) dum,dum,dum,dum,dum,dum,dum,                            &
                      itabcolli1(i,jjj,jjjj,ii,jj,jjj2,jjjj2),                           &
                      itabcolli2(i,jjj,jjjj,ii,jj,jjj2,jjjj2)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
 enddo

 close(unit=10)

!------------------------------------------------------------------------------------------!

! Generate lookup table for rain shape parameter mur
! this is very fast so it can be generated at the start of each run
! make a 150x1 1D lookup table, this is done in parameter
! space of a scaled mean size proportional qr/Nr -- initlamr

 do i = 1,150              ! loop over lookup table values
    initlamr = 1./((real(i)*2.)*1.e-6 + 250.e-6)

! iterate to get mur
! mur-lambda relationship is from Cao et al. (2008), eq. (7)

! start with first guess, mur = 0

    mur = 0.

    do ii = 1,50
       lamr = initlamr*((mur+3.)*(mur+2.)*(mur+1.)*sxth)**thrd

! new estimate for mur based on lambda
! set max lambda in formula for mur to 20 mm-1, so Cao et al.
! formula is not extrapolated beyond Cao et al. data range
       dum = min(20.,lamr*1.e-3)
       mur = max(0.,-0.0201*dum**2+0.902*dum-1.718)

! if lambda is converged within 0.1%, then exit loop
       if (ii.ge.2) then
          if (abs((lamold-lamr)/lamr).lt.0.001) goto 111
       end if

       lamold=lamr

    enddo

111 continue

! assign lookup table values
    mur_table(i) = mur

 enddo

!------------------------------------------------------------------------------------------!
! Generate lookup table for rain fallspeed and ventilation parameters
! the lookup table is two dimensional as a function of number-weighted mean size
! proportional to qr/Nr and shape parameter mur

! loop over
 mur_loop: do ii = 1,9   !range of mur is 0-8

    mur = real(ii-1)  ! values of mu

! loop over number-weighted mean size
    meansize_loop: do jj = 1,300

       if (jj.le.20) then
          dm = (real(jj)*10.-5.)*1.e-6      ! mean size (meter)
       elseif (jj.gt.20) then
          dm = (real(jj-20)*30.+195.)*1.e-6 ! mean size (meter)
       endif

! calculate PSD parameters from dm and mur

       lamr = (mur+1.)/dm

! do numerical integration over PSD

       dum1 = 0. ! numerator, number-weighted fallspeed
       dum2 = 0. ! denominator, number-weighted fallspeed
       dum3 = 0. ! numerator, mass-weighted fallspeed
       dum4 = 0. ! denominator, mass-weighted fallspeed
       dum5 = 0. ! term for ventilation factor in evap
       dd   = 2.

! loop over PSD to numerically integrate number and mass-weighted mean fallspeeds
       psd_loop: do kk = 1,10000

          dia = (real(kk)*dd-dd*0.5)*1.e-6 ! size bin (meter)
          amg = piov6rho*dia**3            ! mass (kg)
          amg = amg*1000.                  ! convert kg to g

! get fallspeed as a function of size, in cm/s
          if (dia*1.e6.le.134.43) then
             vt = 4.5795e5*amg**(2.*thrd)
             goto 101
          endif
          if (dia*1.e6.lt.1511.64) then
             vt = 4.962e3*amg**thrd
             goto 101
          endif
          if (dia*1.e6.lt.3477.84) then
             vt = 1.732e3*amg**sxth
             goto 101
          endif
          vt = 917.
 101      continue

          vt   = vt*1.e-2 ! convert from cm to m
          dum1 = dum1+vt*10.**(mur*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6
          dum2 = dum2+10.**(mur*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6
          dum3 = dum3+vt*10.**((mur+3.)*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6
          dum4 = dum4+10.**((mur+3.)*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6
          dum5 = dum5+(vt*dia)**0.5*10.**((mur+1.)*alog10(dia)+3.*mur)*exp(-lamr*dia)*   &
                 dd*1.e-6

       enddo psd_loop

       vn_table(jj,ii)    = dum1/dum2
       vm_table(jj,ii)    = dum3/dum4
       dum                = alog10(dum5)+(mur+1.)*alog10(lamr)-(3.*mur)
       revap_table(jj,ii) = 10.**(dum)

    enddo meansize_loop

 enddo mur_loop

! end rain lookup table generation
!------------------------------------------------------------------------------------------!

 print*, 'P3_INIT DONE.'
 print*

 END SUBROUTINE p3_init


!==========================================================================================!
SUBROUTINE mp_p3_wrapper_wrf
END SUBROUTINE mp_p3_wrapper_wrf
! ! !   SUBROUTINE mp_p3_wrapper_wrf(th3d,qv3d,qc3d,qr3d,qi13d,qi23d,qg13d,qg23d,qnr3d,qni13d,  &
! ! !                           qni23d,qvolg13d,qvolg23d,pii,p,dz,w,dt_in,itimestep,rainnc,     &
! ! !                           rainncv,sr,snownc,snowncv,n_iceCat,                             &
! ! !                           ids, ide, jds, jde, kds, kde ,                                  &
! ! !                           ims, ime, jms, jme, kms, kme ,                                  &
! ! !                           its, ite, jts, jte, kts, kte ,                                  &
! ! !                           zdbz3d,diag_effc3d,diag_effi3d,diag_vmi3d,di3d,diag_rhopo3d)!,    &
! ! ! !                          diag_vmi3d1,di3d1,diag_rhopo3d1,diag_vmi3d2,di3d2,diag_rhopo3d2)
! ! ! ! ! !
! ! !  !------------------------------------------------------------------------------------------!
! ! !  ! This subroutine is the main WRF interface with the P3 microphysics scheme.  It takes     !
! ! !  ! 3D variables form the driving model and passes 2D slabs (i,k) to the main microphysics   !
! ! !  ! subroutine ('P3_MAIN') over a j-loop.  For each slab, 'P3_MAIN' updates the prognostic   !
! ! !  ! variables (hydrometeor variables, potential temperature, and water vapor).  The wrapper  !
! ! !  ! also updates the accumulated precipitation arrays and then passes back them, the         !
! ! !  ! updated 3D fields, and some diagnostic fields to the driver model.                       !
! ! !  !------------------------------------------------------------------------------------------!
! ! ! ! ! !
! ! !  ! inputs:
! ! ! ! ! !
! ! !  ! th3d --> theta (K)
! ! !  ! qv3d --> vapor mass mixing ratio (kg/kg)
! ! !  ! qc3d --> cloud water mass mixing ratio (kg/kg)
! ! !  ! qi3d --> total ice mixing ratio (kg/kg)
! ! !  ! qr3d --> rain mass mixing ratio (kg/kg)
! ! !  ! qg3d --> rime ice mass mixing ratio (kg/kg)
! ! !  ! qnr3d --> rain number mixing ratio (1/kg)
! ! !  ! qni3d --> ice number mixing ratio (1/kg)
! ! !  ! qvolg3d --> ice rime volume mixing ratio (m^-3 kg^-1)
! ! !  ! pii --> Exner function (nondimensional pressure) (currently not used!)
! ! !  ! p --> pressure (pa)
! ! !  ! dz --> height difference across vertical levels (m)
! ! !  ! w --> vertical air velocity (m/s)
! ! !  ! dt_in --> time step (s)
! ! !  ! itimestep --> integer time step counter
! ! !  ! n_iceCat --> number of ice-phase categories
! ! ! ! ! !
! ! !  ! outputs:
! ! ! ! ! !
! ! !  ! rainnc --> accumulated surface precip (mm)
! ! !  ! rainncv --> one time step accumulated surface precip (mm)
! ! !  ! sr --> ice to liquid surface precip ratio
! ! !  ! snownc --> accumulated surface ice precip (mm)
! ! !  ! snowncv --> one time step accumulated surface ice precip (mm)
! ! !  ! ids...kte --> integer domain/tile bounds
! ! !  ! zdbz3d --> reflectivity (dBZ)
! ! !  ! diag_effc3d --> cloud droplet effective radius (m)
! ! !  ! diag_effi3d --> ice effective radius (m)
! ! !  ! diag_vmi3d --> mean mass weighted ice fallspeed (m/s)
! ! !  ! di3d --> mean mass weighted ice size (m)
! ! !  ! diag_rhopo3d --> mean mass weighted ice density (kg/m3)
! ! ! ! ! !
! ! !   implicit none
! ! ! ! ! !
! ! !   integer, intent(in)            ::   ids, ide, jds, jde, kds, kde ,                      &
! ! !                                       ims, ime, jms, jme, kms, kme ,                      &
! ! !                                       its, ite, jts, jte, kts, kte
! ! !   real, dimension(ims:ime, kms:kme, jms:jme), intent(inout):: th3d,qv3d,qc3d,qr3d,qi13d,   &
! ! !         qg13d,qnr3d,qni13d,qvolg13d,zdbz3d,diag_effc3d,diag_effi3d,diag_vmi3d,di3d,diag_rhopo3d
! ! !   real, dimension(ims:ime, kms:kme, jms:jme), intent(inout), optional :: qi23d,qg23d,qni23d,qvolg23d
! ! !   real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: pii,p,dz,w
! ! !   real, dimension(ims:ime, jms:jme), intent(inout) :: RAINNC,RAINNCV,SR,SNOWNC,SNOWNCV
! ! !   real, intent(in)    :: dt_in
! ! !   integer, intent(in) :: itimestep
! ! !   integer, intent(in) :: n_iceCat
! ! ! ! ! !
! ! !  ! local variables
! ! ! ! ! !
! ! !   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: &
! ! !                 qitot,qirim,nitot,birim,diag_di,diag_vmi,diag_rhopo,diag_effi
! ! !   real, dimension(ims:ime, jms:jme) :: zdbzmax
! ! !   real, dimension(ims:ime, kms:kme, jms:jme) :: diag_vmi3d2,di3d2,diag_rhopo3d2,diag_vmi3d1,di3d1,diag_rhopo3d1
! ! !
! ! ! ! ! !
! ! !   real, dimension(its:ite) :: pcprt_liq,pcprt_sol
! ! !   real                     :: dt
! ! !   integer                  :: I,K,J
! ! !  ! integer, parameter       :: n_iceCat = 2
! ! !   logical, parameter       :: nk_bottom = .false.  !.F. --> nk at model top (as in WRF)
! ! ! ! ! !
! ! !  !----- local variables and parameters:
! ! !   real :: dum1,dum2
! ! ! ! ! !
! ! !   ALLOCATE(qitot(ims:ime, kms:kme,n_iceCat))      ! ice    mixing ratio, mass (total)  kg kg-1
! ! !   ALLOCATE(qirim(ims:ime, kms:kme,n_iceCat))      ! ice    mixing ratio, mass (rime)   kg kg-1
! ! !   ALLOCATE(nitot(ims:ime, kms:kme,n_iceCat))      ! ice    mixing ratio, ratio, number # kg-1
! ! !   ALLOCATE(birim(ims:ime, kms:kme,n_iceCat))      ! ice    mixing ratio, volume m3 kg-1
! ! !   ALLOCATE(diag_di(ims:ime, kms:kme,n_iceCat))    ! mean-mass diameter, ice m
! ! !   ALLOCATE(diag_vmi(ims:ime, kms:kme,n_iceCat))   ! fall speed (mass-weighted), ices    m s-1
! ! !   ALLOCATE(diag_rhopo(ims:ime, kms:kme,n_iceCat)) ! bulk density, ice kg m-3
! ! !   ALLOCATE(diag_effi(ims:ime, kms:kme,n_iceCat))  ! effective radius of ice, for each category
! ! ! ! ! !
! ! !   qitot      = 0.
! ! !   qirim      = 0.
! ! !   nitot      = 0.
! ! !   birim      = 0.
! ! !   diag_di    = 0.
! ! !   diag_vmi   = 0.
! ! !   diag_rhopo = 0.
! ! !   diag_effi  = 0.
! ! ! ! ! !
! ! !  !------------------------------------------------------------------------------------------!
! ! ! ! ! !
! ! !   DT = DT_IN
! ! ! ! ! !
! ! !   do j = jts,jte      ! j loop (north-south)
! ! ! ! ! !
! ! !    !contruct full ice arrays from individual category arrays:
! ! !     qitot(:,:,1) = qi13d(:,:,j)
! ! !     qirim(:,:,1) = qg13d(:,:,j)
! ! !     nitot(:,:,1) = qni13d(:,:,j)
! ! !     birim(:,:,1) = qvolg13d(:,:,j)
! ! ! ! ! !
! ! !     if (n_iceCat .ge. 2) then
! ! !        qitot(:,:,2) = qi23d(:,:,j)
! ! !        qirim(:,:,2) = qg23d(:,:,j)
! ! !        nitot(:,:,2) = qni23d(:,:,j)
! ! !        birim(:,:,2) = qvolg23d(:,:,j)
! ! ! ! ! !
! ! !  !      if (n_iceCat .ge. 3) then
! ! !  !         qitot(:,:,3) = qi33d(:,:,j)
! ! !  !         qirim(:,:,3) = qg33d(:,:,j)
! ! !  !         nitot(:,:,3) = qni33d(:,:,j)
! ! !  !         birim(:,:,3) = qvolg33d(:,:,j)
! ! !  !
! ! !  !         if (n_iceCat == 4) then
! ! !  !            qitot(:,:,4) = qi43d(:,:,j)
! ! !  !            qirim(:,:,4) = qg43d(:,:,j)
! ! !  !            nitot(:,:,4) = qni43d(:,:,j)
! ! !  !            birim(:,:,4) = qvolg43d(:,:,j)
! ! !  !         endif
! ! !  !      endif
! ! !     endif
! ! ! ! ! !
! ! !  !NOTE:  unnesessary passing grid indices (jms, etc.) to P3_MAIN has been removed.  JM
! ! ! ! ! !
! ! !      call P3_MAIN(qc3d(its:ite,kts:kte,j),qr3d(its:ite,kts:kte,j),qnr3d(its:ite,kts:kte,j),                         &
! ! !              th3d(its:ite,kts:kte,j),qv3d(its:ite,kts:kte,j),dt,qitot(its:ite,kts:kte,1:n_iceCat),                  &
! ! !              qirim(its:ite,kts:kte,1:n_iceCat),nitot(its:ite,kts:kte,1:n_iceCat),birim(its:ite,kts:kte,1:n_iceCat), &
! ! !              W(its:ite,kts:kte,j),P(its:ite,kts:kte,j),DZ(its:ite,kts:kte,j),itimestep,                             &
! ! !              pcprt_liq(its:ite),pcprt_sol(its:ite),its,ite,kts,kte,nk_bottom,n_iceCat,                              &
! ! !              zdbz3d(its:ite,kts:kte,j),zdbzmax,diag_effc3d(its:ite,kts:kte,j),                                      &
! ! !              diag_effi(its:ite,kts:kte,1:n_iceCat),diag_vmi(its:ite,kts:kte,1:n_iceCat),                            &
! ! !              diag_di(its:ite,kts:kte,1:n_iceCat),diag_rhopo(its:ite,kts:kte,1:n_iceCat))
! ! !
! ! ! !=====================
! ! !
! ! ! ! ! !
! ! !  ! surface precipitation output
! ! !      pcprt_liq(:) = pcprt_liq(:)*dt*1000.  !convert to mm per time step (from m s-1)
! ! !      pcprt_sol(:) = pcprt_sol(:)*dt*1000.  !convert to mm per time step (from m s-1)
! ! !      RAINNC(its:ite,j)  = RAINNC(its:ite,j) + pcprt_liq(:) + pcprt_sol(:)   !note: this is the TOTAL precipitation
! ! !      RAINNCV(its:ite,j) = pcprt_liq(:)
! ! !      SR(its:ite,j)      = pcprt_sol(:)/(pcprt_liq(:)+1.E-12)
! ! !      SNOWNC(its:ite,j)  = SNOWNC(its:ite,j) + pcprt_sol(:)
! ! !      SNOWNCV(its:ite,j) = pcprt_sol(:)
! ! ! ! ! !
! ! !    !decompose full ice arrays back into individual category arrays:
! ! !     qi13d(:,:,j) = qitot(:,:,1)
! ! !     qg13d(:,:,j) = qirim(:,:,1)
! ! !     qni13d(:,:,j) = nitot(:,:,1)
! ! !     qvolg13d(:,:,j) = birim(:,:,1)
! ! ! ! ! !
! ! !     if (n_iceCat .eq. 1) then
! ! !     diag_vmi3d(:,:,j) = diag_vmi(:,:,1)
! ! !     di3d(:,:,j) = diag_di(:,:,1)
! ! !     diag_rhopo3d(:,:,j) = diag_rhopo(:,:,1)
! ! !     diag_effi3d(:,:,j) = diag_effi(:,:,1)
! ! !     end if
! ! ! ! ! !
! ! !     !............................................................................
! ! !     if (n_iceCat .ge. 2) then
! ! !        qi23d(:,:,j) = qitot(:,:,2)
! ! !        qg23d(:,:,j) = qirim(:,:,2)
! ! !        qni23d(:,:,j) = nitot(:,:,2)
! ! !        qvolg23d(:,:,j) = birim(:,:,2)
! ! !
! ! ! ! ! !
! ! !     do i=its,ite
! ! !        do k=kts,kte
! ! ! ! ! !
! ! !     ! for output fallspeed, size, and density, use mass-weighting of categories
! ! !     ! ****NOTE: this is only for two categories, needs to be modified for more than 2
! ! !           if ((qitot(i,k,1)+qitot(i,k,2)).ge.qsmall) then
! ! !           diag_vmi3d(i,k,j) = (diag_vmi(i,k,1)*qitot(i,k,1)+diag_vmi(i,k,2)*qitot(i,k,2))/(qitot(i,k,1)+qitot(i,k,2))
! ! !           di3d(i,k,j) = (diag_di(i,k,1)*qitot(i,k,1)+diag_di(i,k,2)*qitot(i,k,2))/(qitot(i,k,1)+qitot(i,k,2))
! ! !           diag_rhopo3d(i,k,j) = (diag_rhopo(i,k,1)*qitot(i,k,1)+diag_rhopo(i,k,2)*qitot(i,k,2))/(qitot(i,k,1)+qitot(i,k,2))
! ! !           else  ! set to default values of 0 if ice is not present
! ! !           diag_vmi3d(i,k,j) = 0.
! ! !           di3d(i,k,j) = 0.
! ! !           diag_rhopo3d(i,k,j) = 0.
! ! !           end if
! ! ! ! ! !
! ! !     ! for the combined effective radius, we need to approriately weight by mass and projected area
! ! !           if (qitot(i,k,1).ge.qsmall) then
! ! !           dum1=qitot(i,k,1)/diag_effi(i,k,1)
! ! !           else
! ! !           dum1=0.
! ! !           end if
! ! !           if (qitot(i,k,2).ge.qsmall) then
! ! !           dum2=qitot(i,k,2)/diag_effi(i,k,2)
! ! !           else
! ! !           dum2=0.
! ! !           end if
! ! !           diag_effi3d(i,k,j)=25.e-6  ! set to default 25 microns
! ! !           if (qitot(i,k,1).ge.qsmall.or.qitot(i,k,2).ge.qsmall) then
! ! !           diag_effi3d(i,k,j)=(qitot(i,k,1)+qitot(i,k,2))/(dum1+dum2)
! ! !           end if
! ! !        end do
! ! !     end do
! ! !
! ! ! ! hm 5/4/15, add separate diagnostic output for each category
! ! !           diag_vmi3d1(:,:,j) = diag_vmi(:,:,1)
! ! !           di3d1(:,:,j) = diag_di(:,:,1)
! ! !           diag_rhopo3d1(:,:,j) = diag_rhopo(:,:,1)
! ! !           diag_vmi3d2(:,:,j) = diag_vmi(:,:,2)
! ! !           di3d2(:,:,j) = diag_di(:,:,2)
! ! !           diag_rhopo3d2(:,:,j) = diag_rhopo(:,:,2)
! ! !
! ! ! ! ! !
! ! !  !      if (n_iceCat .ge. 3) then
! ! !  !         qi33d(:,:,j) = qitot(:,:,3)
! ! !  !         qg33d(:,:,j) = qirim(:,:,3)
! ! !  !         qni33d(:,:,j) = nitot(:,:,3)
! ! !  !         qvolg33d(:,:,j) = birim(:,:,3)
! ! !  !
! ! !  !         if (n_iceCat == 4) then
! ! !  !            qi43d(:,:,j) = qitot(:,:,4)
! ! !  !            qg43d(:,:,j) = qirim(:,:,4)
! ! !  !            qni43d(:,:,j) = nitot(:,:,4)
! ! !  !            qvolg43d(:,:,j) = birim(:,:,4)
! ! !  !         endif
! ! !  !      endif
! ! !     endif
! ! ! ! ! !
! ! !   enddo ! j loop
! ! !
! ! !   DEALLOCATE(qitot,qirim,nitot,birim,diag_di,diag_vmi,diag_rhopo,diag_effi)
! ! !
! ! !   END SUBROUTINE mp_p3_wrapper_wrf


!==========================================================================================!
!note: to compile with WRF or kin_1d, comment full subroutine and uncomment the following:
! ! ! SUBROUTINE mp_p3_wrapper_gem
! ! ! END SUBROUTINE mp_p3_wrapper_gem

 SUBROUTINE mp_p3_wrapper_gem(qvap,temp,dt,ww,psfc,sigma,kount,trnch,ni,nk,prt_liq,      &
                              prt_sol,diag_Zet,diag_Zec,diag_effc,diag_effi,diag_vmi,    &
                              diag_di,diag_rhopo,qc,qr,nr,n_iceCat,                      &
                              qitot_1,qirim_1,nitot_1,birim_1,                           &
                              qitot_2,qirim_2,nitot_2,birim_2,                           &
                              qitot_3,qirim_3,nitot_3,birim_3,                           &
                              qitot_4,qirim_4,nitot_4,birim_4)

!------------------------------------------------------------------------------------------!
! This wrapper subroutine is the main GEM interface with the P3 microphysics scheme.  It   !
! prepares some necessary fields (converts temperature to potential temperature, etc.),    !
! passes 2D slabs (i,k) to the main microphysics subroutine ('P3_MAIN') -- which updates   !
! the prognostic variables (hydrometeor variables, temperature, and water vapor) and       !
! computes various diagnostics fields (precipitation rates, reflectivity, etc.) -- and     !
! finally converts the updated potential temperature to temperature.                       !
!------------------------------------------------------------------------------------------!

 implicit none

!----- input/ouput arguments:  ----------------------------------------------------------!

 integer, intent(in)                    :: ni                    ! number of columns in slab
 integer, intent(in)                    :: nk                    ! number of vertical levels
 integer, intent(in)                    :: n_iceCat              ! number of ice categories
 integer, intent(in)                    :: kount                 ! time step counter
 integer, intent(in)                    :: trnch                 ! number of slice
 real, intent(in)                       :: dt                    ! model time step                     s
 real, intent(inout), dimension(ni,nk)  :: qc                    ! cloud mixing ratio, mass            kg kg-1
 real, intent(inout), dimension(ni,nk)  :: qr                    ! rain  mixing ratio, mass            kg kg-1
 real, intent(inout), dimension(ni,nk)  :: nr                    ! rain  mixing ratio, number          #  kg-1
 real, intent(inout), dimension(ni,nk)  :: qitot_1               ! ice   mixing ratio, mass (total)    kg kg-1
 real, intent(inout), dimension(ni,nk)  :: qirim_1               ! ice   mixing ratio, mass (rime)     kg kg-1
 real, intent(inout), dimension(ni,nk)  :: nitot_1               ! ice   mixing ratio, number          #  kg-1
 real, intent(inout), dimension(ni,nk)  :: birim_1               ! ice   mixing ratio, volume          m3 kg-1

 real, intent(inout), dimension(ni,nk), optional  :: qitot_2     ! ice   mixing ratio, mass (total)    kg kg-1
 real, intent(inout), dimension(ni,nk), optional  :: qirim_2     ! ice   mixing ratio, mass (rime)     kg kg-1
 real, intent(inout), dimension(ni,nk), optional  :: nitot_2     ! ice   mixing ratio, number          #  kg-1
 real, intent(inout), dimension(ni,nk), optional  :: birim_2     ! ice   mixing ratio, volume          m3 kg-1

 real, intent(inout), dimension(ni,nk), optional  :: qitot_3     ! ice   mixing ratio, mass (total)    kg kg-1
 real, intent(inout), dimension(ni,nk), optional  :: qirim_3     ! ice   mixing ratio, mass (rime)     kg kg-1
 real, intent(inout), dimension(ni,nk), optional  :: nitot_3     ! ice   mixing ratio, number          #  kg-1
 real, intent(inout), dimension(ni,nk), optional  :: birim_3     ! ice   mixing ratio, volume          m3 kg-1

 real, intent(inout), dimension(ni,nk), optional  :: qitot_4     ! ice   mixing ratio, mass (total)    kg kg-1
 real, intent(inout), dimension(ni,nk), optional  :: qirim_4     ! ice   mixing ratio, mass (rime)     kg kg-1
 real, intent(inout), dimension(ni,nk), optional  :: nitot_4     ! ice   mixing ratio, number          #  kg-1
 real, intent(inout), dimension(ni,nk), optional  :: birim_4     ! ice   mixing ratio, volume          m3 kg-1

 real, intent(inout), dimension(ni,nk)  :: qvap                  ! vapor  mixing ratio, mass           kg kg-1
 real, intent(inout), dimension(ni,nk)  :: temp                  ! temperature                         K
 real, intent(in),    dimension(ni)     :: psfc                  ! surface air pressure                Pa
 real, intent(in),    dimension(ni,nk)  :: sigma                 ! sigma = p(k,:)/psfc(:)
 real, intent(in),    dimension(ni,nk)  :: ww                    ! vertical motion                     m s-1
 real, intent(out),   dimension(ni)     :: prt_liq               ! precipitation rate, liquid          m s-1
 real, intent(out),   dimension(ni)     :: prt_sol               ! precipitation rate, solid           m s-1
 real, intent(out),   dimension(ni,nk)  :: diag_Zet              ! equivalent reflectivity, 3D         dBZ
 real, intent(out),   dimension(ni)     :: diag_Zec              ! equivalent reflectivity, col-max    dBZ
 real, intent(out),   dimension(ni,nk)  :: diag_effc             ! effective radius, cloud             m
 real, intent(out),   dimension(ni,nk)  :: diag_effi             ! effective radius, ice               m
 real, intent(out),   dimension(ni,nk)  :: diag_di               ! mean-mass diameter, ice             m
 real, intent(out),   dimension(ni,nk)  :: diag_vmi              ! fall speed (mass-weighted), ice     m s-1
 real, intent(out),   dimension(ni,nk)  :: diag_rhopo            ! bulk density, ice                   kg m-3

!----------------------------------------------------------------------------------------!

!----- local variables and parameters:
 real, dimension(ni,nk,n_iceCat)  :: qitot      ! ice    mixing ratio, mass (total)       kg kg-1
 real, dimension(ni,nk,n_iceCat)  :: qirim      ! ice    mixing ratio, mass (rime)        kg kg-1
 real, dimension(ni,nk,n_iceCat)  :: nitot      ! ice    mixing ratio, number             #  kg-1
 real, dimension(ni,nk,n_iceCat)  :: birim      ! ice    mixing ratio, volume             m3 kg-1

 real, dimension(ni,nk)  :: theta               ! potential temperature                   K
 real, dimension(ni,nk)  :: pres                ! pressure                                Pa
 real, dimension(ni,nk)  :: rho_air             ! air density                             kg m-3
 real, dimension(ni,nk)  :: DP                  ! difference in pressure between levels   Pa
 real, dimension(ni,nk)  :: DZ                  ! difference in height between levels     m
 logical, parameter      :: nk_BOTTOM = .true.  ! .true. for nk at bottom (GEM)
 integer                 :: k,ktop,kbot,kdir,i_strt,k_strt
!----------------------------------------------------------------------------------------!

   include "thermoconsts.inc"

   i_strt = 1  ! beginning index of slab
   k_strt = 1  ! beginning index of column

  !for nk_bottom = .true. :
   ktop  = 1   ! k index of top level
   kbot  = nk  ! k index of bottom level
   kdir  = -1  ! direction of vertical leveling for 1=bottom, nk=top

  !air pressure:
   do k = kbot,ktop,kdir
      pres(:,k)= psfc(:)*sigma(:,k)
   enddo

  !air density:
   rho_air  = pres/(RGASD*temp)

  !pressure difference between levels:
   DP(:,kbot) = psfc(:)-pres(:,kbot)
   do k = kbot+kdir,ktop,kdir
      DP(:,k) = pres(:,k-kdir)-pres(:,k)
   enddo

  !thickness of layers for sedimentation calculation: (in height coordiates)
   DZ = DP/(rho_air*GRAV)

  !convert to potential temperature:
   theta = temp*(1.e+5/pres)**0.286

  !contruct full ice arrays from individual category arrays:
   qitot(:,:,1) = qitot_1(:,:)
   qirim(:,:,1) = qirim_1(:,:)
   nitot(:,:,1) = nitot_1(:,:)
   birim(:,:,1) = birim_1(:,:)

   if (n_iceCat >= 2) then
      qitot(:,:,2) = qitot_2(:,:)
      qirim(:,:,2) = qirim_2(:,:)
      nitot(:,:,2) = nitot_2(:,:)
      birim(:,:,2) = birim_2(:,:)

      if (n_iceCat >= 3) then
         qitot(:,:,3) = qitot_3(:,:)
         qirim(:,:,3) = qirim_3(:,:)
         nitot(:,:,3) = nitot_3(:,:)
         birim(:,:,3) = birim_3(:,:)

         if (n_iceCat == 4) then
            qitot(:,:,4) = qitot_4(:,:)
            qirim(:,:,4) = qirim_4(:,:)
            nitot(:,:,4) = nitot_4(:,:)
            birim(:,:,4) = birim_4(:,:)
         endif
      endif
   endif

   call p3_main(qc,qr,nr,theta,qvap,dt,qitot,qirim,nitot,birim,ww,pres,DZ,kount,         &
                prt_liq,prt_sol,i_strt,ni,k_strt,nk,nk_bottom,n_iceCat,diag_Zet,         &
                diag_Zec,diag_effc,diag_effi,diag_vmi,diag_di,diag_rhopo)

!-- test: (convert rates to m/s -- still need to find where this "missing" factor is)
 prt_liq = prt_liq*1000.
 prt_sol = prt_sol*1000.
!==

  !decompose full ice arrays back into individual category arrays:
   qitot_1(:,:) = qitot(:,:,1)
   qirim_1(:,:) = qirim(:,:,1)
   nitot_1(:,:) = nitot(:,:,1)
   birim_1(:,:) = birim(:,:,1)

   if (n_iceCat >= 2) then
      qitot_2(:,:) = qitot(:,:,2)
      qirim_2(:,:) = qirim(:,:,2)
      nitot_2(:,:) = nitot(:,:,2)
      birim_2(:,:) = birim(:,:,2)

      if (n_iceCat >= 3) then
         qitot_3(:,:) = qitot(:,:,3)
         qirim_3(:,:) = qirim(:,:,3)
         nitot_3(:,:) = nitot(:,:,3)
         birim_3(:,:) = birim(:,:,3)

         if (n_iceCat == 4) then
            qitot_4(:,:) = qitot(:,:,4)
            qirim_4(:,:) = qirim(:,:,4)
            nitot_4(:,:) = nitot(:,:,4)
            birim_4(:,:) = birim(:,:,4)
         endif
      endif
   endif


  !convert back to temperature:
   temp = theta*(pres*1.e-5)**0.286

 END SUBROUTINE mp_p3_wrapper_gem

!==================================================================================================!
 SUBROUTINE p3_main(qc,qr,nr,th,qv,dt,qitot,qirim,nitot,birim,uzpl,pres,dzq,it,pcprt_liq, &
                    pcprt_sol,its,ite,kts,kte,nk_bottom,n_iceCat,diag_ze,diag_zec,        &
                    diag_effc,diag_effi,diag_vmi,diag_di,diag_rhopo)

!----------------------------------------------------------------------------------------!
!                                                                                        !
! This is the main subroutine for the P3 microphysics scheme.  It is called from the     !
! wrapper subroutine ('MP_P3_WRAPPER') and is passed i,k slabs of all prognostic         !
! variables -- hydrometeor fields, potential temperature, and water vapor mixing ratio.  !
! Microphysical process rates are computed first.  These tendencies are then used to     !
! computed updated values of the prognostic variables.  The hydrometeor variables are    !
! then updated further due to sedimentation.                                             !
!                                                                                        !
! Several diagnostic values are also computed and returned to the wrapper subroutine,    !
! including precipitation rates.                                                         !
!                                                                                        !
!----------------------------------------------------------------------------------------!

 implicit none

!----- arguments:  ----------------------------------------------------------------------!

 integer, intent(in)                                      :: its        ! grid (tile) i-index (start)
 integer, intent(in)                                      :: ite        ! grid (tile) i-index (end)
 integer, intent(in)                                      :: kts        ! grid (tile) k-index (start)
 integer, intent(in)                                      :: kte        ! grid (tile) k-index (end)
 integer, intent(in)                                      :: n_iceCat   ! number of ice-phase categories
 real,    intent(in)                                      :: dt         ! model time step (s)
 logical, intent(in)                                      :: nk_bottom  ! .F. for nk at model top,  .T. for nk at model bottom

 real, intent(inout), dimension(its:ite,kts:kte)          :: qc         ! cloud, mass mixing ratio          kg kg-1
!real, intent(inout), dimension(its:ite,kts:kte)          :: nc         ! cloud, number mixing ratio        #  kg-1 (note: specified in this version)
 real, intent(inout), dimension(its:ite,kts:kte)          :: qr         ! rain,  mass mixing ratio          kg kg-1
 real, intent(inout), dimension(its:ite,kts:kte)          :: nr         ! rain, number mixing ratio         #  kg-1
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: qitot      ! ice, total mass mixing ratio      kg kg-1
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: qirim      ! ice, rime mass mixing ratio       kg kg-1
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: nitot      ! ice, total number mixing ratio    #  kg-1
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: birim      ! ice, rime volume mixing ratio     m3 kg-1

 real, intent(inout), dimension(its:ite,kts:kte)          :: qv         ! water vapor mixing ratio          kg kg-1
 real, intent(inout), dimension(its:ite,kts:kte)          :: th         ! potential temperature             K
 real, intent(in),    dimension(its:ite,kts:kte)          :: uzpl       ! vertical air velocity             m s-1
 real, intent(in),    dimension(its:ite,kts:kte)          :: pres       ! pressure                          Pa
 real, intent(in),    dimension(its:ite,kts:kte)          :: dzq        ! vertical grid spacing             m
 real, intent(out),   dimension(its:ite)                  :: pcprt_liq  ! precipitation rate, liquid        m s-1
 real, intent(out),   dimension(its:ite)                  :: pcprt_sol  ! precipitation rate, solid         m s-1

 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: diag_vmi   ! mass-weighted fall speed, ice     m s-1
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: diag_di    ! mean-mass diameter, ice           m
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: diag_rhopo ! bulk density, ice                 kg m-3
 real, intent(inout), dimension(its:ite,kts:kte,n_iceCat) :: diag_effi  ! effective radius, ice             m
 real, intent(out),   dimension(its:ite,kts:kte)          :: diag_effc  ! effective radius, cloud           m
 real, intent(out),   dimension(its:ite,kts:kte)          :: diag_ze    ! equivalent reflectivity, 3D       dBZ
 real, intent(out),   dimension(its:ite)                  :: diag_zec   ! equivalent reflectivity, col-max  dBZ


!----- local variables and parameters: --------------------------------------------------!

 real, dimension(n_iceCat)                 :: epsi
 real, dimension(its:ite,kts:kte,n_iceCat) :: diam_ice

!** declare local
 real, dimension(its:ite,kts:kte) :: nc   ! cloud, number mixing ratio  # kg-1 (note: specified in this version)
!**

 real, dimension(its:ite,kts:kte) :: mur
 real, dimension(its:ite,kts:kte) :: ssat
 real, dimension(its:ite,kts:kte) :: t ! temperature (K)

! time-varying parameters
 real, dimension(its:ite,kts:kte) :: rho,rhofacr,rhofaci,acn,xxls,xxlv,xlf,qvs,qvi,sup,  &
                                     supi,ss
 real                             :: mu,dv,sc,dqsdt,ab,kap,epsr,epsc,epsi_tot

! variables for cond/evap/sub/dep
 real                             :: xx,aaa,epsilon

! local variables - aerosol/droplet activation
 real                             :: sigvl,aact,alpha,gamm,gg,psi,eta1,eta2,sm1,sm2,     &
                                     smax,uu1,uu2

! local/dummy variables
 real                             :: dum,dum1,dum2,dumt,dumqv,dumqvs,dums,dumqc,ratio,   &
                                     qsat0,udiff

! 2D size distribution and fallspeed parameters
 real, dimension(its:ite,kts:kte) :: lamc
 real, dimension(its:ite,kts:kte) :: lamr
 real, dimension(its:ite,kts:kte) :: n0c
 real, dimension(its:ite,kts:kte) :: n0r
 real, dimension(its:ite,kts:kte) :: pgam
 real, dimension(its:ite,kts:kte) :: effr
 real, dimension(its:ite,kts:kte) :: vtrn
 real, dimension(its:ite,kts:kte) :: vtrm
 real, dimension(its:ite,kts:kte) :: nu
 real, dimension(its:ite,kts:kte) :: vtrnc
 real, dimension(its:ite,kts:kte) :: vtrmc
 real, dimension(its:ite,kts:kte) :: cdist
 real, dimension(its:ite,kts:kte) :: cdist1
 real, dimension(its:ite,kts:kte) :: cdistr

! max/min lambda for cloud droplets and rain
 real    :: lammax,lammin

! parameters for rain fallspeed lookup table
 integer :: dumi
 real    :: dum3,dum4,dum5

! integer indices for ice lookup table
 integer :: i,k,kk,ii,jj

 integer :: it                !time step counter
 integer :: iice              !index for ice category
 integer :: iice_dest         !index for destination ice category for nucleation

! parameters for mur lookup table
 real    :: lamold,rdumii,rdumjj

! local 2D ice variables
 real, dimension(its:ite,kts:kte) :: vtrmi1,vtrnitot

! microphysical process rates based on liquid only
 real :: qrcon   !
 real :: qcacc   !              ADD UNITS TO ALL PROCESS RATE TERMS
 real :: qcaut   !
 real :: ncacc   !
 real :: ncautc  !
 real :: ncslf   !
 real :: nrslf   !
 real :: ncnuc   !
 real :: qccon   !
 real :: qcnuc   !
 real :: qrevp   !
 real :: qcevp   !
 real :: nrevp   ! change in Nr, evaporation
 real :: ncautr  ! change in Nr, autoconversion of cloud water (not in Part 1)

! microphysical process rates based on ice phase
 real, dimension(n_iceCat) :: qccol    ! collection cloud water
 real, dimension(n_iceCat) :: qidep    ! vapor deposition
 real, dimension(n_iceCat) :: qrcol    ! collection rain mass by ice
 real, dimension(n_iceCat) :: qinuc    ! deposition/condensation freezing nuc
 real, dimension(n_iceCat) :: nccol    !
 real, dimension(n_iceCat) :: nrcol    !
 real, dimension(n_iceCat) :: ninuc    !
 real, dimension(n_iceCat) :: qisub    ! sublimation species
 real, dimension(n_iceCat) :: qimlt    ! melting of ice
 real, dimension(n_iceCat) :: nimlt    ! melting of ice
 real, dimension(n_iceCat) :: nisub    !
 real, dimension(n_iceCat) :: nislf    !
 real, dimension(n_iceCat) :: qchetc   ! contact freezing droplets
 real, dimension(n_iceCat) :: qcheti   ! contact freezing droplets
 real, dimension(n_iceCat) :: qrhetc   ! contact freezing rain
 real, dimension(n_iceCat) :: qrheti   ! immersion freezing rain
 real, dimension(n_iceCat) :: nchetc   ! contact freezing droplets
 real, dimension(n_iceCat) :: ncheti   ! immersion freezing droplets
 real, dimension(n_iceCat) :: nrhetc   ! contact freezing rain
 real, dimension(n_iceCat) :: nrheti   ! immersion freezing rain
 real, dimension(n_iceCat) :: nrshdr   ! source for rain number from collision of rain/ice above freezing and shedding
 real, dimension(n_iceCat) :: qcshd    ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
 real, dimension(n_iceCat) :: qcmul    ! change in q, ice multiplication from rime-splitnering of cloud water (not included in the paper)
 real, dimension(n_iceCat) :: qrmul    ! change in q, ice multiplication from rime-splitnering of rain (not included in the paper)
 real, dimension(n_iceCat) :: nimul    ! change in Ni, ice multiplication from rime-splintering (not included in the paper)
 real, dimension(n_iceCat) :: ncshdc   ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
 real, dimension(n_iceCat) :: rhorime_c ! density of rime (from cloud)
!real, dimension(n_iceCat) :: rhorime_r ! density of rime (from rain)
 real, dimension(n_iceCat,n_iceCat) :: nicol   ! change of N due to ice-ice collision between categories
 real, dimension(n_iceCat,n_iceCat) :: qicol   ! change of q due to ice-ice collision between categories
!real, dimension(n_iceCat) :: pracsw1  ! collection of snow mass by rain above freezing, not currently usednot used
!real, dimension(n_iceCat) :: npracsw1 ! collection snow number by rain above freezing, not currently used!not used

 real :: qwgrth
 logical, dimension(n_iceCat) :: wetgrowth

 integer :: j,dumk,dumj,dumii,dumjj

! interpolated quantities from ice lookup table
 real :: f1pr1,f1pr2,f1pr3,f1pr4,f1pr5,f1pr6,f1pr7,f1pr8,f1pr9,f1pr10,f1pr12,f1pr13,     &
         f1pr14,f1pr15,f1pr16,f1pr17,f1pr18

! time-varying parameters for ice microphysics
 real :: dqsidt,abi,dumqvi,dap,nacnt,rhop

! variables for rime density calculation
 real :: v_impact,ri,iTc,D_c,D_r

! real, dimension(its:ite,kts:kte) :: ze_ice,ze_rain  !,ze_cloud
 real :: ze_ice,ze_rain  !,ze_cloud
 real, parameter :: z_min    = 1.e-30
 real, parameter :: zdBZ_min = -99.

! tendencies from sedimentation
 real, dimension(kts:kte) :: qrsten            ! rain SED TEND (KG/KG/S)
 real, dimension(kts:kte) :: nrsten
 real, dimension(kts:kte) :: qcsten            ! cloud SED TEND (KG/KG/S)
 real, dimension(kts:kte) :: ncsten
 real, dimension(kts:kte) :: qisten            ! ice SED TEND (KG/KG/S)
 real, dimension(kts:kte) :: qristen
 real, dimension(kts:kte) :: bgsten
 real, dimension(kts:kte) :: nisten

! local variables for sedimentation calculations
 real, dimension(kts:kte) :: dumqi,dumr,dumfni,dumri,dumbg,dumfnr,dumc,dumfnc,FR,fi,fni, &
                             fnr,fc,fnc,faloutr,falouti,faloutni,faloutnr,faloutri,      &
                             faloutbg,faloutc,faloutnc
 real    :: faltndr,faltndi,faltndni,faltndri,faltndbg,faltndc,faltndnc,faltndnr,rgvm
 integer :: n,nstep
 integer :: ktop,kbot,kdir    ! level indices for generalized direction

! integer variables to skip code calculations
 logical :: log_nucleationPossible,log_hydrometeorsPresent

 integer :: index

! for sedimentation
 logical :: log_qcpresent,log_qrpresent,log_qipresent
 integer :: qcindex,qrindex,qiindex

! local variables
 real   :: dumlr

! varaibles for code optimization
 real, dimension(its:ite,kts:kte)   :: inv_dzq,inv_rho
 real :: inv_nstep,inv_dum,inv_dum3,odt,oxx,oabi,zero,test,test2,test3
 real, parameter :: thrd = 1./3.

 real, parameter :: rhorime_c_dflt  = 400.
!real, parameter :: rhorime_r_dflt  = 400.
 real, parameter :: rho_rimeMin     =  50.
 real, parameter :: rho_rimeMax     = 900.
 real, parameter :: inv_rho_rimeMax = 1./rho_rimeMax

 logical :: log_predictsSsat,log_exitlevel,log_hmossopOn
 real    :: D_new,Q_nuc,N_nuc

! variables for category interaction (collection)
 real, dimension(n_iceCat) :: Eii_fact  !,vmimean
 real                      :: dum1c,dum4c,dum5c,drhop, tmp1,tmp2,tmp3
 integer                   :: dumic,dumiic,dumjjc,catcoll

 real :: deltaD_init

 logical :: ni_add

! end of local variable declarations
!------------------------------------------------------------------------------------------!

! Determine threshold size difference [m] as a function of n_iceCat
! (used for destination category upon ice initiation)
 select case (n_iceCat)
    case (1)
       deltaD_init = 999.    !not used if n_iceCat=1 (but should be defined)
    case (2)
       deltaD_init = 500.e-6
    case (3)
       deltaD_init = 400.e-6
    case (4)
       deltaD_init = 235.e-6
    case (5)
       deltaD_init = 175.e-6
    case (6:)
       deltaD_init = 150.e-6
 end select

! deltaD_init = 250.e-6   !for testing
! deltaD_init = dummy_in   !temporary; passed in from cld1d


  !check if variables inialized in 'P3_INIT' are indeed initialized:
 if (rd /= 287.15) then
    print*, '*** ABORT in s/r P3_MAIN ***'
    print*, '*'
    print*, '* Values supposed to be initialized in P3_INIT are not seen here.'
    print*, '* e.g. rd: ',rd
    print*, '***'
    stop
!   call qqexit(1) !GEM abort
 endif


 !direction of vertical leveling:
 if (nk_bottom) then
   !GEM / kin_1d:
    ktop = kts        !k of top level
    kbot = kte        !k of bottom level
    kdir = -1         !(k: 1=top, nk=bottom)
 else
   !WRF / kin_2d:
    ktop = kte        !k of top level
    kbot = kts        !k of bottom level
    kdir = 1          !(k: 1=bottom, nk=top)
 endif

! set logical for prediction of supersaturation, true or false
 log_predictsSsat = .false.

! calculate inverse model time step
 odt = 1./dt

 i_loop_main: do i = its,ite  ! main i-loop around the entire scheme

    log_hydrometeorsPresent = .false.
    log_nucleationPossible  = .false.

! initialize surface precipitation rate output
    pcprt_liq(i) = 0.
    pcprt_sol(i) = 0.

    diag_zec(i)  = -99.

    k_loop_1: do k = kbot,ktop,kdir

! initialize effective radii to default values
       diag_effc(i,k)   = 10.e-6
       effr(i,k)        = 25.e-6
       diag_effi(i,k,:) = 25.e-6

! initialize output mass-weighted ice fallspeed and mean size
       diag_vmi(i,k,:) = 0.
       diag_di(i,k,:)  = 0.
       diag_rhopo(i,k,:) = 0.

! initialize mean ice size to zero
       diam_ice(i,k,:) = 0.

! initialize rain mu
       mur(i,k)   = 0.

! initialize reflectivity
       diag_ze(i,k) = zdBZ_min

! calculate temperature from theta
       t(i,k) = th(i,k) * (1.e-5*pres(i,k))**(rd*inv_cp)

! calculate some time-varying atmospheric variables
       rho(i,k)     = pres(i,k)/(rd*t(i,k))
       inv_rho(i,k) = 1./rho(i,k)
       xxlv(i,k)    = 3.1484e6-2370.*t(i,k)
       xxls(i,k)    = xxlv(i,k)+0.3337e6
       xlf(i,k)     = xxls(i,k)-xxlv(i,k)
       qvs(i,k)     = ep_2*polysvp1(t(i,k),0)/(pres(i,k)-polysvp1(t(i,k),0))
       qvi(i,k)     = ep_2*polysvp1(t(i,k),1)/(pres(i,k)-polysvp1(t(i,k),1))
       sup(i,k)     = qv(i,k)/qvs(i,k)-1.
       supi(i,k)    = qv(i,k)/qvi(i,k)-1.
       rhofacr(i,k) = (rhosur*inv_rho(i,k))**0.54
       rhofaci(i,k) = (rhosui*inv_rho(i,k))**0.54
       mu           = 1.496e-6*t(i,k)**1.5/(t(i,k)+120.)
       acn(i,k)     = g*rhow/(18.*mu)

! specify droplet concentration
       nc(i,k)      = nccnst*inv_rho(i,k)

! apply mass clipping if conditions are dry and mass mixing ratio is sufficiently small
! (and thus all mass is expected to evaporate/sublimate in one time step)

       do iice = 1,n_iceCat

          if (qitot(i,k,iice).lt.qsmall .or. (qitot(i,k,iice).lt.1.e-8 .and.             &
           supi(i,k).lt.-0.1)) then
             qv(i,k) = qv(i,k)+qitot(i,k,iice)
             th(i,k) = th(i,k)+th(i,k)/t(i,k)*qitot(i,k,iice)*xxls(i,k)*inv_cp
             qitot(i,k,iice)  = 0.
             nitot(i,k,iice)  = 0.
             qirim(i,k,iice)  = 0.
             birim(i,k,iice)  = 0.
          else
             log_hydrometeorsPresent = .true.    ! updated further down
          endif

          if (qitot(i,k,iice).ge.qsmall .and. qitot(i,k,iice).lt.1.e-8 .and.             &
           t(i,k).ge.273.15) then
             qr(i,k) = qr(i,k)+qitot(i,k,iice)
             th(i,k) = th(i,k)+th(i,k)/t(i,k)*qitot(i,k,iice)*xlf(i,k)*inv_cp
             qitot(i,k,iice)  = 0.
             nitot(i,k,iice)  = 0.
             qirim(i,k,iice)  = 0.
             birim(i,k,iice)  = 0.
          endif

       enddo !iice-loop

       if (qc(i,k).lt.qsmall .or. (qc(i,k).lt.1.e-8 .and. sup(i,k).lt.-0.1)) then
          qv(i,k) = qv(i,k)+qc(i,k)
          th(i,k) = th(i,k)+th(i,k)/t(i,k)*qc(i,k)*xxlv(i,k)*inv_cp
          qc(i,k) = 0.
          nc(i,k) = 0.
       else
          log_hydrometeorsPresent = .true.    ! updated further down
       endif

       if (qr(i,k).lt.qsmall .or. (qr(i,k).lt.1.e-8 .and. sup(i,k).lt.-0.1)) then
          qv(i,k) = qv(i,k)+qr(i,k)
          th(i,k) = th(i,k)+th(i,k)/t(i,k)*qr(i,k)*xxls(i,k)*inv_cp
          qr(i,k) = 0.
          nr(i,k) = 0.
       else
          log_hydrometeorsPresent = .true.    ! final update
       endif

! If there is the possibility of nucleation or droplet activation (i.e., if RH
! is relatively high) then calculate microphysical processes even if there
! is no existing condensate
         if ((t(i,k).lt.273.15 .and. supi(i,k).ge.-0.05) .or.                            &
             (t(i,k).ge.273.15 .and. sup(i,k).ge.-0.05 )) log_nucleationPossible = .true.

    enddo k_loop_1

! jump to end of i-loop if log_nucleationPossible=.false.  (i.e. skip everything)
    if (.not. (log_nucleationPossible .or. log_hydrometeorsPresent)) goto 333

    log_hydrometeorsPresent = .false.   ! reset value; used again below

!------------------------------------------------------------------------------------------!
! main k-loop (for processes):

    k_loop_main: do k = kbot,ktop,kdir

       rhorime_c = rhorime_c_dflt

       log_exitlevel = .true.

! if it is dry and there are no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
       if (qc(i,k).ge.qsmall .or. qr(i,k).ge.qsmall) log_exitlevel = .false.

       do iice = 1,n_iceCat
          if (qitot(i,k,iice).ge.qsmall) log_exitlevel = .false.
       enddo !iice-loop

       if (log_exitlevel) then
          if ((t(i,k).lt.273.15 .and. supi(i,k).lt.-0.05) .or.                           &
              (t(i,k).ge.273.15 .and. sup(i,k) .lt.-0.05)) goto 555   !i.e. skip all process rates
       endif

! initialize warm microphysics processes to zero
       qcacc   = 0.;     qrevp   = 0.;     qccon   = 0.
       qcaut   = 0.;     qcevp   = 0.;     qrcon   = 0.
       ncacc   = 0.;     ncnuc   = 0.;     ncslf   = 0.
       ncautc  = 0.;     qcnuc   = 0.;     nrslf   = 0.
       nrevp   = 0.;     ncautr  = 0.

! initialize ice microphysical processes to zero
       qchetc  = 0.;     qisub   = 0.;     nrshdr  = 0.
       qcheti  = 0.;     qrcol   = 0.;     qcshd   = 0.
       qrhetc  = 0.;     qimlt   = 0.;     qcmul   = 0.
       qrheti  = 0.;     qinuc   = 0.;     qrmul   = 0.
       nchetc  = 0.;     nccol   = 0.;     ncshdc  = 0.
       ncheti  = 0.;     nrcol   = 0.;     nislf   = 0.
       nrhetc  = 0.;     ninuc   = 0.;     qidep   = 0.
       nrheti  = 0.;     nisub   = 0.;     nimul   = 0.
       qccol   = 0.;     nimlt   = 0.;

! initialize ice-ice collection processes to zero
       qicol   = 0.
       nicol   = 0.

!!     pracsw1  = 0.
!!     npracsw1 = 0.

       wetgrowth = .false.  ! logical for wetgrowth
!===

!----------------------------------------------------------------------
       predict_supersaturation: if (log_predictsSsat) then

! adjust cloud water and thermodynamics to prognostic supersaturation
! following the method in Grabowski and Morrison (2008)

          dqsdt   = xxlv(i,k)*qvs(i,k)/(rv*t(i,k)*t(i,k))
          ab      = 1.+dqsdt*xxlv(i,k)*inv_cp
          epsilon = (qv(i,k)-qvs(i,k)-ssat(i,k))/ab

! limit adjustment to available water

          epsilon = max(epsilon,-qc(i,k))

! don't adjust upward if subsaturated
! otherwise this could result in positive adjustment
! (spurious generation ofcloud water) in subsaturated conditions

          if (ssat(i,k).lt.0.) epsilon = min(0.,epsilon)

! now do the adjustment
          if (abs(epsilon).ge.1.e-15) then

             qc(i,k)   = qc(i,k)+epsilon
             qv(i,k)   = qv(i,k)-epsilon
             th(i,k)   = th(i,k)+epsilon*th(i,k)/t(i,k)*xxlv(i,k)*inv_cp

! recalculate variables if there was adjustment
! calculate temperature from theta
             t(i,k)    = th(i,k)*(1.e-5*pres(i,k))**(rd*inv_cp)
             qvs(i,k)  = ep_2*polysvp1(t(i,k),0)/(pres(i,k)-polysvp1(t(i,k),0))
             qvi(i,k)  = ep_2*polysvp1(t(i,k),1)/(pres(i,k)-polysvp1(t(i,k),1))
             sup(i,k)  = qv(i,k)/qvs(i,k)-1.
             supi(i,k) = qv(i,k)/qvi(i,k)-1.
             ssat(i,k) = qv(i,k)-qvs(i,k)

          endif

       endif predict_supersaturation
!----------------------------------------------------------------------

! skip micro process calculations except nucleation/acvtivation if there no hydrometeors are present
!       if (it.gt.1 .and. qitot(i,k,iice).lt.qsmall .and. qc(i,k).lt.qsmall .and.        &
!        qr(i,k).lt.qsmall) goto 444   !i.e. skip to nucleation

       log_exitlevel = .true.
       if (qc(i,k).ge.qsmall .or. qr(i,k).ge.qsmall) log_exitlevel = .false.

       do iice = 1,n_iceCat
          if (qitot(i,k,iice).ge.qsmall) log_exitlevel = .false.
       enddo !iice-loop

       if (log_exitlevel) goto 444   !i.e. skip to nucleation


! time/space varying physical variables
       mu        = 1.496e-6*t(i,k)**1.5/(t(i,k)+120.)
       dv        = 8.794e-5*t(i,k)**1.81/pres(i,k)
       sc        = mu/(rho(i,k)*dv)
       dum       = 1./(rv*t(i,k)**2)
       dqsdt     = xxlv(i,k)*qvs(i,k)*dum
       dqsidt    = xxls(i,k)*qvi(i,k)*dum
       ab        = 1.+dqsdt*xxlv(i,k)*inv_cp
       abi       = 1.+dqsidt*xxls(i,k)*inv_cp
       kap       = 1.414e3*mu
       ssat(i,k) = qv(i,k)-qvs(i,k)

! contact freezing currently turned off
!      dum = 7.37*t(i,k)/(288.*10.*pres(i,k))/100.
!      dap = 4.*pi*1.38e-23*t(i,k)*(1.+dum/rin)/(6.*pi*rin*mu)
!      nacnt = exp(-2.80+0.262*(273.15-t(i,k)))*1000.
!      rhofacr(i,k) = (rhosur*inv_rho(i,k))**0.54
!      rhofaci(i,k) = (rhosui*inv_rho(i,k))**0.54
! 'a' parameter for droplet fallspeed calculated from Stokes' law
       acn(i,k)  = g*rhow/(18.*mu)

!..................................................
! get droplet size distribution parameters
       if (qc(i,k).ge.qsmall) then

! set minimum nc to prevent floating point error
          nc(i,k)   = max(nc(i,k),1.e-16)
          pgam(i,k) = 0.0005714*(nc(i,k)*1.e-6*rho(i,k))+0.2714
          pgam(i,k) = 1./(pgam(i,k)**2)-1.
          pgam(i,k) = max(pgam(i,k),2.)
          pgam(i,k) = min(pgam(i,k),15.)  !note: pgam is mu_c

! interpolate for mass distribution spectral shape parameter (for SB warm processes)
          if (iparam.eq.1) then
             dumi    = int(pgam(i,k))
             nu(i,k) = dnu(dumi)+(dnu(dumi+1)-dnu(dumi))*(pgam(i,k)-dumi)
          endif

! calculate lamc
          lamc(i,k) = (cons1*nc(i,k)*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.)/       &
                      qc(i,k))**thrd

! apply lambda limiters
          lammin = (pgam(i,k)+1.)*2.5e+4   ! min: 40 micron mean diameter
          lammax = (pgam(i,k)+1.)*1.e+6    ! max:  1 micron mean diameter

          if (lamc(i,k).lt.lammin) then
             lamc(i,k) = lammin
             nc(i,k)   = 6.*lamc(i,k)**3*qc(i,k)/(pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)* &
                         (pgam(i,k)+1.))
          elseif (lamc(i,k).gt.lammax) then
             lamc(i,k) = lammax
             nc(i,k)   = 6.*lamc(i,k)**3*qc(i,k)/(pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)* &
                         (pgam(i,k)+1.))
          endif

          cdist(i,k)  = nc(i,k)*(pgam(i,k)+1.)/lamc(i,k)
          cdist1(i,k) = nc(i,k)/gamma(pgam(i,k)+1.)

       else

          lamc(i,k)   = 0.
          cdist(i,k)  = 0.
          cdist1(i,k) = 0.

       endif

!..................................................
! get rain size distribution parameters
       if (qr(i,k).ge.qsmall) then

! set minimum value for Nr to prevent taking 0^(1/3) and thus code crash
          nr(i,k) = max(nr(i,k),1.e-16)

! use lookup table to get mur

! find spot in lookup table
! calculate scaled mean size proportional to N/q for lookup table parameter space
! the scaled mean size is identical to lambda for an exponential distribution
          dum  = (cons1*NR(i,k)*6./(QR(i,k)))**thrd
          inv_dum = 1./dum

          if (inv_dum.lt.282.e-6) then
             mur(i,k) = 8.282
          elseif (inv_dum.ge.282.e-6 .and. inv_dum.lt.502.e-6) then
! find integer index corresponding with scaled mean size dum
             rdumii = (inv_dum-250.e-6)*0.5e+6
             rdumii = max(rdumii,1.)
             rdumii = min(rdumii,150.)
             dumii  = int(rdumii)
             dumii  = min(149,dumii)
! interpolate lookup table values
             mur(i,k) = mur_table(dumii)+(mur_table(dumii+1)-mur_table(dumii))*(rdumii-  &
                        real(dumii))
          elseif (inv_dum.ge.502.e-6) then
             mur(i,k) = 0.
          endif

! recalculate lambda based on mur
          lamr(i,k) = (cons1*nr(i,k)*(mur(i,k)+3.)*(mur(i,k)+2)*(mur(i,k)+1.)/           &
                      qr(i,k))**thrd

! set maximum value for rain lambda
          lammax = (mur(i,k)+1.)*1.e5

! set to small value since breakup is explicitly included (mean size 0.8 mm)
          lammin = (mur(i,k)+1.)*1250.

! apply lambda limiters
          if (lamr(i,k).lt.lammin) then
             lamr(i,k) = lammin
             nr(i,k)   = exp(3.*log(lamr(i,k))+log(qr(i,k))+log(gamma(mur(i,k)+1.))-     &
                         log(gamma(mur(i,k)+4.)))/(cons1)
          elseif (lamr(i,k).gt.lammax) then
             lamr(i,k) = lammax
             nr(i,k)   = exp(3.*log(lamr(i,k))+log(QR(i,k))+log(gamma(mur(i,k)+1.))-     &
                         log(gamma(mur(i,k)+4.)))/(cons1)
          endif

          cdistr(i,k) = nr(i,k)/gamma(mur(i,k)+1.)

! NOTE: n0r is calculated as log10(n0r)
          n0r(i,k) = alog10(nr(i,k))+(mur(i,k)+1.)*alog10(lamr(i,k))-                    &
                     alog10(gamma(mur(i,k)+1))

       else

          lamr(i,k)   = 0.
          n0r(i,k)    = 0.
          cdistr(i,k) = 0.

       endif

!------------------------------------------------------------------------------------------!

!JM: Can the following not be combined with the IF block above?
       if (qr(i,k).ge.qsmall) then

! calculate rain evaporation including ventilation and number-
! weighted fallspeed (used for calculating rime density)
! find appropriate place in 2D rain lookup table

! find location in scaled mean size space
          dum1 = (mur(i,k)+1.)/lamr(i,k)
          if (dum1.le.195.e-6) then
             inv_dum3  = 0.1
             rdumii = (dum1*1.e6+5.)*inv_dum3
             rdumii = max(rdumii,1.)
             rdumii = min(rdumii,20.)
             dumii  = int(rdumii)
             dumii  = max(dumii,1)
             dumii  = min(dumii,20)
          elseif (dum1.gt.195.e-6) then
             inv_dum3  = 0.0333333           !i.e. 1/30
             rdumii = (dum1*1.e6-195.)*inv_dum3 + 20.
             rdumii = max(rdumii,20.)
             rdumii = min(rdumii,300.)
             dumii  = int(rdumii)
             dumii  = max(dumii,20)
             dumii  = min(dumii,299)
          endif

! find location in mur space
          rdumjj = mur(i,k)+1.
          rdumjj = max(rdumjj,1.)
          rdumjj = min(rdumjj,10.)
          dumjj  = int(rdumjj)
          dumjj  = max(dumjj,1)
          dumjj  = min(dumjj,9)

! interpolate value at mur
          dum1 = revap_table(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*                 &
                 (revap_table(dumii+1,dumjj)-revap_table(dumii+1,dumjj))

! interoplate value at mur+1
          dum2 = revap_table(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*               &
                 (revap_table(dumii+1,dumjj+1)-revap_table(dumii+1,dumjj+1))

! final interpolation
          dum  = dum1+(rdumjj-real(dumjj))*(dum2-dum1)
          epsr = 2.*pi*cdistr(i,k)*rho(i,k)*DV*(F1R*gamma(mur(i,k)+2.)/(lamr(i,k))+F2R*  &
                 (rho(i,k)/MU)**0.5*SC**0.33333*dum)

! number-weighted fallspeed
! value at mur
!         dum1 = vn_table(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*(vn_table(dumii+1,  &
!                dumjj)-vn_table(dumii+1,dumjj))

! value at mur+1
!         dum2 = vn_table(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*(vn_table(dumii+1,&
!                dumjj+1)-vn_table(dumii+1,dumjj+1))

! final interpolation
!        vtrn(i,k)=dum1+(rdumjj-real(dumjj))* &
!              (dum2-dum1)
!        vtrn(i,k)=vtrn(i,k)*rhofacr(i,k)

       else
          epsr = 0.
!         vtrn(i,k) = 0.
       endif


       if (qc(i,k).ge.qsmall) then
          epsc = 2.*pi*rho(i,k)*dv*cdist(i,k)
       else
          epsc = 0.
       endif

!------------------------------------------------------------------------------------------!
! initialize inverse supersaturation relaxation timescale for combined ice categories
       epsi_tot = 0.

! limit max total ice concentration for combined ice categories to 500 L-1
! if max is exceeded scale each category to preserve ratio of number between categories
       if (sum(nitot(i,k,:)).ge.1.e-20) then
          dum = 500.e+3*inv_rho(i,k)/sum(nitot(i,k,:))
          nitot(i,k,:) = nitot(i,k,:)*min(dum,1.)
       endif


       iice_loop1: do iice = 1,n_iceCat

! ice lookup table for ice micro processes
          qitot_notsmall_1: if (qitot(i,k,iice).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
             nitot(i,k,iice) = max(nitot(i,k,iice),1.e-16)

             !--compute mean-mass ice diameters (estimated; rigorous approach to be implemented later)
             dum2 = 500. !ice density
             diam_ice(i,k,iice) = ((qitot(i,k,iice)*6.)/(nitot(i,k,iice)*dum2*pi))  &
                                  **thrd
             !==

! calculate predicted bulk rime density
             if (birim(i,k,iice).ge.1.e-15) then
                rhop = qirim(i,k,iice)/birim(i,k,iice)
             else
                rhop = 0.
             endif

! limit rho_rimeMin (50) < rhop < rho_rimeMax (900), adjust bg if needed
             if (rhop.lt.rho_rimeMin) then
                rhop     = rho_rimeMin
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif

             if (rhop.gt.rho_rimeMax) then
                rhop     = rho_rimeMax
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif

             if (qirim(i,k,iice).lt.qsmall) then
                birim(i,k,iice) = 0.
             endif

! set upper constraint on qri to ensure qri cannot be > qi
             if (qirim(i,k,iice).gt.qitot(i,k,iice)) then
                qirim(i,k,iice) = qitot(i,k,iice)
                birim(i,k,iice)  = qirim(i,k,iice)/rhop
             endif

! find indices in 4D ice lookup table
!------------------------------------------------------------------------------------------!

! find index for qi (total ice mass mixing ratio)
             dum1 = (alog10(qitot(i,k,iice))+16.)*1.41328
             dumi = int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
             dum1 = min(dum1,real(isize))
             dum1 = max(dum1,1.)
             dumi = max(1,dumi)
             dumi = min(isize-1,dumi)

! find index for Ni (ice number mixing ratio)
             dum2 = (alog10(nitot(i,k,iice))+10.)*1.10731
             dumk = int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
             dum2 = min(dum2,real(jsize))
             dum2 = max(dum2,1.)
             dumk = max(1,dumk)
             dumk = min(jsize-1,dumk)

! find index for scaled mean rain size
! if no rain, then just choose dumj = 1 and don't calculate rain-ice
! collection processes

             if (qr(i,k).ge.qsmall) then
! calculate scaled mean size for consistency with ice lookup table
                dumlr = (qr(i,k)/(pi*rhow*nr(i,k)))**thrd
                dum3  = (alog10(1.*dumlr)+5.)*10.70415
                dumj  = int(dum3)
! set limits
                dum3  = min(dum3,real(rcollsize))
                dum3  = max(dum3,1.)
                dumj  = max(1,dumj)
                dumj  = min(rcollsize-1,dumj)
             else
                dumj  = 1
                dum3  = 1.
             endif

! find index for rime mass fraction
             dum4  = qirim(i,k,iice)/qitot(i,k,iice)*3. + 1.
             dumii = int(dum4)

! set limits
             dum4  = min(dum4,real(rimsize))
             dum4  = max(dum4,1.)
             dumii = max(1,dumii)
             dumii = min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
             if (rhop.le.650.) then
                dum5 = (rhop-50.)*0.005 + 1.
             else
                dum5 =(rhop-650.)*0.004 + 4.
             endif
             dumjj = int(dum5)

! set limits
             dum5  = min(dum5,real(densize))
             dum5  = max(dum5,1.)
             dumjj = max(1,dumjj)
             dumjj = min(densize-1,dumjj)

!------------------------------------------------------------------------------------------!
!
! call to lookup table interpolation subroutines to get process rates
             call access_lookup_table(dumjj,dumii,dumi,dumk,2,dum1,dum2,dum4,dum5,f1pr2)
             call access_lookup_table(dumjj,dumii,dumi,dumk,3,dum1,dum2,dum4,dum5,f1pr3)
             call access_lookup_table(dumjj,dumii,dumi,dumk,4,dum1,dum2,dum4,dum5,f1pr4)
             call access_lookup_table(dumjj,dumii,dumi,dumk,5,dum1,dum2,dum4,dum5,f1pr5)
             call access_lookup_table(dumjj,dumii,dumi,dumk,7,dum1,dum2,dum4,dum5,f1pr9)
             call access_lookup_table(dumjj,dumii,dumi,dumk,8,dum1,dum2,dum4,dum5,f1pr10)
             call access_lookup_table(dumjj,dumii,dumi,dumk,10,dum1,dum2,dum4,dum5,f1pr14)

! ice-rain collection processes
             if (qr(i,k).ge.qsmall) then
                call access_lookup_table_coll(dumjj,dumii,dumj,dumi,dumk,1,dum1,dum2,    &
                                              dum3,dum4,dum5,f1pr7)
                call access_lookup_table_coll(dumjj,dumii,dumj,dumi,dumk,2,dum1,dum2,    &
                                              dum3,dum4,dum5,f1pr8)
             else
                f1pr7 = 0.
                f1pr8 = 0.
             endif

! adjust Ni if needed to make sure mean size is in bounds (i.e., apply lambda limiters)
             nitot(i,k,iice) = min(nitot(i,k,iice),f1pr9)
             nitot(i,k,iice) = max(nitot(i,k,iice),f1pr10)

! Determine additional collection efficiency factor to be applied to ice-ice collection.
! The value of the mass-weighted fall speed (f1pr2) is used to determine if the ice category
! represents graupel, where it is assumed that collection (by or of) will not occur.
! The computed values of qicol and nicol are multipiled by Eii_fact to gradually shut off collection
! if the ice in iice is considered to be graupel (Frim>0.5 and V>1 m/s).  (This approach is ad hoc;
! in the future the collection efficiencies will be modified in the ice-ice collection lookup table.)
! note: f1pr2 is the ice fall speed (without air density correction)
             if (qirim(i,k,iice)>0.) then
                tmp1 = qirim(i,k,iice)/qitot(i,k,iice)   !rime mass fraction
                if (tmp1>0.5 .and. f1pr2>1.) then
                   Eii_fact(iice) = max(0., (1.-(f1pr2-1.)))
                else
                   Eii_fact(iice) = 1.
                endif
             else
                Eii_fact(iice) = 1.
             endif

          endif qitot_notsmall_1

!......................................................................
! begin ice microphysical processes

!......................................................................
! collection between ice categories

! set up loop over category interactions

!        iceice_interaction1:  if (.false.) then
          iceice_interaction1:  if (iice.ge.2) then

             qitot_notsmall: if (qitot(i,k,iice).ge.qsmall) then

                catcoll_loop: do catcoll = 1,iice-1

                   qitotcatcoll_notsmall: if (qitot(i,k,catcoll).ge.qsmall) then

! find which category is collector and which is the collectee, based on
! difference in mean fallspeed

!             if (vmimean(iice).ge.vmimean(catcoll)) then

! fist, calculate collection of catcoll category by iice category

! find index in lookup table for collector category

! find index for qi (total ice mass mixing ratio)
                      dum1 = (alog10(qitot(i,k,iice))+16.5229)*0.816893
                      dumi = int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
                      dum1 = min(dum1,real(iisize))
                      dum1 = max(dum1,1.)
                      dumi = max(1,dumi)
                      dumi = min(iisize-1,dumi)

! find index for Ni (ice number mixing ratio)
                      dum2 = (alog10(nitot(i,k,iice))+10.5229)*0.640057
                      dumk = int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
                      dum2 = min(dum2,real(jjsize))
                      dum2 = max(dum2,1.)
                      dumk = max(1,dumk)
                      dumk = min(jjsize-1,dumk)

! note that the code below for finding rime mass fraction and density index is
! redundant with code for main ice lookup table and can probably be omitted
! for efficiency, for now it is left in

! find index for rime mass fraction
                      dum4  = qirim(i,k,iice)/qitot(i,k,iice)*3. + 1.
                      dumii = int(dum4)

! set limits
                      dum4  = min(dum4,real(rimsize))
                      dum4  = max(dum4,1.)
                      dumii = max(1,dumii)
                      dumii = min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
                      if (rhop.le.650.) then
                         dum5 = (rhop-50.)*0.005 + 1.
                      else
                         dum5 =(rhop-650.)*0.004 + 4.
                      endif
                      dumjj = int(dum5)

! set limits
                      dum5  = min(dum5,real(densize))
                      dum5  = max(dum5,1.)
                      dumjj = max(1,dumjj)
                      dumjj = min(densize-1,dumjj)

! find index in lookup table for collectee category, here 'q' is a scaled q/N
! find index for qi (total ice mass mixing ratio)
                      dum1c = (alog10(qitot(i,k,catcoll)/nitot(i,k,catcoll))+16.)
                      dumic = int(dum1c)

! set limits to make sure the calculated index doesn't exceed range of lookup table
                      dum1c = min(dum1c,real(iisize))
                      dum1c = max(dum1c,1.)
                      dumic = max(1,dumic)
                      dumic = min(iisize-1,dumic)

! find index for rime mass fraction
                      dum4c  = qirim(i,k,catcoll)/qitot(i,k,catcoll)*3. + 1.
                      dumiic = int(dum4c)

! set limits
                      dum4c  = min(dum4c,real(rimsize))
                      dum4c  = max(dum4c,1.)
                      dumiic = max(1,dumiic)
                      dumiic = min(rimsize-1,dumiic)

! calculate predicted bulk rime density
                      if (birim(i,k,catcoll).ge.1.e-15) then
                         drhop = qirim(i,k,catcoll)/birim(i,k,catcoll)
                      else
                         drhop = 0.
                      endif

! find index for bulk rime density
! account for uneven spacing in lookup table for density
                      if (drhop.le.650.) then
                         dum5c = (drhop-50.)*0.005 + 1.
                      else
                         dum5c =(drhop-650.)*0.004 + 4.
                      endif
                      dumjjc = int(dum5c)

! set limits
                      dum5c  = min(dum5c,real(densize))
                      dum5c  = max(dum5c,1.)
                      dumjjc = max(1,dumjjc)
                      dumjjc = min(densize-1,dumjjc)

                      call access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj,  &
                                                     dumi,dumk,1,dum1c,dum4c,dum5c,dum1,    &
                                                     dum2,dum4,dum5,f1pr17)
                      call access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj,  &
                                                     dumi,dumk,2,dum1c,dum4c,dum5c,dum1,    &
                                                     dum2,dum4,dum5,f1pr18)

! we need to multiply by air density, air density fallspeed correction factor, and N of the collectee category
! for process rates nicol and qicol, first index is the collectee, second is
! the collector
                      nicol(catcoll,iice) = f1pr17*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*      &
                                            nitot(i,k,catcoll)
                      qicol(catcoll,iice) = f1pr18*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*      &
                                            nitot(i,k,catcoll)

                      nicol(catcoll,iice) = Eii_fact(iice)*nicol(catcoll,iice)
                      qicol(catcoll,iice) = Eii_fact(iice)*qicol(catcoll,iice)
  !===

!             else if (vmimean(catcoll).gt.vmimean(iice)) then

!............................................................................................
! second, calculate collection of iice category by catcoll category

!......................................................
! find index in lookup table for collector category

! find index for qi (total ice mass mixing ratio)
                      dum1 = (alog10(qitot(i,k,catcoll))+16.5229)*0.816893
                      dumi = int(dum1)

! set limits to make sure the calculated index does not exceed range of lookup table
                      dum1 = min(dum1,real(iisize))
                      dum1 = max(dum1,1.)
                      dumi = max(1,dumi)
                      dumi = min(iisize-1,dumi)

! find index for Ni (ice number mixing ratio)
                      dum2 = (alog10(nitot(i,k,catcoll))+10.5229)*0.640057
                      dumk = int(dum2)

! set limits to make sure the calculated index does not exceed range of lookup table
                      dum2 = min(dum2,real(jjsize))
                      dum2 = max(dum2,1.)
                      dumk = max(1,dumk)
                      dumk = min(jjsize-1,dumk)

! note that the code below for finding rime mass fraction and density index is
! redundant with code for main ice lookup table and can probably be omitted
! for efficiency, for now it is left in

! find index for rime mass fraction
                      dum4  = qirim(i,k,catcoll)/qitot(i,k,catcoll)*3. + 1.
                      dumii = int(dum4)

! set limits
                      dum4  = min(dum4,real(rimsize))
                      dum4  = max(dum4,1.)
                      dumii = max(1,dumii)
                      dumii = min(rimsize-1,dumii)

! calculate predicted bulk rime density
                      if (birim(i,k,catcoll).ge.1.e-15) then
                         drhop = qirim(i,k,catcoll)/birim(i,k,catcoll)
                      else
                         drhop = 0.
                      endif

! find index for bulk rime density
! account for uneven spacing in lookup table for density
                      if (drhop.le.650.) then
                         dum5 = (drhop-50.)*0.005 + 1.
                      else
                         dum5 =(drhop-650.)*0.004 + 4.
                      endif
                      dumjj = int(dum5)

! set limits
                      dum5  = min(dum5,real(densize))
                      dum5  = max(dum5,1.)
                      dumjj = max(1,dumjj)
                      dumjj = min(densize-1,dumjj)

! find index in lookup table for collectee category, here 'q' is a scaled q/N
! find index for qi (total ice mass mixing ratio)
                      dum1c = (alog10(qitot(i,k,iice)/nitot(i,k,iice))+16.)
                      dumic = int(dum1c)

! set limits to make sure the calculated index doesn't exceed range of lookup table
                      dum1c = min(dum1c,real(iisize))
                      dum1c = max(dum1c,1.)
                      dumic = max(1,dumic)
                      dumic = min(iisize-1,dumic)

! find index for rime mass fraction
                      dum4c  = qirim(i,k,iice)/qitot(i,k,iice)*3. + 1.
                      dumiic = int(dum4c)

! set limits
                      dum4c  = min(dum4c,real(rimsize))
                      dum4c  = max(dum4c,1.)
                      dumiic = max(1,dumiic)
                      dumiic = min(rimsize-1,dumiic)

! calculate predicted bulk rime density
                      if (birim(i,k,iice).ge.1.e-15) then
                         drhop = qirim(i,k,iice)/birim(i,k,iice)
                      else
                         drhop = 0.
                      endif

! find index for bulk rime density
! account for uneven spacing in lookup table for density
                      if (drhop.le.650.) then
                         dum5c = (drhop-50.)*0.005 + 1.
                      else
                         dum5c = (drhop-650.)*0.004 + 4.
                      endif
                      dumjjc = int(dum5c)

! set limits
                      dum5c  = min(dum5c,real(densize))
                      dum5c  = max(dum5c,1.)
                      dumjjc = max(1,dumjjc)
                      dumjjc = min(densize-1,dumjjc)

                      call access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj, &
                                                     dumi,dumk,1,dum1c,dum4c,dum5c,dum1,   &
                                                     dum2,dum4,dum5,f1pr17)
                      call access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj, &
                                                     dumi,dumk,2,dum1c,dum4c,dum5c,dum1,   &
                                                     dum2,dum4,dum5,f1pr18)

! we need to multiply by air density, air density fallspeed correction factor, and N of the collectee category
! for process rates nicol and qicol, first index is the collectee, second is
! the collector

                      nicol(iice,catcoll) = f1pr17*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*     &
                                            nitot(i,k,iice)
                      qicol(iice,catcoll) = f1pr18*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*     &
                                            nitot(i,k,iice)

                      nicol(iice,catcoll) = Eii_fact(iice)*nicol(iice,catcoll)
                      qicol(iice,catcoll) = Eii_fact(iice)*qicol(iice,catcoll)

!===

!             end if ! iice vmi versus catcoll vmi, for determining collector versus collectee

                   endif qitotcatcoll_notsmall

                enddo catcoll_loop

             endif qitot_notsmall

          endif iceice_interaction1  ! iice ge 2

!......................................................................
! other ice microphysical processes

! collection of droplets
! here we multiply rates by air density, air density fallspeed correction
! factor, and collection diag_efficiency since these parameters are not
! included in lookup table calculations

! for T < 273.15, assume collected cloud water is instantly frozen
          if (qitot(i,k,iice).ge.qsmall .and. qc(i,k).ge.qsmall .and. t(i,k).le.273.15) then
             qccol(iice) = rhofaci(i,k)*f1pr4*qc(i,k)*eci*rho(i,k)
             nccol(iice) = rhofaci(i,k)*f1pr4*nc(i,k)*eci*rho(i,k)
          endif

! for T > 273.15, assume cloud water is collected and shed as rain drops
          if (qitot(i,k,iice).ge.qsmall .and. qc(i,k).ge.qsmall .and.t(i,k).gt.273.15) then
! sink for cloud water mass and number, note qcshed is source for rain mass
             qcshd(iice) = rhofaci(i,k)*f1pr4*qc(i,k)*eci*rho(i,k)
             nccol(iice) = rhofaci(i,k)*f1pr4*nc(i,k)*eci*rho(i,k)
! source for rain number, assume 1mm drops are shed
             ncshdc(iice) = qcshd(iice)*1.923e+6
          endif

!............................................................
! collection of rain
! here we multiply rates by air density, air density fallspeed correction
! factor, collection diag_efficiency, and n0r since these parameters are not
! included in lookup table calculations

! for T < 273.15, assume all collected rain mass freezes
! note this is a sink for rain mass and number and a source
! for ice mass

          if (qitot(i,k,iice).ge.qsmall .and. qr(i,k).ge.qsmall .and.                    &
           t(i,k).le.273.15) then
!            qrcol(iice)  = f1pr8*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
!            nrcol(iice) = f1pr7*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
! note: f1pr8 and n0r are already calculated as log_10
             qrcol(iice)  = 10.**(f1pr8+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
             nrcol(iice) = 10.**(f1pr7+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
          endif

! for T > 273.15, assume collected rain number is shed as
! 1 mm drops
! note that melting of ice number is scaled to the loss
! rate of ice mass due to melting
! collection of rain above freezing does not impact total rain mass

          if (qitot(i,k,iice).ge.qsmall .and. qr(i,k).ge.qsmall .and. t(i,k).gt.273.15) then
! rain number sink due to collection
             nrcol(iice) = 10.**(f1pr7+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
! rain number source due to shedding = collected rain mass/mass of 1 mm drop
             dum     = 10.**(f1pr8+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
             nrshdr(iice) = dum*1.923e6   ! 1./5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
! note: pracsw1(iice) not currently used
!            pracsw1(iice)=f1pr8*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
!            npracsw1(iice)=f1pr11*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
          endif

!............................................................
! self-collection of ice
! here we multiply rates by collection diag_efficiency, air density,
! and air density correction factor since these are not included
! in the lookup table calculations

          if (qitot(i,k,iice).ge.qsmall) then

            !--- The following code is copied from multicategory collection section.
            !    (eventually, this will have to be formalized)
             if (qirim(i,k,iice)>0.) then
                tmp1 = qirim(i,k,iice)/qitot(i,k,iice)   !rime mass fraction
                if (tmp1>0.5 .and. f1pr2>1.) then
                   Eii_fact(iice) = max(0., (1.-(f1pr2-1.)))
                else
                   Eii_fact(iice) = 1.
                endif
             else
                Eii_fact(iice) = 1.
             endif
             !===

             nislf(iice) = f1pr3*rho(i,k)*eii*rhofaci(i,k)*Eii_fact(iice)
          endif

!............................................................
! melting
! include accelerated melting due to collection of ice mass by rain (pracsw1(iice))
! this logic follows Lin et al. 1983 and others

          if (qitot(i,k,iice).ge.qsmall .and. t(i,k).gt.273.15) then
             qsat0 = 0.622*e0/(pres(i,k)-e0)
!            dum = cpw/xlf(i,k)*(t(i,k)-273.15)*(pracsw1(iice)+qcshd(iice))
! currently enhanced melting from collision is neglected
!            dum = cpw/xlf(i,k)*(t(i,k)-273.15)*(pracsw1(iice))
             dum = 0.
!            qimlt(iice) = (f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)*    &
!                          (t(i,k)-273.15)*2.*pi*kap/xlf(i,k)+dum
! include RH dependence
             qimlt(iice) = (f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)*    &
                           ((t(i,k)-273.15)*kap-rho(i,k)*xxlv(i,k)*dv*(qsat0-qv(i,k)))*  &
                           2.*pi/xlf(i,k)+dum
             qimlt(iice) = max(qimlt(iice),0.)
!            dum1  = (f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)* &
!                    (t(i,k)-273.15)*2.*pi*rho(i,k)*kap/xlf(i,k)
          endif

!............................................................
! calculate wet growth, similar to Musil (1970), JAS
          if (qitot(i,k,iice).ge.qsmall .and. qc(i,k)+qr(i,k).ge.1.e-6 .and.             &
           t(i,k).lt.273.15) then

             qsat0  = 0.622*e0/(pres(i,k)-e0)
             qwgrth = (f1pr5+f1pr14*sc**thrd*(rhofaci(i,k)*rho(i,k)/mu)**0.5)*2.*pi*     &
                       (rho(i,k)*xxlv(i,k)*dv*(qsat0-qv(i,k))-(t(i,k)-273.15)*kap)/      &
                       (xlf(i,k)+cpw*(t(i,k)-273.15))
             qwgrth = max(qwgrth,0.)
! calculate shedding for wet growth
             dum     = max(0.,(qccol(iice)+qrcol(iice))-qwgrth)

             if (dum.ge.1.e-10) then
                nrshdr(iice) = nrshdr(iice)+dum*1.923e6   ! 1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
                if ((qccol(iice)+qrcol(iice)).ge.1.e-10) then
                   dum1        = 1./(qccol(iice)+qrcol(iice))
                   qcshd(iice) = qcshd(iice)+dum*qccol(iice)*dum1
                   qccol(iice) = qccol(iice)-dum*qccol(iice)*dum1
                   qrcol(iice) = qrcol(iice)-dum*qrcol(iice)*dum1
                endif
! densify due to wet growth
                wetgrowth(iice) = .true.
                endif
             endif


! sensitivity, no immersion freezing
!      qcheti = 0.
!      ncheti = 0.
!      qrheti = 0.
!      nrheti = 0.

!.................................................................................
! calculate rime density

! NOTE: Tc (ambient) is assumed for the surface temperature.  Technically,
! we should diagose graupel surface temperature from heat balance equation.
! (but the ambient temperature is a reasonable approximation; tests show
! very little sensitivity to different assumed values, Milbrandt and Morrison 2012).

! Compute rime density: (based on parameterization of Cober and List, 1993 [JAS])
! for simplicty use mass-weighted ice and droplet/rain fallspeeds

!          rhorime_c(iice) = rhorime_c_dflt

          if (qitot(i,k,iice).ge.qsmall .and. t(i,k).lt.273.15) then

! get mass-weighted mean ice fallspeed
             vtrmi1(i,k) = f1pr2*rhofaci(i,k)
             iTc         = 1./min(-0.001,t(i,k)-273.15)

! First for cloud droplets
             if (qc(i,k).ge.qsmall) then
!...................................................................
! droplet fall speed
! use Stokes' formulation
! all droplets in smallest dactegory fallspeed
! this means we can use analytic solution
                vtrmc(i,k) = acn(i,k)*gamma(4.+bcn+pgam(i,k))/(lamc(i,k)**bcn*           &
                             gamma(pgam(i,k)+4.))
! use mass-weighted mean size
                D_c = (pgam(i,k)+4.)/lamc(i,k)
                V_impact = abs(vtrmi1(i,k)-vtrmc(i,k))
                Ri       = -(0.5e+6*D_c)*V_impact*iTc
!               Ri       = max(1.,min(Ri,8.))
                Ri       = max(1.,min(Ri,12.))
                if (Ri.le.8.) then
                   rhorime_c(iice)  = (0.051 + 0.114*Ri - 0.0055*Ri**2)*1000.
                else
! for Ri > 8 assume a linear fit between 8 and 12,
! rhorime = 900 kg m-3 at Ri = 12
! this is somewhat ad-hoc but allows a smoother transition
! in rime density up to wet growth
                   rhorime_c(iice)  = 611.+72.25*(Ri-8.)
                endif
             endif

          endif ! qi > qsmall and T < 273.15

! Next for rain drops
! assume rime density for rain collecting ice is 900 kg/m3
!      if (qr(i,k).ge.qsmall) then
!         D_r       = (mur(i,k)+1.)/lamr(i,k)
!         V_impact  = abs(vtrmi1(i,k)-vtrm(i,k))
!         Ri        = -(0.5e+6*D_r)*V_impact*iTc
!         Ri        = max(1.,min(Ri,8.))
!         rhorime_r(iice) = (0.051 + 0.114*Ri - 0.0055*Ri**2)*1000.
!      else
!         rhorime_r = rhorime_r_dflt
!      endif

       enddo iice_loop1

       log_hmossopOn = (n_iceCat.gt.1)  !default, off for 1-category
!      log_hmossopOn = .true.           !switch to have Hallet-Mossop ON
!      log_hmossopOn = .false.          !switch to have Hallet-Mossop OFF

       if (log_hmossopOn) then

! find ice category for which rime splintering adds ice mass and number
         !--determine destination ice-phase category:
          D_new = 10.e-6 !assumes ice crystals from rime splintering are tiny
          call icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,       &
                                  ni_add,iice_dest)
         !==

          do iice = 1,n_iceCat

! rime-splintering (Hallett-Mossop 1974)
! calculate rime-splintering AFTER conservation to avoid
! unrealistic riming rates
             if (qitot(i,k,iice).ge.qsmall.and. (qccol(iice).gt.0. .or.                  &
              qrcol(iice).gt.0.)) then

                if (t(i,k).gt.270.15) then
                   dum = 0.
                elseif (t(i,k).le.270.15 .and. t(i,k).gt.268.15) then
                   dum = (270.15-t(i,k))*0.5
                elseif (t(i,k).le.268.15.and.t(i,k).ge.265.15) then
                   dum = (t(i,k)-265.15)*thrd
                elseif (t(i,k).lt.265.15) then
                   dum = 0.
                endif

! rime splintering from riming of cloud droplets
!                dum1 = 35.e4*qccol(iice)*dum*1000. ! 1000 is to convert kg to g
!                dum2 = dum1*pi/6.*900.*(10.e-6)**3  ! assume 10 micron splinters

!                qccol(iice) = qccol(iice)-dum2 ! subtract splintering from rime mass transfer

!                if (qccol(iice) .lt. 0.) then
!                   dum2 = qccol(iice)
!                   qccol(iice) = 0.
!                endif

!                qcmul(iice_dest) = qcmul(iice_dest)+dum2
!                if (ni_add) then
!                nimul(iice_dest) = nimul(iice_dest)+dum2/(pi/6.*900.*(10.e-6)**3)
!                end if

! rime splintering from riming of raindrops
                dum1 = 35.e4*qrcol(iice)*dum*1000.! 1000 is to convert kg to g
                dum2 = dum1*pi/6.*900.*(10.e-6)**3 ! assume 10 micron splinters

                qrcol(iice) = qrcol(iice)-dum2 ! subtract splintering from rime mass transfer

                if (qrcol(iice) .lt. 0.) then
                   dum2 = qrcol(iice)
                   qrcol(iice) = 0.
                endif

                qrmul(iice_dest) = qrmul(iice_dest) + dum2
                if (ni_add) then
                   nimul(iice_dest) = nimul(iice_dest) + dum2/(pi/6.*900.*(10.e-6)**3)
                end if

             endif

          enddo !iice-loop

       endif  ! logical for Hallet-Mossop

       do iice = 1,n_iceCat

! reduce number concentration from sublimation or melting
          if (qisub(iice).gt.0. .and. qitot(i,k,iice).ge.qsmall) then
             dum         = -qisub(iice)*dt/qitot(i,k,iice)
             dum         = max(-1.,dum)
             nisub(iice) = dum*nitot(i,k,iice)*odt
          endif
          if (qimlt(iice).gt.0. .and. qitot(i,k,iice).ge.qsmall) then
             dum          = -qimlt(iice)*dt/qitot(i,k,iice)
             dum          = max(-1.,dum)
             nimlt(iice)  = dum*nitot(i,k,iice)*odt
          endif

! calcualte total inverse ice relaxation timescale combined for all ice categories
          if (qitot(i,k,iice).ge.qsmall .and. t(i,k).lt.273.15) then
             epsi(iice) = (f1pr5+f1pr14*sc**thrd*(rhofaci(i,k)*rho(i,k)/mu)**0.5)*2.*pi* &
                          rho(i,k)*dv
             epsi_tot   = epsi_tot+epsi(iice)
          else
             epsi(iice) = 0.
          endif

       enddo  !iice-loop

! end ice microphysical processes

!............................................................
! vapor deposition/sublimation onto ice
! note: capacitance factor already incluced in lookup table

! condensation/evaporation/deposition/sublimation
! use semi-analytic formulation

       if (t(i,k).lt.273.15) then
          oabi = 1./abi
          xx   = epsc+epsr+epsi_tot*(1.+xxls(i,k)*inv_cp*dqsdt)*oabi
       else
          xx   = epsc+epsr
       endif

       dum    = qvs(i,k)*rho(i,k)*g*uzpl(i,k)/(pres(i,k)-polysvp1(t(i,k),0))
       dumqvi = qvi(i,k)

! 'A' term including ice (bergeron process)
! note: qv and T tendencies due to mixing and radiation are
! currently neglected --> assumed to be much smaller than cooling
! due to vertical motion which IS included

       if (t(i,k).lt.273.15) then
          aaa = -dum-dqsdt*(-uzpl(i,k)*g*inv_cp)-(qvs(i,k)-dumqvi)*(1.+xxls(i,k)*inv_cp* &
                 dqsdt)*oabi*epsi_tot
       else
          aaa = -dum-dqsdt*(-uzpl(i,k)*g*inv_cp)
       endif

! set lower bound on xx to prevent division by zero

       xx  = max(1.e-8,xx)
       oxx = 1./xx

       if (qc(i,k).ge.qsmall) then
          qccon = (aaa*epsc*oxx+(ssat(i,k)-aaa*oxx)*odt*epsc*oxx*(1.-dexp(-dble(xx*dt))))/ab
       endif
       if (qr(i,k).ge.qsmall) then
          qrcon = (aaa*epsr*oxx+(ssat(i,k)-aaa*oxx)*odt*epsr*oxx*(1.-dexp(-dble(xx*dt))))/ab
       endif

       do iice = 1,n_iceCat

          if (qitot(i,k,iice).ge.qsmall .and. t(i,k).lt.273.15) then
             qidep(iice) = (aaa*epsi(iice)*oxx+(ssat(i,k)-aaa*oxx)*odt*epsi(iice)*oxx*    &
                          (1.-dexp(-dble(xx*dt))))*oabi+(qvs(i,k)-dumqvi)*epsi(iice)*oabi
          endif

          if (qv(i,k)-qvi(i,k).lt.0.) qidep(iice) = min(0.,qidep(iice))

! for very small water contents, evaporate instantly
          if (supi(i,k).lt.-0.001 .and. qitot(i,k,iice).lt.1.e-12) then
             qidep(iice) = -qitot(i,k,iice)*odt
          endif

          if (qidep(iice).lt.0.) then
             qisub(iice) = -qidep(iice)
             qidep(iice)  = 0.
          endif

       enddo !iice-loop

       if (ssat(i,k).lt.0.) qccon = min(0.,qccon)
       if (ssat(i,k).lt.0.) qrcon = min(0.,qrcon)

! for very small water contents, evaporate instantly
       if (sup(i,k).lt.-0.001 .and. qc(i,k).lt.1.e-12) then
          qccon = -qc(i,k)*odt
       endif
       if (sup(i,k).lt.-0.001 .and. qr(i,k).lt.1.e-12) then
          qrcon = -qr(i,k)*odt
       endif

       if (qccon.lt.0.) then
          qcevp = -qccon
          qccon  = 0.
       endif
       if (qrcon.lt.0.) then
          qrevp = -qrcon
          qrcon  = 0.
       endif

! in wetgrowth regime, add ice deposition to rain condensation and shed
! neglect for now, instead assume deposition onto wet growth hail goes to ice
!       if (wetgrowth(iice)) then
!          if (qidep(iice).gt.0.) then
!             nrshdr(iice) = nrshdr(iice) + qidep(iice)*1.923e6   ! 1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
!             qrcon     = qrcon + qidep(iice)
!             qidep(iice)    = 0.
!          endif
!       endif


!------------------------------------------------------------------------------------------!
! begin liquid microphysical process calculations

!............................................................
! autoconversion

       if (qc(i,k).ge.1.e-8) then

! Seifert and Beheng (2001)
          if (iparam.eq.1) then

             dum   = 1.-qc(i,k)/(qc(i,k)+qr(i,k))
             dum1  = 600.*dum**0.68*(1.-dum**0.68)**3
             qcaut = kc*1.92307e+5*(nu(i,k)+2.)*(nu(i,k)+4.)/(nu(i,k)+1.)**2*(rho(i,k)*  &
                     qc(i,k)*0.001)**4/(rho(i,k)*nc(i,k)/1.e6)**2*(1.+dum1/(1.-dum)**2)* &
                     1000.*inv_rho(i,k)
             ncautc = qcaut*7.6923e+9

! Beheng (1994)
          elseif (iparam.eq.2) then

             if (nc(i,k)*1.e-6*rho(i,k).lt.100.) then
                qcaut = 6.e28*inv_rho(i,k)*pgam(i,k)**(-1.7)*(1.e-6*rho(i,k)*nc(i,k))**  &
                        (-3.3)*(1.e-3*rho(i,k)*qc(i,k))**4.7
             else

! 2D interpolation of tabled logarithmic values
                dum   = 41.46+(nc(i,k)*1.e-6*rho(i,k)-100.)*(37.53-41.46)*0.005
                dum1  = 39.36+(nc(i,k)*1.e-6*rho(i,k)-100.)*(30.72-39.36)*0.005
                qcaut = dum+(pgam(i,k)-5.)*(dum1-dum)*0.1

! 1000/rho is for conversion from g cm-3/s to kg/kg
                qcaut  = exp(qcaut)*(1.e-3*rho(i,k)*qc(i,k))**4.7*1000.*inv_rho(i,k)

             endif

             ncautc = 7.7e9*qcaut

! Khroutdinov and Kogan (2000)
          elseif (iparam.eq.3) then

             dum   = qc(i,k)
             qcaut   = 1350.*dum**2.47*(nc(i,k)*1.e-6*rho(i,k))**(-1.79)
! note: ncautr is change in Nr,
! ncautc is change in Nc
             ncautr = qcaut*cons3
             ncautc  = qcaut*nc(i,k)/qc(i,k)

          endif

          if (qcaut.eq.0.)  ncautc = 0.
          if (ncautc.eq.0.) qcaut  = 0.

       endif
!............................................................
! self collection of droplets
       if (qc(i,k).ge.qsmall) then

          if (iparam.eq.1.) then
             ncslf = -kc*(1.e-3*rho(i,k)*qc(i,k))**2*(nu(i,k)+2.)/(nu(i,k)+1.)*1.e+6*    &
                     inv_rho(i,k)+ncautc
          elseif (iparam.eq.2.) then
             ncslf = -5.5e16*inv_rho(i,k)*pgam(i,k)**(-0.63)*(1.e-3*rho(i,k)*qc(i,k))**2
          elseif (iparam.eq.3) then
             ncslf = 0.
          endif

       endif
!............................................................
! accretion of cloud droplets by rain
       if (qr(i,k).ge.qsmall .and. qc(i,k).ge.qsmall) then

          if (iparam.eq.1) then

! Seifert and Beheng (2001) formulation:
             dum   = 1.-qc(i,k)/(qc(i,k)+qr(i,k))
             dum1  = (dum/(dum+5.e-4))**4
             qcacc = kr*rho(i,k)*0.001*qc(i,k)*qr(i,k)*dum1
             ncacc = qcacc*rho(i,k)*0.001*(nc(i,k)*rho(i,k)*1.e-6)/(qc(i,k)*rho(i,k)*    &
                    0.001)*1.e6*inv_rho(i,k)

          elseif (iparam.eq.2) then
! Beheng (1994) formulation:
             dum   = (qc(i,k)*qr(i,k))
             qcacc = 6.*rho(i,k)*dum
             ncacc = qcacc*rho(i,k)*0.001*(nc(i,k)*rho(i,k)*1.e-6)/(qc(i,k)*rho(i,k)*    &
                    0.001)*1.e6*inv_rho(i,k)

          elseif (iparam.eq.3) then
! Khairoutdinov and Kogan (2000) formulation:
             qcacc = 67.*(qc(i,k)*qr(i,k))**1.15
             ncacc = qcacc*nc(i,k)/qc(i,k)

          endif

          if (qcacc.eq.0.) ncacc = 0.
          if (ncacc.eq.0.) qcacc = 0.
!         ncacc = ncacc*(qcacc.ne.0.)*(-1.)  !optimization for: if (qcacc.eq.0.) ncacc = 0.
!         qcacc = qcacc*(ncacc.ne.0.)*(-1.)  !optimization for: if (ncacc.eq.0.) qcacc  = 0.

       endif
!............................................................
! self-collection and breakup of rain
! add breakup following modified Verlinde and Cotton scheme
       if (qr(i,k).ge.qsmall) then

! include breakup
          dum1 = 300.e-6

! use mass-mean diameter (do this by using
! the old version of lambda w/o mu dependence)
! note there should be a factor of 6^(1/3), but we
! want to keep breakup threshold consistent so 'dum'
! is expressed in terms of lambda rather than mass-mean D
          dum2 = (QR(i,k)/(PI*rhoW*NR(i,k)))**thrd
          if (dum2.lt.dum1) then
             dum = 1.
          else
             dum = 2.-exp(2300.*(dum2-dum1))
          endif

          if (iparam.eq.1.) then
             nrslf = -dum*kr*1.e-3*qr(i,k)*nr(i,k)*rho(i,k)
          elseif (iparam.eq.2.or.iparam.eq.3) then
             nrslf = -dum*5.78*nr(i,k)*qr(i,k)*rho(i,k)
          endif

! add spontaneous breakup (i.e., independent of collection)
! treat as nudging over a 10 sec timescale
          if (dum2.gt.550.e-6) then
             dum   = qr(i,k)*cons4
             dum1  = (dum-nr(i,k))/max(10.,dt)
             nrslf = nrslf+dum1
          endif

       endif

!................................................................
! evaporate/melt number concentration of rain
!!       onrevp(i,k)=0.
       if (qrevp.gt.0. .and. qr(i,k).ge.qsmall) then
          dum   = -qrevp*dt/qr(i,k)
          dum   = max(-1.,dum)
!         dum1  = 0.5
! mu dependence from Seifert 2008, JAS, neglecting size dependence
          dum1  = exp(-0.2*mur(i,k))
          nrevp = dum1*dum*nr(i,k)
!!          onrevp(i,k) = dum1*dum*nr(i,k)*odt
       endif

!............................................................
! contact and immersion freezing droplets

       if (qc(i,k).ge.qsmall .and. t(i,k).le.269.15) then
!         qchetc(iice) = pi*pi/3.*Dap*Nacnt*rhow*cdist1(i,k)*gamma(pgam(i,k)+5.)/lamc(i,k)**4
!         nchetc(iice) = 2.*pi*Dap*Nacnt*cdist1(i,k)*gamma(pgam(i,k)+2.)/lamc(i,k)
! for future: calculate gamma(pgam+4) in one place since its used multiple times
          dum    = (1./lamc(i,k))**3
!         qcheti(iice_dest) = cons6*cdist1(i,k)*gamma(7.+pgam(i,k))*exp(aimm*(273.15-t(i,k)))*dum**2
!         ncheti(iice_dest) = cons5*cdist1(i,k)*gamma(pgam(i,k)+4.)*exp(aimm*(273.15-t(i,k)))*dum
          Q_nuc = cons6*cdist1(i,k)*gamma(7.+pgam(i,k))*exp(aimm*(273.15-t(i,k)))*dum**2
          N_nuc = cons5*cdist1(i,k)*gamma(pgam(i,k)+4.)*exp(aimm*(273.15-t(i,k)))*dum
         !--determine destination ice-phase category:
          dum1      = 900.     !density of new ice
          D_new     = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
          call icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,      &
                                  ni_add,iice_dest)
         !==
          qcheti(iice_dest) = Q_nuc
          if (ni_add) ncheti(iice_dest) = N_nuc
       endif
!............................................................
! immersion freezing of rain
! for future: get rid of log statements below for rain freezing
!        if (qr(i,k).ge.qsmall.and.t(i,k).le.269.15) then
!           qrheti(iice_dest) = cons6*exp(log(cdistr(i,k))+log(gamma(7.+mur(i,k)))-6.*     &
!                    log(lamr(i,k)))*exp(aimm*(273.15-T(i,k)))
!           nrheti(iice_dest) = cons5*exp(log(cdistr(i,k))+log(gamma(mur(i,k)+4.))-3.*     &
!                    log(lamr(i,k)))*exp(aimm*(273.15-T(i,k)))
!        endif
       if (qr(i,k).ge.qsmall.and.t(i,k).le.269.15) then

          Q_nuc = cons6*exp(log(cdistr(i,k))+log(gamma(7.+mur(i,k)))-6.*log(lamr(i,k)))* &
                  exp(aimm*(273.15-T(i,k)))
          N_nuc = cons5*exp(log(cdistr(i,k))+log(gamma(mur(i,k)+4.))-3.*log(lamr(i,k)))* &
                  exp(aimm*(273.15-T(i,k)))
         !--determine destination ice-phase category:
          dum1      = 900.     !density of new ice
          D_new     = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
          call icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,        &
                                  ni_add,iice_dest)
         !==
          qrheti(iice_dest) = Q_nuc
          if (ni_add) nrheti(iice_dest) = N_nuc

       endif

! end liquid microphysical processes

 444   continue

!................................................................
! deposition/condensation-freezing nucleation
! allow ice nucleation if < -5 C and > 0.1% ice supersaturation

       if (t(i,k).lt.258.15 .and. supi(i,k).ge.0.05) then

!         dum = exp(-0.639+0.1296*100.*supi(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
          dum = 0.005*exp(0.304*(273.15-t(i,k)))*1000.*inv_rho(i,k)   !Cooper (1986)
          dum = min(dum,100.e3*inv_rho(i,k))

          N_nuc = max(0.,(dum-sum(nitot(i,k,:)))*odt)

          if (N_nuc.ge.1.e-20) then
             Q_nuc = max(0.,(dum-sum(nitot(i,k,:)))*mi0*odt)
             !--determine destination ice-phase category:
             dum1      = 900.     !density of new ice
             D_new     = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
             call icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,    &
                                     ni_add,iice_dest)
             !==
             qinuc(iice_dest) = Q_nuc
             if (ni_add) then
             ninuc(iice_dest) = N_nuc
             end if
           endif

       endif

!.................................................................
! droplet activation

! for specified Nc, make sure droplets are present if conditions are supersaturated

       if (sup(i,k).gt.1.e-6 .or. it.eq.0) then
          dum   = nccnst*inv_rho(i,k)*cons7-qc(i,k)
          dum   = max(0.,dum)
          dqsdt = xxlv(i,k)*qvs(i,k)/(rv*t(i,k)*t(i,k))
          ab    = 1.+dqsdt*xxlv(i,k)*inv_cp

! don't over-deplete supersaturation
          dum   = min(dum,(qv(i,k)-qvs(i,k))/ab)
          qcnuc  = dum*odt
       endif

       goto 713

       if (sup(i,k).gt.1.e-6 .or. it.eq.0) then
          sigvl = 0.0761-1.55e-4*(t(i,k)-273.15)
          aact  = 2.*mw/(rhow*rr*t(i,k))*sigvl
          sm1   = 2.*inv_bact**0.5*(aact*thrd*inv_rm1)**1.5
          sm2   = 2.*inv_bact**0.5*(aact*thrd*inv_rm2)**1.5
          if (it.eq.0) sup(i,k) = 0.001
          uu1   = 2.*log(sm1/sup(i,k))/(4.242*log(sig1))
!         uu2   = 2.*log(sm2/sup(i,k))/(4.242*log(sig2))
          dum1  = nanew1*0.5*(1.-derf(uu1))
!         dum2  = nanew2*0.5*(1.-derf(uu2))
          dum2  = dum1*inv_rho(i,k)  !convert to kg-1
! make sure this value is not greater than total number of aerosol
          dum2  = min((nanew1+nanew2)*inv_rho(i,k),dum2)
          dum2  = (dum2-nc(i,k))*odt
          dum2  = max(0.,dum2)
          ncnuc = dum2
          if (it.eq.0) then
             qcnuc  = 0.
          else
             qcnuc = ncnuc*cons7
          endif

       endif

 713   continue

!................................................................
! saturation adjustment to get initial cloud water
! This is only called once at the beginning of the simulation
! to remove any supersaturation in the intial conditions
       if (it.eq.0) then
!          dumt   = th(i,k)*(1.e-5*pres(i,k))**(rd*inv_cp)
!          dumqv  = qv(i,k)
! Fix to prevent over-depletion of vapor if conditions are saturated at first
! time step
          dumt   = th(i,k)*(1.e-5*pres(i,k))**(rd*inv_cp)+qcnuc*xxlv(i,k)*inv_cp*dt
          dumqv  = qv(i,k)-qcnuc*dt
          dumqvs = ep_2*polysvp1(dumt,0)/(pres(i,k)-polysvp1(dumt,0))
          dums   = dumqv-dumqvs
          qccon  = dums/(1.+xxlv(i,k)**2*dumqvs/(cp*rv*dumt**2))*odt
          qccon  = max(0.,qccon)
          if (qccon.le.1.e-7) qccon = 0.
       endif

!------------------------------------------------------------------------------------------!
!conservation of water

! cloud water
       dum  = (qcaut+qcacc+qcevp+sum(qchetc)+sum(qcheti)+sum(qcshd)+sum(qccol)+          &
              sum(qcmul))*dt
       dum1 = qc(i,k)+(qccon+qcnuc)*dt
       if (dum.gt.dum1 .and. dum.ge.1.e-20) then
          ratio   = dum1/dum
          qcaut   = qcaut*ratio
          qcacc   = qcacc*ratio
          qcevp   = qcevp*ratio
! NOTE: processes below are multiplied by ratio for each array element
          qchetc  = qchetc*ratio
          qcheti  = qcheti*ratio
          qcshd   = qcshd*ratio
          qccol   = qccol*ratio
          qcmul   = qcmul*ratio
       endif

! rain
       dum  = (qrevp+sum(qrcol)+sum(qrhetc)+sum(qrheti)+sum(qrmul))*dt
       dum1 = qr(i,k)+(qrcon+qcaut+qcacc+sum(qimlt)+sum(qcshd))*dt
       if (dum.gt.dum1 .and. dum.ge.1.e-20) then
          ratio  = dum1/dum
          qrevp   = qrevp*ratio
! NOTE: processes below are multiplied by ratio for each array element
          qrhetc = qrhetc*ratio
          qrheti = qrheti*ratio
          qrcol = qrcol*ratio
          qrmul = qrmul*ratio
       endif


! ice
       do iice = 1,n_iceCat

          dum = (qisub(iice)+qimlt(iice))*dt

! category interaction leading to sink for iice category
          do catcoll = 1,n_iceCat
             dum = dum+qicol(iice,catcoll)*dt
          enddo

          dum1 = qitot(i,k,iice)+(qidep(iice)+qinuc(iice)+qrcol(iice)+qccol(iice)+       &
                 qrhetc(iice)+qrheti(iice)+qchetc(iice)+qcheti(iice)+qcmul(iice)+        &
                 qrmul(iice))*dt

! category interaction leading to source for iice category
          do catcoll = 1,n_iceCat
             dum1 = dum1 + qicol(catcoll,iice)*dt
          enddo

          if (dum.gt.dum1 .and. dum.ge.1.e-20) then
             ratio       = dum1/dum
             qisub(iice) = qisub(iice)*ratio
             qimlt(iice) = qimlt(iice)*ratio
             do catcoll = 1,n_iceCat
                qicol(iice,catcoll) = qicol(iice,catcoll)*ratio
             enddo

          endif

       enddo !iice-loop

       iice_loop2: do iice = 1,n_iceCat

!........................................................................
! update microphysics, thermodynamic variables

! update qc and qr due to ice-phase processes (update for liquid processes occurs below)
          qc(i,k) = qc(i,k) + (-qccol(iice)-qcshd(iice)-qchetc(iice)-qcheti(iice)-       &
                    qcmul(iice))*dt
          qr(i,k) = qr(i,k) + (-qrcol(iice)+qimlt(iice)+qcshd(iice)-qrhetc(iice)-        &
                    qrheti(iice)-qrmul(iice))*dt

! update nc and nr due to ice-phase processes (update for liquid processes occurs below)
! NOTE: nc is currently specified, thus code below is commented out
! number concentration variables
!      nc(i,k) = nc(i,k)+(-nccol(iice))*dt

          if (iparam.eq.1 .or. iparam.eq.2) then
             nr(i,k) = nr(i,k) + (-nrcol(iice)-nimlt(iice)+nrshdr(iice)+ncshdc(iice)-    &
                       nrhetc(iice)-nrheti(iice))*dt
          else
             nr(i,k) = nr(i,k) + (-nrcol(iice)-nimlt(iice)+nrshdr(iice)+ncshdc(iice)-    &
                       nrhetc(iice)-nrheti(iice))*dt
          endif

! sink terms for ice (sublimation and melting)
          if (qitot(i,k,iice).ge.qsmall) then
! add sink terms, assume density stays constant for sink terms
             birim(i,k,iice) = birim(i,k,iice)- ((qisub(iice)+qimlt(iice))/              &
                               qitot(i,k,iice))*dt*birim(i,k,iice)
             qirim(i,k,iice) = qirim(i,k,iice)-((qisub(iice)+qimlt(iice))*               &
                               qirim(i,k,iice)/qitot(i,k,iice))*dt
             qitot(i,k,iice) = qitot(i,k,iice)-(qisub(iice)+qimlt(iice))*dt
          endif

! update total ice and rime ice mass mixing ratios (from sources)
          dum  = (qrcol(iice)+qccol(iice)+qrhetc(iice)+qrheti(iice)+qchetc(iice)+        &
                 qcheti(iice)+qcmul(iice)+qrmul(iice))*dt
          qitot(i,k,iice) = qitot(i,k,iice)+(qidep(iice)+qinuc(iice))*dt+dum
          qirim(i,k,iice) = qirim(i,k,iice)+dum

! update bulk volume mixing ratio
          birim(i,k,iice) = birim(i,k,iice) + (qrcol(iice)*inv_rho_rimeMax +             &
                            qccol(iice)/rhorime_c(iice) +(qrhetc(iice)+qrheti(iice)+     &
                            qchetc(iice)+qcheti(iice)+qcmul(iice)+qrmul(iice))*          &
                            inv_rho_rimeMax)*dt

          nitot(i,k,iice) = nitot(i,k,iice) + (ninuc(iice)+nimlt(iice)+nisub(iice)-      &
                            nislf(iice)+nimul(iice)+nrhetc(iice)+nrheti(iice)+           &
                            nchetc(iice)+ncheti(iice))*dt

! add ice-ice category interaction collection tendencies
! note: nicol is a sink for the collectee category, but NOT a source for collector
          do catcoll =1,n_iceCat

             qitot(i,k,catcoll) = qitot(i,k,catcoll)-qicol(catcoll,iice)*dt
             nitot(i,k,catcoll) = nitot(i,k,catcoll)-nicol(catcoll,iice)*dt
             qitot(i,k,iice)    = qitot(i,k,iice)+qicol(catcoll,iice)*dt

! now modify rime mass and density, assume collection does not modify rime mass
! fraction or density of the collectee, consistent with the assumption that
! these are constant over the PSD
             if (qitot(i,k,catcoll).ge.qsmall) then
! source for collector category
                qirim(i,k,iice) = qirim(i,k,iice)+qicol(catcoll,iice)*dt*                &
                                  qirim(i,k,catcoll)/qitot(i,k,catcoll)
                birim(i,k,iice) = birim(i,k,iice)+qicol(catcoll,iice)*dt*                &
                                  birim(i,k,catcoll)/qitot(i,k,catcoll)
! sink for collectee category
                qirim(i,k,catcoll) = qirim(i,k,catcoll)-qicol(catcoll,iice)*dt*          &
                                     qirim(i,k,catcoll)/qitot(i,k,catcoll)
                birim(i,k,catcoll) = birim(i,k,catcoll)-qicol(catcoll,iice)*dt*          &
                                     birim(i,k,catcoll)/qitot(i,k,catcoll)
             endif

          enddo ! catcoll loop

          if (qirim(i,k,iice) .lt. 0.) then
             qirim(i,k,iice)  = 0.
             birim(i,k,iice)  = 0.
          endif

! densify under wet growth
          if (wetgrowth(iice)) then   ! i.e., if wetgrowth is diagnosed as true
             qirim(i,k,iice) = qitot(i,k,iice)
             birim(i,k,iice) = qirim(i,k,iice)*inv_rho_rimeMax
          endif

! add tendencies for qv and theta for ice processes

          qv(i,k) = qv(i,k) + (-qidep(iice)+qisub(iice)-qinuc(iice))*dt
          th(i,k) = th(i,k) + th(i,k)/t(i,k)*((qidep(iice)-qisub(iice)+qinuc(iice))*     &
                              xxls(i,k)*inv_cp+(qrcol(iice)+qccol(iice)+qchetc(iice)+    &
                              qcheti(iice)+qrhetc(iice)+qrheti(iice)+qcmul(iice)+        &
                              qrmul(iice)-qimlt(iice))*xlf(i,k)*inv_cp)*dt

       enddo iice_loop2

! add tendencies for qc and qr for liquid processes
       qc(i,k) = qc(i,k) + (-qcacc-qcaut+qcnuc+qccon-qcevp)*dt
       qr(i,k) = qr(i,k) + (qcacc+qcaut+qrcon-qrevp)*dt

! add tendencies for qv and theta for liquid processes (note: tendency for ice-phase processes added above)
       qv(i,k) = qv(i,k) + (-qcnuc-qccon-qrcon+qcevp+qrevp)*dt
       th(i,k) = th(i,k) + th(i,k)/t(i,k)*((qcnuc+qccon+qrcon-qcevp-qrevp)*xxlv(i,k)*    &
                           inv_cp)*dt

! add tendencies for nc and nr for liquid processes
!      nc(i,k) = nc(i,k)+(-ncacc-ncautc+ncslf+ncnuc- &
!      nccol(iice)-nchetc-ncheti)*dt

! specify constant nc
       nc(i,k) = nccnst*inv_rho(i,k)
       if (iparam.eq.1 .or. iparam.eq.2) then
          nr(i,k) = nr(i,k) + (0.5*ncautc+nrslf+nrevp)*dt
       else
          nr(i,k) = nr(i,k) + (ncautr+nrslf+nrevp)*dt
       endif

!--- store processes for output  (for diagnostics only)
!      revap(i,k)  = qrevp
!      rcon(i,k)   = qrcon
!      ccon(i,k)   = p1c
!      auto(i,k)   = qcaut
!      acc(i,k)    = qcacc
!      act(i,k)    = ncnuc
!      opre(i,k)   = qrcon
!      opra(i,k)   = qcacc
!      oprc(i,k)   = qcaut
!      onpra(i,k)  = ncacc
!      onprc(i,k)  = ncautc
!      oncagg(i,k) = ncslf
!      onragg(i,k) = nrslf
!      onpccn(i,k) = ncnuc
!      pcc(i,k)    = qccon
!      opccn(i,k)  = qcnuc
!      oncautr(i,k) = ncautr
!===

       if (qc(i,k).lt.qsmall) then
          qv(i,k) = qv(i,k) + qc(i,k)
          th(i,k) = th(i,k) + th(i,k)/t(i,k)*qc(i,k)*xxlv(i,k)*inv_cp
          qc(i,k) = 0.
          nc(i,k) = 0.
       else
          log_hydrometeorsPresent = .true.
       endif

       if (qr(i,k).lt.qsmall) then
          qv(i,k) = qv(i,k) + qr(i,k)
          th(i,k) = th(i,k) + th(i,k)/t(i,k)*qr(i,k)*xxls(i,k)*inv_cp
          qr(i,k) = 0.
          nr(i,k) = 0.
       else
          log_hydrometeorsPresent = .true.
       endif

       do iice = 1,n_iceCat
          if (qitot(i,k,iice).lt.qsmall) then
             qv(i,k)   = qv(i,k) + qitot(i,k,iice)
             th(i,k)   = th(i,k) + th(i,k)/t(i,k)*qitot(i,k,iice)*xxls(i,k)*inv_cp
             qitot(i,k,iice)  = 0.
             nitot(i,k,iice)  = 0.
             qirim(i,k,iice) = 0.
             birim(i,k,iice)  = 0.
          else
             log_hydrometeorsPresent = .true.
          endif
       enddo !iice-loop

 555   continue

    enddo k_loop_main  !k-loop

!   if (ltrue.eq.0) goto 333
    if (.not. log_hydrometeorsPresent) goto 333

!------------------------------------------------------------------------------------------!

! initialize logicals for presence of hydrometeor species to .false.
    log_qcpresent = .false.

    do k = ktop,kbot,-kdir

! define 1/dzq
       inv_dzq(i,k) = 1./dzq(i,k)

!------------------------------------------------------------------------------------------!
! calculate number and mass weighted fallspeeds and find highest k level
! that a given species is present

!..................................................
! get droplet size distribution parameters
       if (qc(i,k).ge.qsmall) then

! set minimum nc to prevent floating point error
          nc(i,k)   = max(nc(i,k),1.e-16)

          pgam(i,k) = 0.0005714*(nc(i,k)*1.e-6*rho(i,k))+0.2714
!         pgam(i,k) = 0.146-5.964e-2*log(nc(i,k)/1.e6*rho(i,k)/2000.)
          pgam(i,k) = 1./(pgam(i,k)**2)-1.
          pgam(i,k) = max(pgam(i,k),2.)
          pgam(i,k) = min(pgam(i,k),15.)

! calculate lamc
          lamc(i,k) = (cons1*nc(i,k)*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.)/       &
                      qc(i,k))**thrd
          lammin = (pgam(i,k)+1.)*2.5e+4   ! min: 40 micron mean diameter
          lammax = (pgam(i,k)+1.)*1.e+6    ! max:  1 micron mean diameter

          if (lamc(i,k).lt.lammin) then
             lamc(i,k) = lammin
             nc(i,k)   = 6.*lamc(i,k)**3*qc(i,k)/(pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)* &
                         (pgam(i,k)+1.))
          elseif (lamc(i,k).gt.lammax) then
             lamc(i,k) = lammax
             nc(i,k)   = 6.*lamc(i,k)**3*qc(i,k)/(pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)* &
                         (pgam(i,k)+1.))
             if (.not. log_qcpresent) then
                qcindex = k
             endif
             log_qcpresent = .true.
          endif

       else
          lamc(i,k) = 0.
       endif

! droplet fall speed
! all droplets in smallest category fallspeed
! this means we can use analytic solution

       if (qc(i,k).ge.qsmall) then
          dum = 1./lamc(i,k)**bcn
!         vtrnc(i,k) =  acn(i,k)*gamma(1.+bcn+pgam(i,k))*dum/(gamma(pgam(i,k)+1.))
          vtrmc(i,k) = acn(i,k)*gamma(4.+bcn+pgam(i,k))*dum/(gamma(pgam(i,k)+4.))
       else
!         vtrnc(i,k) = 0.
          vtrmc(i,k) = 0.
       endif

    enddo ! k-loop

!------------------------------------------------------------------------------------------!
! cloud droplet sedimentation

    if (log_qcpresent) then

       nstep = 1

       do k = qcindex+kdir,kbot,-kdir

          fc(K)  = vtrmc(i,k)
!         fnc(K) = vtrnc(i,k)

          if (kdir.eq.1) then

             if (k.le.qcindex-kdir) then
                if (fc(k).lt.1.E-10) then
                   fc(k)=fc(k+kdir)
                endif
             endif

          elseif (kdir.eq.-1) then

             if (k.ge.qcindex-kdir) then
                if (fc(k).lt.1.e-10) then
                   fc(k) = fc(k+kdir)
                endif
             endif

          endif

! CALCULATE NUMBER OF SPLIT TIME STEPS
          rgvm   = fc(k)
          nstep  = max(int(rgvm*DT*inv_dzq(i,k)+1.),nstep)
          dumc(k)   = qc(i,k)*rho(i,k)
!         dumfnc(k) = nc(i,k)*rho(i,k)
          qcsten(K) = 0.
          ncsten(K) = 0.

       enddo ! k-loop

       inv_nstep = 1./real(nstep)
       if (nstep.ge.100) print*,'CLOUD nstep LARGE:',i,nstep

! calculate sedimentation using first-order upwind method
       tmp1 = 0.
       do n = 1,nstep

          do k = kbot,qcindex,kdir
             faloutc(k) = fc(k)*dumc(k)
!            faloutnc(k) = fnc(k)*dumfnc(k)
          enddo
          tmp1 = tmp1 + faloutc(kbot)  !sum flux at lowest level for averaging over sub-stepping

! top level with hydrometeor present

          k = qcindex
          faltndc   = faloutc(k)*inv_dzq(i,k)
!         faltndnc  = faloutnc(k)*inv_dzq(i,k)
          qcsten(k) = qcsten(k)-faltndc*inv_nstep*inv_rho(i,k)
!         ncsten(k) = ncsten(k)-faltndnc*inv_nstep*inv_rho(i,k)
          dumc(k)   = dumc(k)-faltndc*dt*inv_nstep
!         dumfnc(k) = dumfnc(k)-faltndnc*dt*inv_nstep

! loop from sceond to top level of hydrometeor to surface

          do k = qcindex-kdir,kbot,-kdir
             faltndc   = (faloutc(k+kdir)-faloutc(K))*inv_dzq(i,k)
!            faltndnc  = (faloutnc(k+kdir)-faloutnc(K))*inv_dzq(i,k)
             qcsten(k) = qcsten(k)+faltndc*inv_nstep*inv_rho(i,k)
!            ncsten(k) = ncsten(k)+faltndnc*inv_nstep*inv_rho(i,k)
             dumc(k)   = dumc(k)+faltndc*dt*inv_nstep
!            dumfnc(k) = dumfnc(k)+faltndnc*dt*inv_nstep
          enddo ! k loop

       enddo ! nstep-loop

! update prognostic variables with sedimentation tendencies

       do k = kbot,qcindex,kdir
          qc(i,k) = qc(i,k)+qcsten(k)*dt
!         nc(i,k) = nc(i,k)+ncsten(k)*dt
       enddo

! compute cloud contribution to liquid precipitation rate at surface
       tmp1 = tmp1*inv_nstep           !flux at surface, averaged over sub-step
       pcprt_liq(i) = tmp1*inv_rhow    !convert flux (kg m-2 s-1) to pcp rate (m s-1)

    endif ! log_qcpresent

!------------------------------------------------------------------------------------------!
    log_qrpresent = .false.

    do k = ktop,kbot,-kdir

!..................................................
! get rain size distribution parameters
       if (qr(i,k).ge.qsmall) then

! use lookup table to get mu
! mu-lambda relationship is from Cao et al. (2008), eq. (7)

! find spot in lookup table
! scaled N/q for lookup table parameter space
          nr(i,k) = max(nr(i,k),1.e-16)
          inv_dum = (qr(i,k)/(cons1*nr(i,k)*6.))**thrd

          if (inv_dum.lt.282.e-6) then
             mur(i,k) = 8.282
          elseif (inv_dum.ge.282.e-6 .and. inv_dum.lt.502.e-6) then
! interpolate
             rdumii   = (inv_dum-250.e-6)*1.e6*0.5
             rdumii   = max(rdumii,1.)
             rdumii   = min(rdumii,150.)
             dumii    = int(rdumii)
             dumii    = min(149,dumii)
             mur(i,k) = mur_table(dumii)+(mur_table(dumii+1)-mur_table(dumii))*(rdumii-  &
                        real(dumii))
          elseif (inv_dum.ge.502.e-6) then
             mur(i,k) = 0.
          endif

! recalculate slope based on mur
          lamr(i,k) = (cons1*nr(i,k)*(mur(i,k)+3.)*(mur(i,k)+2)*(mur(i,k)+1.)/(qr(i,k))) &
                      **thrd

! check for slope
          lammax = (mur(i,k)+1.)*1.e5

! set to small value since breakup is explicitly included (mean size 0.8 mm)
          lammin = (mur(i,k)+1.)*1250.

! apply lambda limiters for rain
          if (lamr(i,k).lt.lammin) then
             lamr(i,k) = lammin
             nr(i,k)   = exp(3.*log(lamr(i,k))+log(qr(i,k))+log(gamma(mur(i,k)+1.))-     &
                         log(gamma(mur(i,k)+4.)))/(cons1)
          elseif (lamr(i,k).gt.lammax) then
             lamr(i,k) = lammax
             nr(i,k)   = exp(3.*log(lamr(i,k))+log(qr(i,k))+log(gamma(mur(i,k)+1.))-     &
                         log(gamma(mur(i,k)+4.)))/(cons1)
          endif

          if (.not. log_qrpresent) then
             qrindex = k
          endif
          log_qrpresent = .true.

       else
          lamr(i,k) = 0.
       endif

       if (qr(i,k).ge.qsmall) then
!...................................................................
! read in fall mass and number weighted fall velocity from table

! find location in mean size space
          dum1 = (mur(i,k)+1.)/lamr(i,k)
          if (dum1.le.195.e-6) then
             inv_dum3  = 0.1
             rdumii = (dum1*1.e6+5.)*inv_dum3
             rdumii = max(rdumii,1.)
             rdumii = min(rdumii,20.)
             dumii  = int(rdumii)
             dumii  = max(dumii,1)
             dumii  = min(dumii,20)
          elseif (dum1.gt.195.e-6) then
             inv_dum3  = 0.0333333     ! i.e. 1/30
             rdumii = (dum1*1.e6-195.)*inv_dum3+20.
             rdumii = max(rdumii,20.)
             rdumii = min(rdumii,300.)
             dumii  = int(rdumii)
             dumii  = max(dumii,20)
             dumii  = min(dumii,299)
          endif

! find location in mur space
          rdumjj = mur(i,k)+1.
          rdumjj = max(rdumjj,1.)
          rdumjj = min(rdumjj,10.)
          dumjj  = int(rdumjj)
          dumjj  = max(dumjj,1)
          dumjj  = min(dumjj,9)

! number-weighted fallspeed
! value at mur
          dum1 = vn_table(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*                    &
                 (vn_table(dumii+1,dumjj)-vn_table(dumii+1,dumjj))

! value at mur+1
          dum2 = vn_table(dumii,dumjj+1)+(rdumii-real(dumii))*                           &
                 inv_dum3*(vn_table(dumii+1,dumjj+1)-vn_table(dumii+1,dumjj+1))

          vtrn(i,k) = dum1+(rdumjj-real(dumjj))*(dum2-dum1)

! mass-weighted fallspeed
! value at mur
          dum1 = vm_table(dumii,dumjj)+(rdumii-real(dumii))*inv_dum3*                    &
                 (vm_table(dumii+1,dumjj)-vm_table(dumii+1,dumjj))

! value at mur+1
          dum2 = vm_table(dumii,dumjj+1)+(rdumii-real(dumii))*inv_dum3*                  &
                 (vm_table(dumii+1,dumjj+1)-vm_table(dumii+1,dumjj+1))
          vtrm(i,k) = dum1+(rdumjj-real(dumjj))*(dum2-dum1)
          vtrn(i,k) = vtrn(i,k)*rhofacr(i,k)
          vtrm(i,k) = vtrm(i,k)*rhofacr(i,k)

       else
          vtrn(i,k) = 0.
          vtrm(i,k) = 0.
       endif

    enddo ! k-loop

!------------------------------------------------------------------------------------------!
! rain sedimentation

    if (log_qrpresent) then

       nstep = 1

       do k = qrindex+kdir,kbot,-kdir

          fr(k)  = vtrm(i,k)
          fnr(k) = vtrn(i,k)

          if (kdir.eq.1) then
             if (k.le.qrindex-kdir) then
                if (fr(k).lt.1.e-10) then
                   fr(k) = fr(k+kdir)
                endif
                if (fnr(k).lt.1.e-10) then
                   fnr(k) = fnr(k+kdir)
                endif
             endif

          elseif (kdir.eq.-1) then
             if (k.ge.qrindex-kdir) then
                if (fr(k).lt.1.e-10) then
                   fr(k) = fr(k+kdir)
                endif
                if (fnr(k).lt.1.e-10) then
                   fnr(k) = fnr(k+kdir)
                endif
             endif
          endif

! CALCULATE NUMBER OF SPLIT TIME STEPS
          rgvm      = max(fr(k),fnr(k))
          nstep     = max(int(rgvm*DT*inv_dzq(i,k)+1.),nstep)
          dumr(k)   = qr(i,k)*rho(i,k)
          dumfnr(k) = nr(i,k)*rho(i,k)
          qrsten(k) = 0.
          nrsten(k) = 0.

       enddo ! k loop

       inv_nstep = 1./real(nstep)
       if (nstep .ge. 100) print*,'RAIN nstep LARGE:',i,nstep

!--test:  explicitly calculate pcp rate:
!pcprt_liq(i) = qr(i,kbot)*rho(i,kbot)*vtrm(i,kbot)*1.e-3  !m s-1
!==

! calculate sedimentation using first-order upwind method
       tmp1 = 0.
       do n = 1,nstep

          do k = kbot,qrindex,kdir
             faloutr(k)  = fr(k)*dumr(k)
             faloutnr(k) = fnr(k)*dumfnr(k)
          enddo
          tmp1 = tmp1 + faloutr(kbot)  !sum flux at lowest level for averaging over sub-stepping

! top level with hydrometeor present
          k         = qrindex
          faltndr   = faloutr(k)*inv_dzq(i,k)
          faltndnr  = faloutnr(k)*inv_dzq(i,k)
          qrsten(k) = qrsten(k) - faltndr*inv_nstep*inv_rho(i,k)
          nrsten(k) = nrsten(k) - faltndnr*inv_nstep*inv_rho(i,k)
          dumr(k)   = dumr(k)   - faltndr*dt*inv_nstep
          dumfnr(k) = dumfnr(k) - faltndnr*dt*inv_nstep

! loop from second to top level of hydrometeor to surface
          do k = qrindex-kdir,kbot,-kdir
             faltndr   = (faloutr(k+kdir) - faloutr(K))*inv_dzq(i,k)
             faltndnr  = (faloutnr(k+kdir)-faloutnr(K))*inv_dzq(i,k)
             qrsten(k) = qrsten(k) + faltndr*inv_nstep*inv_rho(i,k)
             nrsten(k) = nrsten(k) + faltndnr*inv_nstep*inv_rho(i,k)
             dumr(k)   = dumr(k)   + faltndr*dt*inv_nstep
             dumfnr(k) = dumfnr(k) + faltndnr*dt*inv_nstep
          enddo ! k loop

!pcprt_liq(i) = pcprt_liq(i) + faloutr(kbot)*dt*inv_nstep   !orig

       enddo ! nstep loop

! update prognostic variables with sedimentation tendencies
       do k = kbot,qrindex,kdir
          qr(i,k) = qr(i,k) + qrsten(k)*dt
          nr(i,k) = nr(i,k) + nrsten(k)*dt
       enddo

! add rain component of liquid precipitation rate at surface
       tmp1 = tmp1*inv_nstep               !flux at surface, averaged over sub-step
       tmp1 = tmp1*inv_rhow                !convert flux (kg m-2 s-1) to pcp rate (m s-1)
       pcprt_liq(i) = pcprt_liq(i) + tmp1  !add pcp rate from cloud and rain

    endif ! log_qrpresent

!------------------------------------------------------------------------------------------!
! limit max total ice concentration for combined ice categories to 500 L-1
! if max is exceeded scale each category to preserve ratio of number between categories
    do k = ktop,kbot,-kdir
       if (sum(nitot(i,k,:)).ge.1.e-20) then
          dum = 500.e+3*inv_rho(i,k)/sum(nitot(i,k,:))
          nitot(i,k,:) = nitot(i,k,:)*min(dum,1.)
       endif
    enddo

    do iice = 1,n_iceCat

       log_qipresent = .false.  !note: this applies to ice category 'iice' only

       do k = ktop,kbot,-kdir

!.......................................................................
! get ice fallspeed for updated variables
          if (qitot(i,k,iice).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
             nitot(i,k,iice) = max(nitot(i,k,iice),1.e-16)

! calculate predicted bulk rime density
             if (birim(i,k,iice).ge.1.e-15) then
                rhop = qirim(i,k,iice)/birim(i,k,iice)
             else
                rhop = 0.
             endif
! limit 50 < rhop < 900, adjust bg if needed
             if (rhop.lt.rho_rimeMin) then
                rhop            = rho_rimeMin
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif
             if (rhop.gt.rho_rimeMax) then
                rhop            = rho_rimeMax
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif
             if (qirim(i,k,iice).lt.qsmall) then
                birim(i,k,iice) = 0.
             endif

! set upper constraint on qri to ensure qri cannot be > qi
             if (qirim(i,k,iice).gt.qitot(i,k,iice)) then
                qirim(i,k,iice) = qitot(i,k,iice)
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif

! find indices in 4D ice lookup table
!------------------------------------------------------------------------------------------!

! find index for qi (total ice mass mixing ratio)
             dum1 = (alog10(qitot(i,k,iice))+16.)*1.41328
             dumi = int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
             dum1 = min(dum1,real(isize))
             dum1 = max(dum1,1.)
             dumi = max(1,dumi)
             dumi = min(isize-1,dumi)

! find index for Ni (ice number mixing ratio)
             dum2 = (alog10(nitot(i,k,iice))+10.)*1.10731
             dumk = int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
             dum2 = min(dum2,real(jsize))
             dum2 = max(dum2,1.)
             dumk = max(1,dumk)
             dumk = min(jsize-1,dumk)

! find index for rime mass fraction
             dum4  = qirim(i,k,iice)/qitot(i,k,iice)*3.+1.
             dumii = int(dum4)

! set limits
             dum4  = min(dum4,real(rimsize))
             dum4  = max(dum4,1.)
             dumii = max(1,dumii)
             dumii = min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
             if (rhop.le.650.) then
                dum5 = (rhop-50.)*0.005 + 1.
             else
                dum5 = (rhop-650.)*0.004 + 4.
             endif
             dumjj = int(dum5)

! set limits
             dum5  = min(dum5,real(densize))
             dum5  = max(dum5,1.)
             dumjj = max(1,dumjj)
             dumjj = min(densize-1,dumjj)

! call subroutine to interpolate ice lookup table values
             call access_lookup_table(dumjj,dumii,dumi,dumk,1,dum1,dum2,dum4,dum5,f1pr1)
             call access_lookup_table(dumjj,dumii,dumi,dumk,2,dum1,dum2,dum4,dum5,f1pr2)
             call access_lookup_table(dumjj,dumii,dumi,dumk,7,dum1,dum2,dum4,dum5,f1pr9)
             call access_lookup_table(dumjj,dumii,dumi,dumk,8,dum1,dum2,dum4,dum5,f1pr10)

!...................................................
! make sure mean ice size is in bounds (i.e., apply lambda limiters)
             nitot(i,k,iice) = min(nitot(i,k,iice),f1pr9)
             nitot(i,k,iice) = max(nitot(i,k,iice),f1pr10)

             if (.not. log_qipresent) then
                qiindex = k
             endif
             log_qipresent = .true.

          endif ! qitot < qsmall

!........................................................................
          if (qitot(i,k,iice).ge.qsmall) then
             vtrnitot(i,k) = f1pr1*rhofaci(i,k)
             vtrmi1(i,k) = f1pr2*rhofaci(i,k)
! output fallspeed, w/o density correction
             diag_vmi(i,k,iice)    = f1pr2
          else
             vtrnitot(i,k) = 0.
             vtrmi1(i,k) = 0.
          endif

       enddo ! k-loop

!------------------------------------------------------------------------------------------!
! ice sedimentation
       if (log_qipresent) then

          nstep = 1

          do k = qiindex+kdir,kbot,-kdir

             fi(k)  = vtrmi1(i,k)
             fni(k) = vtrnitot(i,k)
             if (kdir.eq.1) then
                if (k.le.qiindex-kdir) then
                   if (fi(k).lt.1.e-10) then
                      fi(k) = fi(k+kdir)
                   endif
                   if (fni(k).lt.1.e-10) then
                      fni(k) = fni(k+kdir)
                   endif
                endif

             elseif (kdir.eq.-1) then

                if (k.ge.qiindex-kdir) then
                   if (fi(k).lt.1.e-10) then
                      fi(k) = fi(k+kdir)
                   endif
                   if (fni(k).lt.1.e-10) then
                      fni(k) = fni(k+kdir)
                   endif
                endif

             endif ! kdir

! calculate number of split time steps
             rgvm       = max(fi(k),fni(k))
             nstep      = max(int(rgvm*dt*inv_dzq(i,k)+1.),nstep)
             dumqi(k)   = qitot(i,k,iice)*rho(i,k)
             dumri(k)   = qirim(i,k,iice)*rho(i,k)
             dumbg(k)   = birim(i,k,iice)*rho(i,k)
             dumfni(k)  = nitot(i,k,iice)*rho(i,k)
             qisten(k)  = 0.
             qristen(k) = 0.
             bgsten(k)  = 0.
             nisten(k)  = 0.

          enddo ! k loop

          inv_nstep = 1./real(nstep)
          if (nstep.ge.100) print*,'ICE nstep LARGE:',i,nstep

! calculate sedimentation using first-order upwind method
          tmp1 = 0.
          do n = 1,nstep

             do k = kbot,qiindex+kdir,kdir
                falouti(k)  = fi(k)*dumqi(k)
                faloutni(k) = fni(k)*dumfni(k)
                faloutri(k) = fi(k)*dumri(k)
                faloutbg(k) = fi(k)*dumbg(k)
             enddo
            tmp1 = tmp1 + falouti(kbot)  !sum flux at lowest level for averaging over sub-stepping

! top level with hydrometeor present
             k          = qiindex
             faltndi    = falouti(k)*inv_dzq(i,k)
             faltndri   = faloutri(k)*inv_dzq(i,k)
             faltndbg   = faloutbg(k)*inv_dzq(i,k)
             faltndni   = faloutni(k)*inv_dzq(i,k)
             qisten(k)  = qisten(k)  - faltndi*inv_nstep*inv_rho(i,k)
             qristen(k) = qristen(k) -faltndri*inv_nstep*inv_rho(i,k)
             bgsten(k)  = bgsten(k)  - faltndbg*inv_nstep*inv_rho(i,k)
             nisten(k)  = nisten(k)  - faltndni*inv_nstep*inv_rho(i,k)
             dumqi(k)   = dumqi(k)   - faltndi*dt*inv_nstep
             dumri(k)   = dumri(k)   - faltndri*dt*inv_nstep
             dumbg(k)   = dumbg(k)   - faltndbg*dt*inv_nstep
             dumfni(k)  = dumfni(k)  - faltndni*dt*inv_nstep

! loop from sceond to top level of hydrometeor to surface
             do k = qiindex-kdir,kbot,-kdir
                faltndi    = (falouti(k+kdir)  - falouti(k) )*inv_dzq(i,k)
                faltndri   = (faloutri(k+kdir) - faloutri(k))*inv_dzq(i,k)
                faltndbg   = (faloutbg(k+kdir) - faloutbg(k))*inv_dzq(i,k)
                faltndni   = (faloutni(k+kdir) - faloutni(k))*inv_dzq(i,k)
                qisten(k)  = qisten(k)  + faltndi*inv_nstep*inv_rho(i,k)
                qristen(k) = qristen(k) + faltndri*inv_nstep*inv_rho(i,k)
                bgsten(k)  = bgsten(k)  + faltndbg*inv_nstep*inv_rho(i,k)
                nisten(k)  = nisten(k)  + faltndni*inv_nstep*inv_rho(i,k)
                dumqi(k)   = dumqi(k)   + faltndi*dt*inv_nstep
                dumri(k)   = dumri(k)   + faltndri*dt*inv_nstep
                dumbg(k)   = dumbg(k)   + faltndbg*dt*inv_nstep
                dumfni(k)  = dumfni(k)  + faltndni*dt*inv_nstep
             enddo ! k loop

          enddo ! nstep loop

! update prognostic variables with sedimentation tendencies
          do k = kbot,qiindex,kdir
             qitot(i,k,iice) = qitot(i,k,iice) + qisten(k)*dt
             qirim(i,k,iice) = qirim(i,k,iice) + qristen(k)*dt
             birim(i,k,iice) = birim(i,k,iice) + bgsten(k)*dt
             nitot(i,k,iice) = nitot(i,k,iice) + nisten(k)*dt
          enddo

! add contirubtion from iice to solid precipitation rate at surface
       tmp1 = tmp1*inv_nstep   !flux at surface, averaged over sub-step
       tmp1 = tmp1*inv_rhow    !convert flux (kg m-2 s-1) to pcp rate (m s-1), liquid-equivalent
       pcprt_sol(i) = pcprt_sol(i) + tmp1  !add pcp rate from cloud and rain

       endif ! log_qipresent

    enddo  !iice-loop

!------------------------------------------------------------------------------------------!
! end sedimentation
!------------------------------------------------------------------------------------------!

!.............................
! homogeneous freezing of cloud droplets and rain

    do k = kbot,ktop,kdir

       diam_ice(i,k,:) = 0.

       do iice=1,n_iceCat

          if (qitot(i,k,iice).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
             nitot(i,k,iice) = max(nitot(i,k,iice),1.e-16)

             !--compute mean-mass ice diameters (estimated; rigorous approach to be implemented later)
             dum2 = 500. !ice density
             diam_ice(i,k,iice) = ((qitot(i,k,iice)*6.)/(nitot(i,k,iice)*dum2*pi))**thrd
             !==
       endif

       enddo  !iice-loop

       if (qc(i,k).ge.qsmall .and. t(i,k).lt.233.15) then
          Q_nuc = qc(i,k)
          N_nuc = max(nc(i,k),1.e-16)
         !--determine destination ice-phase category:
          dum1      = 900.     !density of new ice
          D_new     = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
          call icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,      &
                                  ni_add,iice_dest)
         !==
          qirim(i,k,iice_dest) = qirim(i,k,iice_dest) + Q_nuc
          qitot(i,k,iice_dest)  = qitot(i,k,iice_dest)  + Q_nuc
          birim(i,k,iice_dest)  = birim(i,k,iice_dest)  + Q_nuc*inv_rho_rimeMax
          if (ni_add) then
          nitot(i,k,iice_dest)  = nitot(i,k,iice_dest)  + N_nuc
          end if
          th(i,k)   = th(i,k)   + th(i,k)/t(i,k)*Q_nuc*xlf(i,k)*inv_cp
          qc(i,k)   = 0.  != qc(i,k) - Q_nuc
          nc(i,k)   = 0.  != nc(i,k) - N_nuc
       endif

       if (qr(i,k).ge.qsmall .and. t(i,k).lt.233.15) then
          Q_nuc = qr(i,k)
          N_nuc = nr(i,k)
          N_nuc = max(nr(i,k),1.e-16)
         !--determine destination ice-phase category:
          dum1      = 900.     !density of new ice
          D_new     = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
          call icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,      &
                                  ni_add,iice_dest)
         !==
          qirim(i,k,iice_dest) = qirim(i,k,iice_dest) + Q_nuc
          qitot(i,k,iice_dest) = qitot(i,k,iice_dest) + Q_nuc
          birim(i,k,iice_dest) = birim(i,k,iice_dest) + Q_nuc*inv_rho_rimeMax
          if (ni_add) nitot(i,k,iice_dest) = nitot(i,k,iice_dest) + N_nuc
          th(i,k)   = th(i,k)   + th(i,k)/t(i,k)*Q_nuc*xlf(i,k)*inv_cp
          qr(i,k)   = 0.  ! = qr(i,k) - Q_nuc
          nr(i,k)   = 0.  ! = nr(i,k) - N_nuc
       endif

    enddo ! k loop

!.......................................................................
! calculate final size distribution parameters/effective radii
! also final checks to ensure consistency of mass/number

! get size distribution parameters

    do k = kbot,ktop,kdir

       ze_ice  = z_min
       ze_rain = z_min

       if (qc(i,k).lt.qsmall) then
          qv(i,k) = qv(i,k) + qc(i,k)
          th(i,k) = th(i,k) + th(i,k)/t(i,k)*qc(i,k)*xxlv(i,k)*inv_cp
          qc(i,k) = 0.
          nc(i,k) = 0.
       endif

       if (qr(i,k).lt.qsmall) then
          qv(i,k) = qv(i,k) + qr(i,k)
          th(i,k) = th(i,k) + th(i,k)/t(i,k)*qr(i,k)*xxls(i,k)*inv_cp
          qr(i,k) = 0.
          nr(i,k) = 0.
       endif

!..................................................
! get droplet size distribution parameters
       if (qc(i,k).ge.qsmall) then

! set minimum nc to prevent floating point error
          nc(i,k)   = max(nc(i,k),1.e-16)
          pgam(i,k) = 0.0005714*(nc(i,k)*1.e-6*rho(i,k))+0.2714
!         pgam(i,k) = 0.146-5.964e-2*log(nc(i,k)/1.e6*rho(i,k)/2000.)
          pgam(i,k) = 1./(pgam(i,k)**2)-1.
          pgam(i,k) = max(pgam(i,k),2.)
          pgam(i,k) = min(pgam(i,k),15.)

! calculate lamc
          lamc(i,k) = (cons1*nc(i,k)*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.)/       &
                      qc(i,k))**thrd
          lammin    = (pgam(i,k)+1.)*2.5e+4   ! lammin, 40 micron mean diameter
          lammax    = (pgam(i,k)+1.)*1.e+6    ! lammax,  1 micron mean diameter
          if (lamc(i,k).lt.lammin) then
             lamc(i,k) = lammin
             nc(i,k)   = 6.*lamc(i,k)**3*qc(i,k)/(pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)* &
                         (pgam(i,k)+1.))
          elseif (lamc(i,k).gt.lammax) then
             lamc(i,k) = lammax
             nc(i,k)   = 6.*lamc(i,k)**3*qc(i,k)/(pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)* &
                         (pgam(i,k)+1.))
          endif

       else
          lamc(i,k) = 0.
       endif

!..................................................
! get rain size distribution parameters
       if (qr(i,k).ge.qsmall) then

! use lookup table to get mu
! mu-lambda relationship is from Cao et al. (2008), eq. (7)

! find spot in lookup table
! scaled N/q for lookup table parameter space
          nr(i,k) = max(nr(i,k),1.e-16)
          dum     = (cons1*NR(i,k)*6./(QR(i,k)))**thrd
          inv_dum    = 1./dum

          if (inv_dum.lt.282.e-6) then
             mur(i,k) = 8.282
          elseif (inv_dum.ge.282.e-6 .and. inv_dum.lt.502.e-6) then
! interpolate
             rdumii   = (inv_dum-250.e-6)*0.5e+6
             rdumii   = max(rdumii,1.)
             rdumii   = min(rdumii,150.)
             dumii    = int(rdumii)
             dumii    = min(149,dumii)
             mur(i,k) = mur_table(dumii)+(mur_table(dumii+1)-mur_table(dumii))*(rdumii-  &
                        real(dumii))
          elseif (inv_dum.ge.502.e-6) then
             mur(i,k) = 0.
          endif

! recalculate slope based on mur
          lamr(i,k) = (cons1*nr(i,k)*(mur(i,k)+3.)*(mur(i,k)+2)*(mur(i,k)+1.)/           &
                      qr(i,k))**thrd

! check for slope
          lammax = (mur(i,k)+1.)*1.e5

! set to small value since breakup is explicitly included (mean size 0.8 mm)
          lammin = (mur(i,k)+1.)*1250.

! adjust vars
          if (lamr(i,k).lt.lammin) then
             lamr(i,k) = lammin
             NR(i,k)   = exp(3.*log(lamr(i,k))+log(QR(i,k))+log(gamma(mur(i,k)+1.))-     &
                         log(gamma(mur(i,k)+4.)))/(cons1)
          elseif (lamr(i,k).gt.lammax) then
             lamr(i,k) = lammax
             nr(i,k)   = exp(3.*log(lamr(i,k))+log(QR(i,k))+log(gamma(mur(i,k)+1.))-     &
                         log(gamma(mur(i,k)+4.)))/(cons1)
          endif

       else
          lamr(i,k) = 0.
       endif

! calculate effective radius for cloud droplets and rain
       if (qr(i,k).ge.qsmall) then
          effr(i,k) = 1.5/lamr(i,k)*1.e6
       else
          effr(i,k) = 25.
       endif

! use default of 10 microns
!       if (qc(i,k).ge.qsmall) then
!          diag_effc(i,k) = (pgam(i,k)+3.)/lamc(i,k)*0.5
!       else
!          diag_effc(i,k) = 25.e-6
!       endif
!.......................................................................
! size distribution parameters for ice

! limit max total ice concentration for combined ice categories to 500 L-1
! if max is exceeded scale each category to preserve ratio of number between categories

       if (sum(nitot(i,k,:)).ge.1.e-20) then
          dum = 500.e+3*inv_rho(i,k)/sum(nitot(i,k,:))
          nitot(i,k,:) = nitot(i,k,:)*min(dum,1.)
       endif

       do iice = 1,n_iceCat

! calculate total mixing ratio from deposition and riming mixing ratios
          if (qitot(i,k,iice).lt.qsmall) then
             qv(i,k)   = qv(i,k) + qitot(i,k,iice)
             th(i,k)   = th(i,k) + th(i,k)/t(i,k)*qitot(i,k,iice)*xxls(i,k)*inv_cp
             qitot(i,k,iice) = 0.
             nitot(i,k,iice) = 0.
             qirim(i,k,iice) = 0.
             birim(i,k,iice) = 0.
          endif

          if (qitot(i,k,iice).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
             nitot(i,k,iice) = max(nitot(i,k,iice),1.e-16)

! calculate predicted bulk rime density
             if (birim(i,k,iice).ge.1.e-15) then
                rhop = qirim(i,k,iice)/birim(i,k,iice)
             else
                rhop = 0.
             endif

! limit 50 < rhop < 900, adjust bg if needed
             if (rhop.lt.rho_rimeMin) then
                rhop = rho_rimeMin
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif
             if (rhop.gt.rho_rimeMax) then
                rhop = rho_rimeMax
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif
             if (qirim(i,k,iice).lt.qsmall) then
                birim(i,k,iice) = 0.
             endif

! set upper constraint on qri to ensure qri cannot be > qi
             if (qirim(i,k,iice).gt.qitot(i,k,iice)) then
                qirim(i,k,iice) = qitot(i,k,iice)
                birim(i,k,iice) = qirim(i,k,iice)/rhop
             endif

! find indices in 4D ice lookup table
!------------------------------------------------------------------------------------------!

! find index for qi (total ice mass mixing ratio)
             dum1  = (alog10(qitot(i,k,iice))+16.)*1.41328
             dumi  = int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
             dum1  = min(dum1,real(isize))
             dum1  = max(dum1,1.)
             dumi  = max(1,dumi)
             dumi  = min(isize-1,dumi)

! find index for Ni (ice number mixing ratio)
             dum2  = (alog10(nitot(i,k,iice))+10.)*1.10731
             dumk  = int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
             dum2  = min(dum2,real(jsize))
             dum2  = max(dum2,1.)
             dumk  = max(1,dumk)
             dumk  = min(jsize-1,dumk)

! find index for rime mass fraction
             dum4  = qirim(i,k,iice)/qitot(i,k,iice)*3.+1.
             dumii = int(dum4)

! set limits
             dum4  = min(dum4,real(rimsize))
             dum4  = max(dum4,1.)
             dumii = max(1,dumii)
             dumii = min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
             if (rhop.le.650.) then
                dum5  = (rhop-50.)*0.005 + 1.
             else
                dum5  = (rhop-650.)*0.004 + 4.
             endif
             dumjj = int(dum5)

! set limits
             dum5  = min(dum5,real(densize))
             dum5  = max(dum5,1.)
             dumjj = max(1,dumjj)
             dumjj = min(densize-1,dumjj)

! call subroutine to interpolate lookup table values
             call access_lookup_table(dumjj,dumii,dumi,dumk,6,dum1,dum2,dum4,dum5,f1pr6)
             call access_lookup_table(dumjj,dumii,dumi,dumk,7,dum1,dum2,dum4,dum5,f1pr9)
             call access_lookup_table(dumjj,dumii,dumi,dumk,8,dum1,dum2,dum4,dum5,f1pr10)
             call access_lookup_table(dumjj,dumii,dumi,dumk,9,dum1,dum2,dum4,dum5,f1pr13)
             call access_lookup_table(dumjj,dumii,dumi,dumk,11,dum1,dum2,dum4,dum5,f1pr15)
             call access_lookup_table(dumjj,dumii,dumi,dumk,12,dum1,dum2,dum4,dum5,f1pr16)

!......................................................
! make sure number concentration is within bounds
             nitot(i,k,iice) = min(nitot(i,k,iice),f1pr9)
             nitot(i,k,iice) = max(nitot(i,k,iice),f1pr10)

! make sure qri is > qsmall
             if (qirim(i,k,iice).lt.qsmall) then
                qirim(i,k,iice) = 0.
                birim(i,k,iice) = 0.
             endif

             diag_effi(i,k,iice) = f1pr6 ! units are in m

          else
             diag_effi(i,k,iice) = 25.e-6
          endif ! qitot < qsmall

          if (qitot(i,k,iice).ge.qsmall) then
! mean ice size for output
             diag_di(i,k,iice)    = f1pr15
             diag_rhopo(i,k,iice) = f1pr16
             ze_ice = ze_ice + 0.1892*f1pr13   ! sum contribution from each ice category
             ze_ice = max(ze_ice,z_min)
          endif

       enddo !iice-loop

       if (qr(i,k).ge.qsmall) then
!         ze_rain(i,k) = n0r(i,k)*720./lamr(i,k)**3/lamr(i,k)**3/lamr(i,k)
! non-exponential rain
          ze_rain = nr(i,k)*(mur(i,k)+6.)*(mur(i,k)+5.)*(mur(i,k)+4.)*                   &
                    (mur(i,k)+3.)*(mur(i,k)+2.)*(mur(i,k)+1.)/lamr(i,k)**6
          ze_rain = max(ze_rain,z_min)
       endif

!        if (qc(i,k).ge.qsmall) then
!           ze_cloud = 720.*nc(i,k)/lamc(i,k)**6
!           ze_cloud = max(ze_cloud,z_min)
!        endif

! convert to dbz
!      diag_ze(i,k) = 10.*log10((ze_cloud + ze_rain + ze_ice)*1.d18)
       diag_ze(i,k) = 10.*log10((ze_rain + ze_ice)*1.d18)
       diag_ze(i,k) = max(diag_ze(i,k),zdBZ_min)

!....................................................................
! merge ice categories if properties are similar
       do iice = n_iceCat,2,-1

          if (abs(diag_di(i,k,iice)-diag_di(i,k,iice-1)).le.150.e-6   .and.              &
              abs(diag_rhopo(i,k,iice)-diag_rhopo(i,k,iice-1)).le.100.) then

            !merge:
             qitot(i,k,iice-1) = qitot(i,k,iice-1) + qitot(i,k,iice)
             nitot(i,k,iice-1) = nitot(i,k,iice-1) + nitot(i,k,iice)
             qirim(i,k,iice-1) = qirim(i,k,iice-1) + qirim(i,k,iice)
             birim(i,k,iice-1) = birim(i,k,iice-1) + birim(i,k,iice)

             qitot(i,k,iice) = 0.
             nitot(i,k,iice) = 0.
             qirim(i,k,iice) = 0.
             birim(i,k,iice) = 0.

          end if
       end do
!....................................................................

    enddo ! k-loop

    diag_zec(i) = maxval(diag_ze(i,:))

333 continue

! recalculate supersaturation from T and qv
! calculate temperature from theta

    if (log_predictsSsat) then
       do k = kbot,ktop,kdir
          t(i,k)    = th(i,k)*(1.e-5*pres(i,k))**(rd*inv_cp)
          dum       = 0.622*polysvp1(t(i,k),0)/(pres(i,k)-polysvp1(t(i,k),0))
          ssat(i,k) = qv(i,k)-dum
       enddo
    endif

!.....................................................

 enddo i_loop_main

! end main microphysics routine

 return

 END SUBROUTINE p3_main

!==========================================================================================!

 SUBROUTINE access_lookup_table(dumjj,dumii,dumi,dumk,index,dum1,dum2,dum4,dum5,proc)

 implicit none

 real    :: dum1,dum2,dum4,dum5,proc,dproc1,dproc2,iproc1,gproc1,tmp1,tmp2
 integer :: dumjj,dumii,dumi,dumk,index

! get value at current density index

! first interpolate for current rimed fraction index
   dproc1 = itab(dumjj,dumii,dumi,dumk,index)+(dum1-real(dumi))*(itab(dumjj,dumii,       &
            dumi+1,dumk,index)-itab(dumjj,dumii,dumi,dumk,index))

   dproc2 = itab(dumjj,dumii,dumi,dumk+1,index)+(dum1-real(dumi))*(itab(dumjj,dumii,     &
          dumi+1,dumk+1,index)-itab(dumjj,dumii,dumi,dumk+1,index))

   iproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! linearly interpolate to get process rates for rimed fraction index + 1
   dproc1 = itab(dumjj,dumii+1,dumi,dumk,index)+(dum1-real(dumi))*(itab(dumjj,dumii+1,   &
          dumi+1,dumk,index)-itab(dumjj,dumii+1,dumi,dumk,index))

   dproc2 = itab(dumjj,dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))*(itab(dumjj,dumii+1, &
          dumi+1,dumk+1,index)-itab(dumjj,dumii+1,dumi,dumk+1,index))

   gproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)
   tmp1   = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! get value at density index + 1

! first interpolate for current rimed fraction index
   dproc1 = itab(dumjj+1,dumii,dumi,dumk,index)+(dum1-real(dumi))*(itab(dumjj+1,dumii,   &
            dumi+1,dumk,index)-itab(dumjj+1,dumii,dumi,dumk,index))

   dproc2 = itab(dumjj+1,dumii,dumi,dumk+1,index)+(dum1-real(dumi))*(itab(dumjj+1,dumii, &
            dumi+1,dumk+1,index)-itab(dumjj+1,dumii,dumi,dumk+1,index))

   iproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! linearly interpolate to get process rates for rimed fraction index + 1
   dproc1 = itab(dumjj+1,dumii+1,dumi,dumk,index)+(dum1-real(dumi))*(itab(dumjj+1,       &
            dumii+1,dumi+1,dumk,index)-itab(dumjj+1,dumii+1,dumi,dumk,index))

   dproc2 = itab(dumjj+1,dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))*(itab(dumjj+1,     &
            dumii+1,dumi+1,dumk+1,index)-itab(dumjj+1,dumii+1,dumi,dumk+1,index))

   gproc1 = dproc1+(dum2-real(dumk))*(dproc2-dproc1)
   tmp2   = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! get final process rate
   proc   = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

END SUBROUTINE access_lookup_table

!------------------------------------------------------------------------------------------!
SUBROUTINE access_lookup_table_coll(dumjj,dumii,dumj,dumi,dumk,index,dum1,dum2,dum3,     &
                                    dum4,dum5,proc)

 implicit none

 real    :: dum1,dum2,dum3,dum4,dum5,proc,dproc1,dproc2,iproc1,gproc1,tmp1,tmp2,dproc11, &
            dproc12,dproc21,dproc22
 integer :: dumjj,dumii,dumj,dumi,dumk,index


! This subroutine interpolates lookup table values for rain/ice collection processes

! current density index

! current rime fraction index
   dproc11 = itabcoll(dumjj,dumii,dumi,dumk,dumj,index)+(dum1-real(dumi))*               &
             (itabcoll(dumjj,dumii,dumi+1,dumk,dumj,index)-itabcoll(dumjj,dumii,dumi,    &
             dumk,dumj,index))

   dproc21 = itabcoll(dumjj,dumii,dumi,dumk+1,dumj,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj,dumii,dumi+1,dumk+1,dumj,index)-itabcoll(dumjj,dumii,dumi,  &
             dumk+1,dumj,index))

   dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)

   dproc12 = itabcoll(dumjj,dumii,dumi,dumk,dumj+1,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj,dumii,dumi+1,dumk,dumj+1,index)-itabcoll(dumjj,dumii,dumi,  &
             dumk,dumj+1,index))

   dproc22 = itabcoll(dumjj,dumii,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj,dumii,dumi+1,dumk+1,dumj+1,index)-itabcoll(dumjj,dumii,     &
             dumi,dumk+1,dumj+1,index))

   dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)
   iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

! rime fraction index + 1

   dproc11 = itabcoll(dumjj,dumii+1,dumi,dumk,dumj,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj,dumii+1,dumi+1,dumk,dumj,index)-itabcoll(dumjj,dumii+1,     &
                 dumi,dumk,dumj,index))

   dproc21 = itabcoll(dumjj,dumii+1,dumi,dumk+1,dumj,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj,dumii+1,dumi+1,dumk+1,dumj,index)-itabcoll(dumjj,dumii+1,   &
             dumi,dumk+1,dumj,index))

   dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)

   dproc12 = itabcoll(dumjj,dumii+1,dumi,dumk,dumj+1,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj,dumii+1,dumi+1,dumk,dumj+1,index)-itabcoll(dumjj,dumii+1,   &
             dumi,dumk,dumj+1,index))

   dproc22 = itabcoll(dumjj,dumii+1,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*         &
             (itabcoll(dumjj,dumii+1,dumi+1,dumk+1,dumj+1,index)-itabcoll(dumjj,dumii+1, &
             dumi,dumk+1,dumj+1,index))

   dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)

   gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
   tmp1    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! density index + 1

! current rime fraction index
   dproc11 = itabcoll(dumjj+1,dumii,dumi,dumk,dumj,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj+1,dumii,dumi+1,dumk,dumj,index)-itabcoll(dumjj+1,dumii,     &
                 dumi,dumk,dumj,index))

   dproc21 = itabcoll(dumjj+1,dumii,dumi,dumk+1,dumj,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj+1,dumii,dumi+1,dumk+1,dumj,index)-itabcoll(dumjj+1,dumii,   &
             dumi,dumk+1,dumj,index))

   dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)

   dproc12 = itabcoll(dumjj+1,dumii,dumi,dumk,dumj+1,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj+1,dumii,dumi+1,dumk,dumj+1,index)-itabcoll(dumjj+1,dumii,   &
             dumi,dumk,dumj+1,index))

   dproc22 = itabcoll(dumjj+1,dumii,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*         &
                 (itabcoll(dumjj+1,dumii,dumi+1,dumk+1,dumj+1,index)-itabcoll(dumjj+1,   &
                 dumii,dumi,dumk+1,dumj+1,index))

   dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)
   iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

! rime fraction index + 1

   dproc11 = itabcoll(dumjj+1,dumii+1,dumi,dumk,dumj,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj+1,dumii+1,dumi+1,dumk,dumj,index)-itabcoll(dumjj+1,dumii+1, &
             dumi,dumk,dumj,index))

   dproc21 = itabcoll(dumjj+1,dumii+1,dumi,dumk+1,dumj,index)+(dum1-real(dumi))*         &
             (itabcoll(dumjj+1,dumii+1,dumi+1,dumk+1,dumj,index)-itabcoll(dumjj+1,       &
                 dumii+1,dumi,dumk+1,dumj,index))

   dproc1  = dproc11+(dum2-real(dumk))*(dproc21-dproc11)

   dproc12 = itabcoll(dumjj+1,dumii+1,dumi,dumk,dumj+1,index)+(dum1-real(dumi))*         &
             (itabcoll(dumjj+1,dumii+1,dumi+1,dumk,dumj+1,index)-itabcoll(dumjj+1,       &
                 dumii+1,dumi,dumk,dumj+1,index))

   dproc22 = itabcoll(dumjj+1,dumii+1,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))*       &
             (itabcoll(dumjj+1,dumii+1,dumi+1,dumk+1,dumj+1,index)-itabcoll(dumjj+1,     &
                 dumii+1,dumi,dumk+1,dumj+1,index))

   dproc2  = dproc12+(dum2-real(dumk))*(dproc22-dproc12)

   gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
   tmp2    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! interpolate over density to get final values
   proc    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

 END SUBROUTINE access_lookup_table_coll

!------------------------------------------------------------------------------------------!

 SUBROUTINE access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj,dumi,dumk,    &
                                     index,dum1c,dum4c,dum5c,dum1,dum2,dum4,dum5,proc)

 implicit none

 real    :: dum1,dum2,dum3,dum4,dum5,dum1c,dum4c,dum5c,proc,dproc1,dproc2,iproc1,iproc2, &
            gproc1,gproc2,rproc1,rproc2,tmp1,tmp2,dproc11,dproc12
 integer :: dumjj,dumii,dumj,dumi,dumk,index,dumjjc,dumiic,dumic


! This subroutine interpolates lookup table values for rain/ice collection processes

! current density index collectee category

! current rime fraction index for collectee category

! current density index collector category

! current rime fraction index for collector category

  if (index.eq.1) then

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*    &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk,dumii,dumjj)-                     &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj)-                   &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj)-                   &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj)-                 &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj)-                   &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj)-                 &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj)-                 &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)-               &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc1  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!............................................................................................................
!............................................................................................................
! collectee density index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj)-                      &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc2  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!..........................................................................................
! final process rate interpolation over collectee density

   proc    = rproc1+(dum5c-real(dumjjc))*(rproc2-rproc1)

 else if (index.eq.2) then

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk,dumii,dumjj)-                      &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc1  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!............................................................................................................
!............................................................................................................
! collectee density index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj)-                      &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj+1)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc1  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*     &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj+1)-                      &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk,dumii+1,dumjj+1))

   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))*   &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1)-                    &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumk+1,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumk+1,dumii+1,dumjj+1))

   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc2  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!..........................................................................................
! final process rate interpolation over collectee density

   proc    = rproc1+(dum5c-real(dumjjc))*(rproc2-rproc1)

 endif ! index =1 or 2

 END SUBROUTINE access_lookup_table_colli

!==========================================================================================!

 real function polysvp1(T,TYPE)

!-------------------------------------------
!  COMPUTE SATURATION VAPOR PRESSURE
!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!-------------------------------------------

      implicit none

      real    :: DUM,T
      integer :: TYPE

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /           &
        6.11147274, 0.503160820, 0.188439774e-1,           &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7,     &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                    &
        6.11239921, 0.443987641, 0.142986287e-1,           &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7,    &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

!-------------------------------------------

      if (TYPE.EQ.1 .and. T.lt.273.15) then
! ICE

!       Flatau formulation:
         dt       = max(-80.,t-273.16)
         polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+       &
                    a8i*dt)))))))
         polysvp1 = polysvp1*100.

!       Goff-Gratch formulation:
!        POLYSVP1 = 10.**(-9.09718*(273.16/T-1.)-3.56654*                 &
!          log10(273.16/T)+0.876793*(1.-T/273.16)+                        &
!          log10(6.1071))*100.


      elseif (TYPE.EQ.0 .or. T.ge.273.15) then
! LIQUID

!       Flatau formulation:
         dt       = max(-80.,t-273.16)
         polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
         polysvp1 = polysvp1*100.

!       Goff-Gratch formulation:
!        POLYSVP1 = 10.**(-7.90298*(373.16/T-1.)+                         &
!             5.02808*log10(373.16/T)-                                    &
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+                  &
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+                &
!             log10(1013.246))*100.

         endif


 end function polysvp1

!------------------------------------------------------------------------------------------!

 real function gamma(X)
!----------------------------------------------------------------------
! THIS ROUTINE CALCULATES THE gamma FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE gamma
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!----------------------------------------------------------------------
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH gamma(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  gamma(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!----------------------------------------------------------------------
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, log, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      integer :: I,N
      logical :: PARITY
      real ::                                                       &
          CONV,EPS,FACT,HALF,ONE,res,sum,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      real, dimension(7) :: C
      real, dimension(8) :: P
      real, dimension(8) :: Q
      real, parameter    :: constant1 = 0.9189385332046727417803297

!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      data ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      data XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      data P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      data Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      data C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,                 &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,    &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      if (Y.LE.ZERO) then
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        res=Y-Y1
        if (res.NE.ZERO) then
          if(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*res)
          Y=Y+ONE
        else
          res=XINF
          goto 900
        endif
      endif
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      if (Y.LT.EPS) then
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        if (Y.GE.XMININ) then
          res=ONE/Y
        else
          res=XINF
          goto 900
        endif
      elseif (Y.LT.TWELVE) then
        Y1=Y
        if (Y.LT.ONE) then
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        else
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        endif
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        do I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        enddo
        res=XNUM/XDEN+ONE
        if (Y1.LT.Y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          res=res/Y1
        elseif (Y1.GT.Y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          do I=1,N
            res=res*Y
            Y=Y+ONE
          enddo
        endif
      else
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        if (Y.LE.XBIG) then
          YSQ=Y*Y
          sum=C(7)
          do I=1,6
            sum=sum/YSQ+C(I)
          enddo
          sum=sum/Y-Y+constant1
          sum=sum+(Y-HALF)*log(Y)
          res=exp(sum)
        else
          res=XINF
          goto 900
        endif
      endif
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      if (PARITY)res=-res
      if (FACT.NE.ONE)res=FACT/res
  900 gamma=res
      return
! ---------- LAST LINE OF gamma ----------

 end function gamma

!------------------------------------------------------------------------------------------!

 real function DERF(X)

 implicit none

 real :: X
 real, dimension(0 : 64) :: A, B
 real :: W,T,Y
 integer :: K,I
      data A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,                            &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,                            &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0,                            &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,                            &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      data (B(I), I = 0, 12) /                                 &
         -0.00000000029734388465E0,  0.00000000269776334046E0, &
         -0.00000000640788827665E0, -0.00000001667820132100E0, &
         -0.00000021854388148686E0,  0.00000266246030457984E0, &
          0.00001612722157047886E0, -0.00025616361025506629E0, &
          0.00015380842432375365E0,  0.00815533022524927908E0, &
         -0.01402283663896319337E0, -0.19746892495383021487E0, &
          0.71511720328842845913E0 /
      data (B(I), I = 13, 25) /                                &
         -0.00000000001951073787E0, -0.00000000032302692214E0, &
          0.00000000522461866919E0,  0.00000000342940918551E0, &
         -0.00000035772874310272E0,  0.00000019999935792654E0, &
          0.00002687044575042908E0, -0.00011843240273775776E0, &
         -0.00080991728956032271E0,  0.00661062970502241174E0, &
          0.00909530922354827295E0, -0.20160072778491013140E0, &
          0.51169696718727644908E0 /
      data (B(I), I = 26, 38) /                                &
         0.00000000003147682272E0, -0.00000000048465972408E0,  &
         0.00000000063675740242E0,  0.00000003377623323271E0,  &
        -0.00000015451139637086E0, -0.00000203340624738438E0,  &
         0.00001947204525295057E0,  0.00002854147231653228E0,  &
        -0.00101565063152200272E0,  0.00271187003520095655E0,  &
         0.02328095035422810727E0, -0.16725021123116877197E0,  &
         0.32490054966649436974E0 /
      data (B(I), I = 39, 51) /                                &
         0.00000000002319363370E0, -0.00000000006303206648E0,  &
        -0.00000000264888267434E0,  0.00000002050708040581E0,  &
         0.00000011371857327578E0, -0.00000211211337219663E0,  &
         0.00000368797328322935E0,  0.00009823686253424796E0,  &
        -0.00065860243990455368E0, -0.00075285814895230877E0,  &
         0.02585434424202960464E0, -0.11637092784486193258E0,  &
         0.18267336775296612024E0 /
      data (B(I), I = 52, 64) /                                &
        -0.00000000000367789363E0,  0.00000000020876046746E0,  &
        -0.00000000193319027226E0, -0.00000000435953392472E0,  &
         0.00000018006992266137E0, -0.00000078441223763969E0,  &
        -0.00000675407647949153E0,  0.00008428418334440096E0,  &
        -0.00017604388937031815E0, -0.00239729611435071610E0,  &
         0.02064129023876022970E0, -0.06905562880005864105E0,  &
         0.09084526782065478489E0 /
      W = ABS(X)
      if (W .LT. 2.2D0) then
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T +    &
              A(K + 11)) * T + A(K + 12)) * W
      elseif (W .LT. 6.9D0) then
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T +     &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T +     &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T +    &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      else
          Y = 1
      endif
      if (X .LT. 0) Y = -Y
      DERF = Y

 end function DERF

!------------------------------------------------------------------------------------------!

 logical function isnan(arg1)
       real,intent(in) :: arg1
       isnan=( arg1  .ne. arg1 )
       return
 end function isnan

!------------------------------------------------------------------------------------------!

!==========================================================================================!
 subroutine icecat_destination(Qi,Di,D_nuc,deltaD_init,ni_add,iice_dest)

 !--------------------------------------------------------------------------------------!
 ! Returns the index of the destination ice category into which new ice is nucleated.
 !
 ! New ice will be nucleated into the category in which the existing ice is
 ! closest in size to the ice being nucleated.  The exception is that if the
 ! size difference between the nucleated ice and existing ice exceeds a threshold
 ! value for all categories, then ice is initiated into a new category.
 !
 ! D_nuc        = mean diameter of new particles being added to a category
 ! D(i)         = mean diameter of particles in category i
 ! diff(i)      = |D(i) - D_nuc|
 ! deltaD_init  = threshold size difference to consider a new (empty) category
 ! mindiff      = minimum of all diff(i) (for non-empty categories)
 !
 ! POSSIBLE CASES                      DESTINATION CATEGORY
 !---------------                      --------------------
 ! case 1:  all empty                  category 1
 ! case 2:  all full                   category with smallest diff
 ! case 3:  partly full
 !  case 3a:  mindiff <  diff_thrs     category with smallest diff
 !  case 3b:  mindiff >= diff_thrs     first empty category
 !--------------------------------------------------------------------------------------!

 implicit none

! arguments:
 real, intent(in), dimension(:) :: Qi,Di
 real, intent(in)               :: D_nuc,deltaD_init
 integer, intent(out)           :: iice_dest
 logical, intent(out)           :: ni_add

! local variables:
 logical                        :: all_full,all_empty
 integer                        :: i_firstEmptyCategory,iice,i_mindiff,n_cat
 real                           :: mindiff,diff
 real, parameter                :: qsmall_loc = 1.e-14

 !--------------------------------------------------------------------------------------!

 n_cat     = size(Qi)
 ni_add    = .true.
 iice_dest = -99

 if (sum(Qi(:))<qsmall_loc) then

 !case 1:
    iice_dest = 1
    return

 else

    all_full  = .true.
    all_empty = .false.
    mindiff   = 9.e+9
    i_firstEmptyCategory = 0

    do iice = 1,n_cat
       if (Qi(iice) .ge. qsmall_loc) then
          all_empty = .false.
          diff      = abs(Di(iice)-D_nuc)
          if (diff .lt. mindiff) then
             mindiff   = diff
             i_mindiff = iice
          endif
       else
          all_full = .false.
          if (i_firstEmptyCategory.eq.0) i_firstEmptyCategory = iice
       endif
    enddo

    if (all_full) then
 !case 2:
       iice_dest = i_mindiff
       if (mindiff .ge. 100.e-6) ni_add=.false.
       return
    else
       if (mindiff .lt. deltaD_init) then
 !case 3a:
          iice_dest = i_mindiff
          return
       else
 !case 3b:
          iice_dest = i_firstEmptyCategory
          return
       endif
    endif

 endif

 print*, 'ERROR in s/r icecat_destination -- made it to end'
 stop

 end subroutine icecat_destination


!======================================================================================!

END MODULE MP_P3
