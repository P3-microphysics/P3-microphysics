!--------------------------------------------------------------------------------------------------!
! This module is the Predicted Particle Property scheme (P3-ICE). For details                      !
! see Morrison and Milbrandt (2014 [in preparation]).                                              !
!                                                                                                  !
!--------------------------------------------------------------------------------------------------!
! Version:        v17                                                                              !
! Last modified:  2013-11-22                                                                       !
!--------------------------------------------------------------------------------------------------!

MODULE MODULE_MP_P3

   IMPLICIT NONE

   REAL, PARAMETER :: PI = 3.1415926535897932384626434
   REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297

   PUBLIC  ::  MP_P3_WRAPPER
   PUBLIC  ::  POLYSVP1

   PRIVATE :: gamma, derf
   PRIVATE :: PI, SQRTPI
!  PRIVATE :: P3_MAIN   !WRF
   PUBLIC  :: P3_MAIN   !GEM
   private :: isnan

! integer switch for warm rain autoconversion/accretion schemes
   integer, private :: iparam

! ice microphysics lookup table array dimensions
   integer, parameter :: isize=20
   integer, parameter :: jsize=20
   integer, parameter :: densize=5
   integer, parameter :: rimsize=4
   integer, parameter :: rcollsize=30

! number of ice microphysical quantities used from lookup table
   integer, parameter :: tabsize=12

! number of ice-rain collection microphysical quantities used from lookup table
   integer, parameter :: colltabsize=2

   private :: isize,jsize,densize,rimsize,rcollsize,tabsize,colltabsize

!ice lookup table values
      real :: itab(densize,rimsize,isize,jsize,tabsize)

!ice lookup table values for ice-rain collision/collection
      double precision itabcoll(densize,rimsize,isize,jsize,rcollsize,colltabsize)

      private :: itab,itabcoll

! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
! warm rain autoconversion/accretion option only (iparam = 1)
      real dnu(16)
      private :: dnu

! physical constants
   real, private :: rhosur,rhosui,ar,br,f1r,f2r,ecr,qsmall,rhow,kr,kc, &
                    bimm,aimm,rin,mi0,eci,eri,eii,bcn,cpw,e0,cons1, &
                    cons2,cons3,cons4,cons5,cons6,cons7,nccnst

! physical constants
   real, private :: cp,g,rd,rv,ep_2,ocp

! aerosol and droplet activation parameters
   real, private :: mw,osm,vi,epsm,rhoa,map,ma,rr,bact,rm1,sig1,nanew1, &
       f11,f21,rm2,sig2,nanew2,f12,f22

! lookup table values for rain shape parameter mur
     real, private :: mur_table(150)
! lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
     real, private :: vn_table(300,10),vm_table(300,10),revap_table(300,10)

       contains

!==================================================================================================!

SUBROUTINE P3_INIT

!--------------------------------------------------------------------------------------------------!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AND PARAMETERS
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!--------------------------------------------------------------------------------------------------!

 IMPLICIT NONE

 integer i,j,k,ii,jj,kk

! local variables for calculating mur lookup table
 real lamr,mur,lamold,dum,initlamr

! local variables for rain fallspeed/ventilation lookup table
 real dm,dum1,dum2,dum3,dum4,dum5
 real dd,amg,vt,dia,vn,vm

! switch for warm-rain parameterization
! = 1 seifert and beheng 2001
! = 2 beheng 1994
! = 3 khairoutdinov and kogan 2000
 iparam = 3

! droplet concentration (m-3)
 nccnst=250.e6

! parameters for Seifert and Beheng (2001) autoconversion/accretion
 kc = 9.44e9
 kr = 5.78e3

! physical constants
 cp = 1005.
 ocp=1./cp
 g=9.816
 rd=287.15
 rv=461.51
 ep_2=0.622

 rhosur = 100000./(rd*273.15)
 rhosui = 60000./(rd*253.15)
 ar = 841.99667
 br = 0.8
 f1r = 0.78
 f2r = 0.32
 ecr = 1.
 qsmall = 1.e-14
 rhow = 997.
 cpw = 4218.

! Bigg (1953)
!      bimm = 100.
!      aimm = 0.66
! Barklie and Gokhale (1959)
 bimm = 2.
 aimm = 0.65
 rin=0.1e-6
 mi0=4./3.*pi*900.*1.e-18

 eci=1.0
!eci=0.5
 eri=1.
 eii=0.1
 bcn=2.

! saturation pressure at T = 0 C
 e0 = polysvp1(273.15,0)

 cons1=pi/6.*rhow
 cons2=4./3.*pi*rhow
 cons3=1./(cons2*(25.e-6)**3)
 cons4=1./((600.e-6)**3*pi*rhow)
 cons5=pi/6.*bimm
 cons6=pi*pi/36.*rhow*bimm
 cons7=4./3.*pi*rhow*(1.e-6)**3

! aerosol/droplet activation parameters
 mw = 0.018
 osm = 1.
 vi = 3.
 epsm = 0.9
 rhoa = 1777.
 map = 0.132
 ma = 0.0284
 rr = 8.3187
 bact = vi*osm*epsm*mw*rhoa/(map*rhow)

! mode 1
 rm1 = 0.05e-6
 sig1 = 2.0
 nanew1 = 300.e6
 f11 = 0.5*exp(2.5*(log(sig1))**2)
 f21 = 1.+0.25*log(sig1)

! mode 2
 rm2 = 1.3e-6
 sig2 = 2.5
 nanew2 = 0.
 f12 = 0.5*exp(2.5*(log(sig2))**2)
 f22 = 1.+0.25*log(sig2)

! parameters for droplet mass spectral shape, used by Seifert and Beheng (2001)
! warm rain scheme only (iparam = 1)
 dnu(1) = 0.
 dnu(2) = -0.557
 dnu(3) = -0.430
 dnu(4) = -0.307
 dnu(5) = -0.186
 dnu(6) = -0.067
 dnu(7) = 0.050
 dnu(8) = 0.167
 dnu(9) = 0.282
 dnu(10) = 0.397
 dnu(11) = 0.512
 dnu(12) = 0.626
 dnu(13) = 0.739
 dnu(14) = 0.853
 dnu(15) = 0.966
 dnu(16) = 0.966

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in ice microphysics table

!open(unit=10,file='/sawtooth/morrison/kinematic_ice/lookup_ice_prog_de_v13.dat',status='old')
open(unit=10,file='lookup_ice_prog_de_v13.dat',status='old')

 do jj=1,densize
    do ii=1,rimsize
       do i=1,isize
          do k=1,jsize
             read(10,*) dum,dum,dum,dum,itab(jj,ii,i,k,1),itab(jj,ii,i,k,2), &
                  itab(jj,ii,i,k,3),itab(jj,ii,i,k,4), &
                  itab(jj,ii,i,k,5),itab(jj,ii,i,k,6), &
                  itab(jj,ii,i,k,7),itab(jj,ii,i,k,8),dum, &
                  itab(jj,ii,i,k,9),itab(jj,ii,i,k,10),itab(jj,ii,i,k,11),itab(jj,ii,i,k,12)
           end do
        end do
! read in table for ice-rain collection
       do i=1,isize
          do k=1,jsize
             do j=1,rcollsize
                read(10,*) dum,dum,dum,dum,dum,itabcoll(jj,ii,i,k,j,1), &
                 itabcoll(jj,ii,i,k,j,2),dum
                 itabcoll(jj,ii,i,k,j,1)=dlog10(itabcoll(jj,ii,i,k,j,1))
                 itabcoll(jj,ii,i,k,j,2)=dlog10(itabcoll(jj,ii,i,k,j,2))
             end do
          end do
       end do
    end do
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! generate lookup table for rain shape parameter mur
! this is very fast so it can be generated at the start of each run
! make a 150x1 1D lookup table, this is done in parameter
! space of a scaled mean size proportional qr/Nr -- initlamr

 do i=1,150              ! loop over lookup table values
    initlamr = (real(i)*2.)/1.e6+250.e-6
    initlamr = 1./initlamr
! iterate to get mur
! mur-lambda relationship is from Cao et al. (2008), eq. (7)

! start with first guess, mur = 0

    mur=0.

    do ii=1,50
       LAMR = initlamr*((mur+3.)*(mur+2.)*(mur+1.)/6.)**0.33333

! new estimate for mur based on lambda
! set max lambda in formula for mur to 20 mm-1, so Cao et al.
! formula is not extrapolated beyond Cao et al. data range
       dum=min(20.,lamr/1.e3)
       mur=max(0.,-0.0201*dum**2+0.902*dum-1.718)

! if lambda is converged within 0.1%, then exit loop
       if (ii.ge.2) then
          if (abs((lamold-lamr)/lamr).lt.0.001) goto 111
       end if

       lamold=lamr

    end do

111 continue

! assign lookup table values
    mur_table(i)=mur

 end do

!.......................................................................
! generate lookup table for rain fallspeed and ventilation parameters
! the lookup table is two dimensional as a function of number-weighted mean size
! proportional to qr/Nr and shape parameter mur

! loop over mur
 do ii=1,10

    mur=real(ii-1)  ! values of mu

! loop over number-weighted mean size
    do jj=1,300

       if (jj.le.20) then
          dm=(real(jj)*10.-5.)*1.e-6      ! mean size (meter)
       else if (jj.gt.20) then
          dm=(real(jj-20)*30.+195.)*1.e-6 ! mean size (meter)
       end if

! calculate PSD parameters from dm and mur

       lamr=(mur+1)/dm

! do numerical integration over PSD

       dum1=0. ! numerator, number-weighted fallspeed
       dum2=0. ! denominator, number-weighted fallspeed
       dum3=0. ! numerator, mass-weighted fallspeed
       dum4=0. ! denominator, mass-weighted fallspeed
       dum5=0. ! term for ventilation factor in evap

       dd=2.

! loop over PSD to numerically integrate number and mass-weighted mean fallspeeds
       do kk=1,10000
          dia=(real(kk)*dd-dd/2.)*1.e-6 ! size bin (meter)
          amg=pi/6.*997.*dia**3 ! mass (kg)
          amg=amg*1000.  ! convert kg to g

! get fallspeed as a function of size, in cm/s
          if(dia*1.e6.le.134.43) then
             vt=4.5795e5*amg**(2./3.)
             go to 101
          endif
          if(dia*1.e6.lt.1511.64) then
             vt=4.962e3*amg**(1./3.)
             go to 101
          endif
          if(dia*1.e6.lt.3477.84) then
             vt=1.732e3*amg**(1./6.)
             go to 101
          endif
          vt=917.

 101      continue

          vt=vt*1.e-2 ! convert from cm to m

          dum1=dum1+vt*10.**(mur*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6
          dum2=dum2+10.**(mur*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6

          dum3=dum3+vt*10.**((mur+3.)*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6
          dum4=dum4+10.**((mur+3.)*alog10(dia)+4.*mur)*exp(-lamr*dia)*dd*1.e-6

          dum5=dum5+(vt*dia)**0.5*10.**((mur+1.)*alog10(dia)+3.*mur)*exp(-lamr*dia)*dd*1.e-6

       end do ! loop over PSD

       vn_table(jj,ii)=dum1/dum2
       vm_table(jj,ii)=dum3/dum4
       dum=alog10(dum5)+(mur+1.)*alog10(lamr)-(3.*mur)
       revap_table(jj,ii)=10.**(dum)

    end do ! loop over number-weighted mean size

 end do ! loop over mur

! end rain lookup table generation
!..................................................

END SUBROUTINE P3_INIT
!==================================================================================================!

!--------------------------------------------------------------------------------------------------!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE P3 MICROPHYSICS SCHEME.
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 2D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE P3_MAIN).
! 2D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE.
! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES IN SUBROUTINE P3_MAIN
! BEFORE BEING PASSED BACK TO DRIVER MODEL.
!
! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND JASON MILBRANDT (ENVIRONMENT CANADA)
!
! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916
!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE MP_P3_WRAPPER(th3d,qv3d,qc3d,qr3d,qi3d,qg3d,qnr3d,qni3d,qvolg3d,     &
                               pii,p,dz,w,dt_in,itimestep,rainnc,rainncv,sr,    &
                               ids, ide, jds, jde, kds, kde ,                   &
                               ims, ime, jms, jme, kms, kme ,                   &
                               its, ite, jts, jte, kts, kte ,                   &
                               zdbz3d,effc3d,effi3d,vmi3d,di3d,rhopo3d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                                       ims, ime, jms, jme, kms, kme , &
                                       its, ite, jts, jte, kts, kte

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
         th3d,qv3d,qc3d,qr3d,qi3d,qg3d,qnr3d,qni3d,qvolg3d,zdbz3d,effc3d,effi3d,vmi3d,di3d,rhopo3d

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                          pii, p, dz, w

   REAL, INTENT(IN):: dt_in
   INTEGER, INTENT(IN):: itimestep

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          RAINNC, RAINNCV, SR

! local variables

!    REAL, DIMENSION(its:ite,kts:kte) ::                                &
!          th2d,qv2d,p2d,dz2d,w2d,qc2d,qr2d,qdi2d,qri2d,                 &
!          nc2d,nr2d,ni2d,bg2d,z_dbz,murr,vmi,di,rhopo,revap

   REAL, DIMENSION(its:ite) :: PRECPRT, SNOWRT

   INTEGER I,K,J

   REAL DT

   LOGICAL, parameter :: nk_BOTTOM     = .false.  !.F. --> nk at model top (as in WRF)

!..........................................................................

   DT = DT_IN

   do j=jts,jte      ! j loop (north-south)
!
! Transfer 3D arrays into 2D for microphysical calculations
!

! hm , initialize 1d tendency arrays to zero

!   do i=its,ite      ! i loop (east-west)
!      do k=kts,kte   ! k loop (vertical)
!
!          qc2d(i,k)      = qc3d(i,k,j)
!          qr2d(i,k)      = qr3d(i,k,j)
!          qdi2d(i,k)     = qi3d(i,k,j)
!          qri2d(i,k)     = qg3d(i,k,j)
!          nr2d(i,k)      = qnr3d(i,k,j)
!          ni2d(i,k)      = qni3d(i,k,j)
!          bg2d(i,k)      = qvolg3d(i,k,j)
!           th2d(i,k)       = th3d(i,k,j)
!           QV2D(I,k)       = QV3d(i,k,j)
!           P2D(I,k)        = P(i,k,j)
!           DZ2D(I,k)       = DZ(i,k,j)
!           W2D(I,k)        = W(i,k,j)
!
!     end do ! k loop
!     end do ! i loop

! For future:
!       CALL P3_MAIN(qc(its:ite,kts:kte,j), ...  !then remove qc2d

!!!    ...and do it like this:
       CALL P3_MAIN(qc3d(its:ite,kts:kte,j),qr3d(its:ite,kts:kte,j),qnr3d(its:ite,kts:kte,j),th3d(its:ite,kts:kte,j),     &
                    qv3d(its:ite,kts:kte,j),dt,qi3d(its:ite,kts:kte,j),qg3d(its:ite,kts:kte,j),qni3d(its:ite,kts:kte,j),  &
                    qvolg3d(its:ite,kts:kte,j),W(its:ite,kts:kte,j),P(its:ite,kts:kte,j),DZ(its:ite,kts:kte,j),itimestep, &
                    precprt,snowrt,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte,zdbz3d(its:ite,kts:kte,j),             &
                    effc3d(its:ite,kts:kte,j),effi3d(its:ite,kts:kte,j),                                                  &
                    vmi3d(its:ite,kts:kte,j),di3d(its:ite,kts:kte,j),rhopo3d(its:ite,kts:kte,j),                          &
                    nk_BOTTOM)

!!!
!        CALL P3_MAIN(qc2d,qr2d,nr2d,th2d,qv2d,dt,qdi2d,qri2d,ni2d,bg2d,w2d,p2d,dz2d,      &
!                     itimestep,precprt,snowrt,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,    &
!                     kts,kte,z_dbz,vmi,di,rhopo,revap,nk_BOTTOM)

!    do i=its,ite      ! i loop (east-west)
!       do k=kts,kte
! hm, add tendencies to update global variables
! HM, TENDENCIES FOR Q AND N NOW ADDED IN BUPPS MICRO, SO WE
! ONLY NEED TO TRANSFER 2D PROGNOSTIC VARIABLES BACK TO 3D
!          qc3d(i,k,j)      = qc2d(i,k)
!          qr3d(i,k,j)      = qr2d(i,k)
!          qi3d(i,k,j)      = qdi2d(i,k)
!          qg3d(i,k,j)      = qri2d(i,k)
!          qnr3d(i,k,j)     = nr2d(i,k)
!          qni3d(i,k,j)     = ni2d(i,k)
!          qvolg3d(i,k,j)   = bg2d(i,k)
!          qv3d(i,k,j)     = qv2d(i,k)
!          th3d(i,k,j)     = th2d(i,k)
!          zdbz3d(i,k,j)   = z_dbz(i,k)
!          vmi3d(i,k,j)   = vmi(i,k)
!          di3d(i,k,j)   = di(i,k)
!          rhopo3d(i,k,j)   = rhopo(i,k)
!
!       end do ! k loop

! surface precipitation output
      RAINNC(its:ite,j)  = RAINNC(its:ite,j)+PRECPRT(:)
      RAINNCV(its:ite,j) = PRECPRT(:)
      SR(its:ite,j)      = SNOWRT(:)/(PRECPRT(:)+1.E-12)

!       RAINNC(i,j) = RAINNC(I,J)+PRECPRT(i)
!       RAINNCV(i,j) = PRECPRT(i)
!       SR(i,j) = SNOWRT(i)/(PRECPRT(i)+1.E-12)

!   end do ! i loop

!for future:
!   RAINNC(:,j) = RAINNC(:,J)+PRECPRT(:)

   end do ! j loop

END SUBROUTINE MP_P3_WRAPPER

!==================================================================================================!

SUBROUTINE P3_MAIN(qc,qr,nr,th,qv,dt,qi1,qri1,ni1,bg1,uzpl,pres,dzq,it,precprt,snowrt,   &
                   ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte,z_dbz,effc,effi1,     &
                   vmi,di,rhopo,nk_BOTTOM)

!--------------------------------------------------------------------------------------------------!
! This code written by H. Morrison (12/20/12) MMM Division, NCAR.
!--------------------------------------------------------------------------------------------------!

 implicit none

! declarations

 LOGICAL, intent(in)              :: nk_BOTTOM   ! .F. for k=nk at model top (WRF) / .T. for k=nk at model BOTTOM (GEM)

! INPUT/OUTPUT PARAMETERS
 INTEGER, INTENT( IN)  :: IMS,IME, JMS,JME, KMS,KME,          &
                               ITS,ITE, JTS,JTE, KTS,KTE
! 2d input/output variables
 real, dimension(its:ite,kts:kte) :: qc   ! cloud water mixing ratio (kg/kg)
 real, dimension(its:ite,kts:kte) :: qr
 real, dimension(its:ite,kts:kte) :: nc   ! note: specified in this version of scheme
 real, dimension(its:ite,kts:kte) :: nr
 real, dimension(its:ite,kts:kte) :: qi1
 real, dimension(its:ite,kts:kte) :: qri1
 real, dimension(its:ite,kts:kte) :: ni1
 real, dimension(its:ite,kts:kte) :: bg1
 real, dimension(its:ite,kts:kte) :: qv
 real, dimension(its:ite,kts:kte) :: th
 real, dimension(its:ite,kts:kte) :: uzpl
 real, dimension(its:ite,kts:kte) :: pres
 real, dimension(its:ite,kts:kte) :: dzq
 real, dimension(its:ite,kts:kte) :: mur
 real, dimension(its:ite,kts:kte) :: ssat
 real, dimension(its:ite)         :: precprt,snowrt
 real, dimension(its:ite,kts:kte) :: effc
 real, dimension(its:ite,kts:kte) :: effi1

! output microphysical process rates (currently not used!)
 real, dimension(its:ite,kts:kte) :: rcon
 real, dimension(its:ite,kts:kte) :: ccon
 real, dimension(its:ite,kts:kte) :: auto
 real, dimension(its:ite,kts:kte) :: acc
 real, dimension(its:ite,kts:kte) :: act
 real, dimension(its:ite,kts:kte) :: opre
 real, dimension(its:ite,kts:kte) :: opra
 real, dimension(its:ite,kts:kte) :: oprc
 real, dimension(its:ite,kts:kte) :: onpra
 real, dimension(its:ite,kts:kte) :: onprc
 real, dimension(its:ite,kts:kte) :: oncagg
 real, dimension(its:ite,kts:kte) :: onragg
 real, dimension(its:ite,kts:kte) :: onpccn
 real, dimension(its:ite,kts:kte) :: opcc
 real, dimension(its:ite,kts:kte) :: opccn
 real, dimension(its:ite,kts:kte) :: onprc1
 real, dimension(its:ite,kts:kte) :: onsubr

 real dt  ! model time step (sec)

! local variables
 real, dimension(its:ite,kts:kte) :: t ! temperature (K)

! time-varying parameters
 real mu,dv,sc,dqsdt,ab,kap,epsr,epsc,epsi
 real, dimension(its:ite,kts:kte) :: rho,rhofacr,rhofaci,acn,xxls,xxlv,xlf, &
           qvs,qvi,sup,supi,ss

! variables for cond/evap/sub/dep
 real xx,aaa,epsilon

! local variables - aerosol/droplet activation
 real sigvl,aact,alpha,gamm,gg,psi,eta1,eta2,sm1,sm2, &
      smax,uu1,uu2

! local/dummy variables
 real dum,dum1,dum2,dumt,dumqv,dumqvs,dums,dumqc,ratio,qsat0
 real udiff

! warm microphysical process rates
 real pre,pra,prc,npra,nprc,ncagg,nragg,npccn,pcc,pccn, &
         nprc1,pre1,pcc1,nsubr

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
 real lammax,lammin

! parameters for rain fallspeed lookup table
 integer dumi
 real dum3,dum4,dum5

! integer indices for ice lookup table
 integer i,k,kk,ii,jj
 integer it

! parameters for mur lookup table
 real lamold
 real rdumii
 real rdumjj

! local 2D ice variables
 real, dimension(its:ite,kts:kte) :: vtrmi1
 real, dimension(its:ite,kts:kte) :: vtrni1

! ice microphysical process rates

 real psacw1,prd1,pracs1 &
       ,mnuccc,mnuccr,mnucic,mnucir,mnuccd,npsacw1 &
      ,nprd1,npracs1,nnuccc &
      ,nnuccr,nnucic,nnucir,nnuccd,prd11,qmlt1, &
      nmlt1,nsubi1,niagg1, &
      nmult1,pracsw1,npracsw1,nrshed1,ncshed1,qcshed1, &
      qwgrth1

 integer j,dumk,dumj,dumii,dumjj

! interpolated quantities from ice lookup table

 real f1pr1,f1pr2,f1pr3,f1pr4,f1pr5,f1pr6,f1pr7,f1pr8
 real f1pr9,f1pr10,f1pr12,f1pr13,f1pr14,f1pr15,f1pr16

! time-varying parameters for ice microphysics
 real dqsidt,abi,dumqvi
 real dap,nacnt,rhop

! variables for rime density calculation
 real v_impact,ri,rhorime_c,rhorime_r,iTc,D_c,D_r

! 2D reflectivity variables
 real, dimension(its:ite,kts:kte) :: z_dbz,ze_ice,ze_rain

! mass-weighted fallspeed (w/o air density correction) and mean ice size
 real, dimension(its:ite,kts:kte) :: vmi,di,rhopo,prec,revap

! tendencies from sedimentation
 REAL, DIMENSION(KTS:KTE) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
 REAL, DIMENSION(KTS:KTE) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
 REAL, DIMENSION(KTS:KTE) ::  QRISTEN
 REAL, DIMENSION(KTS:KTE) ::  BGSTEN
 REAL, DIMENSION(KTS:KTE) ::  NISTEN
 REAL, DIMENSION(KTS:KTE) ::  NRSTEN
 REAL, DIMENSION(KTS:KTE) ::  QCSTEN
 REAL, DIMENSION(KTS:KTE) ::  NCSTEN

! local variables for sedimentation calculations
 REAL, DIMENSION(KTS:KTE) ::    DUMQI,DUMR,DUMFNI,DUMRI,DUMBG,DUMFNR,DUMC,DUMFNC
 REAL, DIMENSION(KTS:KTE) ::    FR, FI, FNI, FNR, FC, FNC
 REAL, DIMENSION(KTS:KTE) ::   FALOUTR,FALOUTI,FALOUTNI,FALOUTNR, &
                                    FALOUTRI,FALOUTBG,FALOUTC,FALOUTNC
 REAL FALTNDR,FALTNDI,FALTNDNI,FALTNDRI,FALTNDBG,FALTNDC,FALTNDNC,FALTNDNR
 real rgvm
 integer n,nstep
 integer ktop,kbot,kdir    ! level indices for generalized direction

! integer switch to skip code calculations
! for future: use logicals instead
 integer ltrue,htrue

 integer index

! for sedimentation
 logical qcpresent,qrpresent,qipresent
 integer qcindex,qrindex,qiindex

! local variables
 real dumlr,wetgrowth

! varaibles for speeding up code
 real onstep,odum,odum3,odt,odzq(its:ite,kts:kte),orho(its:ite,kts:kte),oxx,oabi

 real zero,test,test2,test3

! end variable declarations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	test2=1.
!	test=1.
!	zero=0.0
!	test=test/zero
!	test2=zero/zero
!	print*,test
!	print*,test2
!	print*,test3
!	print*,'here'
!	stop

 !direction of vertical leveling:
 if (nk_BOTTOM) then
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


! calculate inverse model time step
 odt = 1./dt

 do i = its,ite  ! main loop around the entire scheme (i loop)

    ltrue=0
    htrue=0

! initialize surface precipitation rate output
    precprt(i) = 0.
    snowrt(i) = 0.

!        do k = kts,kte
    do k = kbot,ktop,kdir

! sensitivity
!         qdi1(i,k)=0.
!         qri1(i,k)=1.e-3
!         ni1(i,k)=10.e3
!         qr(i,k)=1.e-3
!         nr(i,k)=0.7e5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize effective radii to default values
       effc(i,k) = 10.e-6
       effr(i,k) = 25.e-6
       effi1(i,k) = 25.e-6

! initialize reflectivity output
       z_dbz(i,k) = -35.

! initialize output mass-weighted ice fallspeed and mean size
       vmi(i,k) = 0.
       di(i,k) = 0.
       rhopo(i,k) = 0.
       prec(i,k) = 0.
       revap(i,k) = 0.

! initialize rain mu
       mur(i,k)   = 0.

! initialize reflectivity
       ze_ice(i,k) = 1.e-22
       ze_rain(i,k) = 1.e-22

! calculate temperature from theta
       dum = (100000./pres(i,k))**(rd*ocp)
       t(i,k) = th(i,k)/dum

! calculate some time-varying atmospheric variables
       rho(i,k) = pres(i,k)/(rd*t(i,k))
       orho(i,k) = 1./rho(i,k)
       xxlv(i,k) = 3.1484e6-2370.*t(i,k)
       xxls(i,k) = xxlv(i,k)+0.3337e6
       xlf(i,k)=xxls(i,k)-xxlv(i,k)

       qvs(i,k) = ep_2*polysvp1(t(i,k),0)/ &
            (pres(i,k)-polysvp1(t(i,k),0))
       qvi(i,k) = ep_2*polysvp1(t(i,k),1)/ &
            (pres(i,k)-polysvp1(t(i,k),1))
       sup(i,k)=qv(i,k)/qvs(i,k)-1.
       supi(i,k)=qv(i,k)/qvi(i,k)-1.

!!write(*,*) 'SUP: ',sup(i,k),t(i,k),qvs(i,k),qvi(i,k),uzpl(i,k)

! specify droplet concentration
       nc(i,k)=nccnst*orho(i,k)

! apply mass clipping if conditions are dry and mass mixing ratio is sufficiently small
! (and thus all mass is expected to evaporate/sublimate in one time step)
       if (qi1(i,k).lt.qsmall.or.(qi1(i,k).lt.1.e-8.and.supi(i,k).lt.0.9)) then
          qv(i,k)=qv(i,k)+qi1(i,k)
          th(i,k)=th(i,k)+th(i,k)/t(i,k)*qi1(i,k)*xxls(i,k)*ocp
          qi1(i,k)=0.
          ni1(i,k)=0.
          qri1(i,k)=0.
          bg1(i,k)=0.
       else
          ltrue=1
       end if

       if (qi1(i,k).ge.qsmall.and.qi1(i,k).lt.1.e-8.and.t(i,k).ge.273.15) then
          qr(i,k)=qr(i,k)+qi1(i,k)
          th(i,k)=th(i,k)+th(i,k)/t(i,k)*qi1(i,k)*xlf(i,k)*ocp
          qi1(i,k)=0.
          ni1(i,k)=0.
          qri1(i,k)=0.
          bg1(i,k)=0.
        end if

        if (qc(i,k).lt.qsmall.or.(qc(i,k).lt.1.e-8.and.sup(i,k).lt.0.9)) then
           qv(i,k) = qv(i,k)+qc(i,k)
           th(i,k) = th(i,k)+th(i,k)/t(i,k)*qc(i,k)*xxlv(i,k)*ocp
           qc(i,k) = 0.
           nc(i,k) = 0.
        else
           ltrue=1
        end if

        if (qr(i,k).lt.qsmall.or.(qr(i,k).lt.1.e-8.and.sup(i,k).lt.0.9)) then
           qv(i,k) = qv(i,k)+qr(i,k)
           th(i,k) = th(i,k)+th(i,k)/t(i,k)*qr(i,k)*xxls(i,k)*ocp
           qr(i,k) = 0.
           nr(i,k) = 0.
        else
           ltrue=1
        end if

! if ltrue is 1 then set htrue is 1 so that microphysical processes ARE calculated
        if (ltrue.eq.1.) htrue=1
! if there is the possibility of nucleation or droplet activation, i.e., if RH
! is relatively high, then calculate microphysical processes even if there
! is no existing condensate
        if (t(i,k).lt.273.15.and.supi(i,k).ge.-0.05) htrue=1
        if (t(i,k).ge.273.15.and.sup(i,k).ge.-0.05) htrue=1

     end do

! jump to end of loop if htrue.eq.0
         if (htrue.eq.0) goto 333

! reset ltrue
         ltrue = 0

!..................................................
! main k loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       do k = kts,kte
        do k = kbot,ktop,kdir

! if it's dry and there are no hydrometeors at this vertical level, skip to end of k loop
        if (it.gt.1.and.qi1(i,k).lt.qsmall.and.qc(i,k).lt.qsmall.and.qr(i,k).lt.qsmall) then
           if (t(i,k).lt.273.15.and.supi(i,k).lt.-0.05) goto 555
           if (t(i,k).ge.273.15.and.sup(i,k).lt.-0.05) goto 555
        end if

! initialize warm microphysics processes to zero
            pra = 0.
            prc = 0.
            npra = 0.
            nprc = 0.
            nprc1 = 0.
            nragg = 0.
            ncagg = 0.
            pre = 0.
            pcc = 0.
            pre1 = 0.
            pcc1 = 0.
            npccn = 0.
            pccn=0.
            nsubr = 0.

! initialize ice microphysical processes to zero
            psacw1=0.  ! collection cloud water species 1
            prd1=0.    ! vapor deposition species 1
            prd11=0.    ! sublimation species 1
            pracs1=0.  ! collection rain mass by snow
            pracsw1=0.  ! collection snow mass by rain
            mnuccc=0.  ! contact freezing droplets
            mnuccr=0.  ! contact freezing rain
            mnucic=0.  ! immersion freezing droplets
            mnucir=0.  ! immersion freezing rain
            mnuccd=0.  ! deposition/condensation freezing nuc
            qmlt1=0.
            npsacw1=0.
            nprd1=0.
            npracs1=0.
            npracsw1=0.
            nrshed1=0.
            ncshed1=0.
            qcshed1=0.
            nnuccc=0.
            nnuccr=0.
            nnucic=0.
            nnucir=0.
            nnuccd=0.
            nsubi1=0.
            nmlt1=0.
            niagg1=0.
            nmult1=0.
            qwgrth1=0.
            wetgrowth=0.

	    goto 711

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! adjust cloud water and thermodynamics to prognostic supersaturation
! following the method in Grabowski and Morrison (2008)

         dqsdt = xxlv(i,k)*qvs(i,k)/(rv*t(i,k)*t(i,k))
         ab = 1.+dqsdt*xxlv(i,k)*ocp
         epsilon = (qv(i,k)-qvs(i,k)-ssat(i,k))/ab

! limit adjustment to available water

         epsilon=max(epsilon,-qc(i,k))

! don't adjust upward if subsaturated
! otherwise this could result in positive adjustment
! (spurious generation ofcloud water) in subsaturated conditions

         if (ssat(i,k).lt.0.) epsilon=min(0.,epsilon)

! now do the adjustment
         if (abs(epsilon).ge.1.e-15) then

         qc(i,k)=qc(i,k)+epsilon
         qv(i,k)=qv(i,k)-epsilon
         th(i,k)=th(i,k)+epsilon*th(i,k)/t(i,k)*xxlv(i,k)*ocp

! recalculate variables if there was adjustment
! calculate temperature from theta
         dum = (100000./pres(i,k))**(rd*ocp)
         t(i,k) = th(i,k)/dum
         qvs(i,k) = ep_2*polysvp1(t(i,k),0)/ &
           (pres(i,k)-polysvp1(t(i,k),0))
         qvi(i,k) = ep_2*polysvp1(t(i,k),1)/ &
           (pres(i,k)-polysvp1(t(i,k),1))
         sup(i,k)=qv(i,k)/qvs(i,k)-1.
         supi(i,k)=qv(i,k)/qvi(i,k)-1.
         ssat(i,k)=qv(i,k)-qvs(i,k)

         end if

 711     continue

! skip micro process calculations except nucleation/acvtivation if there no hydrometeors are present
        if (it.gt.1.and.qi1(i,k).lt.qsmall.and.qc(i,k).lt.qsmall.and.qr(i,k).lt.qsmall) goto 444

! time/space varying physical variables
         mu = 1.496e-6*t(i,k)**1.5/(t(i,k)+120.)
         dv = 8.794e-5*t(i,k)**1.81/pres(i,k)
         sc = mu/(rho(i,k)*dv)
	 dum=1./(rv*t(i,k)**2)
         dqsdt = xxlv(i,k)*qvs(i,k)*dum
         dqsidt = xxls(i,k)*qvi(i,k)*dum

         ab = 1.+dqsdt*xxlv(i,k)*ocp
         abi = 1.+dqsidt*xxls(i,k)*ocp
         kap = 1.414e3*mu

         ssat(i,k)=qv(i,k)-qvs(i,k)

! contact freezing currently turned off
!         dum=7.37*t(i,k)/(288.*10.*pres(i,k))/100.
!         dap=4.*pi*1.38e-23*t(i,k)*(1.+dum/rin)/ &
!                (6.*pi*rin*mu)
!         nacnt=exp(-2.80+0.262*(273.15-t(i,k)))*1000.
         rhofacr(i,k)=(rhosur*orho(i,k))**0.54
         rhofaci(i,k)=(rhosui*orho(i,k))**0.54
! 'a' parameter for droplet fallspeed calculated from Stokes' law
         acn(i,k)=g*rhow/(18.*mu)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!..................................................
! get droplet size distribution parameters
	if (qc(i,k).ge.qsmall) then

! set minimum nc to prevent floating point error
        nc(i,k)=max(nc(i,k),1.e-16)

        pgam(i,k)=0.0005714*(nc(i,k)*1.e-6*rho(i,k))+0.2714
        pgam(i,k)=1./(pgam(i,k)**2)-1.
        pgam(i,k)=max(pgam(i,k),2.)
        pgam(i,k)=min(pgam(i,k),15.)

! interpolate for mass distribution spectral shape parameter (for SB warm processes)

	if (iparam.eq.1) then
        dumi=int(pgam(i,k))
        nu(i,k)=dnu(dumi)+(dnu(dumi+1)-dnu(dumi))* &
               (pgam(i,k)-dumi)
	end if

! calculate lamc

        lamc(i,k)=(cons1*nc(i,k)*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.)/ &
            (qc(i,k)))**0.33333

! apply lambda limiters
! lammin, 40 micron mean diameter
! lammax, 1 mean micron

	lammin = (pgam(i,k)+1.)/40.e-6
	lammax = (pgam(i,k)+1.)/1.e-6

	if (lamc(i,k).lt.lammin) then
	lamc(i,k) = lammin

	nc(i,k) = 6.*lamc(i,k)**3*qc(i,k)/ &
               (pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.))

	else if (lamc(i,k).gt.lammax) then
	lamc(i,k) = lammax

	nc(i,k) = 6.*lamc(i,k)**3*qc(i,k)/ &
               (pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.))


	end if

        cdist(i,k)=nc(i,k)*(pgam(i,k)+1.)/lamc(i,k)
        cdist1(i,k) = nc(i,k)/gamma(pgam(i,k)+1.)
	else
	lamc(i,k) = 0.
        cdist(i,k)=0.
        cdist1(i,k)=0.
	end if

!..................................................
! get rain size distribution parameters
	if (qr(i,k).ge.qsmall) then

! set minimum value for Nr to prevent taking 0^(1/3) and thus code crash
         nr(i,k)=max(nr(i,k),1.e-16)

! use lookup table to get mur

! find spot in lookup table
! calculate scaled mean size proportional to N/q for lookup table parameter space
! the scaled mean size is identical to lambda for an exponential distribution
         dum = (cons1*NR(i,k)*6./(QR(i,k)))**0.33333
	 odum=1./dum

         if (odum.lt.282.e-6) then
            mur(i,k) = 8.282
         else if (odum.ge.282.e-6.and.odum.lt.502.e-6) then
! find integer index corresponding with scaled mean size dum
            rdumii = (odum-250.e-6)*1.e6/2.
            rdumii = max(rdumii,1.)
            rdumii = min(rdumii,150.)
            dumii = int(rdumii)
            dumii = min(149,dumii)
! interpolate lookup table values
            mur(i,k) = mur_table(dumii)+(mur_table(dumii+1)-mur_table(dumii))* &
                       (rdumii-real(dumii))
         else if (odum.ge.502.e-6) then
            mur(i,k)=0.
         end if

! recalculate lambda based on mur
         LAMR(i,k) = (cons1*nr(i,k)*(mur(i,k)+3.)*(mur(i,k)+2)*(mur(i,k)+1.)/   &
                 (qr(i,k)))**0.33333

! set maximum value for rain lambda
 	 lammax = (mur(i,k)+1.)*1.e5

! set to small value since breakup is explicitly included (mean size 0.8 mm)
	 lammin = (mur(i,k)+1.)*1250.

! apply lambda limiters
        IF (lamr(i,k).LT.LAMMIN) THEN

        LAMR(i,k) = LAMMIN

        NR(i,k) = EXP(3.*LOG(LAMR(i,k))+LOG(QR(i,k))+              &
                LOG(GAMMA(MUR(i,k)+1.))-LOG(GAMMA(MUR(i,k)+4.)))/(cons1)

        ELSE IF (LAMR(i,k).GT.LAMMAX) THEN

        LAMR(i,k) = LAMMAX

        NR(i,k) = EXP(3.*LOG(LAMR(i,k))+LOG(QR(i,k))+              &
                LOG(GAMMA(MUR(i,k)+1.))-LOG(GAMMA(MUR(i,k)+4.)))/(cons1)

        END IF

        CDISTR(i,k) = NR(i,k)/GAMMA(MUR(i,k)+1.)

! NOTE: n0r is calculated as log10(n0r)
        n0r(i,k) = alog10(nr(i,k))+(mur(i,k)+1.)*alog10(lamr(i,k)) &
             -alog10(gamma(mur(i,k)+1))

	else
	lamr(i,k) = 0.
	n0r(i,k) = 0.
        cdistr(i,k)=0.
	end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!................................................................................

      if (qr(i,k).ge.qsmall) then

! calculate rain evaporation including ventilation and number-
! weighted fallspeed (used for calculating rime density)
! find appropriate place in 2D rain lookup table

! find location in scaled mean size space
        dum1=(mur(i,k)+1.)/lamr(i,k)
        if (dum1.le.195.e-6) then
           rdumii=(dum1*1.e6+5.)/10.
           rdumii=max(rdumii,1.)
           rdumii=min(rdumii,20.)
           dumii=int(rdumii)
           dumii=max(dumii,1)
           dumii=min(dumii,20)
           dum3=10.
        else if (dum1.gt.195.e-6) then
           rdumii=(dum1*1.e6-195.)/30.+20.
           rdumii=max(rdumii,20.)
           rdumii=min(rdumii,300.)
           dumii=int(rdumii)
           dumii=max(dumii,20)
           dumii=min(dumii,299)
           dum3=30.
        end if

! find location in mur space
        rdumjj=mur(i,k)+1.
        rdumjj=max(rdumjj,1.)
        rdumjj=min(rdumjj,10.)
        dumjj=int(rdumjj)
        dumjj=max(dumjj,1)
        dumjj=min(dumjj,9)

! interpolate value at mur
        dum1=revap_table(dumii,dumjj)+(rdumii-real(dumii))/dum3*  &
        (revap_table(dumii+1,dumjj)-revap_table(dumii+1,dumjj))

! interoplate value at mur+1
        dum2=revap_table(dumii,dumjj+1)+(rdumii-real(dumii))/dum3*  &
        (revap_table(dumii+1,dumjj+1)-revap_table(dumii+1,dumjj+1))

! final interpolation
        dum=dum1+(rdumjj-real(dumjj))* &
              (dum2-dum1)

        EPSR = 2.*PI*CDISTR(i,k)*RHO(i,k)*DV*                           &
                   (F1R*GAMMA(MUR(i,k)+2.)/(LAMR(i,k))+                       &
                    F2R*(RHO(i,k)/MU)**0.5*                      &
                    SC**0.33333*dum)

! number-weighted fallspeed
! value at mur
!        dum1=vn_table(dumii,dumjj)+(rdumii-real(dumii))/dum3*  &
!        (vn_table(dumii+1,dumjj)-vn_table(dumii+1,dumjj))

! value at mur+1
!        dum2=vn_table(dumii,dumjj+1)+(rdumii-real(dumii))/dum3*  &
!        (vn_table(dumii+1,dumjj+1)-vn_table(dumii+1,dumjj+1))

! final interpolation
!        vtrn(i,k)=dum1+(rdumjj-real(dumjj))* &
!              (dum2-dum1)
!        vtrn(i,k)=vtrn(i,k)*rhofacr(i,k)

         else
            epsr = 0.
!            vtrn(i,k) = 0.
         end if

         if (qc(i,k).ge.qsmall) then
            epsc=2.*pi*rho(i,k)*dv*cdist(i,k)
         else
            epsc=0.
         end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.....................................................................

! ice lookup table for ice micro processes

         if (qi1(i,k).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
         ni1(i,k)=max(ni1(i,k),1.e-16)

! calculate predicted bulk rime density
         if (bg1(i,k).ge.1.e-15) then
         rhop=qri1(i,k)/bg1(i,k)
         else
         rhop=0.
         end if
! limit 50 < rhop < 900, adjust bg if needed
         if (rhop.lt.50.) then
            rhop=50.
            bg1(i,k)=qri1(i,k)/rhop
         end if
         if (rhop.gt.900.) then
            rhop=900.
            bg1(i,k)=qri1(i,k)/rhop
         end if
         if (qri1(i,k).lt.qsmall) then
            bg1(i,k)=0.
         end if

! set upper constraint on qri to ensure qri cannot be > qi
	 if (qri1(i,k).gt.qi1(i,k)) then
         qri1(i,k)=qi1(i,k)
	 bg1(i,k)=qri1(i,k)/rhop
	 end if

! find indices in 4D ice lookup table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

! find index for qi (total ice mass mixing ratio)
         dum1 = (alog10(qi1(i,k))+16.)/0.70757
         dumi=int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
         dum1=min(dum1,real(isize))
         dum1=max(dum1,1.)
         dumi=max(1,dumi)
         dumi=min(isize-1,dumi)

! find index for Ni (ice number mixing ratio)
         dum2 = (alog10(ni1(i,k))+10.)/0.90309
         dumk=int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
         dum2=min(dum2,real(jsize))
         dum2=max(dum2,1.)
         dumk=max(1,dumk)
         dumk=min(jsize-1,dumk)

! find index for scaled mean rain size
! if no rain, then just choose dumj = 1 and don't calculate rain-ice
! collection processes

         if (qr(i,k).ge.qsmall) then
! calculate scaled mean size for consistency with ice lookup table
         dumlr=(pi*rhow*nr(i,k)/qr(i,k))**0.33333
         dum3 = (alog10(1./dumlr)+5.)/0.0934217
         dumj=int(dum3)

! set limits
         dum3=min(dum3,real(rcollsize))
         dum3=max(dum3,1.)
         dumj=max(1,dumj)
         dumj=min(rcollsize-1,dumj)
         else
            dumj=1
	    dum3=1.
         end if

! find index for rime mass fraction
         dum4 = qri1(i,k)/qi1(i,k)*3.+1.
         dumii=int(dum4)

! set limits
         dum4=min(dum4,real(rimsize))
         dum4=max(dum4,1.)
         dumii=max(1,dumii)
         dumii=min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
         if (rhop.le.650.) then
         dum5 = (rhop-50.)/200.+1.
         else
         dum5=(rhop-650.)/250.+4.
         end if
         dumjj=int(dum5)

! set limits
         dum5=min(dum5,real(densize))
         dum5=max(dum5,1.)
         dumjj=max(1,dumjj)
         dumjj=min(densize-1,dumjj)

!cccccccccccccccccccccccccccccccccccccccccccc
!
! call to lookup table interpolation subroutines to get process rates

         index=2
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr2)
         index=3
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr3)
         index=4
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr4)
         index=5
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr5)
         index=7
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr9)
         index=8
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr10)
         index=10
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr14)

! ice-rain collection processes
	 if (qr(i,k).ge.qsmall) then
         index=1
         call access_lookup_table_coll(dumjj,dumii,dumj,dumi,dumk,index, &
                   dum1,dum2,dum3,dum4,dum5,f1pr7)
         index=2
         call access_lookup_table_coll(dumjj,dumii,dumj,dumi,dumk,index, &
                   dum1,dum2,dum3,dum4,dum5,f1pr8)
	 else
	 f1pr7=0.
	 f1pr8=0.
	 end if

! adjust Ni if needed to make sure mean size is in bounds (i.e., apply lambda limiters)
            ni1(i,k)=min(ni1(i,k),f1pr9)
            ni1(i,k)=max(ni1(i,k),f1pr10)

! limit max ice concentration to 500 L-1
            ni1(i,k)=min(ni1(i,k),500.e3*orho(i,k))

         end if   ! qi1 > qsmall

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!............................................................
!............................................................
! begin microphysical process calculations

! vapor deposition/sublimation onto ice
! note: capacitance factor already incluced in lookup table

         if (qi1(i,k).ge.qsmall.and.t(i,k).lt.273.15) then
           epsi=(f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)* &
           2.*pi*rho(i,k)*dv
	 else
	   epsi=0.
	 end if

! condensation/evaporation/deposition/sublimation
! use semi-analytic formulation

         if (t(i,k).lt.273.15) then
         oabi=1./abi
         xx=epsc+epsr+epsi*(1.+xxls(i,k)*ocp*dqsdt)*oabi
         else
         xx=epsc+epsr
         end if

         dum = qvs(i,k)*rho(i,k)*g*uzpl(i,k)/(pres(i,k)-polysvp1(t(i,k),0))

! 'A' term including ice (bergeron process)
! note: qv and T tendencies due to mixing and radiation are
! currently neglected --> assumed to be much smaller than cooling
! due to vertical motion which IS included


         if (t(i,k).lt.273.15) then
         aaa=-dum-dqsdt*(-uzpl(i,k)*g*ocp)- &
         (qvs(i,k)-qvi(i,k))*(1.+xxls(i,k)*ocp*dqsdt)*oabi*epsi
         else
         aaa=-dum-dqsdt*(-uzpl(i,k)*g*ocp)
         end if

! set lower bound on xx to prevent division by zero

         xx=max(1.e-8,xx)
	 oxx=1./xx

	 if (qc(i,k).ge.qsmall) then
         pcc=(aaa*epsc*oxx+(ssat(i,k)-aaa*oxx)*odt*epsc*oxx*(1.-exp(-xx*dt)))/ab
	 end if
	 if (qr(i,k).ge.qsmall) then
         pre=(aaa*epsr*oxx+(ssat(i,k)-aaa*oxx)*odt*epsr*oxx*(1.-exp(-xx*dt)))/ab
	 end if
         if (qi1(i,k).ge.qsmall.and.t(i,k).lt.273.15) then
         prd1=(aaa*epsi*oxx+(ssat(i,k)-aaa*oxx)*odt*epsi*oxx* &
         (1.-exp(-xx*dt)))*oabi+(qvs(i,k)-qvi(i,k))*epsi*oabi
         end if

         if (ssat(i,k).lt.0.) pcc=min(0.,pcc)
         if (ssat(i,k).lt.0.) pre=min(0.,pre)
         if (qv(i,k)-qvi(i,k).lt.0.) prd1=min(0.,prd1)

! for very small water contents, evaporate instantly
         if (supi(i,k).lt.-0.001.and.qi1(i,k).lt.1.e-12) then
         prd1=-qi1(i,k)*odt
         end if

         if (prd1.lt.0.) then
            prd11 = -prd1
            prd1 = 0.
         end if

! for very small water contents, evaporate instantly
         if (sup(i,k).lt.-0.001.and.qc(i,k).lt.1.e-12) then
         pcc=-qc(i,k)*odt
         end if
         if (sup(i,k).lt.-0.001.and.qr(i,k).lt.1.e-12) then
         pre=-qr(i,k)*odt
         end if

         if (pcc.lt.0.) then
            pcc1 = -pcc
            pcc = 0.
         end if
         if (pre.lt.0.) then
            pre1 = -pre
            pre = 0.
         end if
!............................................................
! autoconversion

	if (qc(i,k).ge.1.e-8) then

! seifert and beheng (2001)

        if (iparam.eq.1) then

        dum = 1.-qc(i,k)/(qc(i,k)+qr(i,k))
        dum1 = 600.*dum**0.68*(1.-dum**0.68)**3

        prc = kc/(20.*2.6e-7)* &
       (nu(i,k)+2.)*(nu(i,k)+4.)/(nu(i,k)+1.)**2* &
      (rho(i,k)*qc(i,k)/1000.)**4/(rho(i,k)*nc(i,k)/1.e6)**2* &
      (1.+dum1/(1.-dum)**2)*1000.*orho(i,k)

        nprc = prc*2./2.6e-7*1000.

! beheng (1994)

        else if (iparam.eq.2) then

        if (nc(i,k)/1.e6*rho(i,k).lt.100.) then

	prc = 6.e28*orho(i,k)*pgam(i,k)**(-1.7)* &
                (1.e-6*rho(i,k)*nc(i,k))**(-3.3)* &
                (1.e-3*rho(i,k)*qc(i,k))**4.7
        else

! 2D interpolation of tabled logarithmic values

           dum = 41.46+ &
      (nc(i,k)/1.e6*rho(i,k)-100.)*(37.53-41.46)/200.
           dum1 = 39.36+ &
      (nc(i,k)/1.e6*rho(i,k)-100.)*(30.72-39.36)/200.

           prc = dum+(pgam(i,k)-5.)*(dum1-dum)/10.

! 1000/rho is for conversion from g cm-3/s to kg/kg

           prc = exp(prc)*(1.e-3*rho(i,k)*qc(i,k))**4.7* &
                 1000.*orho(i,k)

        end if

	nprc = 7.7e9*prc

! Khroutdinov and Kogan (2000)

        else if (iparam.eq.3) then

           dum=qc(i,k)

        prc = 1350.*dum**2.47* &
           (nc(i,k)*1.e-6*rho(i,k))**(-1.79)
! note: nprc1 is change in Nr,
! nprc is change in Nc

        nprc1 = prc*cons3
        nprc = prc*nc(i,k)/qc(i,k)

        end if

        if (prc.eq.0.) nprc = 0.
        if (nprc.eq.0.) prc = 0.

	end if
!............................................................
! self collection of droplets

	if (qc(i,k).ge.qsmall) then

        if (iparam.eq.1.) then
        ncagg = -kc*(1.e-3*rho(i,k)*qc(i,k))**2* &
        (nu(i,k)+2.)/(nu(i,k)+1.)*1.e6*orho(i,k)+nprc


        else if (iparam.eq.2.) then
	ncagg = -5.5e16*orho(i,k)*pgam(i,k)**(-0.63)* &
                  (1.e-3*rho(i,k)*qc(i,k))**2
        else if (iparam.eq.3) then
           ncagg = 0.
        end if
	end if
!............................................................
! accretion of cloud droplets by rain

	if (qr(i,k).ge.qsmall .and. qc(i,k).ge.qsmall) then

        if (iparam.eq.1) then
! seifert and beheng (2001) formulation

           dum = 1.-qc(i,k)/(qc(i,k)+qr(i,k))
           dum1 = (dum/(dum+5.e-4))**4
           pra = kr*rho(i,k)*0.001*qc(i,k)*qr(i,k)*dum1
           npra = pra*rho(i,k)*0.001*(nc(i,k)*rho(i,k)*1.e-6)/ &
           (qc(i,k)*rho(i,k)*0.001)*1.e6*orho(i,k)

        else if (iparam.eq.2) then

! beheng, 1994 formula
           dum=(qc(i,k)*qr(i,k))

           pra = 6.*rho(i,k)*dum
           npra = pra*rho(i,k)/1000.*(nc(i,k)*rho(i,k)/1.e6)/ &
           (qc(i,k)*rho(i,k)/1000.)*1.e6*orho(i,k)

        else if (iparam.eq.3) then
           dum=(qc(i,k)*qr(i,k))
           pra = 67.*(dum)**1.15
           npra = pra*nc(i,k)/qc(i,k)
        end if

        if (pra.eq.0.) npra = 0.
        if (npra.eq.0.) pra = 0.

	end if
!............................................................
! self-collection and breakup of rain
! add breakup following modified Verlinde and Cotton scheme

	if (qr(i,k).ge.qsmall) then

! include breakup
           dum1=350.e-6

! use mass-mean diameter (do this by using
! the old version of lambda w/o mu dependence)
! note there should be a factor of 6^(1/3), but we
! want to keep breakup threshold consistent so 'dum'
! is expressed in terms of lambda rather than mass-mean D

           dum2 = (QR(i,k)/(PI*RHOW*NR(i,k)))**0.33333
           if (dum2.lt.dum1) then
           dum=1.
           else if (dum2.ge.dum1) then
           dum=2.-exp(2300.*(dum2-dum1))
           end if

           if (iparam.eq.1.) then
           nragg = -dum*kr*1.e-3*qr(i,k)*nr(i,k)*rho(i,k)
           else if (iparam.eq.2.or.iparam.eq.3) then
   	   nragg = -dum*5.78*nr(i,k)*qr(i,k)*rho(i,k)
           end if

! add spontaneous breakup (i.e., independent of collection)
! treat as nudging over a 10 sec timescale
           if (dum2.gt.600.e-6) then
              dum=qr(i,k)*cons4
              dum1=(dum-nr(i,k))/max(10.,dt)
              nragg=nragg+dum1
           end if

	end if

!......................................................................
! ice processes

! collection of droplets
! here we multiply rates by air density, air density fallspeed correction
! factor, and collection efficiency since these parameters are not
! included in lookup table calculations

! for T < 273.15, assume collected cloud water is instantly frozen

        if (qi1(i,k).ge.qsmall.and.qc(i,k).ge.qsmall.and. &
        t(i,k).le.273.15) then
           psacw1=rhofaci(i,k)*f1pr4*qc(i,k)*eci*rho(i,k)
           npsacw1=rhofaci(i,k)*f1pr4*nc(i,k)*eci*rho(i,k)
        end if

! for T > 273.15, assume cloud water is collected and shed
! as rain drops

        if (qi1(i,k).ge.qsmall.and.qc(i,k).ge.qsmall.and. &
        t(i,k).gt.273.15) then
! sink for cloud water mass and number, note qcshed is source for rain mass
           qcshed1=rhofaci(i,k)*f1pr4*qc(i,k)*eci*rho(i,k)
           npsacw1=rhofaci(i,k)*f1pr4*nc(i,k)*eci*rho(i,k)
! source for rain number, assume 1mm drops are shed
           ncshed1=qcshed1*1.923e6
        end if

!............................................................
! collection of rain
! here we multiply rates by air density, air density fallspeed correction
! factor, collection efficiency, and n0r since these parameters are not
! included in lookup table calculations

! for T < 273.15, assume all collected rain mass freezes
! note this is a sink for rain mass and number and a source
! for ice mass

        if (qi1(i,k).ge.qsmall.and.qr(i,k).ge.qsmall.and. &
        t(i,k).le.273.15) then
!           pracs1=f1pr8*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
!           npracs1=f1pr7*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
! note: f1pr8 and n0r are already calculated as log_10
           pracs1=10.**(f1pr8+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
           npracs1=10.**(f1pr7+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
        end if

! for T > 273.15, assume collected rain number is shed as
! 1 mm drops
! note that melting of ice number is scaled to the loss
! rate of ice mass due to melting
! collection of rain above freezing does not impact total rain mass

        if (qi1(i,k).ge.qsmall.and.qr(i,k).ge.qsmall.and. &
        t(i,k).gt.273.15) then

! rain number sink due to collection
           npracs1=10.**(f1pr7+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
! rain number source due to shedding = collected rain mass/mass of 1 mm drop
           dum=10.**(f1pr8+n0r(i,k))*rho(i,k)*rhofaci(i,k)*eri
           nrshed1=dum*1.923e6   ! 1./5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop

! pracsw1 not currently used!!!!!!!
!           pracsw1=f1pr8*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
!           npracsw1=f1pr11*n0r(i,k)*rho(i,k)*rhofaci(i,k)*eri
        end if

!............................................................
! self-collection of ice
! here we multiply rates by collection efficiency, air density,
! and air density correction factor since these are not included
! in the lookup table calculations

        if (qi1(i,k).ge.qsmall) then
           niagg1=f1pr3*rho(i,k)*eii*rhofaci(i,k)
        end if

!............................................................
! melting
! include accelerated melting due to collection of ice mass by rain (pracsw1)
! this logic follows Lin et al. 1983 and others

        if (qi1(i,k).ge.qsmall.and.t(i,k).gt.273.15) then
           qsat0=0.622*e0/(pres(i,k)-e0)
!           dum=cpw/xlf(i,k)*(t(i,k)-273.15)*(pracsw1+qcshed1)
! currently enhanced melting from collision is neglected
!           dum=cpw/xlf(i,k)*(t(i,k)-273.15)*(pracsw1)
           dum=0.
!           qmlt1=(f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)* &
!                 (t(i,k)-273.15)*2.*pi*kap/xlf(i,k)+dum
! include RH dependence
           qmlt1=(f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)* &
                 ((t(i,k)-273.15)*kap-rho(i,k)*xxlv(i,k)*dv*(qsat0-qv(i,k)))*2.*pi/xlf(i,k)+dum
           qmlt1=max(qmlt1,0.)
!           dum1=(f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)* &
!                 (t(i,k)-273.15)*2.*pi*rho(i,k)*kap/xlf(i,k)
        end if

!............................................................
! calculate wet growth
! similar to Musil (1970), JAS

        if (qi1(i,k).ge.qsmall.and.qc(i,k)+qr(i,k).ge.1.e-6.and. &
                 t(i,k).lt.273.15) then
           qsat0=0.622*e0/(pres(i,k)-e0)
           qwgrth1=(f1pr5+f1pr14*sc**0.3333*(rhofaci(i,k)*rho(i,k)/mu)**0.5)* &
                    2.*pi*(rho(i,k)*xxlv(i,k)*dv* &
                    (qsat0-qv(i,k))-(t(i,k)-273.15)*kap)/(xlf(i,k)+cpw*(t(i,k)-273.15))
           qwgrth1=max(qwgrth1,0.)
! calculate shedding for wet growth
           dum=max(0.,(psacw1+pracs1)-qwgrth1)

	   if (dum.ge.1.e-10) then

           nrshed1=nrshed1+dum*1.923e6   ! 1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
           if ((psacw1+pracs1).ge.1.e-10) then
           dum1=1./(psacw1+pracs1)
           qcshed1=qcshed1+dum*psacw1*dum1
           psacw1=psacw1-dum*psacw1*dum1
           pracs1=pracs1-dum*pracs1*dum1
	   end if

! densify due to wet growth
           wetgrowth=1.

           end if

        end if

!............................................................
! contact and immersion freezing droplets

        if (qc(i,k).ge.qsmall.and.t(i,k).le.269.15) then
!           mnuccc=pi*pi/3.*Dap*Nacnt*rhow*cdist1(i,k)* &
!                   gamma(pgam(i,k)+5.)/lamc(i,k)**4
!           nnuccc=2.*pi*Dap*Nacnt*cdist1(i,k)* &
!                    gamma(pgam(i,k)+2.)/ &
!                    lamc(i,k)
! for future: calculate gamma(pgam+4) in one place since its used multiple times

	   dum=(1./lamc(i,k))**3
           mnucic=cons6* &
                  cdist1(i,k)*gamma(7.+pgam(i,k))* &
                   exp(aimm*(273.15-t(i,k)))*dum*dum

           nnucic=cons5*cdist1(i,k)*gamma(pgam(i,k)+4.)* &
        exp(aimm*(273.15-t(i,k)))*dum
        end if
!............................................................
! immersion freezing of rain
! for future: get rid of log statements below for rain freezing

        if (qr(i,k).ge.qsmall.and.t(i,k).le.269.15) then
           mnucir = cons6*                   &
                  EXP(LOG(CDISTR(i,k))+LOG(GAMMA(7.+MUR(i,k)))-6.*LOG(LAMR(i,k)))*             &
                  EXP(AIMM*(273.15-T(i,k)))

           nnucir =                                  &
            cons5*EXP(LOG(CDISTR(i,k))+LOG(GAMMA(MUR(i,k)+4.))-3.*LOG(LAMR(i,k)))              &
            *EXP(AIMM*(273.15-T(i,k)))

        end if

! sensitivity, no immersion freezing
!        mnucic=0.
!        nnucic=0.
!        mnucir=0.
!        nnucir=0.

!............................................................
! ***sensitivity no ice

!        mnuccc=0.
!        nnuccc=0.
!        mnucic=0.
!        nnucic=0.
!        mnucir=0.
!        nnucir=0.
!        mnuccd=0.
!        nnuccd=0.

!................................................................
! evaporate/melt number concentration

        onsubr(i,k)=0.
        if (pre1.gt.0..and.qr(i,k).ge.qsmall) then
           dum=-pre1*dt/qr(i,k)
           dum=max(-1.,dum)
           dum1=0.5
! sensitivity*******
!           dum1=exp(-0.2*mur(i,k))
           nsubr=dum1*dum*nr(i,k)
           onsubr(i,k)=dum1*dum*nr(i,k)*odt
        end if
        if (prd11.gt.0..and.qi1(i,k).ge.qsmall) then
           dum=-prd11*dt/qi1(i,k)
           dum=max(-1.,dum)
           nsubi1=dum*ni1(i,k)*odt
        end if
        if (qmlt1.gt.0..and.qi1(i,k).ge.qsmall) then
           dum=-qmlt1*dt/qi1(i,k)
           dum=max(-1.,dum)
           nmlt1=dum*ni1(i,k)*odt
        end if

 444    continue

!................................................................
! deposition/condensation-freezing nucleation
! allow ice nucleation if < -5 C and > 0.1% ice supersaturation

!        if (t(i,k).lt.268.15.and.supi(i,k).ge.0.05) then
        if (t(i,k).lt.258.15.and.supi(i,k).ge.0.05) then
!           dum=exp(-0.639+0.1296*100.*supi(i,k))*1000.*orho(i,k)
! replace w/ cooper curve
           dum=0.005*exp(0.304*(273.15-t(i,k)))*1000.*orho(i,k)
           dum=min(dum,100.e3*orho(i,k))
           if (ni1(i,k).lt.dum) then
              mnuccd=(dum-ni1(i,k))*mi0*odt
              nnuccd=(dum-ni1(i,k))*odt
           end if
        end if

!.................................................................
! droplet activation

! for specified Nc, make sure droplets are present if conditions are supersaturated

        if (sup(i,k).gt.1.e-6.or.it.eq.0) then
	   dum=nccnst*orho(i,k)*cons7-qc(i,k)
	   dum=max(0.,dum)
           dqsdt = xxlv(i,k)*qvs(i,k)/(rv*t(i,k)*t(i,k))
           ab = 1.+dqsdt*xxlv(i,k)*ocp

! don't over-deplete supersaturation
           dum=min(dum,(qv(i,k)-qvs(i,k))/ab)
	   pccn=dum*odt
	end if

        goto 713

        if (sup(i,k).gt.1.e-6.or.it.eq.0) then
           sigvl = 0.0761-1.55e-4*(t(i,k)-273.15)
           aact = 2.*mw/(rhow*rr*t(i,k))*sigvl
           sm1 = 2./bact**0.5*(aact/(3.*rm1))**1.5
           sm2 = 2./bact**0.5*(aact/(3.*rm2))**1.5
           if (it.eq.0) sup(i,k) = 0.001
           uu1 = 2.*log(sm1/sup(i,k))/(4.242*log(sig1))
!           uu2 = 2.*log(sm2/sup(i,k))/(4.242*log(sig2))
           dum1 = nanew1/2.*(1.-derf(uu1))
!           dum2 = nanew2/2.*(1.-derf(uu2))
           dum2 = dum1*orho(i,k)  !convert to kg-1
! make sure this value isn't greater than total number of aerosol
	   dum2 = min((nanew1+nanew2)*orho(i,k),dum2)
	   dum2 = (dum2-nc(i,k))*odt
	   dum2 = max(0.,dum2)
	   npccn = dum2
           if (it.eq.0) then
              pccn=0.
           else
              pccn=npccn*cons7
           end if

        end if

 713    continue

!................................................................
! saturation adjustment to get initial cloud water
! This is only called once at the beginning of the simulation
! to remove any supersaturation in the intial conditions

        if (it.eq.0) then
        dum = (100000./pres(i,k))**(rd*ocp)
        dumt = th(i,k)/dum
        dumqv = qv(i,k)
        dumqvs = ep_2*polysvp1(dumt,0)/ &
         (pres(i,k)-polysvp1(dumt,0))
        dums = dumqv-dumqvs
        pcc = dums/(1.+xxlv(i,k)**2*dumqvs/(cp*rv*dumt**2))*odt
        pcc=max(0.,pcc)
        if (pcc.le.1.e-7) pcc=0.
        end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!.................................................................
! conservation of water
!.................................................................

! cloud water
        dum = (prc+pra+psacw1+pcc1+mnuccc+mnucic+qcshed1)*dt
        dum1 = qc(i,k)+(pcc+pccn)*dt
        if (dum.gt.dum1.and.dum.ge.1.e-20) then
           ratio = dum1/dum
           prc = prc*ratio
           pra = pra*ratio
           pcc1 = pcc1*ratio
           psacw1 = psacw1*ratio
           mnuccc = mnuccc*ratio
           mnucic = mnucic*ratio
           qcshed1 = qcshed1*ratio
        end if

! rain
        dum = (pracs1+mnuccr+mnucir+pre1)*dt
        dum1 = qr(i,k)+(pre+prc+pra+qmlt1+qcshed1)*dt
        if (dum.gt.dum1.and.dum.ge.1.e-20) then
           ratio = dum1/dum
           pracs1 = pracs1*ratio
           mnuccr = mnuccr*ratio
           mnucir = mnucir*ratio
           pre1 = pre1*ratio
        end if

! ice
        dum = (prd11+qmlt1)*dt
        dum1 = qi1(i,k)+(prd1+mnuccd+ &
         pracs1+psacw1)*dt
        if (dum.gt.dum1.and.dum.ge.1.e-20) then
           ratio = dum1/dum
           prd11 = prd11*ratio
           qmlt1 = qmlt1*ratio
        end if

! rime-splintering (Hallett-Mossop 1974)
! calculate rime-splintering AFTER conservation to avoid
! unrealistic riming rates
!        if (qi1(i,k).ge.qsmall.and. &
!            (psacw1.gt.0..or.pracs1.gt.0.)) then
!           if (t(i,k).gt.270.15) then
!              dum=0.
!           else if (t(i,k).le.270.15.and.t(i,k).gt.268.15) then
!              dum=(270.15-t(i,k))/2.
!           else if (t(i,k).le.268.15.and.t(i,k).ge.265.15) then
!              dum=(t(i,k)-265.15)/3.
!           else if (t(i,k).lt.265.15) then
!              dum=0.
!           end if
!           nmult1=35.e4*(pracs1+psacw1)*dum*1000. ! 1000 is to convert kg to g
!        end if

!.................................................................................
! calculate rime density

! NOTE: Tc (ambient) is assumed for the surface temperature.  Technically,
! we should diagose graupel surface temperature from heat balance equation.
! (but the ambient temperature is a reasonable approximation; tests show
! very little sensitivity to different assumed values, Milbrandt and Morrison 2012).

! Compute rime density: (based on parameterization of Cober and List, 1993 [JAS])
! for simplicty use mass-weighted ice and droplet/rain fallspeeds

! initialize to 400 kg/m3
      rhorime_c=400.

      if (qi1(i,k).ge.qsmall.and.t(i,k).lt.273.15) then

! get mass-weighted mean ice fallspeed
        vtrmi1(i,k)=f1pr2*rhofaci(i,k)

        iTc   = 1./min(-0.001,t(i,k)-273.15)

! First for cloud droplets
          if (qc(i,k).ge.qsmall) then

!...................................................................
! droplet fall speed
! use Stokes' formulation
! all droplets in smallest dactegory fallspeed
! this means we can use analytic solution
          vtrmc(i,k) = acn(i,k)*gamma(4.+bcn+pgam(i,k))/ &
             (lamc(i,k)**bcn*gamma(pgam(i,k)+4.))

! use mass-weighted mean size
          D_c = (pgam(i,k)+4.)/lamc(i,k)
          V_impact  = abs(vtrmi1(i,k)-vtrmc(i,k))
          Ri        = -(0.5e+6*D_c)*V_impact*iTc
          Ri        = max(1.,min(Ri,8.))
          rhorime_c  = (0.051 + 0.114*Ri - 0.0055*Ri*Ri)*1000.
          end if

	  end if ! qi > qsmall and T < 273.15

! Next for rain drops
! assume rime density for rain collecting ice is 900 kg/m3
!          if (qr(i,k).ge.qsmall) then
!          D_r = (mur(i,k)+1.)/lamr(i,k)
!          V_impact  = abs(vtrmi1(i,k)-vtrm(i,k))
!          Ri        = -(0.5e+6*D_r)*V_impact*iTc
!          Ri        = max(1.,min(Ri,8.))
!          rhorime_r  = (0.051 + 0.114*Ri - 0.0055*Ri*Ri)*1000.
!          else
!             rhorime_r = 400.
!          end if

!........................................................................
! update microphysics, thermodynamic variables

        qc(i,k) = qc(i,k)+(-pra-prc+pccn+pcc-pcc1- &
      mnuccc-mnucic-psacw1-qcshed1)*dt
        qr(i,k) = qr(i,k)+(pra+prc+pre-pre1-pracs1+ &
      qmlt1-mnuccr-mnucir+qcshed1)*dt

! sink terms for ice (sublimation and melting)

        if (qi1(i,k).ge.qsmall) then

! add sink terms, assume density stays constant for sink terms
           bg1(i,k)=bg1(i,k)- &
            ((prd11+qmlt1)/qi1(i,k))*dt*bg1(i,k)
           qri1(i,k)=qri1(i,k)- &
            ((prd11+qmlt1)*qri1(i,k)/qi1(i,k))*dt
           qi1(i,k)=qi1(i,k)- &
           (prd11+qmlt1)*dt

        end if

! total ice mass mixing ratio growth rate
        dum=(pracs1+psacw1+mnuccr+mnucir+mnuccc+mnucic)*dt
        qi1(i,k)=qi1(i,k)+(prd1+mnuccd)*dt+dum

! riming growth rate
        qri1(i,k)=qri1(i,k)+dum

! for freezing of liquid drops, assume 900 kg m-3

!        if (pracs1.gt.1.e-20)print*,'pracs',rhorime_r
!        if (psacw1.gt.1.e-20)print*,'psacw',rhorime_c

! update bulk volume mixing ratio
! this assumes that vapor growth of graupel doesn't change the particle density
        bg1(i,k)=bg1(i,k)+(pracs1/900.+ &
          psacw1/rhorime_c+ &
          (mnuccr+mnucir+mnuccc+mnucic)/900.)*dt

        if (qri1(i,k).lt.0.) then
           qri1(i,k)=0.
           bg1(i,k)=0.
        end if

! densify under wet growth
        if (wetgrowth.gt.0.5) then
           bg1(i,k)=qri1(i,k)/900.
        end if

        qv(i,k) = qv(i,k)+(-pccn-pcc-pre-prd1 &
           +pcc1+pre1+prd11-mnuccd)*dt
        th(i,k) = th(i,k)+th(i,k)/t(i,k)*( &
            (pccn+pcc+pre-pcc1-pre1)*xxlv(i,k)*ocp &
        +(prd1-prd11+mnuccd)*xxls(i,k)*ocp &
        +(pracs1+psacw1+mnuccc+mnucic+ &
      mnuccr+mnucir-qmlt1)*xlf(i,k)*ocp)*dt

! number concentration variables

!        nc(i,k) = nc(i,k)+(-npra-nprc+ncagg+npccn- &
!      npsacw1-nnuccc-nnucic)*dt

! specify constant nc
        nc(i,k)=nccnst*orho(i,k)
        if (iparam.eq.1.or.iparam.eq.2) then
        nr(i,k) = nr(i,k)+(0.5*nprc+nragg &
      -npracs1-nnuccr-nnucir-nmlt1+nsubr+nrshed1+ncshed1)*dt
        else
        nr(i,k) = nr(i,k)+(nprc1+nragg &
      -npracs1-nnuccr-nnucir-nmlt1+nsubr+nrshed1+ncshed1)*dt
        end if

        ni1(i,k)=ni1(i,k)+(nnuccd+nmlt1 &
       +nsubi1-niagg1+nmult1+nnuccr+nnucir+ &
        nnuccc+nnucic)*dt

! store processes for output

        revap(i,k)=pre1

!        rcon(i,k)=pre
!        ccon(i,k)=pcc
!        auto(i,k)=prc
!        acc(i,k)=pra
!        act(i,k)=npccn
!        opre(i,k)=pre
!        opra(i,k)=pra
!        oprc(i,k)=prc
!        onpra(i,k)=npra
!        onprc(i,k)=nprc
!        oncagg(i,k)=ncagg
!        onragg(i,k)=nragg
!        onpccn(i,k)=npccn
!        opcc(i,k)=pcc
!        opccn(i,k)=pccn
!        onprc1(i,k)=nprc1

            if (qc(i,k).lt.qsmall) then
               qv(i,k) = qv(i,k)+qc(i,k)
               th(i,k) = th(i,k)+th(i,k)/t(i,k)*qc(i,k)*xxlv(i,k)*ocp
               qc(i,k) = 0.
               nc(i,k) = 0.
            else
               ltrue=1
            end if

            if (qr(i,k).lt.qsmall) then
               qv(i,k) = qv(i,k)+qr(i,k)
               th(i,k) = th(i,k)+th(i,k)/t(i,k)*qr(i,k)*xxls(i,k)*ocp
               qr(i,k) = 0.
               nr(i,k) = 0.
            else
               ltrue=1
            end if

            if (qi1(i,k).lt.qsmall) then
               qv(i,k)=qv(i,k)+qi1(i,k)
               th(i,k)=th(i,k)+th(i,k)/t(i,k)*qi1(i,k)*xxls(i,k)*ocp
               qi1(i,k)=0.
               ni1(i,k)=0.
               qri1(i,k)=0.
               bg1(i,k)=0.
            else
               ltrue=1
            end if

 555        continue

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! end main k loop

      end do

      if (ltrue.eq.0) goto 333

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! initialize logicals for presence of hydrometeor species to .false.
  qcpresent=.false.
  qrpresent=.false.
  qipresent=.false.

!     do k = kte,kts,-1
      do k = ktop,kbot,-kdir

! define 1/dzq
       odzq(i,k)=1./dzq(i,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate number and mass weighted fallspeeds and find highest k level
! that a given species is present

!..................................................
! get droplet size distribution parameters

	if (qc(i,k).ge.qsmall) then

! set minimum nc to prevent floating point error
        nc(i,k)=max(nc(i,k),1.e-16)

        pgam(i,k)=0.0005714*(nc(i,k)*1.e-6*rho(i,k))+0.2714
!         pgam(i,k)=0.146-5.964e-2*log(nc(i,k)/1.e6*rho(i,k)/2000.)
        pgam(i,k)=1./(pgam(i,k)**2)-1.
        pgam(i,k)=max(pgam(i,k),2.)
        pgam(i,k)=min(pgam(i,k),15.)

! calculate lamc

        lamc(i,k)=(cons1*nc(i,k)*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.)/ &
            (qc(i,k)))**0.33333

! lammin, 40 micron mean diameter
! lammax, 1 mean micron

	lammin = (pgam(i,k)+1.)/40.e-6
	lammax = (pgam(i,k)+1.)/1.e-6

	if (lamc(i,k).lt.lammin) then
	lamc(i,k) = lammin

	nc(i,k) = 6.*lamc(i,k)**3*qc(i,k)/ &
               (pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.))

	else if (lamc(i,k).gt.lammax) then
	lamc(i,k) = lammax

	nc(i,k) = 6.*lamc(i,k)**3*qc(i,k)/ &
               (pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.))

	if (.not. qcpresent) then
	qcindex=k
	end if
	qcpresent=.true.

	end if

	else
	lamc(i,k) = 0.
	end if

!..................................................
! get rain size distribution parameters

	if (qr(i,k).ge.qsmall) then

! use lookup table to get mu
! mu-lambda relationship is from Cao et al. (2008), eq. (7)

! find spot in lookup table
! scaled N/q for lookup table parameter space
         nr(i,k)=max(nr(i,k),1.e-16)
         dum = (cons1*NR(i,k)*6./   &
                 (QR(i,k)))**0.33333
         odum=1./dum

         if (odum.lt.282.e-6) then
            mur(i,k) = 8.282
         else if (odum.ge.282.e-6.and.odum.lt.502.e-6) then
! interpolate
            rdumii = (odum-250.e-6)*1.e6/2.
            rdumii = max(rdumii,1.)
            rdumii = min(rdumii,150.)
            dumii = int(rdumii)
            dumii = min(149,dumii)
            mur(i,k) = mur_table(dumii)+(mur_table(dumii+1)-mur_table(dumii))* &
                       (rdumii-real(dumii))
         else if (odum.ge.502.e-6) then
            mur(i,k)=0.
         end if

! recalculate slope based on mur
         LAMR(i,k) = (cons1*nr(i,k)*(mur(i,k)+3.)*(mur(i,k)+2)*(mur(i,k)+1.)/   &
                 (qr(i,k)))**0.33333

! check for slope
        lammax = (mur(i,k)+1.)*1.e5

! set to small value since breakup is explicitly included (mean size 0.8 mm)
        lammin = (mur(i,k)+1.)*1250.

! apply lambda limiters for rain
        IF (lamr(i,k).LT.LAMMIN) THEN

        LAMR(i,k) = LAMMIN

        NR(i,k) = EXP(3.*LOG(LAMR(i,k))+LOG(QR(i,k))+              &
                LOG(GAMMA(MUR(i,k)+1.))-LOG(GAMMA(MUR(i,k)+4.)))/(cons1)

        ELSE IF (LAMR(i,k).GT.LAMMAX) THEN

        LAMR(i,k) = LAMMAX

        NR(i,k) = EXP(3.*LOG(LAMR(i,k))+LOG(QR(i,k))+              &
                LOG(GAMMA(MUR(i,k)+1.))-LOG(GAMMA(MUR(i,k)+4.)))/(cons1)

        END IF

	if (.not. qrpresent) then
	qrindex=k
	end if
	qrpresent=.true.

        else
        lamr(i,k) = 0.
        end if

      if (qr(i,k).ge.qsmall) then
!...................................................................
! read in fall mass and number weighted fall velocity from table

! find location in mean size space
        dum1=(mur(i,k)+1.)/lamr(i,k)
        if (dum1.le.195.e-6) then
           rdumii=(dum1*1.e6+5.)/10.
           rdumii=max(rdumii,1.)
           rdumii=min(rdumii,20.)
           dumii=int(rdumii)
           dumii=max(dumii,1)
           dumii=min(dumii,20)
           dum3=10.
        else if (dum1.gt.195.e-6) then
           rdumii=(dum1*1.e6-195.)/30.+20.
           rdumii=max(rdumii,20.)
           rdumii=min(rdumii,300.)
           dumii=int(rdumii)
           dumii=max(dumii,20)
           dumii=min(dumii,299)
           dum3=30.
        end if

! find location in mur space

        rdumjj=mur(i,k)+1.
        rdumjj=max(rdumjj,1.)
        rdumjj=min(rdumjj,10.)
        dumjj=int(rdumjj)
        dumjj=max(dumjj,1)
        dumjj=min(dumjj,9)

	odum3=1./dum3

! number-weighted fallspeed
! value at mur
        dum1=vn_table(dumii,dumjj)+(rdumii-real(dumii))*odum3*  &
        (vn_table(dumii+1,dumjj)-vn_table(dumii+1,dumjj))

! value at mur+1
        dum2=vn_table(dumii,dumjj+1)+(rdumii-real(dumii))*odum3*  &
        (vn_table(dumii+1,dumjj+1)-vn_table(dumii+1,dumjj+1))

        vtrn(i,k)=dum1+(rdumjj-real(dumjj))* &
              (dum2-dum1)

! mass-weighted fallspeed
! value at mur
        dum1=vm_table(dumii,dumjj)+(rdumii-real(dumii))*odum3*  &
        (vm_table(dumii+1,dumjj)-vm_table(dumii+1,dumjj))

! value at mur+1
        dum2=vm_table(dumii,dumjj+1)+(rdumii-real(dumii))*odum3*  &
        (vm_table(dumii+1,dumjj+1)-vm_table(dumii+1,dumjj+1))

        vtrm(i,k)=dum1+(rdumjj-real(dumjj))* &
              (dum2-dum1)

        vtrn(i,k)=vtrn(i,k)*rhofacr(i,k)
        vtrm(i,k)=vtrm(i,k)*rhofacr(i,k)

      else
         vtrn(i,k) = 0.
         vtrm(i,k) = 0.
      end if

! droplet fall speed
! all droplets in smallest dactegory fallspeed
! this means we can use analytic solution

        if (qc(i,k).ge.qsmall) then
	dum=1./lamc(i,k)**bcn
!        vtrnc(i,k) =  acn(i,k)*gamma(1.+bcn+pgam(i,k))*dum/ &
!               (gamma(pgam(i,k)+1.))
        vtrmc(i,k) = acn(i,k)*gamma(4.+bcn+pgam(i,k))*dum/ &
             (gamma(pgam(i,k)+4.))
        else
!           vtrnc(i,k)=0.
           vtrmc(i,k)=0.
        end if

!.......................................................................
! get ice fallspeed for updated variables

         if (qi1(i,k).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
         ni1(i,k)=max(ni1(i,k),1.e-16)

! calculate predicted bulk rime density
         if (bg1(i,k).ge.1.e-15) then
         rhop=qri1(i,k)/bg1(i,k)
         else
         rhop=0.
         end if
! limit 50 < rhop < 900, adjust bg if needed
         if (rhop.lt.50.) then
            rhop=50.
            bg1(i,k)=qri1(i,k)/rhop
         end if
         if (rhop.gt.900.) then
            rhop=900.
            bg1(i,k)=qri1(i,k)/rhop
         end if
         if (qri1(i,k).lt.qsmall) then
            bg1(i,k)=0.
         end if

! set upper constraint on qri to ensure qri cannot be > qi
	 if (qri1(i,k).gt.qi1(i,k)) then
         qri1(i,k)=qi1(i,k)
	 bg1(i,k)=qri1(i,k)/rhop
	 end if

! find indices in 4D ice lookup table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

! find index for qi (total ice mass mixing ratio)
         dum1 = (alog10(qi1(i,k))+16.)/0.70757
         dumi=int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
         dum1=min(dum1,real(isize))
         dum1=max(dum1,1.)
         dumi=max(1,dumi)
         dumi=min(isize-1,dumi)

! find index for Ni (ice number mixing ratio)
         dum2 = (alog10(ni1(i,k))+10.)/0.90309
         dumk=int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
         dum2=min(dum2,real(jsize))
         dum2=max(dum2,1.)
         dumk=max(1,dumk)
         dumk=min(jsize-1,dumk)

! find index for rime mass fraction
         dum4 = qri1(i,k)/qi1(i,k)*3.+1.
         dumii=int(dum4)

! set limits
         dum4=min(dum4,real(rimsize))
         dum4=max(dum4,1.)
         dumii=max(1,dumii)
         dumii=min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
         if (rhop.le.650.) then
         dum5 = (rhop-50.)/200.+1.
         else
         dum5=(rhop-650.)/250.+4.
         end if
         dumjj=int(dum5)

! set limits
         dum5=min(dum5,real(densize))
         dum5=max(dum5,1.)
         dumjj=max(1,dumjj)
         dumjj=min(densize-1,dumjj)

! call subroutine to interpolate ice lookup table values

         index=1
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr1)
         index=2
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr2)
         index=7
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr9)
         index=8
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr10)

!...................................................

! make sure mean ice size is in bounds (i.e., apply lambda limiters)
            ni1(i,k)=min(ni1(i,k),f1pr9)
            ni1(i,k)=max(ni1(i,k),f1pr10)

! limit max ice concentration to 500 L-1
            ni1(i,k)=min(ni1(i,k),500.e3*orho(i,k))

  	   if (.not. qipresent) then
	   qiindex=k
	   end if
	   qipresent=.true.

         end if ! qi1 < qsmall

!........................................................................

      if (qi1(i,k).ge.qsmall) then
         vtrni1(i,k)=f1pr1*rhofaci(i,k)
         vtrmi1(i,k)=f1pr2*rhofaci(i,k)
! output fallspeed, w/o density correction
	 vmi(i,k)=f1pr2
      else
         vtrni1(i,k)=0.
         vtrmi1(i,k)=0.
      end if

      end do ! main k loop

!.....................................................
! cloud droplet sedimentation

     if (qcpresent) then

      nstep = 1
!     do k=qcindex,kts,-1
      do k=qcindex,kbot,-kdir
        FC(K) = vtrmc(i,k)
!        FNC(K) = vtrnc(i,k)
        IF (K.LE.qcindex-kdir) THEN
        IF (FC(K).LT.1.E-10) THEN
        FC(K)=FC(K+kdir)
        END IF
!        IF (FNC(K).LT.1.E-10) THEN
!        FNC(K)=FNC(K+kdir)
!        END IF
        END IF

! CALCULATE NUMBER OF SPLIT TIME STEPS
      RGVM = FC(K)
      NSTEP = MAX(INT(RGVM*DT*odzq(i,k)+1.),NSTEP)
      onstep=1./real(nstep)

      DUMC(k) = qc(i,k)*RHO(i,K)
!      DUMFNC(k) = nc(i,k)*RHO(i,K)
      QCSTEN(K) = 0.
      NCSTEN(K) = 0.

      end do ! k loop

! calculate sedimentation using first-order upwind method

      DO N = 1,NSTEP

!     do k = kts,qcindex
      do k = kbot,qcindex,kdir
      FALOUTC(K) = FC(K)*DUMC(K)
!      FALOUTNC(K) = FNC(K)*DUMFNC(K)
      end do

! top level with hydrometeor present

      k=qcindex
      FALTNDC = FALOUTC(K)*ODZQ(i,k)
!      FALTNDNC = FALOUTNC(K)*ODZQ(i,k)
      QCSTEN(K) = QCSTEN(K)-FALTNDC*ONSTEP*ORHO(i,k)
!      NCSTEN(K) = NCSTEN(K)-FALTNDNC*ONSTEP*ORHO(i,k)
      DUMC(K) = DUMC(K)-FALTNDC*DT*ONSTEP
!      DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT*ONSTEP

! loop from sceond to top level of hydrometeor to surface

!     DO K = qcindex-1,KTS,-1
      DO K = qcindex-kdir,kbot,-kdir
      FALTNDC = (FALOUTC(K+kdir)-FALOUTC(K))*ODZQ(i,K)
!      FALTNDNC = (FALOUTNC(K+kdir)-FALOUTNC(K))*ODZQ(i,K)
      QCSTEN(K) = QCSTEN(K)+FALTNDC*ONSTEP*ORHO(i,k)
!      NCSTEN(K) = NCSTEN(K)+FALTNDNC*ONSTEP*ORHO(i,k)
      DUMC(K) = DUMC(K)+FALTNDC*DT*ONSTEP
!      DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT*ONSTEP
       prec(i,k)=prec(i,k)+(FALOUTC(k))*ONSTEP*3600. ! convert to mm/h
      end do ! k loop

!       PRECPRT(i) = PRECPRT(i)+(FALOUTC(KTS))*DT*ONSTEP
        PRECPRT(i) = PRECPRT(i)+(FALOUTC(kbot))*DT*ONSTEP

      END DO ! nstep loop

! update prognostic variables with sedimentation tendencies
!     do k=kts,qcindex
      do k=kbot,qcindex,kdir
      qc(i,k)=qc(i,k)+QCSTEN(K)*dt
!      nc(i,k)=nc(i,k)+NCSTEN(K)*dt
      end do

      end if ! qcpresent

!.....................................................
! rain sedimentation

     if (qrpresent) then

      nstep = 1
!     do k=qrindex,kts,-1
      do k=qrindex,kbot,-kdir

        FR(K) = vtrm(i,k)
        FNR(K) = vtrn(i,k)
        IF (K.LE.qrindex-kdir) THEN
        IF (FR(K).LT.1.E-10) THEN
        FR(K)=FR(K+kdir)
        END IF
        IF (FNR(K).LT.1.E-10) THEN
        FNR(K)=FNR(K+kdir)
        END IF
        END IF

! CALCULATE NUMBER OF SPLIT TIME STEPS
      RGVM = MAX(FR(K),FNR(K))
      NSTEP = MAX(INT(RGVM*DT*ODZQ(i,K)+1.),NSTEP)
      onstep=1./real(nstep)

      DUMR(k) = qr(i,k)*RHO(i,K)
      DUMFNR(k) = nr(i,k)*RHO(i,K)
      QRSTEN(K) = 0.
      NRSTEN(K) = 0.

      end do ! k loop

! calculate sedimentation using first-order upwind method

      DO N = 1,NSTEP

!     do k = kts,qrindex
      do k = kbot,qrindex,kdir
      FALOUTR(K) = FR(K)*DUMR(K)
      FALOUTNR(K) = FNR(K)*DUMFNR(K)
      end do

! top level with hydrometeor present

      k=qrindex
      FALTNDR = FALOUTR(K)*ODZQ(i,k)
      FALTNDNR = FALOUTNR(K)*ODZQ(i,k)
      QRSTEN(K) = QRSTEN(K)-FALTNDR*ONSTEP*ORHO(i,k)
      NRSTEN(K) = NRSTEN(K)-FALTNDNR*ONSTEP*ORHO(i,k)
      DUMR(K) = DUMR(K)-FALTNDR*DT*ONSTEP
      DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT*ONSTEP

! loop from second to top level of hydrometeor to surface

!     DO K = qrindex-1,KTS,-1
      DO K = qrindex-kdir,kbot,-kdir
      FALTNDR   = (FALOUTR(K+kdir)-FALOUTR(K))*ODZQ(i,K)
      FALTNDNR  = (FALOUTNR(K+kdir)-FALOUTNR(K))*ODZQ(i,K)
      QRSTEN(K) = QRSTEN(K)+FALTNDR*ONSTEP*ORHO(i,k)
      NRSTEN(K) = NRSTEN(K)+FALTNDNR*ONSTEP*ORHO(i,k)
      DUMR(K)   = DUMR(K)+FALTNDR*DT*ONSTEP
      DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT*ONSTEP
      prec(i,k)=prec(i,k)+(FALOUTR(k))*ONSTEP*3600. ! convert to mm/h
      end do ! k loop

!       PRECPRT(i) = PRECPRT(i)+(FALOUTR(KTS))*DT*ONSTEP
        PRECPRT(i) = PRECPRT(i)+(FALOUTR(kbot))*DT*ONSTEP

      END DO ! nstep loop

! update prognostic variables with sedimentation tendencies
!     do k=kts,qrindex
      do k=kbot,qrindex,kdir
      qr(i,k)=qr(i,k)+QRSTEN(K)*dt
      nr(i,k)=nr(i,k)+NRSTEN(K)*dt
      end do

      end if ! qrpresent

!.....................................................
! ice sedimentation

     if (qipresent) then

      nstep = 1
!     do k=qiindex,kts,-1
      do k=qiindex,kbot,-kdir
        FI(K) = vtrmi1(i,k)
        FNI(K) = vtrni1(i,k)
        IF (K.LE.qiindex-kdir) THEN
        IF (FI(K).LT.1.E-10) THEN
        FI(K)=FI(K+kdir)
        END IF
        IF (FNI(K).LT.1.E-10) THEN
        FNI(K)=FNI(K+kdir)
        END IF
        END IF

! CALCULATE NUMBER OF SPLIT TIME STEPS
      RGVM = MAX(FI(K),FNI(K))
      NSTEP = MAX(INT(RGVM*DT*ODZQ(i,K)+1.),NSTEP)
      onstep=1./real(nstep)

      DUMQI(k) = qi1(i,k)*RHO(i,K)
      DUMRI(k) = qri1(i,k)*RHO(i,K)
      DUMBG(k) = bg1(i,k)*RHO(i,K)
      DUMFNI(k) = ni1(i,k)*RHO(i,K)

      QISTEN(K) = 0.
      QRISTEN(K) = 0.
      BGSTEN(K) = 0.
      NISTEN(K) = 0.

      end do ! k loop

! calculate sedimentation using first-order upwind method

      DO N = 1,NSTEP

!     do k = kts,qiindex
      do k = kbot,qiindex,kdir
      FALOUTI(K) = FI(K)*DUMQI(K)
      FALOUTNI(K) = FNI(K)*DUMFNI(K)
      FALOUTRI(K) = FI(K)*DUMRI(K)
      FALOUTBG(K) = FI(K)*DUMBG(K)
      end do

! top level with hydrometeor present

      k=qiindex

      FALTNDI = FALOUTI(K)*ODZQ(i,k)
      FALTNDRI = FALOUTRI(K)*ODZQ(i,k)
      FALTNDBG = FALOUTBG(K)*ODZQ(i,k)
      FALTNDNI = FALOUTNI(K)*ODZQ(i,k)
      QISTEN(K) = QISTEN(K)-FALTNDI*ONSTEP*ORHO(i,k)
      QRISTEN(K) = QRISTEN(K)-FALTNDRI*ONSTEP*ORHO(i,k)
      BGSTEN(K) = BGSTEN(K)-FALTNDBG*ONSTEP*ORHO(i,k)
      NISTEN(K) = NISTEN(K)-FALTNDNI*ONSTEP*ORHO(i,k)
      DUMQI(K) = DUMQI(K)-FALTNDI*DT*ONSTEP
      DUMRI(K) = DUMRI(K)-FALTNDRI*DT*ONSTEP
      DUMBG(K) = DUMBG(K)-FALTNDBG*DT*ONSTEP
      DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT*ONSTEP

! loop from sceond to top level of hydrometeor to surface

!     DO K = qiindex-1,KTS,-1
      DO K = qiindex-kdir,kbot,-kdir

      FALTNDI    = (FALOUTI(K+kdir)-FALOUTI(K))*ODZQ(i,K)
      FALTNDRI   = (FALOUTRI(K+kdir)-FALOUTRI(K))*ODZQ(i,K)
      FALTNDBG   = (FALOUTBG(K+kdir)-FALOUTBG(K))*ODZQ(i,K)
      FALTNDNI   = (FALOUTNI(K+kdir)-FALOUTNI(K))*ODZQ(i,K)
      QISTEN(K)  = QISTEN(K)+FALTNDI*ONSTEP*ORHO(i,k)
      QRISTEN(K) = QRISTEN(K)+FALTNDRI*ONSTEP*ORHO(i,k)
      BGSTEN(K)  = BGSTEN(K)+FALTNDBG*ONSTEP*ORHO(i,k)
      NISTEN(K)  = NISTEN(K)+FALTNDNI*ONSTEP*ORHO(i,k)
      DUMQI(K)   = DUMQI(K)+FALTNDI*DT*ONSTEP
      DUMRI(K)   = DUMRI(K)+FALTNDRI*DT*ONSTEP
      DUMBG(K)   = DUMBG(K)+FALTNDBG*DT*ONSTEP
      DUMFNI(K)  = DUMFNI(K)+FALTNDNI*DT*ONSTEP
      prec(i,k)=prec(i,k)+(FALOUTI(k))*ONSTEP*3600. ! convert to mm/h

      end do ! k loop

!       PRECPRT(i) = PRECPRT(i)+(FALOUTI(KTS)+FALOUTRI(KTS))*DT*ONSTEP
!       SNOWRT(i) = SNOWRT(i)+(FALOUTI(KTS)+FALOUTRI(KTS))*DT*ONSTEP
        PRECPRT(i) = PRECPRT(i) + FALOUTI(kbot)*DT*ONSTEP
        SNOWRT(i)  = SNOWRT(i)  + FALOUTI(kbot)*DT*ONSTEP

      END DO ! nstep loop

! update prognostic variables with sedimentation tendencies
!     do k=kts,qiindex
      do k=kbot,qiindex,kdir
      qi1(i,k)=qi1(i,k)+QISTEN(K)*dt
      qri1(i,k)=qri1(i,k)+QRISTEN(K)*dt
      bg1(i,k)=bg1(i,k)+BGSTEN(K)*dt
      ni1(i,k)=ni1(i,k)+NISTEN(K)*dt
      end do

      end if ! qipresent

!end sedimentation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!.............................
! homogeneous freezing of cloud droplets and rain
!     do k=kts,kte
      do k=kbot,ktop,kdir

      if (qc(i,k).ge.qsmall.and.t(i,k).lt.233.15) then
         qri1(i,k)=qri1(i,k)+qc(i,k)
	 qi1(i,k)=qi1(i,k)+qc(i,k)
         bg1(i,k)=bg1(i,k)+qc(i,k)/900.
         ni1(i,k)=ni1(i,k)+nc(i,k)
         th(i,k)=th(i,k)+th(i,k)/t(i,k)*qc(i,k)*xlf(i,k)*ocp
         qc(i,k)=0.
         nc(i,k)=0.
      end if

      if (qr(i,k).ge.qsmall.and.t(i,k).lt.233.15) then
         qri1(i,k)=qri1(i,k)+qr(i,k)
	 qi1(i,k)=qi1(i,k)+qr(i,k)
         bg1(i,k)=bg1(i,k)+qr(i,k)/900.
         ni1(i,k)=ni1(i,k)+nr(i,k)
         th(i,k)=th(i,k)+th(i,k)/t(i,k)*qr(i,k)*xlf(i,k)*ocp
         qr(i,k)=0.
         nr(i,k)=0.
      end if

      end do ! k loop

!.......................................................................
! calculate final size distribution parameters/effective radii
! also final checks to ensure consistency of mass/number

! get size distribution parameters

!       do k=kts,kte
        do k=kbot,ktop,kdir

            if (qc(i,k).lt.qsmall) then
               qv(i,k) = qv(i,k)+qc(i,k)
               th(i,k) = th(i,k)+th(i,k)/t(i,k)*qc(i,k)*xxlv(i,k)*ocp
               qc(i,k) = 0.
               nc(i,k) = 0.
            end if

            if (qr(i,k).lt.qsmall) then
               qv(i,k) = qv(i,k)+qr(i,k)
               th(i,k) = th(i,k)+th(i,k)/t(i,k)*qr(i,k)*xxls(i,k)*ocp
               qr(i,k) = 0.
               nr(i,k) = 0.
            end if

!..................................................
! get droplet size distribution parameters

	if (qc(i,k).ge.qsmall) then

! set minimum nc to prevent floating point error
        nc(i,k)=max(nc(i,k),1.e-16)

        pgam(i,k)=0.0005714*(nc(i,k)*1.e-6*rho(i,k))+0.2714
!         pgam(i,k)=0.146-5.964e-2*log(nc(i,k)/1.e6*rho(i,k)/2000.)
        pgam(i,k)=1./(pgam(i,k)**2)-1.
        pgam(i,k)=max(pgam(i,k),2.)
        pgam(i,k)=min(pgam(i,k),15.)

! calculate lamc

        lamc(i,k)=(cons1*nc(i,k)*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.)/ &
            (qc(i,k)))**0.33333

! lammin, 40 micron mean diameter
! lammax, 1 mean micron

	lammin = (pgam(i,k)+1.)/40.e-6
	lammax = (pgam(i,k)+1.)/1.e-6

	if (lamc(i,k).lt.lammin) then
	lamc(i,k) = lammin

	nc(i,k) = 6.*lamc(i,k)**3*qc(i,k)/ &
               (pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.))

	else if (lamc(i,k).gt.lammax) then
	lamc(i,k) = lammax

	nc(i,k) = 6.*lamc(i,k)**3*qc(i,k)/ &
               (pi*rhow*(pgam(i,k)+3.)*(pgam(i,k)+2.)*(pgam(i,k)+1.))

	end if

	else
	lamc(i,k) = 0.
	end if

!..................................................
! get rain size distribution parameters

	if (qr(i,k).ge.qsmall) then

! use lookup table to get mu
! mu-lambda relationship is from Cao et al. (2008), eq. (7)

! find spot in lookup table
! scaled N/q for lookup table parameter space
         nr(i,k)=max(nr(i,k),1.e-16)
         dum = (cons1*NR(i,k)*6./   &
                 (QR(i,k)))**0.33333
         odum=1./dum

         if (odum.lt.282.e-6) then
            mur(i,k) = 8.282
         else if (odum.ge.282.e-6.and.odum.lt.502.e-6) then
! interpolate
            rdumii = (odum-250.e-6)*1.e6/2.
            rdumii = max(rdumii,1.)
            rdumii = min(rdumii,150.)
            dumii = int(rdumii)
            dumii = min(149,dumii)
            mur(i,k) = mur_table(dumii)+(mur_table(dumii+1)-mur_table(dumii))* &
                       (rdumii-real(dumii))
         else if (odum.ge.502.e-6) then
            mur(i,k)=0.
         end if

! recalculate slope based on mur
         LAMR(i,k) = (cons1*nr(i,k)*(mur(i,k)+3.)*(mur(i,k)+2)*(mur(i,k)+1.)/   &
                 (qr(i,k)))**0.33333

! check for slope

        lammax = (mur(i,k)+1.)*1.e5

! set to small value since breakup is explicitly included (mean size 0.8 mm)
        lammin = (mur(i,k)+1.)*1250.

! adjust vars

        IF (lamr(i,k).LT.LAMMIN) THEN

        LAMR(i,k) = LAMMIN

        NR(i,k) = EXP(3.*LOG(LAMR(i,k))+LOG(QR(i,k))+              &
                LOG(GAMMA(MUR(i,k)+1.))-LOG(GAMMA(MUR(i,k)+4.)))/(cons1)

        ELSE IF (LAMR(i,k).GT.LAMMAX) THEN

        LAMR(i,k) = LAMMAX

        NR(i,k) = EXP(3.*LOG(LAMR(i,k))+LOG(QR(i,k))+              &
                LOG(GAMMA(MUR(i,k)+1.))-LOG(GAMMA(MUR(i,k)+4.)))/(cons1)

        END IF

	else
	lamr(i,k) = 0.
	end if

! calculate effective radius for cloud droplets and rain

	if (qr(i,k).ge.qsmall) then
	   effr(i,k) = 1.5/lamr(i,k)*1.e6
	else
	   effr(i,k) = 25.
	end if

	if (qc(i,k).ge.qsmall) then
	effc(i,k) = (pgam(i,k)+3.)/lamc(i,k)*0.5
	else
	effc(i,k) = 25.e-6
	end if

!.......................................................................
! size distribution parameters for ice

! calculate total mixing ratio from deposition and riming mixing ratios

         if (qi1(i,k).lt.qsmall) then
            qv(i,k)=qv(i,k)+qi1(i,k)
            th(i,k)=th(i,k)+th(i,k)/t(i,k)*qi1(i,k)*xxls(i,k)*ocp
            qi1(i,k)=0.
            ni1(i,k)=0.
            qri1(i,k)=0.
            bg1(i,k)=0.
         end if

         if (qi1(i,k).ge.qsmall) then

! set lower limit on ni to prevent taking log of # < 0
         ni1(i,k)=max(ni1(i,k),1.e-16)

! calculate predicted bulk rime density
         if (bg1(i,k).ge.1.e-15) then
         rhop=qri1(i,k)/bg1(i,k)
         else
         rhop=0.
         end if
! limit 50 < rhop < 900, adjust bg if needed
         if (rhop.lt.50.) then
            rhop=50.
            bg1(i,k)=qri1(i,k)/rhop
         end if
         if (rhop.gt.900.) then
            rhop=900.
            bg1(i,k)=qri1(i,k)/rhop
         end if
         if (qri1(i,k).lt.qsmall) then
            bg1(i,k)=0.
         end if

! set upper constraint on qri to ensure qri cannot be > qi
         if (qri1(i,k).gt.qi1(i,k)) then
         qri1(i,k)=qi1(i,k)
	 bg1(i,k)=qri1(i,k)/rhop
	 end if

! find indices in 4D ice lookup table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

! find index for qi (total ice mass mixing ratio)
         dum1 = (alog10(qi1(i,k))+16.)/0.70757
         dumi=int(dum1)

! set limits to make sure the calculated index doesn't exceed range of lookup table
         dum1=min(dum1,real(isize))
         dum1=max(dum1,1.)
         dumi=max(1,dumi)
         dumi=min(isize-1,dumi)

! find index for Ni (ice number mixing ratio)
         dum2 = (alog10(ni1(i,k))+10.)/0.90309
         dumk=int(dum2)

! set limits to make sure the calculated index doesn't exceed range of lookup table
         dum2=min(dum2,real(jsize))
         dum2=max(dum2,1.)
         dumk=max(1,dumk)
         dumk=min(jsize-1,dumk)

! find index for rime mass fraction
         dum4 = qri1(i,k)/qi1(i,k)*3.+1.
         dumii=int(dum4)

! set limits
         dum4=min(dum4,real(rimsize))
         dum4=max(dum4,1.)
         dumii=max(1,dumii)
         dumii=min(rimsize-1,dumii)

! find index for bulk rime density
! account for uneven spacing in lookup table for density
         if (rhop.le.650.) then
         dum5 = (rhop-50.)/200.+1.
         else
         dum5=(rhop-650.)/250.+4.
         end if
         dumjj=int(dum5)

! set limits
         dum5=min(dum5,real(densize))
         dum5=max(dum5,1.)
         dumjj=max(1,dumjj)
         dumjj=min(densize-1,dumjj)

!cccccccccccccccccccccccccccccccccccccccccccc
!
! call subroutine to interpolate lookup table values

         index=6
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr6)
         index=7
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr9)
         index=8
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr10)
         index=9
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr13)
         index=11
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr15)
         index=12
         call access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,f1pr16)

!......................................................

! make sure number concentration is within bounds
            ni1(i,k)=min(ni1(i,k),f1pr9)
            ni1(i,k)=max(ni1(i,k),f1pr10)

! make sure qri is > qsmall
            if (qri1(i,k).lt.qsmall) then
            qri1(i,k)=0.
	    bg1(i,k)=0.
	    end if

! limit max ice concentration to 500 L-1
            ni1(i,k)=min(ni1(i,k),500.e3*orho(i,k))

            effi1(i,k)=f1pr6 ! units are in m

         else

            effi1(i,k)=25.e-6

         end if ! qi1 < qsmall

         if (qi1(i,k).ge.qsmall) then

! mean ice size for output
         di(i,k)=f1pr15
         rhopo(i,k)=f1pr16

         ze_ice(i,k) = (0.176/0.93)*f1pr13
         ze_ice(i,k)=max(ze_ice(i,k),1.e-22)
         end if
         if (qr(i,k).ge.qsmall) then
!         ze_rain(i,k) = n0r(i,k)*720./lamr(i,k)**3/lamr(i,k)**3/lamr(i,k)
! non-exponential rain
         ze_rain(i,k) = nr(i,k)*(mur(i,k)+6.)*(mur(i,k)+5.)*(mur(i,k)+4.)* &
         (mur(i,k)+3.)*(mur(i,k)+2.)*(mur(i,k)+1.)/lamr(i,k)**6

         ze_rain(i,k) = max(ze_rain(i,k),1.e-22)
         end if

! convert to dbz
         z_dbz(i,k) = 10.*log10((ze_rain(i,k)+ze_ice(i,k))*1.d18)

         end do ! k loop

 333     continue

! recalculate supersaturation from T and qv
! calculate temperature from theta

!         do k=kts,kte
!         do k=kbot,ktop,kdir
!         dum = (100000./pres(i,k))**(rd*ocp)
!         t(i,k) = th(i,k)/dum

!         dum=0.622*polysvp1(t(i,k),0)/(pres(i,k)-polysvp1(t(i,k),0))
!	 ssat(i,k)=qv(i,k)-dum
!	 end do

!.....................................................

      end do ! main i loop

! end main microphysics routine

        return

        end subroutine P3_MAIN

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine access_lookup_table(dumjj,dumii,dumi,dumk,index, &
                   dum1,dum2,dum4,dum5,proc)

	implicit none

	real dum1,dum2,dum4,dum5,proc
	real dproc1,dproc2,iproc1,gproc1,tmp1,tmp2
	integer dumjj,dumii,dumi,dumk,index

! get value at current density index

! first interpolate for current rimed fraction index

       dproc1 =itab(dumjj,dumii,dumi,dumk,index)+(dum1-real(dumi))* &
       (itab(dumjj,dumii,dumi+1,dumk,index)-itab(dumjj,dumii,dumi,dumk,index))

       dproc2 =itab(dumjj,dumii,dumi,dumk+1,index)+(dum1-real(dumi))* &
      (itab(dumjj,dumii,dumi+1,dumk+1,index)-itab(dumjj,dumii,dumi,dumk+1,index))

       iproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! linearly interpolate to get process rates for rimed fraction index + 1

       dproc1 =itab(dumjj,dumii+1,dumi,dumk,index)+(dum1-real(dumi))* &
       (itab(dumjj,dumii+1,dumi+1,dumk,index)- &
      itab(dumjj,dumii+1,dumi,dumk,index))

       dproc2 =itab(dumjj,dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))* &
     (itab(dumjj,dumii+1,dumi+1,dumk+1,index)- &
      itab(dumjj,dumii+1,dumi,dumk+1,index))

       gproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)

       tmp1=iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! get value at density index + 1

! first interpolate for current rimed fraction index

       dproc1 =itab(dumjj+1,dumii,dumi,dumk,index)+(dum1-real(dumi))* &
      (itab(dumjj+1,dumii,dumi+1,dumk,index)- &
      itab(dumjj+1,dumii,dumi,dumk,index))

       dproc2 =itab(dumjj+1,dumii,dumi,dumk+1,index)+ &
      (dum1-real(dumi))*(itab(dumjj+1,dumii,dumi+1,dumk+1,index)- &
      itab(dumjj+1,dumii,dumi,dumk+1,index))

       iproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)

! linearly interpolate to get process rates for rimed fraction index + 1

       dproc1 =itab(dumjj+1,dumii+1,dumi,dumk,index)+(dum1-real(dumi))* &
       (itab(dumjj+1,dumii+1,dumi+1,dumk,index)- &
       itab(dumjj+1,dumii+1,dumi,dumk,index))

       dproc2 =itab(dumjj+1,dumii+1,dumi,dumk+1,index)+(dum1-real(dumi))* &
      (itab(dumjj+1,dumii+1,dumi+1,dumk+1,index)- &
       itab(dumjj+1,dumii+1,dumi,dumk+1,index))

       gproc1=dproc1+(dum2-real(dumk))*(dproc2-dproc1)

       tmp2=iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! get final process rate

       proc=tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

       end subroutine access_lookup_table

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine access_lookup_table_coll(dumjj,dumii,dumj,dumi,dumk,index, &
                   dum1,dum2,dum3,dum4,dum5,proc)

	implicit none

	real dum1,dum2,dum3,dum4,dum5,proc
	real dproc1,dproc2,iproc1,gproc1,tmp1,tmp2,dproc11,dproc12,dproc21,dproc22
	integer dumjj,dumii,dumj,dumi,dumk,index


! This subroutine interpolates lookup table values for rain/ice collection processes

! current density index

! current rime fraction index

       dproc11 =itabcoll(dumjj,dumii,dumi,dumk,dumj,index)+(dum1-real(dumi))* &
         (itabcoll(dumjj,dumii,dumi+1,dumk,dumj,index)- &
           itabcoll(dumjj,dumii,dumi,dumk,dumj,index))

       dproc21 =itabcoll(dumjj,dumii,dumi,dumk+1,dumj,index)+(dum1-real(dumi))* &
             (itabcoll(dumjj,dumii,dumi+1,dumk+1,dumj,index)- &
         itabcoll(dumjj,dumii,dumi,dumk+1,dumj,index))

       dproc1=dproc11+(dum2-real(dumk))*(dproc21-dproc11)

       dproc12 =itabcoll(dumjj,dumii,dumi,dumk,dumj+1,index)+(dum1-real(dumi))* &
         (itabcoll(dumjj,dumii,dumi+1,dumk,dumj+1,index)- &
         itabcoll(dumjj,dumii,dumi,dumk,dumj+1,index))

       dproc22 =itabcoll(dumjj,dumii,dumi,dumk+1,dumj+1,index)+(dum1-real(dumi))* &
             (itabcoll(dumjj,dumii,dumi+1,dumk+1,dumj+1,index)- &
          itabcoll(dumjj,dumii,dumi,dumk+1,dumj+1,index))

       dproc2=dproc12+(dum2-real(dumk))*(dproc22-dproc12)

       iproc1=dproc1+(dum3-real(dumj))*(dproc2-dproc1)

! rime fraction index + 1

       dproc11 =itabcoll(dumjj,dumii+1,dumi,dumk,dumj,index)+ &
        (dum1-real(dumi))*(itabcoll(dumjj,dumii+1,dumi+1,dumk,dumj,index)- &
         itabcoll(dumjj,dumii+1,dumi,dumk,dumj,index))

       dproc21 =itabcoll(dumjj,dumii+1,dumi,dumk+1,dumj,index)+ &
        (dum1-real(dumi))*(itabcoll(dumjj,dumii+1,dumi+1,dumk+1,dumj,index)- &
         itabcoll(dumjj,dumii+1,dumi,dumk+1,dumj,index))

       dproc1=dproc11+(dum2-real(dumk))*(dproc21-dproc11)

       dproc12 =itabcoll(dumjj,dumii+1,dumi,dumk,dumj+1,index)+ &
         (dum1-real(dumi))*(itabcoll(dumjj,dumii+1,dumi+1,dumk,dumj+1,index)- &
          itabcoll(dumjj,dumii+1,dumi,dumk,dumj+1,index))

       dproc22 =itabcoll(dumjj,dumii+1,dumi,dumk+1,dumj+1,index)+ &
         (dum1-real(dumi))*(itabcoll(dumjj,dumii+1,dumi+1,dumk+1,dumj+1,index)- &
         itabcoll(dumjj,dumii+1,dumi,dumk+1,dumj+1,index))

       dproc2=dproc12+(dum2-real(dumk))*(dproc22-dproc12)

       gproc1=dproc1+(dum3-real(dumj))*(dproc2-dproc1)

       tmp1=iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! density index + 1

! current rime fraction index

       dproc11 =itabcoll(dumjj+1,dumii,dumi,dumk,dumj,index)+ &
      (dum1-real(dumi))*(itabcoll(dumjj+1,dumii,dumi+1,dumk,dumj,index)- &
           itabcoll(dumjj+1,dumii,dumi,dumk,dumj,index))

       dproc21 =itabcoll(dumjj+1,dumii,dumi,dumk+1,dumj,index)+ &
      (dum1-real(dumi))*(itabcoll(dumjj+1,dumii,dumi+1,dumk+1,dumj,index)- &
         itabcoll(dumjj+1,dumii,dumi,dumk+1,dumj,index))

       dproc1=dproc11+(dum2-real(dumk))*(dproc21-dproc11)

       dproc12 =itabcoll(dumjj+1,dumii,dumi,dumk,dumj+1,index)+ &
      (dum1-real(dumi))*(itabcoll(dumjj+1,dumii,dumi+1,dumk,dumj+1,index)- &
         itabcoll(dumjj+1,dumii,dumi,dumk,dumj+1,index))

       dproc22 =itabcoll(dumjj+1,dumii,dumi,dumk+1,dumj+1,index)+ &
      (dum1-real(dumi))*(itabcoll(dumjj+1,dumii,dumi+1,dumk+1,dumj+1,index)- &
          itabcoll(dumjj+1,dumii,dumi,dumk+1,dumj+1,index))

       dproc2=dproc12+(dum2-real(dumk))*(dproc22-dproc12)

       iproc1=dproc1+(dum3-real(dumj))*(dproc2-dproc1)

! rime fraction index + 1

       dproc11 =itabcoll(dumjj+1,dumii+1,dumi,dumk,dumj,index)+ &
      (dum1-real(dumi))*(itabcoll(dumjj+1,dumii+1,dumi+1,dumk,dumj,index)- &
         itabcoll(dumjj+1,dumii+1,dumi,dumk,dumj,index))

       dproc21 =itabcoll(dumjj+1,dumii+1,dumi,dumk+1,dumj,index)+ &
        (dum1-real(dumi))*(itabcoll(dumjj+1,dumii+1,dumi+1,dumk+1,dumj,index)- &
         itabcoll(dumjj+1,dumii+1,dumi,dumk+1,dumj,index))

       dproc1=dproc11+(dum2-real(dumk))*(dproc21-dproc11)

       dproc12 =itabcoll(dumjj+1,dumii+1,dumi,dumk,dumj+1,index)+ &
         (dum1-real(dumi))*(itabcoll(dumjj+1,dumii+1,dumi+1,dumk,dumj+1,index)- &
          itabcoll(dumjj+1,dumii+1,dumi,dumk,dumj+1,index))

       dproc22 =itabcoll(dumjj+1,dumii+1,dumi,dumk+1,dumj+1,index)+ &
         (dum1-real(dumi))*(itabcoll(dumjj+1,dumii+1,dumi+1,dumk+1,dumj+1,index)- &
         itabcoll(dumjj+1,dumii+1,dumi,dumk+1,dumj+1,index))

       dproc2=dproc12+(dum2-real(dumk))*(dproc22-dproc12)

       gproc1=dproc1+(dum3-real(dumj))*(dproc2-dproc1)


       tmp2=iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! interpolate over density to get final values

       proc=tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

       end subroutine access_lookup_table_coll

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL FUNCTION POLYSVP1 (T,TYPE)

!-------------------------------------------

!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

      IMPLICIT NONE

      REAL DUM
      REAL T
      INTEGER TYPE

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE.EQ.1.and.T.lt.273.15) THEN

!         POLYSVP1 = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
!          LOG10(6.1071))*100.


      dt = max(-80.,t-273.16)
      polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
      polysvp1 = polysvp1*100.

      END IF

! LIQUID

      IF (TYPE.EQ.0.or.T.gt.273.15) THEN

       dt = max(-80.,t-273.16)
       polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp1 = polysvp1*100.

!         POLYSVP1 = 10.**(-7.90298*(373.16/T-1.)+                        &
!             5.02808*LOG10(373.16/T)-									&
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
!             LOG10(1013.246))*100.

         END IF


      END FUNCTION POLYSVP1

!------------------------------------------------------------------------------

      REAL FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
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
!*******************************************************************
!*******************************************************************
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
!     INT, DBLE, EXP, LOG, REAL, SIN
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
      INTEGER I,N
      LOGICAL PARITY
      REAL                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,				   &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA


      REAL FUNCTION DERF(X)
      IMPLICIT NONE
      REAL X
      REAL, DIMENSION(0 : 64) :: A, B
      REAL W,T,Y
      INTEGER K,I
      DATA A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,	                          &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,							  &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0, 							  &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,							  &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      DATA (B(I), I = 0, 12) /                                  &
         -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
         -0.00000000640788827665E0, -0.00000001667820132100E0,  &
         -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
          0.00001612722157047886E0, -0.00025616361025506629E0, 	&
          0.00015380842432375365E0,  0.00815533022524927908E0, 	&
         -0.01402283663896319337E0, -0.19746892495383021487E0,  &
          0.71511720328842845913E0 /
      DATA (B(I), I = 13, 25) /                                 &
         -0.00000000001951073787E0, -0.00000000032302692214E0,  &
          0.00000000522461866919E0,  0.00000000342940918551E0, 	&
         -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
          0.00002687044575042908E0, -0.00011843240273775776E0, 	&
         -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
          0.00909530922354827295E0, -0.20160072778491013140E0, 	&
          0.51169696718727644908E0 /
      DATA (B(I), I = 26, 38) /                                 &
         0.00000000003147682272E0, -0.00000000048465972408E0,   &
         0.00000000063675740242E0,  0.00000003377623323271E0, 	&
        -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
         0.00001947204525295057E0,  0.00002854147231653228E0, 	&
        -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
         0.02328095035422810727E0, -0.16725021123116877197E0, 	&
         0.32490054966649436974E0 /
      DATA (B(I), I = 39, 51) /                                 &
         0.00000000002319363370E0, -0.00000000006303206648E0,   &
        -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
         0.00000011371857327578E0, -0.00000211211337219663E0, 	&
         0.00000368797328322935E0,  0.00009823686253424796E0, 	&
        -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
         0.02585434424202960464E0, -0.11637092784486193258E0, 	&
         0.18267336775296612024E0 /
      DATA (B(I), I = 52, 64) /                                 &
        -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
        -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
         0.00000018006992266137E0, -0.00000078441223763969E0, 	&
        -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
        -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
         0.02064129023876022970E0, -0.06905562880005864105E0,   &
         0.09084526782065478489E0 /
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
              A(K + 11)) * T + A(K + 12)) * W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF (X .LT. 0) Y = -Y
      DERF = Y
      END FUNCTION DERF

!+---+-----------------------------------------------------------------+

!+---+-----------------------------------------------------------------+

       logical function isnan(arg1)
       real,intent(IN) :: arg1
       isnan=( arg1  .ne. arg1 )
       return
       end function isnan

END MODULE MODULE_MP_P3
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
