PROGRAM create_p3_lookuptable_2

!______________________________________________________________________________________
!
! This program creates the lookup table 'p3_lookupTable_2.dat' for inter-category ice-ice
! interactions for the multi-ice-category configuration of the P3 microphysics scheme.
!
!--------------------------------------------------------------------------------------
! Writes table version (parameter 'version' below; p3_init checks it).
! The code itself is versioned with the P3 package (vX.Y.Z of microphy_p3.f90).
!______________________________________________________________________________________
!
! Liquid fraction.  The table holds n_agg and q_agg for four combinations of the liquid
! fractions of category 1 (the COLLECTEE, table dims 1-3) and category 2 (the COLLECTOR,
! dims 4-6).  p3_init (liqfrac=.true. branch) reads 14 columns: 6 index columns, then
!
!   columns   (Fl1,Fl2)   read into
!    7,  8     (0,0)      itabcolli001, itabcolli002
!    9, 10     (0,1)      itabcolli011, itabcolli012
!   11, 12     (1,0)      itabcolli101, itabcolli102
!   13, 14     (1,1)      itabcolli111, itabcolli112
!
! i.e. table names are itabcolli[Fl1][Fl2][quantity]; proc_from_LUT_ii blends Fl2 with the
! collector liquid fraction and Fl1 with the collectee liquid fraction.  One run of this
! program computes one (Fl1,Fl2) combination and writes the 6 index columns plus its
! 2 data columns; the driver script checks that the four runs' index columns agree row by
! row and joins them in the order above.
!______________________________________________________________________________________
!
! To generate the table, use the driver script, which runs the blocks in parallel and
! checks the assembled table:
!
!    src/make_luts.sh lt2 [-j N]
!
! Each run of the program computes one block.  Compile it together with the module of
! shared procedures, module first, with double-precision default reals, e.g.
!
!    ifx -r8 -o create_p3_lookupTable_2 p3_lookupTable_common.f90 create_p3_lookupTable_2.f90
!    gfortran -fdefault-real-8 -fdefault-double-8 -o create_p3_lookupTable_2 \
!             p3_lookupTable_common.f90 create_p3_lookupTable_2.f90
!
! and run it as
!
!    ./create_p3_lookupTable_2 Fl1 Fl2 i1         -> 'lt2_<Fl1><Fl2>_<i1>.dat'
!
! with Fl1, Fl2 = 0 or 1 and i1 = 1..25 (the collectee normalized-Q index).  The header is
! written by block (0,0,1) only.
!______________________________________________________________________________________


 use p3_lookupTable_common

 implicit none

 character(len=16), parameter :: version = '6.2'

 real :: Fl1                   ! liquid fraction of cat 1 (0 or 1, from the command line)
 real :: Fl2                   ! liquid fraction of cat 2 (0 or 1, from the command line)

 real    :: pi,g,p,t,rho,mu,pgam,ds,cs,bas,aas,dcrit,eii
 integer :: k,ii,jj,kk,dumii,j2

 integer, parameter :: n_rhor  =  5
 integer, parameter :: n_Fr    =  4
 integer, parameter :: n_Qnorm = 25

 integer, parameter :: num_bins1 =  1000   ! number of bins for numerical integration of fall speeds and total N
!integer, parameter :: num_bins2 =  9000   ! number of bins for solving PSD parameters  [based on Dm_max = 2000.]
 integer, parameter :: num_bins2 = 11000   ! number of bins for solving PSD parameters  [based on Dm_max = 400000.]
 real,    parameter :: dd        = 20.e-6  ! width of size bin (m)

 integer            :: i_rhor,i_rhor1,i_rhor2  ! indices for rho_rime                 [1 .. n_rhor]
 integer            :: i_Fr,i_Fr1,i_Fr2        ! indices for rime-mass-fraction loop  [1 .. n_Fr]
 integer            :: i,i1,i2                 ! indices for normalized (by N) Q loop [1 .. n_Qnorm]  (i is i_Qnorm in LT1)

 real :: N,q,qdum,dum1,dum2,cs1,ds1,lam,n0,lamf,qerror,del0,c0,c1,c2,sum1,sum2,          &
         sum3,sum4,xx,a0,b0,a1,b1,dum,bas1,aas1,aas2,bas2,d1,d2,delu,lamold,             &
         cap,lamr,dia,amg,dv,n0dum,sum5,sum6,sum7,sum8,dg,cg,bag,aag,dcritg,dcrits,      &
         dcritr,csr,dsr,duml,dum3,rhodep,cgpold,m1,m2,m3,dt,mur,initlamr,lamv,Fr,        &
         rdumii,lammin,lammax,cs2,ds2,intgrR1,intgrR2,intgrR3,intgrR4,q_agg,n_agg,area1, &
         area2,mass1,cs5,rhom1,rhom2,intgrR5

 real, parameter :: mu_i_min =  0.   ! limits of the diagnostic mu_i
 real, parameter :: mu_i_max = 20.   !note: for orig 2-mom, mu_i_max=6.

 real :: Q_normalized   ! function
 real :: lambdai        ! function

 real, parameter :: Dm_max1 = 5000.e-6   ! max. mean ice [m] size for lambda limiter
 real, parameter :: Dm_max2 = 20000.e-6  ! max. mean ice [m] size for lambda limiter
!real, parameter :: Dm_max =  2000.e-6   ! max. mean ice [m] size for lambda limiter
 real, parameter :: Dm_min =     2.e-6   ! min. mean ice [m] size for lambda limiter

 real, parameter :: thrd = 1./3.
 real, parameter :: sxth = 1./6.

!real, dimension(n_Qnorm,n_Fr,n_rhor)    :: qsave,nsave,qon1
 real, dimension(n_Qnorm,n_Fr,n_rhor)    :: qsave
 real, dimension(n_rhor,n_Fr)            :: dcrits1,dcritr1,csr1,dsr1,dcrits2,dcritr2,csr2,dsr2
 real, dimension(n_rhor,n_Fr,n_Qnorm)    :: n01,mu_i1,lam1,n02,mu_i2,lam2
 real, dimension(n_rhor)                 :: cgp,crp
 real, dimension(n_rhor,n_Fr)            :: cgp1,cgp2
 real, dimension(n_Fr)                   :: Fr_arr
 real, dimension(n_rhor,n_Fr,num_bins1)  :: fall1,fall2
 real, dimension(num_bins1)              :: num1,num2

 logical, dimension(n_rhor,n_Fr,n_Qnorm) :: log_lamIsMax

 character(len=2014) :: filename

! command line: one block, Fl1 Fl2 i1
 integer              :: i1_1,i1_2,n_arg,ios1,ios2
 character(len=32)    :: arg

!-------------------------------------------------------------------------------

! the table is computed in double precision: refuse a build without promoted reals
 if (kind(1.) /= 8 .or. kind(1.d0) /= 8) then
    print*, 'create_p3_lookupTable_2: default real and double precision must both be 8 bytes'
    print*, '  (ifx/ifort: -r8;  gfortran: -fdefault-real-8 -fdefault-double-8)'
    stop 1
 endif

 n_arg = command_argument_count()
 ios1 = 1;  ios2 = 1;  Fl1 = -1.;  Fl2 = -1.
 if (n_arg == 3) then
    call get_command_argument(1,arg);  read(arg,*,iostat=ios1) Fl1
    call get_command_argument(2,arg);  read(arg,*,iostat=ios2) Fl2
 endif
 if (ios1/=0 .or. ios2/=0 .or. .not.(Fl1==0. .or. Fl1==1.) .or. .not.(Fl2==0. .or. Fl2==1.)) then
    print*, 'usage: create_p3_lookupTable_2 Fl1 Fl2 i1    (Fl1, Fl2 = 0 or 1; i1 = 1..',n_Qnorm,')'
    stop 1
 endif
 call get_command_argument(3,arg);  read(arg,*) i1_1
 if (i1_1<1 .or. i1_1>n_Qnorm) then
    print*, 'create_p3_lookupTable_2: i1 out of range 1..',n_Qnorm;  stop 1
 endif
 i1_2 = i1_1                ! the loop below runs once: one block per process

! set constants and parameters

! assume 600 hPa, 253 K for p and T for fallspeed calcs (for reference air density)
 pi  = acos(-1.)
 g   = 9.861                           ! gravity
 p   = 60000.                          ! air pressure (pa)
 t   = 253.15                          ! temp (K)
 rho = p/(287.15*t)                    ! air density (kg m-3)
 mu  = 1.496E-6*t**1.5/(t+120.)/rho    ! viscosity of air
 dv  = 8.794E-5*t**1.81/p              ! diffusivity of water vapor in air
 dt  = 10.

! parameters for surface roughness of ice particle
! see mitchell and heymsfield 2005
 del0 = 5.83
 c0   = 0.6
 c1   = 4./(del0**2*c0**0.5)
 c2   = del0**2/4.

!--- specified mass-dimension relationship (cgs units) for unrimed crystals:
! ms = cs*D^ds
!
! for graupel:
! mg = cg*D^dg     no longer used, since bulk volume is predicted

! Heymsfield et al. 2006
!      ds=1.75
!      cs=0.0040157+6.06e-5*(-20.)
! sector-like branches (P1b)
!      ds=2.02
!      cs=0.00142
! bullet-rosette
!     ds=2.45
!      cs=0.00739
! side planes
!      ds=2.3
!      cs=0.00419
! radiating assemblages of plates (mitchell et al. 1990)
!      ds=2.1
!      cs=0.00239
! aggreagtes of side planes, bullets, etc. (Mitchell 1996)
!      ds=2.1
!      cs=0.0028
! Brown and Francis (1995)
 ds = 1.9
! cs = 0.01855 ! original, based on assumption of Dmax
 cs = 0.0121 ! scaled value based on assumtion of Dmean from Hogan et al. 2012, JAMC

! note: if using brown and francis, already in mks units!!!!!
! uncomment line below if using other snow m-D relationships
!      cs=cs*100.**ds/1000.  ! convert from cgs units to mks
!===

! applicable for prognostic graupel density
!  note:  cg is not constant, due to variable density
 dg = 3.


!--- projected area-diam relationship (mks units) for unrimed crystals:
!       note: projected area = aas*D^bas
! sector-like branches (P1b)
!      bas = 1.97
!      aas = 0.55*100.**bas/(100.**2)
! bullet-rosettes
!      bas = 1.57
!      aas = 0.0869*100.**bas/(100.**2)
! graupel (values for hail)
!      bag=2.0
!      aag=0.625*100.**bag/(100.**2)
! aggreagtes of side planes, bullets, etc.
 bas = 1.88
 aas = 0.2285*100.**bas/(100.**2)
!===

!--- projected area-diam relationship (mks units) for graupel:
!      (assumed spheres)
!       note: projected area = aag*D^bag
 aag = pi*0.25
 bag = 2.
!===

! calculate critical diameter separating small spherical ice from crystalline ice
! "Dth" in Morrison and Grabowski 2008

! !  !open file to write to look-up table (which gets used by P3 scheme)
! !  open(unit=1,file='./p3_lookup_table_2.dat-v4',status='unknown')

!.........................................................

!dcrit = (pi/(6.*cs)*0.9)**(1./(ds-3.))
 dcrit = (pi/(6.*cs)*900.)**(1./(ds-3.))
!dcrit=dcrit/100.  ! convert from cm to m

!.........................................................
! main loop over graupel density

! 1D array for RIME density (not ice/graupel density)
 crp(1) =  50.*pi*sxth
 crp(2) = 250.*pi*sxth
 crp(3) = 450.*pi*sxth
 crp(4) = 650.*pi*sxth
 crp(5) = 900.*pi*sxth

! array for rime fraction, Fr
 Fr_arr(1) = 0.
 Fr_arr(2) = 0.333
 Fr_arr(3) = 0.667
 Fr_arr(4) = 1.

!...........................................................................................

! parameters for category 1

! main loop over graupel density

 i_rhor_loop_1: do i_rhor = 1,n_rhor

!------------------------------------------------------------------------
! main loops around N, q, Fr for lookup tables

! find threshold with rimed mass added

! loop over rimed mass fraction (4 points)

    i_Fr_loop_1: do i_Fr = 1,n_Fr   ! loop for rime mass fraction, Fr

      ! Rime mass fraction for the lookup table (specific values in model are interpolated between points)
       Fr = Fr_arr(i_Fr)

! ! ! calculate critical dimension separate graupel and nonspherical ice
! ! ! "Dgr" in morrison and grabowski (2008)
! !
! ! !   dcrits = (cs/crp(i_rhor))**(1./(dg-ds)) ! calculated below for variable graupel density
! ! !   dcrits = dcrits/100.  ! convert from cm to m
! !
! ! ! check to make sure projected area at dcrit not greater than than of solid sphere
! ! ! stop and give warning message if that is the case
! !
! ! !    if (pi/4.*dcrit**2.lt.aas*dcrit**bas) then
! ! !       print*,'STOP, area > area of solid ice sphere, unrimed'
! ! !       stop
! ! !    endif
! ! !    if (pi/4.*dcrits1(i_rhor,i_Fr)**2.lt.aag*dcrits1(i_rhor,i_Fr)**bag) then
! ! !       print*,'STOP, area > area of solid ice sphere, graupel'
! ! !       stop
! ! !    endif
! !
! ! !      cg=cg*100.**dg/1000.  ! convert from cgs units to mks
! !
! ! !      print*,cg,dg
! ! !      stop
! ! !      do jj=1,100
! ! !         dd=real(jj)*30.e-6
! ! !         write(6,'5e15.5')dd,aas*dd**bas,pi/4.*dd**2,
! ! !     1      cs*dd**ds,pi*sxth*917.*dd**3
! ! !      end do

! calculate mass-dimension relationship for partially-rimed crystals
! msr = csr*D^dsr
! formula from morrison grabowski 2008

! dcritr is critical size separating graupel from partially-rime crystal
! same as "Dcr" in morrison and grabowski 2008

! first guess, set cgp=crp
       cgp1(i_rhor,i_Fr) = crp(i_rhor)

! case of no riming (Fr = 0%), then we need to set dcrits and dcritr to arbitrary large values

       if (i_Fr.eq.1) then

          dcrits1(i_rhor,i_Fr) = 1.e+6
          dcritr1(i_rhor,i_Fr) = dcrits1(i_rhor,i_Fr)
          csr1(i_rhor,i_Fr)    = cs
          dsr1(i_rhor,i_Fr)    = ds
! case of partial riming (Fr between 0 and 100%)

       elseif (i_Fr.eq.2.or.i_Fr.eq.3) then

          do
             dcrits1(i_rhor,i_Fr) = (cs/cgp1(i_rhor,i_Fr))**(1./(dg-ds))
             dcritr1(i_rhor,i_Fr) = ((1.+Fr/(1.-Fr))*cs/cgp1(i_rhor,i_Fr))**(1./(dg-ds))
             csr1(i_rhor,i_Fr)    = cs*(1.+Fr/(1.-Fr))
             dsr1(i_rhor,i_Fr)    = ds
! get mean density of vapor deposition/aggregation grown ice
             rhodep = 1./(dcritr1(i_rhor,i_Fr)-dcrits1(i_rhor,i_Fr))*6.*cs/(pi*(ds-2.))* &
                      (dcritr1(i_rhor,i_Fr)**(ds-2.)-dcrits1(i_rhor,i_Fr)**(ds-2.))
! get graupel density as rime mass fraction weighted rime density plus
! density of vapor deposition/aggregation grown ice
             cgpold = cgp1(i_rhor,i_Fr)
             cgp1(i_rhor,i_Fr) = crp(i_rhor)*Fr+rhodep*(1.-Fr)*pi*sxth
             if (abs((cgp1(i_rhor,i_Fr)-cgpold)/cgp1(i_rhor,i_Fr)).lt.0.01) goto 115
          enddo

 115  continue

! case of complete riming (Fr=100%)
       else

! set threshold size for pure graupel arbitrary large
          dcrits1(i_rhor,i_Fr) = (cs/cgp1(i_rhor,i_Fr))**(1./(dg-ds))
          dcritr1(i_rhor,i_Fr) = 1.e6
          csr1(i_rhor,i_Fr)    = cgp1(i_rhor,i_Fr)
          dsr1(i_rhor,i_Fr)    = dg

       endif  !if i_Fr.eq.1

!---------------------------------------------------------------------------------------
! set up particle fallspeed arrays
! fallspeed is a function of mass dimension and projected area dimension relationships
! following mitchell and heymsfield (2005), jas

! set up array of particle fallspeed to make computationally efficient
!.........................................................
! ****
!  note: this part could be incorporated into the longer (every 2 micron) loop
! ****
       jj_loop_1: do jj = 1,num_bins1

          d1 = real(jj)*dd - 0.5*dd   !particle size [m]

          if (d1.le.dcrit) then
             cs1  = pi*sxth*900.
             ds1  = 3.
             bas1 = 2.
             aas1 = pi/4.
          else if (d1.gt.dcrit.and.d1.le.dcrits1(i_rhor,i_Fr)) then
             cs1  = cs
             ds1  = ds
             bas1 = bas
             aas1 = aas
          else if (d1.gt.dcrits1(i_rhor,i_Fr).and.d1.le.dcritr1(i_rhor,i_Fr)) then
             cs1  = cgp1(i_rhor,i_Fr)
             ds1  = dg
             bas1 = bag
             aas1 = aag
          else if (d1.gt.dcritr1(i_rhor,i_Fr)) then
             cs1  = csr1(i_rhor,i_Fr)
             ds1  = dsr1(i_rhor,i_Fr)
             if (i_Fr.eq.1) then
                aas1 = aas
                bas1 = bas
             else
! for area, keep bas1 constant, but modify aas1 according to rimed fraction
                bas1 = bas
                dum1 = aas*d1**bas
                dum2 = aag*d1**bag
!               dum3 = (1.-Fr)*dum1+Fr*dum2
                m1   = cs1*d1**ds1
                m2   = cs*d1**ds
                m3   = cgp1(i_rhor,i_Fr)*d1**dg
                dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)  !linearly interpolate based on particle mass
                aas1 = dum3/(d1**bas)
             endif
          endif

! correction for turbulence
!            if (d1.lt.500.e-6) then
          a0 = 0.
          b0 = 0.
!            else
!               a0=1.7e-3
!               b0=0.8
!            end if

! fall speed for ice
! best number
          xx = 2.*cs1*g*rho*d1**(ds1+2.-bas1)/(aas1*(mu*rho)**2)
! drag terms
          b1 = c1*xx**0.5/(2.*((1.+c1*xx**0.5)**0.5-1.)*(1.+c1*xx**0.5)**0.5)-a0*b0*xx** &
               b0/(c2*((1.+c1*xx**0.5)**0.5-1.)**2)
          a1 = (c2*((1.+c1*xx**0.5)**0.5-1.)**2-a0*xx**b0)/xx**b1
! velocity in terms of drag terms


     !------------------------------------
     ! fall speed for rain particle

          dia = d1  ! diameter m
          amg = pi*sxth*997.*dia**3 ! mass [kg]
          amg = amg*1000.           ! convert kg to g

          if (dia.le.134.43e-6) then
             dum2 = 4.5795e5*amg**(2.*thrd)
             goto 101
          endif

          if(dia.lt.1511.64e-6) then
             dum2 = 4.962e3*amg**thrd
            goto 101
          endif

          if(dia.lt.3477.84e-6) then
             dum2 = 1.732e3*amg**sxth
             goto 101
          endif

          dum2 = 917.

101       continue

          fall1(i_rhor,i_Fr,jj) = (1.-Fl1)*(a1*mu**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.))+Fl1*dum2*1.e-2

       enddo jj_loop_1

!---------------------------------------------------------------------------------
! main loops around normalized q for lookup table
!
! q = normalized ice mass mixing ratio = q/N, units are kg^-1

       i_Qnorm_loop_1: do i = 1,n_Qnorm

          q = Q_normalized(i)

          print*,'i,Fr,i_rhor ',i,i_Fr,i_rhor
          print*,'q* ',q

          qerror = 1.e+20   !initialize qerror to arbitrarily large value

!.....................................................................................
! find parameters for gamma distribution

! size distribution for ice is assumed to be
! N(D) = n0 * D^pgam * exp(-lam*D)

! for the given q and N, we need to find n0, pgam, and lam

! approach for finding lambda:
! cycle through a range of lambda, find closest lambda to produce correct q

! start with lam, range of lam from 100 to 5 x 10^6 is large enough to
! cover full range over mean size from 2 to 5000 micron

          ii_loop: do ii = 1,num_bins2

             lam1(i_rhor,i_Fr,i) = lambdai(ii)
             rhom1 = (1.-Fl1)*cgp1(i_rhor,i_Fr)+Fl1*1000.*pi*sxth
             mu_i1(i_rhor,i_Fr,i) = diagnostic_mui(mu_i_min,mu_i_max,lam1(i_rhor,i_Fr,i),q,rhom1,Fr,pi)

! set min,max lam corresponding to Dm_max,Dm_min:
             dum = Dm_max1+Dm_max2*Fr**2.
             lam1(i_rhor,i_Fr,i) = max(lam1(i_rhor,i_Fr,i),(mu_i1(i_rhor,i_Fr,i)+1.)/dum)
             lam1(i_rhor,i_Fr,i) = min(lam1(i_rhor,i_Fr,i),(mu_i1(i_rhor,i_Fr,i)+1.)/Dm_min)
! get normalized n0 = n0/N
             n01(i_rhor,i_Fr,i) = lam1(i_rhor,i_Fr,i)**(mu_i1(i_rhor,i_Fr,i)+1.)/(gamma(mu_i1(i_rhor,i_Fr,i)+1.))

! calculate integral for each of the 4 parts of the size distribution
! check difference with respect to q

! set up m-D relationship for solid ice with D < Dcrit
             cs1  = pi*sxth*900.
             ds1  = 3.
             cs5  = pi*sxth*1000.

             call intgrl_section_Fl(lam1(i_rhor,i_Fr,i),mu_i1(i_rhor,i_Fr,i), ds1,ds,dg,          &
                                    dsr1(i_rhor,i_Fr),dcrit,dcrits1(i_rhor,i_Fr),                 &
                                    dcritr1(i_rhor,i_Fr),intgrR1,intgrR2,intgrR3,intgrR4,intgrR5)

! intgrR1 is integral from 0 to dcrit (solid ice)
! intgrR2 is integral from dcrit to dcrits (snow)
! intgrR3 is integral from dcrits to dcritr (graupel)
! intgrR4 is integral from dcritr to inf (rimed snow)

! sum of the integrals from the 4 regions of the size distribution
             qdum = n01(i_rhor,i_Fr,i)*( (1.-Fl1)*(cs1*intgrR1 + cs*intgrR2 + cgp1(i_rhor,i_Fr)*intgrR3 +   &
                     csr1(i_rhor,i_Fr)*intgrR4) + Fl1*cs5*intgrR5)

             if (ii.eq.1) then
                qerror = abs(q-qdum)
                lamf   = lam1(i_rhor,i_Fr,i)
             endif

! find lam with smallest difference between q and estimate of q, assign to lamf
             if (abs(q-qdum).lt.qerror) then
                lamf   = lam1(i_rhor,i_Fr,i)
                qerror = abs(q-qdum)
             endif

          enddo ii_loop

! check and print relative error in q to make sure it is not too large
! note: large error is possible if size bounds are exceeded!!!!!!!!!!

          print*,'qerror (%)',qerror/q*100.

! find n0 based on final lam value
! set final lamf to 'lam' variable
! this is the value of lam with the smallest qerror
          lam1(i_rhor,i_Fr,i) = lamf
! recalculate mu_i based on final lam
          mu_i1(i_rhor,i_Fr,i) = diagnostic_mui(mu_i_min,mu_i_max,lam1(i_rhor,i_Fr,i),q,rhom1,Fr,pi)

!            n0 = N*lam**(pgam+1.)/(gamma(pgam+1.))

! find n0 from lam and q
! this is done instead of finding n0 from lam and N, since N
! may need to be adjusted to constrain mean size within reasonable bounds

          call intgrl_section_Fl(lam1(i_rhor,i_Fr,i),mu_i1(i_rhor,i_Fr,i), ds1,ds,dg,             &
                                    dsr1(i_rhor,i_Fr),dcrit,dcrits1(i_rhor,i_Fr),                 &
                                    dcritr1(i_rhor,i_Fr),intgrR1,intgrR2,intgrR3,intgrR4,intgrR5)

! normalized n0
          cs5  = pi*sxth*1000.
          n01(i_rhor,i_Fr,i) = q/( (1.-Fl1)*(cs1*intgrR1 + cs*intgrR2 + cgp1(i_rhor,i_Fr)*intgrR3 +       &
                               csr1(i_rhor,i_Fr)*intgrR4) + Fl1*cs5*intgrR5)

          print*,'lam,N0:',lam1(i_rhor,i_Fr,i),n01(i_rhor,i_Fr,i)
          print*,'mu_i:',mu_i1(i_rhor,i_Fr,i)
          print*,'mean size:',(mu_i1(i_rhor,i_Fr,i)+1.)/lam1(i_rhor,i_Fr,i)

! At this point, we have solve for all of the size distribution parameters

! NOTE: In the code it is assumed that mean size and number have already been
! adjusted, so that mean size will fall within allowed bounds. Thus, we do
! not apply a lambda limiter here.

       enddo i_Qnorm_loop_1

    enddo i_Fr_loop_1
 enddo i_rhor_loop_1

!--------------------------------------------------------------------
!.....................................................................................
! now calculate parameters for category 2

 i_rhor_loop_2: do i_rhor = 1,n_rhor

    i_Fr_loop_2: do i_Fr = 1,n_Fr

       Fr = Fr_arr(i_Fr)

! calculate mass-dimension relationship for partially-rimed crystals
! msr = csr*D^dsr
! formula from morrison grabowski 2008

! dcritr is critical size separating graupel from partially-rime crystal
! same as "Dcr" in morrison and grabowski 2008

! first guess, set cgp=crp
       cgp2(i_rhor,i_Fr) = crp(i_rhor)

! case of no riming (Fr = 0%), then we need to set dcrits and dcritr to arbitrary large values

       if (i_Fr.eq.1) then

          dcrits2(i_rhor,i_Fr) = 1.e+6
          dcritr2(i_rhor,i_Fr) = dcrits2(i_rhor,i_Fr)
          csr2(i_rhor,i_Fr)    = cs
          dsr2(i_rhor,i_Fr)    = ds

! case of partial riming (Fr between 0 and 100%)
       elseif (i_Fr.eq.2.or.i_Fr.eq.3) then

          do
             dcrits2(i_rhor,i_Fr) = (cs/cgp2(i_rhor,i_Fr))**(1./(dg-ds))
             dcritr2(i_rhor,i_Fr) = ((1.+Fr/(1.-Fr))*cs/cgp2(i_rhor,i_Fr))**(1./(dg-ds))
             csr2(i_rhor,i_Fr)    = cs*(1.+Fr/(1.-Fr))
             dsr2(i_rhor,i_Fr)    = ds
! get mean density of vapor deposition/aggregation grown ice
             rhodep = 1./(dcritr2(i_rhor,i_Fr)-dcrits2(i_rhor,i_Fr))*6.*cs/(pi*(ds-2.))* &
                      (dcritr2(i_rhor,i_Fr)**(ds-2.)-dcrits2(i_rhor,i_Fr)**(ds-2.))
! get graupel density as rime mass fraction weighted rime density plus
! density of vapor deposition/aggregation grown ice
             cgpold = cgp2(i_rhor,i_Fr)
             cgp2(i_rhor,i_Fr) = crp(i_rhor)*Fr+rhodep*(1.-Fr)*pi*sxth
             if (abs((cgp2(i_rhor,i_Fr)-cgpold)/cgp2(i_rhor,i_Fr)).lt.0.01) goto 116
          enddo

 116  continue

! case of complete riming (Fr=100%)
       else

! set threshold size for pure graupel arbitrary large
          dcrits2(i_rhor,i_Fr) = (cs/cgp2(i_rhor,i_Fr))**(1./(dg-ds))
          dcritr2(i_rhor,i_Fr) = 1.e+6
          csr2(i_rhor,i_Fr)    = cgp2(i_rhor,i_Fr)
          dsr2(i_rhor,i_Fr)    = dg

       endif

!---------------------------------------------------------------------------------------
! set up particle fallspeed arrays
! fallspeed is a function of mass dimension and projected area dimension relationships
! following mitchell and heymsfield (2005), jas

! set up array of particle fallspeed to make computationally efficient

       jj_loop_2: do jj = 1,num_bins1

          d1 = real(jj)*dd - 0.5*dd   !particle size [m]

          if (d1.le.dcrit) then
             cs1  = pi*sxth*900.
             ds1  = 3.
             bas1 = 2.
             aas1 = pi/4.
          else if (d1.gt.dcrit.and.d1.le.dcrits2(i_rhor,i_Fr)) then
             cs1  = cs
             ds1  = ds
             bas1 = bas
             aas1 = aas
          else if (d1.gt.dcrits2(i_rhor,i_Fr).and.d1.le.dcritr2(i_rhor,i_Fr)) then
             cs1  = cgp2(i_rhor,i_Fr)
             ds1  = dg
             bas1 = bag
             aas1 = aag
          else if (d1.gt.dcritr2(i_rhor,i_Fr)) then
             cs1  = csr2(i_rhor,i_Fr)
             ds1  = dsr2(i_rhor,i_Fr)
             if (i_Fr.eq.1) then
                aas1 = aas
                bas1 = bas
             else
! for area, keep bas1 constant, but modify aas1 according to rimed fraction
                bas1 = bas
                dum1 = aas*d1**bas
                dum2 = aag*d1**bag
!               dum3 = (1.-Fr)*dum1+Fr*dum2
                m1   = cs1*d1**ds1
                m2   = cs*d1**ds
                m3   = cgp2(i_rhor,i_Fr)*d1**dg
                dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)  !linearly interpolate based on particle mass
                aas1 = dum3/(d1**bas)
             endif
          endif

! correction for turbulence
!            if (d1.lt.500.e-6) then
          a0 = 0.
          b0 = 0.
!            else
!               a0=1.7e-3
!               b0=0.8
!            end if

! fall speed for ice
! best number
          xx = 2.*cs1*g*rho*d1**(ds1+2.-bas1)/(aas1*(mu*rho)**2)

! drag terms
          b1 = c1*xx**0.5/(2.*((1.+c1*xx**0.5)**0.5-1.)*(1.+c1*xx**0.5)**0.5)-a0*b0*xx** &
               b0/(c2*((1.+c1*xx**0.5)**0.5-1.)**2)

          a1 = (c2*((1.+c1*xx**0.5)**0.5-1.)**2-a0*xx**b0)/xx**b1

     !------------------------------------
     ! fall speed for rain particle

          dia = d1  ! diameter m
          amg = pi*sxth*997.*dia**3 ! mass [kg]
          amg = amg*1000.           ! convert kg to g

          if (dia.le.134.43e-6) then
             dum2 = 4.5795e5*amg**(2.*thrd)
             goto 102
          endif

          if(dia.lt.1511.64e-6) then
             dum2 = 4.962e3*amg**thrd
            goto 102
          endif

          if(dia.lt.3477.84e-6) then
             dum2 = 1.732e3*amg**sxth
             goto 102
          endif

          dum2 = 917.

102       continue

! velocity in terms of drag terms
          fall2(i_rhor,i_Fr,jj) = (1.-Fl2)*(a1*mu**(1.-2.*b1)*(2.*cs1*g/(rho*aas1))**b1*d1**(b1*(ds1-bas1+2.)-1.))+Fl2*dum2*1.e-2

!---------------------------------------------------------------
       enddo jj_loop_2


!---------------------------------------------------------------------------------

! q = total ice mixing ratio (vapor dep. plus rime mixing ratios), normalized by N

       i_Qnorm_loop_2: do i = 1,n_Qnorm

          lamold = 0.

          q = Q_normalized(i)

          print*,'&&&&&&&&&&&i_rhor',i_rhor
          print*,'***************',i,k
          print*,'Fr',Fr
          print*,'q,N',q,N

          qerror = 1.e+20   ! initialize qerror to arbitrarily large value

!.....................................................................................
! find parameters for gamma distribution

! size distribution for ice is assumed to be
! N(D) = n0 * D^pgam * exp(-lam*D)

! for the given q and N, we need to find n0, pgam, and lam

! approach for finding lambda:
! cycle through a range of lambda, find closest lambda to produce correct q

! start with lam, range of lam from 100 to 5 x 10^6 is large enough to
! cover full range over mean size from 2 to 5000 micron

          ii_loop_2: do ii = 1,num_bins2

             lam2(i_rhor,i_Fr,i) = lambdai(ii)
             rhom2 = (1.-Fl2)*cgp2(i_rhor,i_Fr)+Fl2*1000.*pi*sxth
             mu_i2(i_rhor,i_Fr,i) = diagnostic_mui(mu_i_min,mu_i_max,lam2(i_rhor,i_Fr,i),q,rhom2,Fr,pi)

            !set min,max lam corresponding to Dm_max, Dm_min
             dum = Dm_max1+Dm_max2*Fr**2.
             lam2(i_rhor,i_Fr,i) = max(lam2(i_rhor,i_Fr,i),(mu_i2(i_rhor,i_Fr,i)+1.)/dum)
             lam2(i_rhor,i_Fr,i) = min(lam2(i_rhor,i_Fr,i),(mu_i2(i_rhor,i_Fr,i)+1.)/Dm_min)

            !get n0, note this is normalized
             n02(i_rhor,i_Fr,i) = lam2(i_rhor,i_Fr,i)**(mu_i2(i_rhor,i_Fr,i)+1.)/ &
                   (gamma(mu_i2(i_rhor,i_Fr,i)+1.))

! calculate integral for each of the 4 parts of the size distribution
! check difference with respect to q

! set up m-D relationship for solid ice with D < Dcrit
             cs1  = pi*sxth*900.
             ds1  = 3.
             cs5  = pi*sxth*1000.

             call intgrl_section_Fl(lam2(i_rhor,i_Fr,i),mu_i2(i_rhor,i_Fr,i), ds1,ds,dg,          &
                                    dsr2(i_rhor,i_Fr),dcrit,dcrits2(i_rhor,i_Fr),                 &
                                    dcritr2(i_rhor,i_Fr),intgrR1,intgrR2,intgrR3,intgrR4,intgrR5)

! sum of the integrals from the 4 regions of the size distribution
             qdum = n02(i_rhor,i_Fr,i)*( (1.-Fl2)*(cs1*intgrR1 + cs*intgrR2 + cgp2(i_rhor,i_Fr)*intgrR3 +  &
                    csr2(i_rhor,i_Fr)*intgrR4) + Fl2*cs5*intgrR5)

             if (ii.eq.1) then
                qerror = abs(q-qdum)
                lamf   = lam2(i_rhor,i_Fr,i)
             endif

! find lam with smallest difference between q and estimate of q, assign to lamf
             if (abs(q-qdum).lt.qerror) then
                lamf   = lam2(i_rhor,i_Fr,i)
                qerror = abs(q-qdum)
             endif

          enddo ii_loop_2

! check and print relative error in q to make sure it is not too large
! note: large error is possible if size bounds are exceeded!!!!!!!!!!

          print*,'qerror (%)',qerror/q*100.

! find n0 based on final lam value
! set final lamf to 'lam' variable
! this is the value of lam with the smallest qerror
          lam2(i_rhor,i_Fr,i) = lamf

! recalculate mu based on final lam
!         mu_i2(i_rhor,i_Fr,i) = diagnostic_mui(log_diagmu_orig,lam2(i_rhor,i_Fr,i),q,cgp2(i_rhor,i_Fr),Fr,pi)
          mu_i2(i_rhor,i_Fr,i) = diagnostic_mui(mu_i_min,mu_i_max,lam2(i_rhor,i_Fr,i),q,rhom2,Fr,pi)

!            n0 = N*lam**(pgam+1.)/(gamma(pgam+1.))

! find n0 from lam and q
! this is done instead of finding n0 from lam and N, since N
! may need to be adjusted to constrain mean size within reasonable bounds

          call intgrl_section_Fl(lam2(i_rhor,i_Fr,i),mu_i2(i_rhor,i_Fr,i), ds1,ds,dg,              &
                                    dsr2(i_rhor,i_Fr),dcrit,dcrits2(i_rhor,i_Fr),                  &
                                    dcritr2(i_rhor,i_Fr),intgrR1,intgrR2,intgrR3,intgrR4,intgrR5)

          cs5  = pi*sxth*1000.
          n02(i_rhor,i_Fr,i) = q/( (1.-Fl2)*(cs1*intgrR1 + cs*intgrR2 + cgp2(i_rhor,i_Fr)*intgrR3 +        &  ! n0 is normalized
                               csr2(i_rhor,i_Fr)*intgrR4) + Fl2*cs5*intgrR5)
          print*,'lam,N0:',lam2(i_rhor,i_Fr,i),n02(i_rhor,i_Fr,i)
          print*,'pgam:',mu_i2(i_rhor,i_Fr,i)

          qsave(i,i_Fr,i_rhor) = q      ! q is normalized (Q/N)

          log_lamIsMax(i_rhor,i_Fr,i) = abs(lam2(i_rhor,i_Fr,i)-lamold) .lt. 1.e-8
          lamold = lam2(i_rhor,i_Fr,i)

       enddo i_Qnorm_loop_2

    enddo i_Fr_loop_2
 enddo i_rhor_loop_2


! At this point, we have solve for all of the size distribution parameters

! NOTE: In the code it is assumed that mean size and number have already been
! adjusted, so that mean size will fall within allowed bounds. Thus, we do
! not apply a lambda limiter here.

!--------------------------------------------------------------------

!--------------------------------------------------------------------

!.....................................................................................
! begin category process interaction calculations for the lookup table

!.....................................................................................
! collection of category 1 by category 2
!.....................................................................................

 write (filename, "(a4,2i1,a1,i2.2,a4)") "lt2_",nint(Fl1),nint(Fl2),"_",i1_1,".dat"
 open(unit=1, file=filename, status='unknown')

 ! header: written once, by the (0,0) run's first block
 if (i1_1==1 .and. Fl1==0. .and. Fl2==0.) then
    write(1,*) 'LOOKUP_TABLE_2-version:  ',trim(version)
    write(1,'(a)') ''          ! an empty line, whatever the compiler
 endif

 Qnorm_loop_3: do i1 = i1_1,i1_2
   do i_Fr1 = 1,n_Fr
     do i_rhor1 = 1,n_rhor
       do i2 = 1,n_Qnorm
         do i_Fr2 = 1,n_Fr
           do i_rhor2 = 1,n_rhor

                lamIsMax: if (.not. log_lamIsMax(i_rhor2,i_Fr2,i2)) then

                sum1 = 0.
                sum2 = 0.

              !set up binned distribution of ice from categories 1 and 2 (distributions normalized by N)
                do jj = 1,num_bins1
                   d1 = real(jj)*dd - 0.5*dd
                   num1(jj) = n01(i_rhor1,i_Fr1,i1)*d1**mu_i1(i_rhor1,i_Fr1,i1)*         &
                              exp(-lam1(i_rhor1,i_Fr1,i1)*d1)*dd
                   num2(jj) = n02(i_rhor2,i_Fr2,i2)*d1**mu_i2(i_rhor2,i_Fr2,i2)*         &
                              exp(-lam2(i_rhor2,i_Fr2,i2)*d1)*dd
                enddo

! loop over size distribution
! note: collection of ice within the same bin is neglected

! loop over particle 1
                jj_loop_3: do jj = num_bins1,1,-1

                   d1 = real(jj)*dd - 0.5*dd   !particle size [m]

                   if (d1.le.dcrit) then
                      cs1  = pi*sxth*900.
                      ds1  = 3.
                      bas1 = 2.
                      aas1 = pi*0.25
                   elseif (d1.gt.dcrit.and.d1.le.dcrits1(i_rhor1,i_Fr1)) then
                      cs1  = cs
                      ds1  = ds
                      bas1 = bas
                      aas1 = aas
                   else if (d1.gt.dcrits1(i_rhor1,i_Fr1).and.d1.le.dcritr1(i_rhor1,i_Fr1)) then
                      cs1  = cgp1(i_rhor1,i_Fr1)
                      ds1  = dg
                      bas1 = bag
                      aas1 = aag
                   else if (d1.gt.dcritr1(i_rhor1,i_Fr1)) then
                      cs1 = csr1(i_rhor1,i_Fr1)
                      ds1 = dsr1(i_rhor1,i_Fr1)
                      if (i_Fr1.eq.1) then
                         aas1 = aas
                         bas1 = bas
                      else
! for area, keep bas1 constant, but modify aas1 according to rimed fraction
                         bas1 = bas
                         dum1 = aas*d1**bas
                         dum2 = aag*d1**bag
                         m1   = cs1*d1**ds1
                         m2   = cs*d1**ds
                         m3   = cgp1(i_rhor1,i_Fr1)*d1**dg
! linearly interpolate based on particle mass
                         dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
                         aas1 = dum3/(d1**bas)
                      endif
                   endif

                   mass1 = (1.-Fl1)*cs1*d1**ds1+Fl1*pi*sxth*1000.*d1**3.
                   area1 = (1.-Fl1)*aas1*d1**bas1+Fl1*pi/4.*d1**2.

! loop over particle 2
                   kk_loop: do kk = num_bins1,1,-1

                      d2 = real(kk)*dd - 0.5*dd   !particle size [m]

! parameters for particle 2
                      if (d2.le.dcrit) then
!                        cs2  = pi*sxth*900.
!                        ds2  = 3.
                         bas2 = 2.
                         aas2 = pi*0.25
                      elseif (d2.gt.dcrit.and.d2.le.dcrits2(i_rhor2,i_Fr2)) then
!                        cs2  = cs
!                        ds2  = ds
                         bas2 = bas
                         aas2 = aas
                      else if (d2.gt.dcrits2(i_rhor2,i_Fr2).and.d2.le.dcritr2(i_rhor2,i_Fr2)) then
!                        cs2  = cgp2(i_rhor2,i_Fr2)
!                        ds2  = dg
                         bas2 = bag
                         aas2 = aag
                      else if (d2.gt.dcritr2(i_rhor2,i_Fr2)) then
                         cs2 = csr2(i_rhor2,i_Fr2)
                         ds2 = dsr2(i_rhor2,i_Fr2)
                         if (i_Fr2.eq.1) then
                            aas2 = aas
                            bas2 = bas
                         else
! for area, keep bas1 constant, but modify aas1 according to rimed fraction
                            bas2 = bas
                            dum1 = aas*d2**bas
                            dum2 = aag*d2**bag
                            m1   = cs2*d2**ds2
                            m2   = cs*d2**ds
                            m3   = cgp2(i_rhor2,i_Fr2)*d2**dg
! linearly interpolate based on particle mass
                            dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
                            aas2 = dum3/(d2**bas)
                         endif
                      endif

                    area2 = (1.-Fl2)*aas2*d2**bas2+Fl2*pi/4.*d2**2.

! absolute value, differential fallspeed
!                   delu = abs(fall2(i_rhor2,i_Fr2,kk)-fall1(i_rhor1,i_Fr1,jj))

! calculate collection of category 1 by category 2, which occurs
! in fallspeed of particle in category 2 is greater than category 1

                      if (fall2(i_rhor2,i_Fr2,kk).gt.fall1(i_rhor1,i_Fr1,jj)) then

                         delu = fall2(i_rhor2,i_Fr2,kk)-fall1(i_rhor1,i_Fr1,jj)

! note: in micro code we have to multiply by air density
! correction factor for fallspeed, and collection efficiency

! sum for integral

! sum1 = # of collision pairs
! the assumption is that each collision pair reduces crystal
! number mixing ratio by 1 kg^-1 s^-1 per kg/m^3 of air (this is
! why we need to multiply by air density, to get units of
! 1/kg^-1 s^-1)
! NOTE: For consideration of particle depletion, air density is assumed to be 1 kg m-3
! This problem could be avoided by using number concentration instead of number mixing ratio
! for the lookup table calculations, and then not multipling process rate by air density
! in the P3 code... TO BE FIXED IN THE FUTURE

!                   sum1 = sum1+min((sqrt(area1)+sqrt(area2))**2*delu*num(jj)*num(kk),   &
!                          num(kk)/dt)

! set collection efficiency
!                  eii = 0.1

! accretion of number
                         !sum1 = sum1 + (sqrt(aas1*d1**bas1) + sqrt(aas2*d2**bas2))**2*delu*num1(jj)*num2(kk)
                         sum1 = sum1 + (sqrt(area1)+sqrt(area2))**2*delu*num1(jj)*num2(kk)
! accretion of mass
                         !sum2 = sum2 + cs1*d1**ds1*(sqrt(aas1*d1**bas1) + sqrt(aas2*d2**bas2))**2*delu*num1(jj)*num2(kk)
                         sum2 = sum2 + mass1*(sqrt(area1)+sqrt(area2))**2*delu*num1(jj)*num2(kk)

! remove collected particles from distribution over time period dt, update num1
!  note -- dt is time scale for removal, not necessarily the model time step

                      endif ! fall2(i_rhor2,i_Fr2,kk) > fall1(i_rhor1,i_Fr1,jj)

                   enddo kk_loop
                enddo jj_loop_3

! save for output
                n_agg = sum1
                q_agg = sum2

             else

                print*,'&&&&&&, skip'
                n_agg = -999.
                q_agg = -999.

             endif lamIsMax

             n_agg = dim(n_agg, 1.e-50)
             q_agg = dim(q_agg, 1.e-50)

            !write(6,'(a5,6i5,2e15.5)') 'index:',i1,i_Fr1,i_rhor1,i2,i_Fr2,i_rhor2,n_agg,q_agg
             ! every (Fl1,Fl2) run writes the index columns, so the driver can check the
             ! four runs line up row by row before joining their data columns
             write(1,'(6i5,2e15.6)') i1,i_Fr1,i_rhor1,i2,i_Fr2,i_rhor2,n_agg,q_agg

           enddo   ! i_rhor2 loop
         enddo   ! i_Fr2 loop
       enddo   ! i2 loop
     enddo   ! i_rhor1 loop
   enddo   ! i_Fr1
 enddo Qnorm_loop_3  ! i1 loop  (Qnorm)

 close(1)

END PROGRAM create_p3_lookuptable_2
!______________________________________________________________________________________

 real function Q_normalized(i_qnorm)

 !-----------------
 ! Computes the normalized Q (qitot/nitot) based on a given index (used in the Qnorm loops)
 !-----------------

 implicit none

!arguments:
 integer :: i_qnorm

!Q_normalized = 261.7**((i_qnorm+5)*0.2)*1.e-18     ! v4, v5.0 (range of mean mass diameter from ~ 1 micron to 1 cm)   (for n_Qnorm = 25)
!Q_normalized = 800.**((i_qnorm+10)*0.1)*1.e-18     ! based on LT1-v5.2 (range for Dm_max = 400000.)                   (for n_Qnorm = 50)
 Q_normalized = 800.**(0.2*(i_qnorm+5))*1.e-18      ! (range for Dm_max = 400000.)                                     (for n_Qnorm = 25)

 !--- from LT1:  [TO BE DELTED]
 !       !q = 261.7**((i_Qnorm+10)*0.1)*1.e-18     ! old (strict) lambda limiter
 !         q = 800.**((i_Qnorm+10)*0.1)*1.e-18     ! new lambda limiter               (for n_Qnorm = 50)
 !===

 return

 end function Q_normalized


!______________________________________________________________________________________

 real function lambdai(i)

 !-----------------
 ! Computes the lambda_i for a given index (used in loops to solve for PSD parameters)
 !-----------------

 implicit none

!arguments:
 integer, intent(in) :: i

!lambdai = 1.0013**i*100.   !used with Dm_max =   2000.
 lambdai = 1.0013**i*10.    !used with Dm_max = 400000.

 return

 end function lambdai


