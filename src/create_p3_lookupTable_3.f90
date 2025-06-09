PROGRAM create_p3_lookuptable_3

!______________________________________________________________________________________
!
! This program creates the lookup tables (that are combined into the single file
! 'p3_lookupTable_3.dat') for compute mui with triple-moment ice
!
!--------------------------------------------------------------------------------------
! Version:       1.0
! Last modified: 2025 May
!______________________________________________________________________________________

!______________________________________________________________________________________
!
! To generate 'p3_lookupTable_1.dat' using this code, do the following steps :
!              gfortran -fdefault-real-8 create_p3_lookupTable_3.f90)
!______________________________________________________________________________________

!--------------------------------------------------------------------------------------------

 implicit none

 !-----
 character(len=20), parameter :: version   = '1.4'
 logical, parameter           :: log_3momI = .true.    !switch to create table for 2momI (.false.) or 3momI (.true.)
 !-----

 integer            :: i_Znorm         ! index for normalized (by Q) Z (passed in through script; [1 .. n_Znorm])
 integer            :: i_rhor          ! index for rho_rime (passed in through script; [1 .. n_rhor])
 integer            :: i_Fr            ! index for rime-mass-fraction loop      [1 .. n_Fr]
 integer            :: i_Fl            ! index for liquid-mass-fraction loop    [1 .. n_Fl]
 integer            :: i_Qnorm         ! index for normalized (by N) Q loop     [1 .. n_Qnorm]

! NOTE: n_Znorm (number of i_Znorm values) is currently equal to 80.

 integer, parameter :: n_Znorm   = 80  ! number of indices for i_Znorm loop           (1nd "loop")  [not used in parallelized version]
 integer, parameter :: n_rhor    =  5  ! number of indices for i_rhor  loop           (2nd "loop")
 integer, parameter :: n_Fr      =  4  ! number of indices for i_Fr    loop           (3rd loop)
 integer, parameter :: n_Fl      =  4  ! number of indices for i_Fl    loop           (4th loop)
 integer, parameter :: n_Qnorm   = 50  ! number of indices for i_Qnorm loop           (5th loop)

 real, parameter :: mu_i_min = 0.
 real, parameter :: mu_i_max = 20.

 integer, parameter :: num_int_bins      = 40000 ! number of bins for numerical integration of ice processes

 integer :: i,ii,iii,jj,kk,kkk,dumii,i_iter,n_iter_psdSolve

 real :: N,q,qdum,dum1,dum2,cs1,ds1,lam,n0,lamf,qerror,del0,c0,c1,c2,dd,ddd,sum1,sum2,   &
         sum3,sum4,xx,a0,b0,a1,b1,dum,bas1,aas1,aas2,bas2,gammq,gamma,d1,d2,delu,lamold, &
         cap,lamr,dia,amg,dv,n0dum,sum5,sum6,sum7,sum8,dg,cg,bag,aag,dcritg,dcrits,      &
         dcritr,Fr,csr,dsr,duml,dum3,rhodep,cgpold,m1,m2,m3,dt,mu_r,initlamr,lamv,       &
         rdumii,lammin,lammax,pi,g,p,t,rho,mu,mu_i,ds,cs,bas,aas,dcrit,mu_dum,gdum,      &
         Z_value,sum9,mom3,mom6,intgrR1,intgrR2,intgrR3,intgrR4,dum4,cs2,ds2,mur_constant, &
         boltzman,meanpath,Daw,Dai,wcc,icc,Re,diffin,sc,st,aval,st2,Effw,Effi,eiaw,eiai,  &
         mom6_old,mom3_old,zqerror,zqerror_old

! New parameters with liquid fraction
 real :: area,area1,area2,mass,fac1,fac2,dumfac1,dumfac2,dumfac12,dumfac22,capm,gg,      &
         lamd,mu_id,n0d,qid,Zid,cs5,dum5,rhom,intgrR5,Fl

! function to compute mu for triple moment
 real :: compute_mu_3moment,compute_mu_3moment_2,compute_mu_3moment_1

! function to return diagnostic value of shape paramter, mu_i (mu_id with i_Fl)
 real :: diagnostic_mui

! function to return diagnostic value of shape paramter, mu_i
 real :: diagnostic_mui_Fl

! outputs from lookup table (i.e. "mui" read by access_lookup_table in s/r p3_main)
 real, dimension(n_Znorm,n_rhor,n_Qnorm,n_Fr,n_Fl) :: mu_i_save,rhomm,momen3,momen3_save

! for M6 rates
 real, dimension(n_rhor)            :: cgp,crp

 real, parameter                    :: Dm_max1 =  5000.e-6   ! max. mean ice [m] size for lambda limiter
 real, parameter                    :: Dm_max2 = 20000.e-6   ! max. mean ice [m] size for lambda limiter
 real, parameter                    :: Dm_min  =     2.e-6   ! min. mean ice [m] size for lambda limiter

 real, parameter                    :: thrd = 1./3.
 real, parameter                    :: sxth = 1./6.
 real, parameter                    :: cutoff = 1.e-90

 character(len=2014)                :: filename

!===   end of variable declaration ===

 pi  = acos(-1.)

! set constants and parameters
 dd   =  2.e-6 ! bin width for numerical integration of ice processes (units of m)

!--- specified mass-dimension relationship (cgs units) for unrimed crystals:

! ms = cs*D^ds
!
! for graupel:
! mg = cg*D^dg     no longer used, since bulk volume is predicted
!===

!---- Choice of m-D parameters for large unrimed ice:

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

!-- ** note: if using one of the above (.i.e. not brown and francis, which is already in mks units),
!           then uncomment line below to convert from cgs units to mks
!      cs=cs*100.**ds/1000.
!==

! Brown and Francis (1995)
 ds = 1.9
!cs = 0.01855 ! original (pre v2.3), based on assumption of Dmax
 cs = 0.0121 ! scaled value based on assumtion of Dmean from Hogan et al. 2012, JAMC

!====

! specify m-D parameter for fully rimed ice
!  note:  cg is not constant, due to variable density
 dg = 3.

 dcrit = (pi/(6.*cs)*900.)**(1./(ds-3.))

!.........................................................

! alpha parameter of m-D for rimed ice
 crp(1) =  50.*pi*sxth
 crp(2) = 250.*pi*sxth
 crp(3) = 450.*pi*sxth
 crp(4) = 650.*pi*sxth
 crp(5) = 900.*pi*sxth
!------------------------------------------------------------------------

! open file to write to lookup table:
 if (log_3momI) then
    write (filename, "(A17)") "lookupTable_3.dat"
    filename = trim(filename)
    open(unit=1, file=filename, status='unknown')
 endif

!--
! The values of i_Znorm (and possibly i_rhor) are "passed in" for parallelized version of code for 3-moment.
! The values of i_rhor are "passed in" for parallelized version of code for 2-moment.
! Thus, the loops 'i_Znorm_loop' and 'i_rhor_loop' are commented out accordingingly.
!
!i_Znorm_loop: do i_Znorm = 1,n_Znorm   !normally commented (kept to illustrate the structure (and to run in serial)
!   i_rhor_loop: do i_rhor = 1,n_rhor    !COMMENT OUT FOR PARALLELIZATION OF THIS LOOP (2-MOMENT ONLY)
     i_Fr_loop_1: do i_Fr = 1,n_Fr      !COMMENT OUT FOR PARALLELIZATION OF THIS LOOP (2-MOMENT ONLY)


! 3-moment-ice only:
! compute Z value from input Z index whose value is "passed in" through the script
     Z_value = 2.1**(i_Znorm)*1.e-23 ! range from 2x10^(-23) to 600 using 80 values
     !Z_value = 2.1e-05 ! for testing

       ! write header to first file:
       if (log_3momI .and. i_Znorm==1 .and. i_rhor==1 .and. i_Fr==1) then
          write(1,*) 'LOOKUP_TABLE_3-version:  ',trim(version),'-3momI'
          write(1,*)
       endif

!-- these lines to be replaced by Fr(i_Fr) initialization outside of loops
!  OR:  replace with: Fr = 1./float(n_Fr-1)
       if (i_Fr.eq.1) Fr = 0.
       if (i_Fr.eq.2) Fr = 0.333
       if (i_Fr.eq.3) Fr = 0.667
       if (i_Fr.eq.4) Fr = 1.
!==

       i_Fl_loop_1: do i_Fl = 1,n_Fl   !  loop for liquid mass fraction, Fl

          if (i_Fl.eq.1) Fl = 0.
          if (i_Fl.eq.2) Fl = 0.333
          if (i_Fl.eq.3) Fl = 0.667
          if (i_Fl.eq.4) Fl = 1.

          !print*, Z_value,i_Fr,i_Fl
! calculate mass-dimension relationship for partially-rimed crystals
! msr = csr*D^dsr
! formula from P3 Part 1 (JAS)

! dcritr is critical size separating fully-rimed from partially-rime ice

       cgp(i_rhor) = crp(i_rhor)  ! first guess

       if (i_Fr.eq.1) then   ! case of no riming (Fr = 0), then we need to set dcrits and dcritr to arbitrary large values

          dcrits = 1.e+6
          dcritr = dcrits
          csr    = cs
          dsr    = ds

       elseif (i_Fr.eq.2.or.i_Fr.eq.3) then  ! case of partial riming (Fr between 0 and 1)

          do
             dcrits = (cs/cgp(i_rhor))**(1./(dg-ds))
             dcritr = ((1.+Fr/(1.-Fr))*cs/cgp(i_rhor))**(1./(dg-ds))
             csr    = cs*(1.+Fr/(1.-Fr))
             dsr    = ds
           ! get mean density of vapor deposition/aggregation grown ice
             rhodep = 1./(dcritr-dcrits)*6.*cs/(pi*(ds-2.))*(dcritr**(ds-2.)-dcrits**(ds-2.))
           ! get density of fully-rimed ice as rime mass fraction weighted rime density plus
           ! density of vapor deposition/aggregation grown ice
             cgpold      = cgp(i_rhor)
             cgp(i_rhor) = crp(i_rhor)*Fr+rhodep*(1.-Fr)*pi*sxth
             if (abs((cgp(i_rhor)-cgpold)/cgp(i_rhor)).lt.0.01) goto 115
          enddo
115       continue

       else  ! case of complete riming (Fr=1.0)

        ! set threshold size between partially-rimed and fully-rimed ice as arbitrary large
          dcrits = (cs/cgp(i_rhor))**(1./(dg-ds))
          dcritr = 1.e+6       ! here is the "arbitrary large"
          csr    = cgp(i_rhor)
          dsr    = dg

       endif

!---------------------------------------------------------------------------------

     ! loop around normalized Q (Qnorm)
       i_Qnorm_loop: do i_Qnorm = 1,n_Qnorm

       ! lookup table values of normalized Qnorm
       ! (range of mean mass diameter from ~ 1 micron to x cm)

         !q = 261.7**((i_Qnorm+10)*0.1)*1.e-18    ! old (strict) lambda limiter
          q = 800.**((i_Qnorm+10)*0.1)*1.e-18     ! new lambda limiter

!--uncomment to test and print proposed values of qovn
!         print*,i_Qnorm,(6./(pi*500.)*q)**0.3333
!      enddo
!      stop
!==

! test values
!  N = 5.e+3
!  q = 0.01e-3

        ! initialize qerror to arbitrarily large value:
          qerror = 1.e+20

!.....................................................................................
! Find parameters for gamma distribution

! The size distribution of wet ice N(D_p) for all other processes
! Note D_p is full particle diameter (containing liquid and dry ice)
! Variable will be call mu_i, lam, n0 and are found using q (qi,tot)
! See Cholette et al. 2019 for details
! Note that mu_id is needed to compute diagnostic mu_i

   ! size distribution for ice is assumed to be
   ! N(D_p) = n0 * D_p^mu_i * exp(-lam*D_p)

   ! for the given q and N, we need to find n0, mu_i, and lam

   ! approach for finding lambda:
   ! cycle through a range of lambda, find closest lambda to produce correct q

   ! start with lam, range of lam from 100 to 1 x 10^7 is large enough to
   ! cover full range over mean size from approximately 1 micron to x cm

   ! compute mean density assuming rho_dry is cgp(i_rhor)
   ! rhomdry = cgp(i_rhor) (for 2momI only)
    rhom = (1.-Fl)*cgp(i_rhor)+Fl*1000.*pi*sxth
  !print*,q,rhom,Fl,cgp(i_rhor)
   ! assign provisional values for mom3 (first guess for mom3)
   ! NOTE: these are normalized: mom3 = M3/M0, mom6 = M6/M3 (M3 = 3rd moment, etc.)
    mom3 = q/rhom     !note: cgp is pi/6*(mean_density), computed above
   ! update normalized mom6 based on the updated ratio of normalized mom3 and normalized Q
   ! (we want mom6 normalized by mom3 not q)
    dum = mom3/q
    mom6 = Z_value/dum
    n_iter_psdSolve = 0

    do

          !print*,Z_value,q,rhom,mom3,mom6
          ! first estimate of mu_i
          mu_i = compute_mu_3moment_1(mom3,mom6,mu_i_max)
          mu_i = max(mu_i,mu_i_min)  ! make sure mu_i >= 0 (otherwise size dist is infinity at D = 0)
          mu_i = min(mu_i,mu_i_max)  ! set upper limit

          ii_loop_1: do ii = 1,11000 ! this range of ii to calculate lambda chosen by trial and error for the given lambda limiter values

           ! lam = 1.0013**ii*100.   ! old (strict) lambda_i limiter
             lam = 1.0013**ii*10.    ! new lambda_i limiter

           ! for lambda limiter:
             dum = Dm_max1+Dm_max2*Fr**2.
             lam = max(lam,(mu_i+1.)/dum)     ! set min lam corresponding to mean size of x
             lam = min(lam,(mu_i+1.)/Dm_min)  ! set max lam corresponding to mean size of Dm_min (2 micron)

           ! normalized n0:
             n0 = lam**(mu_i+1.)/(gamma(mu_i+1.))

           ! calculate integral for each of the 4 parts of the size distribution
           ! check difference with respect to Qnorm

           ! set up m-D relationship for solid ice with D < Dcrit:
             cs1  = pi*sxth*900.
             ds1  = 3.
             cs5  = pi*sxth*1000.

             call intgrl_section_Fl(lam,mu_i, ds1,ds,dg,dsr, dcrit,dcrits,dcritr,intgrR1,intgrR2,intgrR3,intgrR4,intgrR5)
           ! intgrR1 is integral from 0 to dcrit       (solid ice)
           ! intgrR2 is integral from dcrit to dcrits  (unrimed large ice)
           ! intgrR3 is integral from dcrits to dcritr (fully rimed ice)
           ! intgrR4 is integral from dcritr to inf    (partially rimed)

           ! sum of the integrals from the 4 regions of the size distribution:
             qdum = n0*((1.-Fl)*(cs1*intgrR1 + cs*intgrR2 + cgp(i_rhor)*intgrR3 + csr*intgrR4)+Fl*cs5*intgrR5)

             if (ii.eq.1) then
                qerror = abs(q-qdum)
                lamf   = lam
             endif

           ! find lam with smallest difference between Qnorm and estimate of Qnorm, assign to lamf
             if (abs(q-qdum).lt.qerror) then
                lamf   = lam
                qerror = abs(q-qdum)
             endif

          enddo ii_loop_1

        ! check and print relative error in q to make sure it is not too large
        ! note: large error is possible if size bounds are exceeded!!!!!!!!!!
        ! print*,'qerror (%)',qerror/q*100.

        ! find n0 based on final lam value
        ! set final lamf to 'lam' variable
        ! this is the value of lam with the smallest qerror
          lam = lamf

        ! n0 = N*lam**(mu_i+1.)/(gamma(mu_i+1.))

        ! find n0 from lam and Qnorm:
        !   (this is done instead of finding n0 from lam and N, since N;
        !    may need to be adjusted to constrain mean size within reasonable bounds)

          call intgrl_section_Fl(lam,mu_i, ds1,ds,dg,dsr, dcrit,dcrits,dcritr,intgrR1,intgrR2,intgrR3,intgrR4,intgrR5)
          n0   = q/((1.-Fl)*(cs1*intgrR1 + cs*intgrR2 + cgp(i_rhor)*intgrR3 + csr*intgrR4)+Fl*cs5*intgrR5)

        ! print*,'lam,N0,mu:',lam,n0,mu_i

        ! calculate normalized mom3 directly from PSD parameters (3-moment-ice only)
          !mom3_old = mom3
          mom3 = n0*gamma(4.+mu_i)/lam**(4.+mu_i)
        ! update normalized mom6 based on the updated ratio of normalized mom3 and normalized Q
        ! (we want mom6 normalized by mom3 not q)
          dum  = mom3/q
          mom6_old = mom6
          mom6 = Z_value/dum
          zqerror_old = zqerror
          zqerror = abs((mom6-mom6_old)/(mom6_old))
          if (abs((mom6-mom6_old)/(mom6_old)).lt.1e-16 .or. abs((zqerror-zqerror_old)/zqerror_old).lt.0.01) goto 116
          n_iter_psdSolve=n_iter_psdSolve+1
          if (n_iter_psdSolve.gt.50) goto 116

       enddo

116 continue

       mu_i_save(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl) = mu_i
       momen3_save(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl) = mom3

!.....................................................................................
! mass-weighted density
!.....................................................................................

! assume conditions for t and p as assumed above (giving rhos), then in microphysics scheme
! multiply by density correction factor (rhos/rho)^0.54, from Heymsfield et al. 2006

! fallspeed formulation from Mitchell and Heymsfield 2005

           ! initialize for numerical integration
          sum1 = 0.
          sum2 = 0.

        ! numerically integrate over size distribution
          ii_loop_2: do ii = 1,num_int_bins

             dum = real(ii)*dd - 0.5*dd   ! particle size

            !assign mass-size parameters (depending on size at ii)
             if (dum.le.dcrit) then
                ds1 = 3.
                cs1 = pi*sxth*900.
             else if (dum.gt.dcrit.and.dum.le.dcrits) then
                ds1 = ds
                cs1 = cs
             elseif (dum.gt.dcrits.and.dum.le.dcritr) then
                ds1 = dg
                cs1 = cgp(i_rhor)
             elseif (dum.gt.dcritr) then
                ds1 = dsr
                cs1 = csr
             endif

        ! These processes assume the liquid and the dry components of ice (add Fl)
        ! See Cholette et al. 2019 for details
             mass = (1.-Fl)*cs1*dum**ds1+Fl*pi*sxth*1000.*dum**3.

            !denominator of mass-weighted V:
             sum2 = sum2+mass*dum**(mu_i)*exp(-lam*dum)*dd

            !numerator of mass-weighted density:
            ! particle density is defined as mass divided by volume of sphere with same D
             sum1 = sum1+mass**2/(pi*sxth*dum**3)*dum**mu_i*exp(-lam*dum)*dd

          enddo ii_loop_2

        ! save mean fallspeeds for lookup table:

          rhomm(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl) = sum1/sum2
          momen3(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl) = 6.*q/(pi*rhomm(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl))

       print*,i_Znorm,i_rhor,i_Fr,i_Fl,i_Qnorm,mu_i_save(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl), &
       rhomm(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl),n_iter_psdSolve,zqerror
       print*,i_Znorm,i_rhor,i_Fr,i_Fl,i_Qnorm,momen3_save(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl), &
       momen3(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl)

       enddo i_Qnorm_loop


         !-- ice table
         i_Qnorm_loop_2:  do i_Qnorm = 1,n_Qnorm
             
          rhomm(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl)     = dim( rhomm(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl),     cutoff)

             write(1,'(5i5,4e15.5)')                             &
                         i_Znorm,i_rhor,i_Fr,i_Fl,i_Qnorm,        &
                         mu_i_save(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl), &
                         rhomm(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl), &
                         momen3_save(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl), &
                         momen3(i_Znorm,i_rhor,i_Qnorm,i_Fr,i_Fl)


         enddo i_Qnorm_loop_2

        enddo i_Fl_loop_1
      enddo i_Fr_loop_1
!    enddo i_rhor_loop
!enddo i_Znorm_loop
!==

 close(1)

END PROGRAM create_p3_lookuptable_3
!______________________________________________________________________________________

! Incomplete gamma function
! from Numerical Recipes in Fortran 77: The Art of
! Scientific Computing

      function gammq(a,x)

      real a,gammq,x

! USES gcf,gser
! Returns the incomplete gamma function Q(a,x) = 1-P(a,x)

      real gammcf,gammser,gln
!     if (x.lt.0..or.a.le.0) pause 'bad argument in gammq'
      if (x.lt.0..or.a.le.0) print*, 'bad argument in gammq'
      if (x.lt.a+1.) then
         call gser(gamser,a,x,gln)
         gammq=1.-gamser
      else
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      end if
      return
      end

!______________________________________________________________________________________

      subroutine gser(gamser,a,x,gln)
      integer itmax
      real a,gamser,gln,x,eps
      parameter(itmax=100,eps=3.e-7)
      integer n
      real ap,del,sum,gamma
      gln = log(gamma(a))
      if (x.le.0.) then
!        if (x.lt.0.) pause 'x < 0 in gser'
         if (x.lt.0.) print*, 'x < 0 in gser'
         gamser = 0.
         return
      end if
      ap=a
      sum=1./a
      del=sum
      do n=1,itmax
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if (abs(del).lt.abs(sum)*eps) goto 1
      end do
!     pause 'a too large, itmax too small in gser'
      print*, 'a too large, itmax too small in gser'
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      end

!______________________________________________________________________________________

      subroutine gcf(gammcf,a,x,gln)
      integer itmax
      real a,gammcf,gln,x,eps,fpmin
      parameter(itmax=100,eps=3.e-7,fpmin=1.e-30)
      integer i
      real an,b,c,d,del,h,gamma
      gln=log(gamma(a))
      b=x+1.-a
      c=1./fpmin
      d=1./b
      h=d
      do i=1,itmax
         an=-i*(i-a)
         b=b+2.
         d=an*d+b
         if(abs(d).lt.fpmin) d=fpmin
         c=b+an/c
         if(abs(c).lt.fpmin) c=fpmin
         d=1./d
         del=d*c
         h = h*del
         if(abs(del-1.).lt.eps)goto 1
      end do
!     pause 'a too large, itmax too small in gcf'
      print*, 'a too large, itmax too small in gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h
      return
      end

!______________________________________________________________________________________

 real function compute_mu_3moment(mom3,mom6,mu_max)

 !--------------------------------------------------------------------------
 ! Computes mu as a function of G(mu), where
 !
 ! G(mu)= N*Z/Q^2 = [(6+mu)(5+mu)(4+mu)]/[(3+mu)(2+mu)(1+mu)]
 !
 ! 2018-08-08
 !--------------------------------------------------------------------------

 implicit none

! Arguments passed:
 real, intent(in) :: mom3,mom6 !normalized moments
 real, intent(in) :: mu_max    !maximum allowable value of mu

! Local variables:
 real             :: mu   ! shape parameter in gamma distribution
 real             :: a1,g1,g2,G

! calculate G from normalized moments
    G = mom6/mom3

!----------------------------------------------------------!
! !Solve alpha numerically: (brute-force)
!      mu= 0.
!      g2= 999.
!      do i=0,4000
!         a1= i*0.01
!         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
!         if(abs(g-g1)<abs(g-g2)) then
!            mu = a1
!            g2= g1
!         endif
!      enddo
!----------------------------------------------------------!

!Piecewise-polynomial approximation of G(mu) to solve for mu:
  if (G >= 20.) then
    mu = 0.
  else
    g2 = G**2
    if (G<20.  .and.G>=13.31) mu = 3.3638e-3*g2 - 1.7152e-1*G + 2.0857e+0
    if (G<13.31.and.G>=7.123) mu = 1.5900e-2*g2 - 4.8202e-1*G + 4.0108e+0
    if (G<7.123.and.G>=4.200) mu = 1.0730e-1*g2 - 1.7481e+0*G + 8.4246e+0
    if (G<4.200.and.G>=2.946) mu = 5.9070e-1*g2 - 5.7918e+0*G + 1.6919e+1
    if (G<2.946.and.G>=1.793) mu = 4.3966e+0*g2 - 2.6659e+1*G + 4.5477e+1
    if (G<1.793.and.G>=1.405) mu = 4.7552e+1*g2 - 1.7958e+2*G + 1.8126e+2
    if (G<1.405.and.G>=1.230) mu = 3.0889e+2*g2 - 9.0854e+2*G + 6.8995e+2
    if (G<1.230) mu = mu_max
  endif

  compute_mu_3moment = mu

 end function compute_mu_3moment

 !______________________________________________________________________________________

 real function compute_mu_3moment_1(mom3,mom6,mu_max)

 !--------------------------------------------------------------------------
 ! Computes mu as a function of G(mu), where
 !
 ! G(mu)= N*Z/Q^2 = [(6+mu)(5+mu)(4+mu)]/[(3+mu)(2+mu)(1+mu)]
 !
 ! 2018-08-08
 !--------------------------------------------------------------------------

 implicit none

! Arguments passed:
 real, intent(in) :: mom3,mom6 !normalized moments
 real, intent(in) :: mu_max    !maximum allowable value of mu

! Local variables:
 real             :: mu   ! shape parameter in gamma distribution
 real             :: a1,g1,g2,G
 integer          :: i

! calculate G from normalized moments
    G = mom6/mom3

!----------------------------------------------------------!
 !Solve alpha numerically: (brute-force)
      mu= 0.
      g2= 999.
      do i=0,4000
         a1= i*0.01
         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
         if(abs(g-g1)<abs(g-g2)) then
            mu = a1
            g2= g1
         endif
      enddo
!----------------------------------------------------------!
  mu = min(max(mu,0.),mu_max)

  compute_mu_3moment_1 = mu

 end function compute_mu_3moment_1

 !==========================================================================================!
 real function compute_mu_3moment_2(mom3,mom6,mu_max)

 !--------------------------------------------------------------------------
 ! Computes mu as a function of moments 0, 3, and 6 of the size distribution
 ! represented by N(D) = No*D^mu*e(-lambda*D).
 !
 ! * solution is done using an analytic cubic root *
 !
 ! For piecewise polynomial approximation solution, use 'compute_mu_3moment_1'
 ! (This is coded as seperate subroutines, rather than a single function with an option,
 ! to avoid a IF/THEN block since this is used in loops.)
 !
 ! note: moment 3 is not equal to the mass mixing ratio (due to variable density)
 !--------------------------------------------------------------------------

 implicit none

! arguments:
! real, intent(in) :: mom0    !0th moment
 real, intent(in) :: mom3    !3th moment  (note, normalized)
 real, intent(in) :: mom6    !6th moment  (note, normalized)
 real, intent(in) :: mu_max  !maximum allowable value of mu

! local:
 real             :: mu      !shape parameter in gamma distribution
 double precision :: G       !function of mu (see comments above)
 double precision :: g2,x1,x2,x3
 real, parameter  :: eps_m3 = 1.e-20

 real :: dum,c1,c2,c3,Q,R,aa,bb

     G = mom6/mom3
     ! set minimum on G, below this the analytic solution breaks down
     G = max(1.3, G)

    !analytic cubic root solution:
     dum = 1./(1.-G)
     c1  = (15.-6.*G)*dum
     c2  = (74.-11.*G)*dum
     c3  = (120.-6.*G)*dum
     Q   = (c1**2-3.*c2)/9.
     R   = (2.*c1**3-9.*c1*c2+27.*c3)/54.

     ! NOTE: R is always < 0, thus we take the following:

     aa = (abs(R)+sqrt(R**2-Q**3))**(1./3.)
     bb = Q/aa

     mu = aa+bb-c1*(1./3.)
     mu = min(max(mu,0.),mu_max)

     compute_mu_3moment_2 = mu

 end function compute_mu_3moment_2

!______________________________________________________________________________________


subroutine intgrl_section_Fl(lam,mu, d1,d2,d3,d4, Dcrit1,Dcrit2,Dcrit3,    &
                          intsec_1,intsec_2,intsec_3,intsec_4,intsec_5)
!-----------------
! Computes and returns partial integrals (partial moments) of ice PSD.
!-----------------

implicit none

!Arguments:
real, intent(in)  :: lam,mu, d1,d2,d3,d4, Dcrit1,Dcrit2,Dcrit3
real, intent(out) :: intsec_1,intsec_2,intsec_3,intsec_4,intsec_5

!Local:
real :: dum,gammq
!-----------------

!Region I -- integral from 0 to Dcrit1  (small spherical ice)
intsec_1 = lam**(-d1-mu-1.)*gamma(mu+d1+1.)*(1.-gammq(mu+d1+1.,Dcrit1*lam))

!Region II -- integral from Dcrit1 to Dcrit2  (non-spherical unrimed ice)
intsec_2 = lam**(-d2-mu-1.)*gamma(mu+d2+1.)*(gammq(mu+d2+1.,Dcrit1*lam))
dum      = lam**(-d2-mu-1.)*gamma(mu+d2+1.)*(gammq(mu+d2+1.,Dcrit2*lam))
intsec_2 = intsec_2-dum

!Region III -- integral from Dcrit2 to Dcrit3  (fully rimed spherical ice)
intsec_3 = lam**(-d3-mu-1.)*gamma(mu+d3+1.)*(gammq(mu+d3+1.,Dcrit2*lam))
dum      = lam**(-d3-mu-1.)*gamma(mu+d3+1.)*(gammq(mu+d3+1.,Dcrit3*lam))
intsec_3 = intsec_3-dum

!Region IV -- integral from Dcrit3 to infinity  (partially rimed ice)
intsec_4 = lam**(-d4-mu-1.)*gamma(mu+d4+1.)*(gammq(mu+d4+1.,Dcrit3*lam))

!Region V -- integral from 0 to infinity  (ice completely metled)
!because d1=3.
intsec_5 = lam**(-d1-mu-1.)*gamma(mu+d1+1.)

return

end subroutine intgrl_section_Fl

