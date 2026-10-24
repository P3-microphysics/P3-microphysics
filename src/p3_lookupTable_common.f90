MODULE p3_lookupTable_common

!______________________________________________________________________________________
!
! Procedures shared by the P3 lookup-table generators create_p3_lookupTable_[123].f90.
!
! Compile this file first, with the same flags as the generator (double-precision default
! reals), e.g.
!
!    ifx -r8 p3_lookupTable_common.f90 create_p3_lookupTable_1.f90
!    gfortran -fdefault-real-8 -fdefault-double-8 p3_lookupTable_common.f90 create_p3_lookupTable_1.f90
!
! The driver src/make_luts.sh does this.
!______________________________________________________________________________________

 implicit none

 private
 public :: gammq, diagnostic_mui, intgrl_section, intgrl_section_Fl

 contains

!______________________________________________________________________________________

! Incomplete gamma function
! from Numerical Recipes in Fortran 77: The Art of
! Scientific Computing

      function gammq(a,x)

      real a,gammq,x

! USES gcf,gser
! Returns the incomplete gamma function Q(a,x) = 1-P(a,x)

      real gammcf,gamser,gln
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
      end function gammq

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
      end subroutine gser

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
      end subroutine gcf

!______________________________________________________________________________________

 real function diagnostic_mui(mu_i_min,mu_i_max,lam,q,cgp,Fr,pi)

!----------------------------------------------------------!
! Compute mu_i diagnostically.
!----------------------------------------------------------!

 implicit none

!Arguments:
 real    :: mu_i_min,mu_i_max,lam,q,Fr,pi
 real    :: cgp    ! coefficient of the spherical-ice mass-size relation m = cgp*D^3

! Local variables:
 real, parameter :: Di_thres = 0.2 !diameter threshold [mm]
!real, parameter :: Di_thres = 0.6 !diameter threshold [mm]
 real            :: mu_i,dum1,dum2,dum3


!-- diagnostic mu_i, original formulation: (from Heymsfield, 2003)
!  mu_i = 0.076*(lam/100.)**0.8-2.   ! /100 is to convert m-1 to cm-1
!  mu_i = max(mu_i,0.)
!  mu_i = min(mu_i,6.)

!-- formulation based on 3-moment results (see 2021 JAS article)
 dum1 = (q/cgp)**(1./3)*1000.              ! estimated Dmvd [mm], assuming spherical
 if (dum1<=Di_thres) then
    !diagnostic mu_i, original formulation: (from Heymsfield, 2003)
    mu_i = 0.076*(lam*0.01)**0.8-2.        ! /100 is to convert m-1 to cm-1
    mu_i = min(mu_i,6.)
 else
    dum2 = (6./pi)*cgp                     ! mean density (total)
    dum3 = max(1., 1.+0.00842*(dum2-400.)) ! adjustment factor for density
    mu_i = 0.25*(dum1-Di_thres)*dum3*Fr
   !mu_i = 4.*(dum1-Di_thres)*dum3*Fr
 endif
 mu_i = max(mu_i,mu_i_min)  ! make sure mu_i >= 0, otherwise size dist is infinity at D = 0
 mu_i = min(mu_i,mu_i_max)

 diagnostic_mui = mu_i

 end function diagnostic_mui

!______________________________________________________________________________________

 subroutine intgrl_section(lam,mu, d1,d2,d3,d4, Dcrit1,Dcrit2,Dcrit3,    &
                           intsec_1,intsec_2,intsec_3,intsec_4)
 !-----------------
 ! Computes and returns partial integrals (partial moments) of ice PSD.
 !-----------------

 implicit none

!Arguments:
 real, intent(in)  :: lam,mu, d1,d2,d3,d4, Dcrit1,Dcrit2,Dcrit3
 real, intent(out) :: intsec_1,intsec_2,intsec_3,intsec_4

!Local:
 real :: dum
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

 return

 end subroutine intgrl_section

!______________________________________________________________________________________

 subroutine intgrl_section_Fl(lam,mu, d1,d2,d3,d4, Dcrit1,Dcrit2,Dcrit3,    &
                              intsec_1,intsec_2,intsec_3,intsec_4,intsec_5)
 !-----------------
 ! Computes and returns partial integrals (partial moments) of ice PSD; as intgrl_section,
 ! plus the integral over the whole PSD of completely melted ice (intsec_5).
 !-----------------

 implicit none

!Arguments:
 real, intent(in)  :: lam,mu, d1,d2,d3,d4, Dcrit1,Dcrit2,Dcrit3
 real, intent(out) :: intsec_1,intsec_2,intsec_3,intsec_4,intsec_5

!Local:
 real :: dum
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

 !Region V -- integral from 0 to infinity  (ice completely melted)
 !because d1=3.
 intsec_5 = lam**(-d1-mu-1.)*gamma(mu+d1+1.)

 return

 end subroutine intgrl_section_Fl

!______________________________________________________________________________________

END MODULE p3_lookupTable_common
