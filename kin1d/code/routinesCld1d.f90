MODULE SUBS_CLD1D

!=============================================================================!
!  The following subroutines are used by the driving subroutine cld1d.ftn90.
!  (with option for 'rpnsl.f')
!
!  Last modified:  2015-04-13
!=============================================================================!

CONTAINS

     SUBROUTINE advec(x,w,H,H0,Hcld,nk,nkcld,dt)

!     Performs Euler FIT-BIS advection step.
!
!     Alternatively, performs semi-Lagrangian advection step by first setting up
!     subdomain (in-"cloud", where w>0) and then calls RPNSL (which appears
!     to require that the entire domain be passed).
!     Note: Subroutine 'rpnsl' does not work on AIX (Azur) due to DIVBYZERO errors.

      IMPLICIT NONE

!  PASSING PARAMETERS:
      integer, INTENT(IN)    :: nk, nkcld
      real,    INTENT(IN)    :: H, H0, Hcld, dt
      !real,    INTENT(INOUT) :: H0
      real, DIMENSION(nk), INTENT(IN)    :: w
      real, DIMENSION(nk), INTENT(INOUT) :: x

!  LOCAL PARAMETERS:
      integer :: k, k1, SCHEME
      real, DIMENSION(nkcld) :: xcld, wcld
      real  :: dx

      SCHEME = 1              !1=Euler; 2=semi-Lagrangian
      IF (SCHEME == 1) THEN

   !Euler FIT-BIS advection:
        dx = (H-H0)/(nk-1)
        do k= 1,nk-1
           if ((w(k)*dt/dx)>1.) then
              print*, 'CFL VIOLATION in advec;  Co: ',w(k)*dt/dx
              STOP
           endif
           x(k) = x(k) -w(k)*dt/(1.*dx) * (x(k)-x(k+1))
        enddo

       ELSE

   !RPN Semi-Lagrangian advection:
     !  !direction of w>0 ==>         + .. w=0 ..
     !
     !   sfc                        cldTop     domTop
     !   k=nk                                  k=1      !IN MAIN CLD1d
     !    +----+----+----+--  ...    -+----+----+
     !   k=1                        k=nkcld             !TO ADVECTION ROUTINE
     !
       !flip arrays (for "in-cld" part)
        do k=1,nk
           k1= k-nk+nkcld
           if (k1>0) then
           wcld(k1) = w(k)
           xcld(k1) = x(k)
           endif
        enddo
!        call rpnsl(xcld,wcld,dt,Hcld,H0,nkcld)   !Semi-Lagrangian advection
       !flip back
        do k= nk-nkcld+1,nk
            k1= k-(nk-nkcld)
            x(k)= xcld(k1)
        enddo

      ENDIF

      END SUBROUTINE advec

!=============================================================================!

      SUBROUTINE NewWprof2(w,w1,wmax,Htop,DIV,COMP,rho,dt,z,nk)

!     Htop = height of cloud top (where w goes to zero)

      IMPLICIT NONE

!  PASSING PARAMETERS:
      integer, INTENT(IN) :: nk
      real, INTENT(IN) :: wmax, Htop, dt
      real, DIMENSION(nk),   INTENT(IN)  :: rho, z
      real, DIMENSION(nk),   INTENT(OUT) :: w1, DIV, COMP
      real, DIMENSION(1,nk), INTENT(OUT) :: w

!  LOCAL PARAMETERS:
      integer         :: k, H
      real            :: Co
      real, parameter :: pi= 3.14159265


      Co= wmax*dt/(Htop/nk)
      if (Co>1.5) print*, '** WARNING: Courant number:',Co
      if (Co>5. ) stop
! NOTE: The above warnings for the Courant number are because too rapid adiabatic
! 	cooling in one time step results in problems in the microphyscis scheme
!	(NANs and crashes)

      DO k=1,nk
         IF (z(k)<Htop) THEN
            w(1,k)= wmax*sin((z(k)-z(nk))/(Htop-z(nk))*pi)
            DIV(k)= pi*wmax/(Htop-z(nk))*cos(pi*(z(k)-z(nk))/(Htop-z(nk)))
         ELSE
            w(1,k) = 0
            DIV(k) = 0
         ENDIF
         w(1,k)= max(0., w(1,k))
         w1(k) = w(1,k)
      ENDDO


!  Define COMP(k) [where TERM3 = -q*COMP; see notes 2000-10-04]
      COMP(1) = w(1,1)/rho(1)*(rho(1)-rho(2))/(z(1)-z(2))             !top
      do k = 2,nk-1                                                   !interior
        if (z(k)<Htop) then
          COMP(k)= w(1,k)/rho(k)*(rho(k-1)-rho(k+1))/(z(k-1)-z(k+1))
        else
          COMP(k)= 0.
        endif
      enddo
      COMP(nk) = w(1,nk)/rho(nk)*(rho(nk-1)-rho(nk))/(z(nk-1)-z(nk))  !bottom

      RETURN
      END SUBROUTINE NewWprof2

!=============================================================================!

      SUBROUTINE IntMass(QxInt,Qx,rho,z,ni,nk)

!  Integrates density-weighted Qx(k) from k=1 to nk:
!  (i.e. calculates vapor/hydrometeor density/content)

      IMPLICIT NONE

!  PASSING PARAMETERS:
      integer, INTENT(IN) :: ni, nk
      real, DIMENSION(ni,nk), INTENT(IN) :: Qx
      real, DIMENSION(nk), INTENT(IN) :: rho, z
      real, INTENT(OUT) :: QxInt

!  LOCAL PARAMETERS:
      integer :: k
      real, parameter :: sixth = 1./6.

      QxInt = 0.
      IF (mod(nk,2)==0) THEN
! 2-PT TRAPEZOID:
        do k = 1,nk-1
          QxInt = QxInt + 0.5*(z(k)-z(k+1))*(rho(k)*Qx(1,k)+rho(k+1)*Qx(1,k+1))
        enddo
      ELSE

! 3-PT SIMPSON:
        do k = 1,nk-2,2
          QxInt = QxInt + (z(k)-z(k+2))*sixth *(rho(k)*Qx(1,k)+4*rho(k+1)*Qx(1,k+1)+rho(k+2)*Qx(1,k+2))
        enddo

      ENDIF

      RETURN
      END SUBROUTINE IntMass

!=============================================================================!

      SUBROUTINE IntN(NxInt,Nx,z,ni,nk)

!  Integrates total number concentration from k=1 to nk:
!   (i.e. calculates the total number [per m2] of the column)

      IMPLICIT NONE

!  PASSING PARAMETERS:
      integer, INTENT(IN) :: ni, nk
      real, DIMENSION(ni,nk), INTENT(IN) :: Nx
      real, DIMENSION(nk), INTENT(IN) :: z
      real, INTENT(OUT) :: NxInt

!  LOCAL PARAMETERS:
      integer :: k


      NxInt = 0
      IF (mod(nk,2)==0) THEN

! 2-PT TRAPEZOID:
        do k = 1,nk-1
          NxInt = NxInt + 0.5*(z(k)-z(k+1))*(Nx(1,k)+Nx(1,k+1))
        enddo

      ELSE

! 3-PT SIMPSON:
        do k = 1,nk-2,2
          NxInt = NxInt + (z(k)-z(k+2))/6*(Nx(1,k)+4*Nx(1,k+1)+Nx(1,k+2))
        enddo

      ENDIF

      RETURN
      END SUBROUTINE IntN

!======================================================================!

 real FUNCTION gamma(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  real, intent(IN) :: xx

! LOCAL PARAMETERS:
  integer  :: j
  real*8   :: ser,stp,tmp,x,y,cof(6),gammadp


  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=dble(xx)
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammadp=tmp+log(stp*ser/x)
  gammadp= exp(gammadp)

  gamma  = sngl(gammadp)

 END FUNCTION gamma
!======================================================================!
! 1D semi-lagrangian with conservation constraint

! Authors A.Qaddouri and jean cote
! revised by Sylvie Gravel


!--------------------------------------------------------------------------!
! Modifications:  (Jason Milbrandt)
!
!  - some recoding to prevent DIV_BY_ZERO
!  - Note:  Bug correction pertaining to calculation of 'coefmoy' may
!           not be done correctly.  'coefmoy' was previously initialized
!           with expression where the denominator (betasum) could be equal
!           to zero; the new coding (2004-10-14) may not be correct, but
!           it no longer allows for DIVBYZERO.
!--------------------------------------------------------------------------!

       subroutine rpnsl(ptmp,u,tau,H,H0,dim)

      implicit none

!--------------------------------------*
!  VARIABLE   DESCRIPTION
!  --------   -----------
!  ptmp       array of advected variable
!  u          velocity array
!  tau        timestep
!  H          upper x-position [m]
!  H0         lower x-position [m]
!  dim        number of grid points
!--------------------------------------*

      integer dim,dimm,maxvois
      parameter (dimm=100)
      parameter (maxvois=2)
      real*4 x(dimm),u(dimm),p(dimm),pL(dimm),pH(dimm),dep(dimm)
      real*4 ptmp(dimm)
      real*4  H0,H,Ai,k1,tau,dx,ph0,alpha(dimm),t
      integer Niter, Nalpha,ip(dimm),initer,iter
      integer i ,j,imax,ita,i1,flag(dimm),flagsum
      real*4 umax,pi,xdis,XX(maxvois)
      real*4  x_x0,x_x1,x_x2,x_x3,d,dd
      real*4 sumetoil,cetoil,beta(dimm),coef(dimm)
      real*4 sum,sum2,sum1,surplus,sumbeta,coefmoy
      real*4 pth,pth0,pvrai(dimm),pmax,pmin
      real*4 minimum,maximum,sum10
      real*4 pH1(dimm),pL1(dimm),p1(dimm),maxp
      real*4  pH_pL, coefmax(dimm), qplus,qminus



! physical parameters

      pi=3.14159265358979

!     Ai=5.0
      Ai=10.0
!       if ((H-H0).eq.0.) then
!       	 print*,'HERE 01'
!       	  STOP !JAM
! 	endif
      k1= Ai*pi/(H-H0)

! time-step

!     tau=30.

! iterations ( number of time step)
!     Niter=60
      Niter=1

! iterations for trajectore calculation (<=2)
      Nalpha=2

! geometrical parameters:
! 	if((dim-1).eq.0) then
! 	  print*,'HERE 02',
!       	  STOP !JAM
! 	endif

      dx= (H-H0)/(dim-1)
      do i= 1,dim
        x(i) =H0+(i-1)*dx
      enddo

!  Initialization:
      do  i=1,dim
        alpha(i)= dx/5.
        ip(i)   = 0
        pL(i)   = 0.0
        pH(i)   = 0.0
        beta(i) =0.0
      enddo

!  Speed profile

       umax= 0.
       do i=1,dim
!Flip variable
         p(dim-i+1) = ptmp(i)
        if( u(i) .gt.  umax) then
          umax = u(i)
          imax=i
        endif
      enddo


!  field profile

!     do i = 1,dim
!      p(i)= 1.e-08*(x(i)-H0)*(H-x(i))*
!    $               cos(pi*(x(imax)-x(i))/(H-H0))**2
!       write(*,'(i5,e20.4)')i,p(i)
!     enddo

!  Invariant is not sum of p(i)*dx but sum of p(i)/u(i)*dx
       sum10 =0.0
       do i=2,dim-1

!          if(u(i).eq.0.)  then
!       	 print*,'HERE 03'
!       	  STOP !JAM
! 	endif

         sum10 = sum10+dx*p(i)/u(i)
       enddo

!  real solution and  real invariant after t=Niter*tau
!  tau is a time step

       t= Niter*tau
       dep(1)= H0
       do i=2,dim

!        if ((H-H0).eq.0.)  then
!       	 print*,'HERE 04'
!       	  STOP !JAM
! 	endif

         pth = pi*((x(i)-H0)/(H-H0))
         pth0= 2.* atan(exp(- k1*t) * tan(pth/2.))
! analytical departure point
!c
!        dep(i)= (pth0*(H-H0))/pi + H0
       enddo

!      do i=1,dim
!       pvrai(i)=   1.e-08* (dep(i)-H0)*(H-dep(i))
!    $            * cos(pi*(x(imax)-dep(i))/(H-H0))**2
!      enddo

!      sum2 =0.0
!      do i=2,dim-1
!        sum2 = sum2+dx*pvrai(i)/u(i)
!      enddo


!  Start semi-Lagrangian:

      do iter=1,Niter

!  Calculation of trajectory:
      do ita= 1,Nalpha
!       do i= 1, dim
        do i= 2, dim      !JAM
          xdis= alpha(i)

! 	  if(dx.eq.0.)  then
!       	   print*,'HERE 05'
!       	   STOP !JAM
!           endif

          ip(i)=i-int(xdis/dx)
          d=alpha(i)-(i-ip(i))*dx
          alpha(i)=tau/2.*((d/dx)*u(ip(i)-1)+((dx-d)/dx)*u(ip(i)))
        enddo
      enddo

      do i=2,dim
        dep(i)= x(i)-2*alpha(i)
      enddo


!  -- Interpolation --

!  Cubic interpolation:

      do i1=2,dim

! 	if (dx.eq.0.)  then
!       	 print*,'HERE 06'
!       	  STOP !JAM
! 	endif

        xdis= x(i1)-dep(i1)
        ip(i1)=i1-int(xdis/dx)
        d=x(ip(i1))-dep(i1)
        if((ip(i1) .eq. 2) .or. (ip(i1) .eq. dim) ) then
          pH(i1)=(d/dx)*p(ip(i1)-1)+((dx-d)/dx)*p(ip(i1))
        else

       x_x0=(2*dx-d)
       x_x1=(dx-d)
       x_x2= -d
       x_x3=-(dx+d)
       pH(i1)= x_x1*x_x2*x_x3/(- dx*2.*dx*3.*dx)*p(ip(i1)-2)    &
                  +x_x0*x_x2*x_x3/(dx*dx*2*dx)*p(ip(i1)-1)      &
                  +x_x0*x_x1*x_x3/(- 2.*dx*dx*dx)*p(ip(i1))     &
                  +x_x0*x_x1*x_x2/(3.*dx*2.*dx*dx)*p(ip(i1)+1)
       endif
      enddo

! end of cubic intepolation

! Linear interpolation:

      do i1= 2 ,dim
        xdis= x(i1)-dep(i1)
        ip(i1)=i1-int(xdis/dx)
        d=x(ip(i1))-dep(i1)
        pL(i1)=(d/dx)*p(ip(i1)-1)+((dx-d)/dx)*p(ip(i1))
      enddo

!  end of linear interpolation


!  Calculate coef : p= coef* pH+(1-coef)*pL

      coef(1)=1.
      do i=2 ,dim
        xdis= x(i)-dep(i)
        ip(i)=i-int(xdis/dx)
        XX(1)=p(ip(i)-1)
        XX(2)=p(ip(i))
!---orig:
!       pmax= maximum(XX,maxvois)
!       pmin= minimum(XX,maxvois)
! note: the above code has compiler issues...
       print*, '** ABORT in s/r RPNSL **'
       print*, ' search:  pmax= maximum(XX,maxvois)'
       stop
!===

! limiting procedure
       qplus= pmax-pL(i)
       qminus= pmin-pL(i)
       pH_pL= pH(i)-pL(i)

!        if(pH_pL.eq.0.)  then
!       	 print*,'HERE 07'
!       	  STOP !JAM
! 	endif

       if(pH_pL.gt.0.0) coefmax(i)=min(1.,qplus/pH_pL)
       if(pH_pL.lt.0.0) coefmax(i)=min(1.,qminus/pH_pL)
       if(pH_pL.eq.0.0)  coefmax(i)=0.0

!! JAM: (changed to prevent DIV_BY_ZERO)
       if (pH_pL.eq.0.0) then
          coefmax(i)=0.0
       else
          if(pH_pL.gt.0.0) coefmax(i)=min(1.,qplus/pH_pL)
          if(pH_pL.lt.0.0) coefmax(i)=min(1.,qminus/pH_pL)
       endif
!!

      enddo

! calculate beta

      sumetoil =0.0
      do i=2,dim-1

!        if(u(i).eq.0.)  then
!       	 print*,'HERE 08'
!       	  STOP !JAM
! 	endif

       beta(i)=(pH(i)-pL(i))/u(i)*dx
       sumetoil= sumetoil+coefmax(i)*beta(i)
      enddo
      dd=0.0
      do i=2,dim-1

!        if(u(i).eq.0.)  then
!       	 print*,'HERE 09'
!       	  STOP !JAM
! 	endif

       dd= dd+pL(i)/u(i)*dx
      enddo
      cetoil=sum10-dd

      do i=2,dim-1
        coef(i)=coefmax(i)
      enddo

      if(sumetoil .ne. cetoil) then
       if(cetoil .lt. sumetoil) then
        do i=1,dim
          beta(i)= -1.*beta(i)
        enddo
        cetoil = -1.*cetoil
       endif

      do i=2,dim-1
       if(beta(i) .le. 0.0) then
        coef(i)=coefmax(i)
        flag(i)=1
       else
        coef(i)=0.0
        flag(i)=0
       endif
      enddo

      flagsum =0
      do i=2,dim-1
        flagsum= flagsum+ flag(i)
      enddo

 51   continue

      initer =0
      surplus=0.0
      sumbeta=0.0
      flagsum=0

      do i=2,dim-1
       if(flag(i).eq.1) then
        surplus= surplus+coef(i)*beta(i)
       else
        sumbeta= sumbeta+ beta(i)
       endif
      enddo

      surplus= cetoil-surplus

!       if (sumbeta.eq.0.)  then
!         print*, 'HERE 10  ',dim !JAM
! 	STOP
!       endif

!      coefmoy= surplus/sumbeta

!! JAM: (changed to prevent DIV_BY_ZERO)
      if (sumbeta.ne.0.) then
         coefmoy= surplus/sumbeta
      else
         coefmoy= 0.      !TEST,  What should it be set it to?
      endif
!!

      do i=2,dim-1
       if( flag(i).eq.0) then
        if ( coefmoy .lt. coefmax(i) ) then
         coef(i)=coefmoy
        else
         coef(i)=coefmax(i)
         flag(i)=1
        endif
       endif
      enddo

      flagsum=0
      do i=2,dim-1
       flagsum= flagsum+ flag(i)
      enddo
      if((flagsum.ne.(dim-2)).and.(initer.le.5)) goto 51
      endif

! conservation field result;
      do i=2,dim-1
        p(i)=(coef(i)*pH(i)+(1.-coef(i))*pL(i))
      enddo

!  Reflip:
 	do i = 1, dim
   	  ptmp(dim-i+1) = p(i)
 	enddo


!  end iteration

      enddo


! calculated invariant

      sum1=0.0
      do i=2,dim-1

!       	if(u(i).eq.0.)  then
!       	 print*,'HERE 11'
!       	  STOP !JAM
!	 	endif

        sum1 = sum1+dx*p(i)/u(i)
      enddo

!     call  cmpar8
!    %        ( p, pvrai,dim,'result', 'sexact')


      return
      end subroutine rpnsl


      real*4 function maximum(XX,dim)
      integer dim
      real*4 XX(dim),g
        g= XX(1)
        do i=1,dim
          if (XX(i).gt.g) g=XX(i)
        enddo
        maximum=g
       return
       end function maximum


      real*4 function minimum(XX,dim)
        integer dim
        real*4 XX(dim),g
        g = XX(1)
        do i=1,dim
          if (XX(i).lt.g) g=XX(i)
        enddo
        minimum=g
       return
       end function minimum


      subroutine cmpar8( f2, f1,ni,var2, var1)
      implicit none
      integer maxi
      parameter (maxi=20000)
      real*4  wkd(maxi),zero
      parameter (zero=0.0)
      integer ni
      real*4 f1(ni),f2(ni)
      real*4 sumfd,ssqfd,sumf1,ssqf1,sumf2,ssqf2,fmax,fmin
      integer i,ijma,ijmi,jma,ima,jmi,imi
      real*4 mxnrm
      character*6 var1,var2
      sumfd = zero
      ssqfd = zero
      sumf1 = zero
      ssqf1 = zero
      sumf2 = zero
      ssqf2 = zero
         do i=1,ni
            wkd(i)=f2(i)-f1(i)
            sumfd = sumfd + wkd(i)
            ssqfd = ssqfd + wkd(i)
            sumf1 = sumf1 + f1(i)
            ssqf1 = ssqf1 + f1(i)
            sumf2 = sumf2 + f2(i)
            ssqf2 = ssqf2 + f2(i)
         enddo
!
! 	if(ni.eq.0.) then
! 		print*, 'HERE 12' !JAM
! 		STOP
! 	endif

      sumfd = sumfd/ni
      ssqfd = sqrt( abs(ssqfd)/ni )
      sumf1 = sumf1/ni
      ssqf1 = sqrt( abs(ssqf1)/ni )
      sumf2 = sumf2/ni
      ssqf2 = sqrt( abs(ssqf2)/ni )
!
      fmax = wkd(1)
      do i=2,ni
         if ( wkd(i) .gt. fmax ) then
            ijma = i
            fmax = wkd(i)
         endif
      enddo
      ijmi = 1
      fmin = wkd(1)
      do i=2,ni
         if ( wkd(i) .lt. fmin ) then
            ijmi = i
            fmin = wkd(i)
         endif
      enddo
      jma = ( ijma - 1 )/ni
      ima = ijma - jma * ni
      jma = jma + 1
      jmi = ( ijmi - 1 )/ni
      imi = ijmi - jmi * ni
      jmi = jmi + 1
!
   10 mxnrm = max( mxnrm, abs( fmin ), abs( fmax ) )
!
      return
      end subroutine cmpar8


!======================================================================!

      subroutine vertint2b (frp,frxp,frfin,frxin,nis,ni)

      implicit none

!-----------------------------------------------------------
!	VARIABLE	MEANING
!	========	=======
!	frp		interpolated value at desired level
!	frxp		desired level (or abscissa value)
!	frfin(ni)	array of uninterpolated funcion values
!	frxin(ni)	array of values of levels (abscissa)
!	nis		number of desired interpolated levels
!	ni		number of uninterpolated levels & function values
!-----------------------------------------------------------

      integer nis,ni,ind,ni2
      parameter (ni2 = 300)
      real frfin(ni),frxin(ni)
      real frf(ni2),frx(ni2),temp1(ni2),temp2(ni2)
      real frp,frxp
!
!OBJECT
!    This program determines the value frp at point frxp by cubic
!    interpolation of the function frf. this function is known at ni
!    points and to each point corresponds the coordinate frx.
!
!    This subroutine handles any point located outside the grid or
!    located between the boundary and the first interior point.
!
!METHOD
!
!EXTERNALS
!
!AUTHOR    Andre Robert                             Feb  1980
!
!HISTORY
!
!     Rene Laprise                                  Mar  1989
!           - cosmetiques
!     Yves Chartier et Michel Desgagne          Oct/Nov  1992
!           - implicit none
!           - nis,njs,nks,ni,nj,nk
!           - structural documentation
!           - comdeck
!           - working vectors memory allocation
!           - doctor2
!           - in lining
!     Jason Milbrandt				  Sept	2000
!           - optional descending order of frx
!

      integer i,pnia,pnnim
      real prda,prdb,prxd,prsaf,prsbf,prsad,prsbd
!-----------------------------------------------------------


! ++++++ JAM
      IF (frxin(ni).lt.frxin(1)) THEN		!assumed condition for descending levels
!       Reverse order of arrays:
         do ind = 1,ni
            temp1(ind) = frfin(ind)
            temp2(ind) = frxin(ind)
         enddo
         do ind = 1,ni
            frf(ind) = temp1(ni+1-ind)
            frx(ind) = temp2(ni+1-ind)
         enddo
      ELSE
         do ind = 1,ni
            frf(ind) = frfin(ind)
            frx(ind) = frxin(ind)
         enddo
      ENDIF
! ++++++


      pnnim= ni-1

      if (frxp.le.frx(1)) then
!
!        * persistance au dessus de la premiere couche.
!
         frp=frf(1)
         goto 99
!
      elseif (frxp.ge.frx(ni)) then
!
!        * persistance en dessous de la derniere couche.
!
         frp=frf(ni)
         goto 99
!
      elseif (frxp.le.frx(2)) then
!
!        * interpolation lineaire dans la premiere couche.
!
         prxd=(frxp-frx(1))/(frx(2)-frx(1))
         frp=(1.0-prxd)*frf(1)+prxd*frf(2)
         goto 99
      elseif (frxp.ge.frx(pnnim)) then
!
!        * interpolation lineaire dans la derniere couche.
!
         prxd=(frxp-frx(pnnim))/(frx(ni)-frx(pnnim))
         frp=(1.0-prxd)*frf(pnnim)+prxd*frf(ni)
         goto 99
!
      else
!
!        * interpolation cubique entre la base de la plus haute couche
!       * et le sommet de la plus basse couche.
!
!        * determine the location of the nearest point to the left.
!
         do 1 i=2,pnnim
            if (frxp.ge.frx(i)) pnia=i
 1       continue
!
!        * the calling program must ensure that frxp is between frx(2)
!        * and frx(pnnim) also,the coordinates frx must be stored in a
!        * monotonically increasing order.
!
!        * the next step consists in preparing the information required
!        * by the interpolation subroutine. the derivative is
!        * calculated at points a and b
!        * it must be noted that it requires a realless value
!        * of frxp.
!
         prxd=(frxp-frx(pnia))/(frx(pnia+1)-frx(pnia))
         prda=((frf(pnia+1)-frf(pnia-1))/(frx(pnia+1)-frx(pnia-1)))*  &
              (frx(pnia+1)-frx(pnia))
         prdb=((frf(pnia+2)-frf(pnia))/(frx(pnia+2)-frx(pnia)))*      &
              (frx(pnia+1)-frx(pnia))
!
!        * fits a cubic to the values frf(pnia) and frf(pnia+1) of the
!        * function frf at points pnia and pnia+1. this cubic also fits
!        * the derivatives prda and prdb of the function at both points.
!        * the value of frp is calculated at point frx and returned to
!        * the calling program.
!
!        * the value frx must vary from 0 to 1 from point a to
!        * point b. the four cubic splines prsaf,prsbf,prsad and prsbd
!        * are used for the interpolation.
!
         if (prxd.ge.0.0.and.prxd.le.1.0) go to 2
         write (6,600)
         stop
!
 2       continue
!
         prsaf=(1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
         prsbf=(3.0-2.0*prxd)*prxd*prxd
         prsad=prxd*(1.0-prxd)*(1.0-prxd)
         prsbd=(1.0-prxd)*prxd*prxd
!
!        * resultat
!
         frp=frf(pnia)*prsaf+frf(pnia+1)*prsbf+prda*prsad-prdb*prsbd
!
      endif
!


 99   continue
      return
!---------------------------------------------------------------------
  600 format (3x,'THE ARGUMENT XD GIVEN TO SUBROUTINE VERTINT MUST NOT be less than zero or greater than one.')
      end subroutine vertint2b
!======================================================================!

!=============================================================================!
END MODULE SUBS_CLD1D
