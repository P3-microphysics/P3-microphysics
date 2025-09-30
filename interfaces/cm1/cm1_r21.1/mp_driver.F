  MODULE mp_driver_module

             !------  Microphysics Driver ------!

  implicit none

  private
  public :: mp_driver

  CONTAINS

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine mp_driver(nstep,dt,qbudget,asq,bsq,xh,ruh,xf,ruf,yh,rvh,yf,rvf,  &
                         mh,rmh,c1,c2,zh,mf,rmf,zf,rain,prate,pi0,th0,rho0,prs0,qv0,   &
                         rho,prs,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,  &
                         w3d,ppi,pp3d,ppten,sten,tha,th3d,thten,qa,q3d,qten,  &
                         p3a,p3o,dum2d1,dum2d2,dum2d3,dum2d4,dum2d5,  &
                         effc,effi,effs,effr,effg,effis,    &
                         tdiag,qdiag,out2d,out3d,tmp2d,tmp1d,  &
                         dowriteout,dorad,dotdwrite,doazimwrite,dorestart,  &
                         getdbz,getvt,dotbud,doqbud)

    use input
    use constants
    use index_module
    use misclibs
    use kessler_module
    use goddard_module, only : goddard,satadj_ice,consat2
    use module_mp_thompson , only : mp_gt_driver
    use lfoice_module, only : lfo_ice_drive,lfoice_init
    use module_mp_graupel , only : mp_graupel
    use module_mp_nssl_2mom, only : nssl_2mom_driver
    use microphy_p3, only : mp_p3_wrapper_wrf
    use module_mp_jensen_ishmael , only : mp_jensen_ishmael
    use mp_wsm6 , only : refl10cm_wsm6,mp_wsm6_init
    use mmm_physics_wrapper_module, only : mp_wsm6_wrapper,mp_wsm6_effectRad_wrapper

    implicit none

    integer, intent(in) :: nstep
    real, intent(inout) :: dt
    double precision, intent(inout), dimension(nbudget) :: qbudget
    double precision, intent(inout), dimension(numq) :: asq,bsq
    real, intent(in), dimension(ib:ie) :: xh,ruh
    real, intent(in), dimension(ib:ie+1) :: xf,ruf
    real, intent(in), dimension(jb:je) :: yh,rvh
    real, intent(in), dimension(jb:je+1) :: yf,rvf
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh,c1,c2,zh
    real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf,zf
    real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain
    real, intent(inout), dimension(ib:ie,jb:je) :: prate
    real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,rho0,prs0,qv0
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
    real, intent(in),    dimension(ib:ie,jb:je,kb:ke+1) :: w3d
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten,sten
    real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten
    real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
    real, intent(inout), dimension(ibp3:iep3,kbp3:kep3,1,np3a) :: p3a
    real, intent(inout), dimension(ibp3:iep3,jbp3:jep3,kbp3:kep3,np3o) :: p3o
      double precision, intent(inout), dimension(ib:ie,jb:je) :: dum2d1,dum2d2,dum2d3,dum2d4,dum2d5
    real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: effc,effi,effs,effr,effg,effis
    real, intent(inout) , dimension(ibdt:iedt,jbdt:jedt,kbdt:kedt,ntdiag) :: tdiag
    real, intent(inout) , dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
    real, intent(inout) , dimension(ib2d:ie2d,jb2d:je2d,nout2d) :: out2d
    real, intent(inout) , dimension(ib3d:ie3d,jb3d:je3d,kb3d:ke3d,nout3d) :: out3d
    real, intent(inout), dimension(ibtmp:ietmp,kbtmp:ketmp  ,ntmp2d ) :: tmp2d
    real, intent(inout), dimension(ibtmp:ietmp,ntmp1d) :: tmp1d
    logical, intent(in) :: dowriteout,dorad,dotdwrite,doazimwrite,dorestart
    logical, intent(inout) :: getdbz,getvt
    logical, intent(in) :: dotbud,doqbud

    !........

    integer :: i,j,k,n
    integer :: has_reqc,has_reqi,has_reqs,do_radar_ref
    real :: rdt

    ! for wsm6
    character(len=40) :: errmsg
    integer :: errflg
    real, dimension(kbtmp:ketmp) :: p1d,t1d,qv1d,qr1d,qs1d,qg1d,dbz1d,zp

    !........

        if( stopit ) getdbz = .true.

        if( getdbz .and. qd_dbz.le.0 )then
          getdbz = .false.
        endif

        IF( efall.eq.1 ) getvt = .true.
        IF( dowriteout .and. output_fallvel.eq.1 ) getvt = .true.

        ! dum1 = T
        ! dum3 = appropriate rho for budget calculations
        ! store copy of T in thten array:

        !$omp parallel do default(shared)  &
        !$omp private(i,j,k,n)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
          thten(i,j,k)=dum1(i,j,k)
          qten(i,j,k,nqv)=q3d(i,j,k,nqv)
        enddo
        enddo
        enddo

          if( getdbz .and. qd_dbz.ge.1 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,n)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              qdiag(i,j,k,qd_dbz)=0.0
            enddo
            enddo
            enddo
          endif

          IF(axisymm.eq.0)THEN
            ! for Cartesian grid:
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,n)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=rho(i,j,k)
            enddo
            enddo
            enddo
          ELSE
            ! for axisymmetric grid:
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,n)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum3(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
            enddo
            enddo
            enddo
          ENDIF

          if( dotbud .and. td_mp.ge.1 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,n)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              tdiag(i,j,k,td_mp) = th3d(i,j,k)
              dum5(i,j,k) = pi0(i,j,k)+pp3d(i,j,k)
            enddo
            enddo
            enddo
          endif
          if( doqbud .and. qd_mp.ge.1 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k,n)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              qdiag(i,j,k,qd_mp) = q3d(i,j,k,nqv)
            enddo
            enddo
            enddo
          endif


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  NOTES:
!
!           dum1   is   T
!           dum3   is   rho for budget calculations
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Kessler scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ifptype:  &
        IF(ptype.eq.1)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call k_fallout(rho,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          if(efall.ge.1)then
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              if( abs(dt*dum4(i,j,k)).ge.tsmall )then
                dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
                prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
            if( dotbud .and. td_efall.ge.1 )then
              !$omp parallel do default(shared)  &
              !$omp private(i,j,k)
              do k=1,nk-1
              do j=1,nj
              do i=1,ni
                tdiag(i,j,k,td_efall) = dum4(i,j,k)
              enddo
              enddo
              enddo
            endif
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call kessler(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,dum1,   &
                       rho,dum3,pp3d,th3d,prs,                            &
                       q3d(ib,jb,kb,nqv),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3))
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Goddard LFO scheme   cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.2)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call goddard(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,             &
                       rho,dum3,prs,pp3d,th3d,                            &
     q3d(ib,jb,kb,1), q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),qten(ib,jb,kb,3),   &
     q3d(ib,jb,kb,4),qten(ib,jb,kb,4),q3d(ib,jb,kb,5),qten(ib,jb,kb,5),   &
     q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          call satadj_ice(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,     &
                          rho,dum3,pp3d,prs,th3d,                     &
              q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),   &
              q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6))
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(efall.ge.1)then
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              if( abs(dt*dum4(i,j,k)).ge.tsmall )then
                dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
                prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
            if( dotbud .and. td_efall.ge.1 )then
              !$omp parallel do default(shared)  &
              !$omp private(i,j,k)
              do k=1,nk-1
              do j=1,nj
              do i=1,ni
                tdiag(i,j,k,td_efall) = dum4(i,j,k)
              enddo
              enddo
              enddo
            endif
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                       q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                       q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                       q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(getdbz .and. qd_dbz.ge.1 ) call calcdbz(rho,q3d(ib,jb,kb,3),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6),  &
                                  qdiag(ibdq,jbdq,kbdq,qd_dbz))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Thompson scheme   ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.3)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! cm1r17:  to make things easier to understand, use same arrays
            !          that are used for morrison code:
            ! dum1 = T  (this should have been calculated already)
            ! dum2 = pi (nondimensional pressure)
            ! dum4 = dz
            ! thten = copy of T  (this should have been calculated already)
            dum2(i,j,k)=pi0(i,j,k)+pp3d(i,j,k)
            dum4(i,j,k)=dz*rmh(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            has_reqc = 1
            has_reqi = 1
            has_reqs = 1
          else
            has_reqc = 0
            has_reqi = 0
            has_reqs = 0
          endif
          do_radar_ref = 1

          call   mp_gt_driver(qv=q3d(ib,jb,kb,1), &
                              qc=q3d(ib,jb,kb,2), &
                              qr=q3d(ib,jb,kb,3), &
                              qi=q3d(ib,jb,kb,4), &
                              qs=q3d(ib,jb,kb,5), &
                              qg=q3d(ib,jb,kb,6), &
                              ni=q3d(ib,jb,kb,7), &
                              nr=q3d(ib,jb,kb,8), &
                  t3d=dum1, pii=dum2, p=prs, w=w3d, dz=dum4, dt_in=dt, itimestep=nstep,                &
                  RAINNC=rain(ib,jb,1), RAINNCV=dum5(ib,jb,1), SR=dum5(ib,jb,2),                       &
                  refl_10cm=qdiag(ibdq,jbdq,kbdq,qd_dbz), diagflag=getdbz, do_radar_ref=do_radar_ref,  &
                  re_cloud=dum6, re_ice=dum7, re_snow=dum8,                  &
                  has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,   &
                  nrain=nrain,dx=dx, dy=dy, cm1dz=dz,                        &
                  tcond=qbudget(1),tevac=qbudget(2),                         &
                  tevar=qbudget(5),train=qbudget(6),                         &
                  ruh=ruh,rvh=rvh,rmh=rmh,rr=dum3,rain=rain,prate=prate,     &
                  ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d, &
                  nout3d=nout3d,out3d=out3d,eqtset=eqtset,getvt=getvt,       &
                  vtc=qten(ib,jb,kb,nqc),                                    &
                  vtr=qten(ib,jb,kb,nqr),                                    &
                  vti=qten(ib,jb,kb,nqi),                                    &
                  vts=qten(ib,jb,kb,nqs),                                    &
                  vtg=qten(ib,jb,kb,nqg),                                    &
                  ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,   &
                  ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,     &
                  its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk )

          ! Get final values for th3d,pp3d,prs:
          ! Note:  dum1 stores temperature, thten stores old temperature:
          IF( eqtset.eq.2 )THEN
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                prs(i,j,k)=rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))*dum1(i,j,k)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          ELSE
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
                rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
              endif
            enddo
            enddo
            enddo
          ENDIF
          if( radopt.eq.2 )then
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              effc(i,j,k) = dum6(i,j,k)
              effi(i,j,k) = dum7(i,j,k)
              effs(i,j,k) = dum8(i,j,k)
            enddo
            enddo
            enddo
          endif

          IF( getvt )THEN
          IF( dowriteout .or. dotdwrite .or. doazimwrite .or. dorestart )THEN
          if( qd_vtc.ge.1 .or. qd_vtr.ge.1 .or. qd_vts.ge.1 .or. qd_vtg.ge.1 .or. qd_vti.ge.1 )then
            if(myid.eq.0) print *,'  Getting fall velocities for Thompson MP scheme '
              if( qd_vtc.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vtc) = qten(i,j,k,nqc)
                enddo
                enddo
                enddo
              endif
              if( qd_vtr.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vtr) = qten(i,j,k,nqr)
                enddo
                enddo
                enddo
              endif
              if( qd_vts.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vts) = qten(i,j,k,nqs)
                enddo
                enddo
                enddo
              endif
              if( qd_vtg.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vtg) = qten(i,j,k,nqg)
                enddo
                enddo
                enddo
              endif
              if( qd_vti.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vti) = qten(i,j,k,nqi)
                enddo
                enddo
                enddo
              endif
          endif
          ENDIF
          ENDIF

          if(timestats.ge.1) time_microphy=time_microphy+mytime()

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   GSR LFO scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.4)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call lfo_ice_drive(dt, mf, pi0, prs0, pp3d, prs, th0, th3d,    &
                             qv0, rho0, q3d, qten, dum1)
          do n=2,numq
            call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                         q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
          enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Morrison scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.5)THEN

!          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
!          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
!          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
!          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
!          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
!          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call pdefqtest(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1),dum2d1,dum2d2,dum2d3,dum2d4,dum2d5)
          call pdefqtest( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2),dum2d1,dum2d2,dum2d3,dum2d4,dum2d5)
          call pdefqtest( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3),dum2d1,dum2d2,dum2d3,dum2d4,dum2d5)
          call pdefqtest( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4),dum2d1,dum2d2,dum2d3,dum2d4,dum2d5)
          call pdefqtest( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5),dum2d1,dum2d2,dum2d3,dum2d4,dum2d5)
          call pdefqtest( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6),dum2d1,dum2d2,dum2d3,dum2d4,dum2d5)
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dum1 = T  (this should have been calculated already)
            ! dum4 = dz
            ! thten = copy of T  (this should have been calculated already)
            dum4(i,j,k)=dz*rmh(i,j,k)
          enddo
          enddo
          enddo
          IF(numq.eq.11)THEN
            ! ppten stores ncc:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              ppten(i,j,k) = q3d(i,j,k,11)
            enddo
            enddo
            enddo
          ENDIF
          ! cm1r17:  get fall velocities (store in qten array)
          call MP_GRAUPEL(nstep,dum1,dum5,                            &
                          q3d(ib,jb,kb, 1),q3d(ib,jb,kb, 2),q3d(ib,jb,kb, 3), &
                          q3d(ib,jb,kb, 4),q3d(ib,jb,kb, 5),q3d(ib,jb,kb, 6), &
                          q3d(ib,jb,kb, 7),q3d(ib,jb,kb, 8),q3d(ib,jb,kb, 9), &
                          q3d(ib,jb,kb,10),ppten,                             &
                               prs,rho,dt,dum4,w3d,rain,prate,                &
                          effc,effi,effs,effr,effg,effis,                     &
                          qbudget(1),qbudget(2),qbudget(5),qbudget(6),        &
                          ruh,rvh,rmh,dum3,getdbz,                            &
                          qten(ib,jb,kb,nqc),qten(ib,jb,kb,nqr),qten(ib,jb,kb,nqi),  &
                          qten(ib,jb,kb,nqs),qten(ib,jb,kb,nqg),getvt,dorad,  &
                          dotbud,doqbud,tdiag,qdiag,out3d)
          IF(numq.eq.11)THEN
            ! ppten stores ncc:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              q3d(i,j,k,11) = ppten(i,j,k)
            enddo
            enddo
            enddo
          ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          IF(efall.eq.1)THEN
            ! dum1 = T
            ! dum2 = cvm
            ! dum4 = T tendency
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqc),qten(ib,jb,kb,nqc))
            call getefall(0,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqr),qten(ib,jb,kb,nqr))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqi),qten(ib,jb,kb,nqi))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqs),qten(ib,jb,kb,nqs))
            call getefall(0,cpi,mf,dum1,dum2,dum4,q3d(ib,jb,kb,nqg),qten(ib,jb,kb,nqg))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
            enddo
            enddo
            enddo
            if( dotbud .and. td_efall.ge.1 )then
              !$omp parallel do default(shared)  &
              !$omp private(i,j,k)
              do k=1,nk-1
              do j=1,nj
              do i=1,ni
                tdiag(i,j,k,td_efall) = dum4(i,j,k)
              enddo
              enddo
              enddo
            endif
          ENDIF
          ! Get final values for th3d,pp3d,prs:
          ! Note:  dum1 stores temperature, thten stores old temperature:
          IF( eqtset.eq.2 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                prs(i,j,k)=rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))*dum1(i,j,k)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
          ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
                rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
              endif
            enddo
            enddo
            enddo
          ENDIF

          IF( getvt )THEN
          IF( dowriteout .or. dotdwrite .or. doazimwrite .or. dorestart )THEN
          if( qd_vtc.ge.1 .or. qd_vtr.ge.1 .or. qd_vts.ge.1 .or. qd_vtg.ge.1 .or. qd_vti.ge.1 )then
            if(myid.eq.0) print *,'  Getting fall velocities for Morrison MP scheme '
              if( qd_vtc.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vtc) = qten(i,j,k,nqc)
                enddo
                enddo
                enddo
              endif
              if( qd_vtr.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vtr) = qten(i,j,k,nqr)
                enddo
                enddo
                enddo
              endif
              if( qd_vts.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vts) = qten(i,j,k,nqs)
                enddo
                enddo
                enddo
              endif
              if( qd_vtg.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vtg) = qten(i,j,k,nqg)
                enddo
                enddo
                enddo
              endif
              if( qd_vti.ge.1 )then
                !$omp parallel do default(shared)  &
                !$omp private(i,j,k)
                do k=1,nk
                do j=1,nj
                do i=1,ni
                  qdiag(i,j,k,qd_vti) = qten(i,j,k,nqi)
                enddo
                enddo
                enddo
              endif
          endif
          ENDIF
          ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   RE87-type scheme   cccccccccccccccccccccccccccccccccccccccccccccc
!ccc   also:  reversible moist thermo. if v_t = 0   cccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.6)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(q3d(i,j,k,2).gt.0.001)then
              qten(i,j,k,2) = v_t
            else
              qten(i,j,k,2) = 0.0
            endif
          enddo
          enddo
          enddo
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          if( efall.ge.1 .and. v_t.gt.1.0e-6 )then
            call getcvm(dum2,q3d)
            call getefall(1,cpl,mf,dum1,dum2,dum4,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk-1
            do j=1,nj
            do i=1,ni
              if( abs(dt*dum4(i,j,k)).ge.tsmall )then
                dum1(i,j,k) = dum1(i,j,k) + dt*dum4(i,j,k)
                prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
                pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
              endif
            enddo
            enddo
            enddo
            if( dotbud .and. td_efall.ge.1 )then
              !$omp parallel do default(shared)  &
              !$omp private(i,j,k)
              do k=1,nk-1
              do j=1,nj
              do i=1,ni
                tdiag(i,j,k,td_efall) = dum4(i,j,k)
              enddo
              enddo
              enddo
            endif
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,prate,dum3,rho,   &
                       q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d)

          IF( v_t.lt.(-0.0001) )THEN
            ! pseudoadiabatic approach of Bryan and Rotunno (2009, JAS, pg 3046)
            !$omp parallel do default(shared)  &
            !$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              q3d(i,j,k,nqc) = min( q3d(i,j,k,nqc) , 0.0001 )
            enddo
            enddo
            enddo
          ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Ziegler/Mansell (NSSL) two-moment scheme
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF( ptype.ge.26 .and. ptype.le.28 )THEN

          DO i = 1,numq
            call pdefq(0.0,asq(i),ruh,rvh,rmh,rho,q3d(ib,jb,kb,i))
          ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k) = pi0(i,j,k)+pp3d(i,j,k)
            dum2(i,j,k) = dz*rmh(i,j,k)
            dum4(i,j,k) = th0(i,j,k)+th3d(i,j,k)
            ! store old theta in thten array:
            thten(i,j,k)=dum4(i,j,k)
          enddo
          enddo
          enddo

          if( radopt.eq.2 )then
            has_reqc = 1
            has_reqi = 1
            has_reqs = 1
          else
            has_reqc = 0
            has_reqi = 0
            has_reqs = 0
          endif

          IF ( ptype .eq. 26 ) THEN  ! graupel only
             call nssl_2mom_driver(                          &
                               th  = dum4,                   &
                               qv  = q3d(ib,jb,kb, lv),      &
                               qc  = q3d(ib,jb,kb, lc),      &
                               qr  = q3d(ib,jb,kb, lr),      &
                               qi  = q3d(ib,jb,kb, li),      &
                               qs  = q3d(ib,jb,kb, ls),      &
                               qh  = q3d(ib,jb,kb, lgh),     &
                               cn  = q3d(ib,jb,kb, lccn),    &
                               cna = q3d(ib,jb,kb, lccna),   &
                               cni = q3d(ib,jb,kb, lcina),   &
                               ccw = q3d(ib,jb,kb, lnc),     &
                               crw = q3d(ib,jb,kb, lnr),     &
                               cci = q3d(ib,jb,kb, lni),     &
                               csw = q3d(ib,jb,kb, lns),     &
                               chw = q3d(ib,jb,kb, lnh),     &
                               zrw = q3d(ib,jb,kb, max(1,lzr)), &
                               zhw = q3d(ib,jb,kb, max(1,lzh)), &
                               vhw = q3d(ib,jb,kb, lvh),     &
                               pii = dum1,                   &
                               p   =  prs,                   &
                               w   =  w3d,                   &
                               dn  =  rho,                   &
                               dz  =  dum2,                  &
                               dtp = dt,                     &
                               itimestep = nstep,            &
                              RAIN = rain,                   &
                              nrain = nrain,                 &
                              prate = prate,                 &
                              dbz = qdiag(ibdq,jbdq,kbdq,qd_dbz), &
                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
                              tcond = qbudget(1),            &
                              tevac = qbudget(2),            &
                              tevar = qbudget(5),            &
                              train = qbudget(6),            &
                              rr    = dum3,                  &
                  re_cloud=effc, re_ice=effi, re_snow=effs,                  &
                  has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,   &
                              dorad = dorad,                 &
                              diagflag = getdbz,             &
                  ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d, &
                  nout3d=nout3d,out3d=out3d,                                 &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
         ELSEIF ( ptype .eq. 27 ) THEN
             call nssl_2mom_driver(                          &
                               th  = dum4,                   &
                               qv  = q3d(ib,jb,kb, lv),      &
                               qc  = q3d(ib,jb,kb, lc),      &
                               qr  = q3d(ib,jb,kb, lr),      &
                               qi  = q3d(ib,jb,kb, li),      &
                               qs  = q3d(ib,jb,kb, ls),      &
                               qh  = q3d(ib,jb,kb, lgh),     &
                               qhl = q3d(ib,jb,kb, lhl),     &
                               cn  = q3d(ib,jb,kb, lccn),    &
                               cna = q3d(ib,jb,kb, lccna),   &
                               cni = q3d(ib,jb,kb, lcina),   &
                               ccw = q3d(ib,jb,kb, lnc),     &
                               crw = q3d(ib,jb,kb, lnr),     &
                               cci = q3d(ib,jb,kb, lni),     &
                               csw = q3d(ib,jb,kb, lns),     &
                               chw = q3d(ib,jb,kb, lnh),     &
                               chl = q3d(ib,jb,kb, lnhl),    &
                               zrw = q3d(ib,jb,kb, max(1,lzr)),     &
                               zhw = q3d(ib,jb,kb, max(1,lzh)),     &
                               zhl = q3d(ib,jb,kb, max(1,lzhl)),    &
                               vhw = q3d(ib,jb,kb, lvh),     &
                               vhl = q3d(ib,jb,kb, lvhl),    &
                               pii = dum1,                   &
                               p   =  prs,                   &
                               w   =  w3d,                   &
                               dn  =  rho,                   &
                               dz  =  dum2,                  &
                               dtp = dt,                     &
                               itimestep = nstep,            &
                              RAIN = rain,                   &
                              nrain = nrain,                 &
                              prate = prate,                 &
                              dbz = qdiag(ibdq,jbdq,kbdq,qd_dbz), &
                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
                              tcond = qbudget(1),            &
                              tevac = qbudget(2),            &
                              tevar = qbudget(5),            &
                              train = qbudget(6),            &
                              rr    = dum3,                  &
                  re_cloud=effc, re_ice=effi, re_snow=effs,                  &
                  has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,   &
                              dorad = dorad,                 &
                              diagflag = getdbz,             &
                  ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d, &
                  nout3d=nout3d,out3d=out3d,                                 &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)
          ELSEIF ( ptype .eq. 28 ) THEN  ! single moment
             call nssl_2mom_driver(                          &
                               th  = dum4,                   &
                               qv  = q3d(ib,jb,kb, 1),       &
                               qc  = q3d(ib,jb,kb, 2),       &
                               qr  = q3d(ib,jb,kb, 3),       &
                               qi  = q3d(ib,jb,kb, 4),       &
                               qs  = q3d(ib,jb,kb, 5),       &
                               qh  = q3d(ib,jb,kb, 6),       &
                               pii = dum1,                   &
                               p   =  prs,                   &
                               w   =  w3d,                   &
                               dn  =  rho,                   &
                               dz  =  dum2,                  &
                               dtp = dt,                     &
                               itimestep = nstep,            &
                              RAIN = rain,                   &
                              nrain = nrain,                 &
                              prate = prate,                 &
                              dbz = qdiag(ibdq,jbdq,kbdq,qd_dbz), &
                              ruh = ruh, rvh = rvh, rmh = rmh, &
                              dx = dx, dy = dy,              &
                              tcond = qbudget(1),            &
                              tevac = qbudget(2),            &
                              tevar = qbudget(5),            &
                              train = qbudget(6),            &
                              rr    = dum3,                  &
                              dorad = .false.,                 &
                              diagflag = getdbz,                  &
                  ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d, &
                  nout3d=nout3d,out3d=out3d,                                 &
                              ims = ib ,ime = ie , jms = jb ,jme = je, kms = kb,kme = ke,  &
                              its = 1 ,ite = ni, jts = 1,jte = nj, kts = 1,kte = nk)

          ENDIF

        IF(eqtset.eq.2)THEN
          ! for mass conservation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
              prs(i,j,k) = rho(i,j,k)*rd*dum4(i,j,k)*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
              pp3d(i,j,k) = (prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
              th3d(i,j,k) = dum4(i,j,k)*dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
            endif
          enddo
          enddo
          enddo
        ELSE
          ! traditional thermodynamics:  p,pi remain unchanged
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
              th3d(i,j,k)= dum4(i,j,k) - th0(i,j,k)
              rho(i,j,k)=prs(i,j,k)/(rd*dum4(i,j,k)*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
            endif
          enddo
          enddo
          enddo
        ENDIF

          if(timestats.ge.1) time_microphy=time_microphy+mytime()


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc     P3 scheme      cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-----------------------------------------------------------------------------------
!  ptype  configuration
!  -----  ----------------------------------------------
!   50    1-mom-cld,  2-mom-rain,  2-mom-ice,  1-ice-category
!   51    2-mom-cld,  2-mom-rain,  3-mom-ice,  1-ice-category,  liq-frac
!   52    2-mom-cld,  2-mom-rain,  3-mom-ice,  2-ice-category,  liq-frac
!   53    2-mom-cld,  2-mom-rain,  3-mom-ice,  3-ice-category,  liq-frac
!   54    2-mom-cld,  2-mom-rain,  3-mom-ice,  4-ice-category,  liq-frac
!-----------------------------------------------------------------------------------

     ELSEIF( ptype.eq.50 .or. ptype.eq.51 .or. ptype.eq.52 .or. ptype.eq.53 .or. ptype.eq.54)THEN

        rdt = 1./dt
        p3a(:,:,:,:) = 0.

        p3_jloop: do j=1,nj

           p3_config: if (ptype.eq.50) then
!---------------------------------------------------------------------------!
!   P3 config:  1-mom-cld, 2-mom-rain, 2-mom-ice, 1-ice-category            !
!---------------------------------------------------------------------------!

              do k=1,nk
              do i=1,ni
                 p3a(i,k,1, 1) = th0(i,j,k)+th3d(i,j,k)
                 p3a(i,k,1, 2) = q3d(i,j,k,nqv)
                 thten(i,j,k)    = p3a(i,k,1,1)
                 qten(i,j,k,nqv) = p3a(i,k,1,2)          ! qv (new)
                 p3a(i,k,1, 3) = q3d(i,j,k,nqc)          ! qc
                 p3a(i,k,1, 4) = q3d(i,j,k,nqr)          ! qr
                 p3a(i,k,1, 5) = q3d(i,j,k,nnr)          ! nr
                 p3a(i,k,1, 6) = th0(i,j,k)+tha(i,j,k)   ! th (old)
                 p3a(i,k,1, 7) = qa(i,j,k,nqv)           ! qv (old)
                 p3a(i,k,1, 8) = pi0(i,j,k)+pp3d(i,j,k)
                 p3a(i,k,1, 9) = prs(i,j,k)
                 p3a(i,k,1,10) = dz*rmh(i,j,k)
                 p3a(i,k,1,11) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
!p3a( ..12) for rainnc
!p3a( ..13) for snownc
!14,15,16 for diag_zdbz, _effc, _effi_ave

                 p3a(i,k,1,17) = q3d(i,j,k,nqi)          ! qit1
                 p3a(i,k,1,18) = q3d(i,j,k,nni1)         ! nit1
                 p3a(i,k,1,19) = q3d(i,j,k,nqir1)        ! qir1
                 p3a(i,k,1,20) = q3d(i,j,k,nbi1)         ! qib1
!21,22,23,24 for four diag_ (cat 1)
!25,... for out2d and out3d diagnostic arrays
              enddo
              enddo

              call mp_p3_wrapper_wrf( th     = p3a(:,:,1,1),                   &
                               qv            = p3a(:,:,1,2),                   &
                               qc            = p3a(:,:,1,3),                   &
                               qr            = p3a(:,:,1,4),                   &
                               qnr           = p3a(:,:,1,5),                   &
                               th_old        = p3a(:,:,1,6),                   &
                               qv_old        = p3a(:,:,1,7),                   &
                               pii           = p3a(:,:,1,8),                   &
                               p             = p3a(:,:,1,9),                   &
                               dz            = p3a(:,:,1,10),                  &
                               w             = p3a(:,:,1,11),                  &
                               dt            = dt,                             &
                               itimestep     = nstep,                          &
                               rainnc        = p3a(:,1,1,12),                  &
                               snownc        = p3a(:,1,1,13),                  &
                              !snownc        = dum1(:,:,4),                    &
                               rainncv       = dum1(:,:,2),                    &
                               sr            = dum1(:,:,3),                    &
                               snowncv       = dum1(:,:,5),                    &
                               ids=1, ide=ni, jds=1, jde=1, kds=1, kde=nk ,    &
                               ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nk ,    &
                               its=1, ite=ni, jts=1, jte=1, kts=1, kte=nk ,    &
                               diag_zdbz     = p3a(:,:,1,14),                  &
                               diag_effc     = p3a(:,:,1,15),                  &
                               diag_effi_ave = p3a(:,:,1,16),                  &
                  n_iceCat      = 1,                              &

                               qit_1         = p3a(:,:,1,17),                  &
                               qni_1         = p3a(:,:,1,18),                  &
                               qir_1         = p3a(:,:,1,19),                  &
                               qib_1         = p3a(:,:,1,20),                  &
                               diag_vmi_1    = p3a(:,:,1,21),                  &
                               diag_dmi_1    = p3a(:,:,1,22),                  &
                               diag_rhoi_1   = p3a(:,:,1,23),                  &
                               diag_dhmax_1  = p3a(:,:,1,24),                  &

                               diag2d_01     = p3a(:,1,1,25),                  &
                               diag2d_02     = p3a(:,1,1,26),                  &
                               diag3d_01     = p3a(:,:,1,27),                  &
                               diag3d_02     = p3a(:,:,1,28),                  &
                               diag3d_03     = p3a(:,:,1,29) )

              do k=1,nk
              do i=1,ni
                 if ( abs(p3a(i,k,1,1)-thten(i,j,k)).ge.tsmall .or.            &
                      abs(p3a(i,k,1,2)-qten(i,j,k,nqv)).ge.qsmall ) then
                    pp3d(i,j,k) = (rho(i,j,k)*(rd+rv*p3a(i,k,1,2))*p3a(i,k,1,1)*rp00)**rddcv
                    prs(i,j,k)  = p00*(pp3d(i,j,k)**cpdrd)
                    pp3d(i,j,k) = pp3d(i,j,k)-pi0(i,j,k)
                    th3d(i,j,k) = th3d(i,j,k)+(p3a(i,k,1,1)-thten(i,j,k))
                 endif
                 q3d(i,j,k,nqv)   = p3a(i,k,1, 2)    ! qv
                 q3d(i,j,k,nqc)   = p3a(i,k,1, 3)    ! qc
                 q3d(i,j,k,nqr)   = p3a(i,k,1, 4)    ! qr
                 q3d(i,j,k,nnr)   = p3a(i,k,1, 5)    ! nr
                 q3d(i,j,k,nqi)   = p3a(i,k,1,17)    ! qit_1
                 q3d(i,j,k,nni1)  = p3a(i,k,1,18)    ! nit_1
                 q3d(i,j,k,nqir1) = p3a(i,k,1,19)    ! qir_1
                 q3d(i,j,k,nbi1)  = p3a(i,k,1,20)    ! qib_1
              enddo
              enddo

              if ( dowriteout .or. getdbz ) then
                 do k=1,nk
                 do i=1,ni
                    p3o(i,j,k, 1) = p3a(i,k,1,21)   ! vmi_1
                    p3o(i,j,k, 2) = p3a(i,k,1,22)   ! dmi_1
                    p3o(i,j,k, 3) = p3a(i,k,1,23)   ! rhoi_1
                    p3o(i,j,k, 4) = p3a(i,k,1,24)   ! dhmax_1
                 enddo
                 enddo
              endif

          ! user-specified diagnostic output
              do i=1,ni
                 out2d(i,j,1) = p3a(i,1,1,25)
                 out2d(i,j,2) = p3a(i,1,1,26)
              enddo
              do k=1,nk
              do i=1,ni
                 out3d(i,j,k,1) = p3a(i,k,1,27)
                 out3d(i,j,k,2) = p3a(i,k,1,28)
                 out3d(i,j,k,3) = p3a(i,k,1,29)
              enddo
              enddo


           elseif (ptype.eq.51) then
!---------------------------------------------------------------------------!
!   P3 config:  2-mom-cld, 2-mom-rain, 3-mom-ice, 1-ice-category, liq-frac  !
!---------------------------------------------------------------------------!

              do k=1,nk
              do i=1,ni
                 p3a(i,k,1, 1) = th0(i,j,k)+th3d(i,j,k)
                 p3a(i,k,1, 2) = q3d(i,j,k,nqv)
                 thten(i,j,k)    = p3a(i,k,1,1)
                 qten(i,j,k,nqv) = p3a(i,k,1,2)          ! qv (new)
                 p3a(i,k,1, 3) = q3d(i,j,k,nqc)          ! qc
                 p3a(i,k,1, 4) = q3d(i,j,k,nqr)          ! qr
                 p3a(i,k,1, 5) = q3d(i,j,k,nnr)          ! nr
                 p3a(i,k,1, 6) = th0(i,j,k)+tha(i,j,k)   ! th (old)
                 p3a(i,k,1, 7) = qa(i,j,k,nqv)           ! qv (old)
                 p3a(i,k,1, 8) = pi0(i,j,k)+pp3d(i,j,k)
                 p3a(i,k,1, 9) = prs(i,j,k)
                 p3a(i,k,1,10) = dz*rmh(i,j,k)
                 p3a(i,k,1,11) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
!p3a( ..12) for rainnc
!p3a( ..13) for snownc
!14,15,16 for diag_zdbz, _effc, _effi_ave

                 p3a(i,k,1,17) = q3d(i,j,k,nqi)          ! qit1
                 p3a(i,k,1,18) = q3d(i,j,k,nni1)         ! nit1
                 p3a(i,k,1,19) = q3d(i,j,k,nqir1)        ! qir1
                 p3a(i,k,1,20) = q3d(i,j,k,nbi1)         ! qib1
!21,22,23,24 for four diag_ (cat 1)

                 p3a(i,k,1,25) = q3d(i,j,k,nnc)          ! nc
                 p3a(i,k,1,26) = q3d(i,j,k,nzi1)         ! zit1
                 p3a(i,k,1,27) = q3d(i,j,k,nqil1)        ! qil1
!28,... for out2d and out3d diagnostic arrays
              enddo
              enddo

              call mp_p3_wrapper_wrf( th     = p3a(:,:,1,1),                   &
                               qv            = p3a(:,:,1,2),                   &
                               qc            = p3a(:,:,1,3),                   &
                               qr            = p3a(:,:,1,4),                   &
                               qnr           = p3a(:,:,1,5),                   &
                               th_old        = p3a(:,:,1,6),                   &
                               qv_old        = p3a(:,:,1,7),                   &
                               pii           = p3a(:,:,1,8),                   &
                               p             = p3a(:,:,1,9),                   &
                               dz            = p3a(:,:,1,10),                  &
                               w             = p3a(:,:,1,11),                  &
                               dt            = dt,                             &
                               itimestep     = nstep,                          &
                               rainnc        = p3a(:,1,1,12),                  &
                               snownc        = p3a(:,1,1,13),                  &
                              !snownc        = dum1(:,:,4),                    &
                               rainncv       = dum1(:,:,2),                    &
                               sr            = dum1(:,:,3),                    &
                               snowncv       = dum1(:,:,5),                    &
                               ids=1, ide=ni, jds=1, jde=1, kds=1, kde=nk ,    &
                               ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nk ,    &
                               its=1, ite=ni, jts=1, jte=1, kts=1, kte=nk ,    &
                               diag_zdbz     = p3a(:,:,1,14),                  &
                               diag_effc     = p3a(:,:,1,15),                  &
                               diag_effi_ave = p3a(:,:,1,16),                  &
                  n_iceCat      = 1,                              &

                               qit_1         = p3a(:,:,1,17),                  &
                               qni_1         = p3a(:,:,1,18),                  &
                               qir_1         = p3a(:,:,1,19),                  &
                               qib_1         = p3a(:,:,1,20),                  &
                               diag_vmi_1    = p3a(:,:,1,21),                  &
                               diag_dmi_1    = p3a(:,:,1,22),                  &
                               diag_rhoi_1   = p3a(:,:,1,23),                  &
                               diag_dhmax_1  = p3a(:,:,1,24),                  &

                               nc            = p3a(:,:,1,25),                  &
                               qzi_1         = p3a(:,:,1,26),                  &
                               qli_1         = p3a(:,:,1,27),                  &

                               diag2d_01     = p3a(:,1,1,28),                  &
                               diag2d_02     = p3a(:,1,1,29),                  &
                               diag3d_01     = p3a(:,:,1,30),                  &
                               diag3d_02     = p3a(:,:,1,31),                  &
                               diag3d_03     = p3a(:,:,1,32) )

              do k=1,nk
              do i=1,ni
                 if ( abs(p3a(i,k,1,1)-thten(i,j,k)).ge.tsmall .or.            &
                      abs(p3a(i,k,1,2)-qten(i,j,k,nqv)).ge.qsmall ) then
                    pp3d(i,j,k) = (rho(i,j,k)*(rd+rv*p3a(i,k,1,2))*p3a(i,k,1,1)*rp00)**rddcv
                    prs(i,j,k)  = p00*(pp3d(i,j,k)**cpdrd)
                    pp3d(i,j,k) = pp3d(i,j,k)-pi0(i,j,k)
                    th3d(i,j,k) = th3d(i,j,k)+(p3a(i,k,1,1)-thten(i,j,k))
                 endif
                 q3d(i,j,k,nqv)   = p3a(i,k,1, 2)    ! qv
                 q3d(i,j,k,nqc)   = p3a(i,k,1, 3)    ! qc
                 q3d(i,j,k,nqr)   = p3a(i,k,1, 4)    ! qr
                 q3d(i,j,k,nnr)   = p3a(i,k,1, 5)    ! nr
                 q3d(i,j,k,nqi)   = p3a(i,k,1,17)    ! qit_1
                 q3d(i,j,k,nni1)  = p3a(i,k,1,18)    ! nit_1
                 q3d(i,j,k,nqir1) = p3a(i,k,1,19)    ! qir_1
                 q3d(i,j,k,nbi1)  = p3a(i,k,1,20)    ! qib_1

                 q3d(i,j,k,nnc)   = p3a(i,k,1,25)    ! nc
                 q3d(i,j,k,nzi1)  = p3a(i,k,1,26)    ! zit_1
                 q3d(i,j,k,nqil1) = p3a(i,k,1,27)    ! qil_1
              enddo
              enddo

              if ( dowriteout .or. getdbz ) then
                 do k=1,nk
                 do i=1,ni
                    p3o(i,j,k, 1) = p3a(i,k,1,21)   ! vmi_1
                    p3o(i,j,k, 2) = p3a(i,k,1,22)   ! dmi_1
                    p3o(i,j,k, 3) = p3a(i,k,1,23)   ! rhoi_1
                    p3o(i,j,k, 4) = p3a(i,k,1,24)   ! dhmax_1
                 enddo
                 enddo
              endif

          ! user-specified diagnostic output
              do i=1,ni
                 out2d(i,j,1) = p3a(i,1,1,28)
                 out2d(i,j,2) = p3a(i,1,1,29)
              enddo
              do k=1,nk
              do i=1,ni
                 out3d(i,j,k,1) = p3a(i,k,1,30)
                 out3d(i,j,k,2) = p3a(i,k,1,31)
                 out3d(i,j,k,3) = p3a(i,k,1,32)
              enddo
              enddo

           elseif (ptype.eq.52) then
!---------------------------------------------------------------------------!
!   P3 config:  2-mom-cld, 2-mom-rain, 3-mom-ice, 2-ice-category, liq-frac  !
!---------------------------------------------------------------------------!

              do k=1,nk
              do i=1,ni
                 p3a(i,k,1, 1) = th0(i,j,k)+th3d(i,j,k)
                 p3a(i,k,1, 2) = q3d(i,j,k,nqv)
                 thten(i,j,k)    = p3a(i,k,1,1)
                 qten(i,j,k,nqv) = p3a(i,k,1,2)          ! qv (new)
                 p3a(i,k,1, 3) = q3d(i,j,k,nqc)          ! qc
                 p3a(i,k,1, 4) = q3d(i,j,k,nqr)          ! qr
                 p3a(i,k,1, 5) = q3d(i,j,k,nnr)          ! nr
                 p3a(i,k,1, 6) = th0(i,j,k)+tha(i,j,k)   ! th (old)
                 p3a(i,k,1, 7) = qa(i,j,k,nqv)           ! qv (old)
                 p3a(i,k,1, 8) = pi0(i,j,k)+pp3d(i,j,k)
                 p3a(i,k,1, 9) = prs(i,j,k)
                 p3a(i,k,1,10) = dz*rmh(i,j,k)
                 p3a(i,k,1,11) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
!p3a( ..12) for rainnc
!p3a( ..13) for snownc
!14,15,16 for diag_zdbz, _effc, _effi_ave

                 p3a(i,k,1,17) = q3d(i,j,k,nqi)          ! qit1
                 p3a(i,k,1,18) = q3d(i,j,k,nni1)         ! nit1
                 p3a(i,k,1,19) = q3d(i,j,k,nqir1)        ! qir1
                 p3a(i,k,1,20) = q3d(i,j,k,nbi1)         ! qib1
!21,22,23,24 for four diag_ (cat 1)

                 p3a(i,k,1,25) = q3d(i,j,k,nnc)          ! nc
                 p3a(i,k,1,26) = q3d(i,j,k,nzi1)         ! zit1
                 p3a(i,k,1,27) = q3d(i,j,k,nqil1)        ! qil1

                 p3a(i,k,1,28) = q3d(i,j,k,nqi2)         ! qit2
                 p3a(i,k,1,29) = q3d(i,j,k,nni2)         ! nii2
                 p3a(i,k,1,30) = q3d(i,j,k,nqir2)        ! qir2
                 p3a(i,k,1,31) = q3d(i,j,k,nbi2)         ! qib2
                 p3a(i,k,1,32) = q3d(i,j,k,nzi2)         ! zit2
                 p3a(i,k,1,33) = q3d(i,j,k,nqil2)        ! qil2
!34,35,36,37 for four diag_ (cat 2)
!38,... for out2d and out3d diagnostic arrays
              enddo
              enddo

              call mp_p3_wrapper_wrf( th     = p3a(:,:,1,1),                   &
                               qv            = p3a(:,:,1,2),                   &
                               qc            = p3a(:,:,1,3),                   &
                               qr            = p3a(:,:,1,4),                   &
                               qnr           = p3a(:,:,1,5),                   &
                               th_old        = p3a(:,:,1,6),                   &
                               qv_old        = p3a(:,:,1,7),                   &
                               pii           = p3a(:,:,1,8),                   &
                               p             = p3a(:,:,1,9),                   &
                               dz            = p3a(:,:,1,10),                  &
                               w             = p3a(:,:,1,11),                  &
                               dt            = dt,                             &
                               itimestep     = nstep,                          &
                               rainnc        = p3a(:,1,1,12),                  &
                               snownc        = p3a(:,1,1,13),                  &
                              !snownc        = dum1(:,:,4),                    &
                               rainncv       = dum1(:,:,2),                    &
                               sr            = dum1(:,:,3),                    &
                               snowncv       = dum1(:,:,5),                    &
                               ids=1, ide=ni, jds=1, jde=1, kds=1, kde=nk ,    &
                               ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nk ,    &
                               its=1, ite=ni, jts=1, jte=1, kts=1, kte=nk ,    &
                               diag_zdbz     = p3a(:,:,1,14),                  &
                               diag_effc     = p3a(:,:,1,15),                  &
                               diag_effi_ave = p3a(:,:,1,16),                  &
                  n_iceCat      = 2,                              &

                               qit_1         = p3a(:,:,1,17),                  &
                               qni_1         = p3a(:,:,1,18),                  &
                               qir_1         = p3a(:,:,1,19),                  &
                               qib_1         = p3a(:,:,1,20),                  &
                               diag_vmi_1    = p3a(:,:,1,21),                  &
                               diag_dmi_1    = p3a(:,:,1,22),                  &
                               diag_rhoi_1   = p3a(:,:,1,23),                  &
                               diag_dhmax_1  = p3a(:,:,1,24),                  &

                               nc            = p3a(:,:,1,25),                  &
                               qzi_1         = p3a(:,:,1,26),                  &
                               qli_1         = p3a(:,:,1,27),                  &

                               qit_2         = p3a(:,:,1,28),                  &
                               qni_2         = p3a(:,:,1,29),                  &
                               qir_2         = p3a(:,:,1,30),                  &
                               qib_2         = p3a(:,:,1,31),                  &
                               qzi_2         = p3a(:,:,1,32),                  &
                               qli_2         = p3a(:,:,1,33),                  &
                               diag_vmi_2    = p3a(:,:,1,34),                  &
                               diag_dmi_2    = p3a(:,:,1,35),                  &
                               diag_rhoi_2   = p3a(:,:,1,36),                  &
                               diag_dhmax_2  = p3a(:,:,1,37),                  &

                               diag2d_01     = p3a(:,1,1,38),                  &
                               diag2d_02     = p3a(:,1,1,39),                  &
                               diag3d_01     = p3a(:,:,1,40),                  &
                               diag3d_02     = p3a(:,:,1,41),                  &
                               diag3d_03     = p3a(:,:,1,42) )

              do k=1,nk
              do i=1,ni
                 if ( abs(p3a(i,k,1,1)-thten(i,j,k)).ge.tsmall .or.            &
                      abs(p3a(i,k,1,2)-qten(i,j,k,nqv)).ge.qsmall ) then
                    pp3d(i,j,k) = (rho(i,j,k)*(rd+rv*p3a(i,k,1,2))*p3a(i,k,1,1)*rp00)**rddcv
                    prs(i,j,k)  = p00*(pp3d(i,j,k)**cpdrd)
                    pp3d(i,j,k) = pp3d(i,j,k)-pi0(i,j,k)
                    th3d(i,j,k) = th3d(i,j,k)+(p3a(i,k,1,1)-thten(i,j,k))
                 endif
                 q3d(i,j,k,nqv)   = p3a(i,k,1, 2)    ! qv
                 q3d(i,j,k,nqc)   = p3a(i,k,1, 3)    ! qc
                 q3d(i,j,k,nqr)   = p3a(i,k,1, 4)    ! qr
                 q3d(i,j,k,nnr)   = p3a(i,k,1, 5)    ! nr
                 q3d(i,j,k,nqi)   = p3a(i,k,1,17)    ! qit_1
                 q3d(i,j,k,nni1)  = p3a(i,k,1,18)    ! nit_1
                 q3d(i,j,k,nqir1) = p3a(i,k,1,19)    ! qir_1
                 q3d(i,j,k,nbi1)  = p3a(i,k,1,20)    ! qib_1

                 q3d(i,j,k,nnc)   = p3a(i,k,1,25)    ! nc
                 q3d(i,j,k,nzi1)  = p3a(i,k,1,26)    ! zit_1
                 q3d(i,j,k,nqil1) = p3a(i,k,1,27)    ! qil_1

                 q3d(i,j,k,nqi2)  = p3a(i,k,1,28)    ! qit_2
                 q3d(i,j,k,nni2)  = p3a(i,k,1,29)    ! nit_2
                 q3d(i,j,k,nqir2) = p3a(i,k,1,30)    ! qir_2
                 q3d(i,j,k,nbi2)  = p3a(i,k,1,31)    ! bir_2
                 q3d(i,j,k,nzi2)  = p3a(i,k,1,32)    ! zit_2
                 q3d(i,j,k,nqil2) = p3a(i,k,1,33)    ! qil_2
              enddo
              enddo

              if ( dowriteout .or. getdbz ) then
                 do k=1,nk
                 do i=1,ni
                    p3o(i,j,k, 1) = p3a(i,k,1,21)   ! vmi_1
                    p3o(i,j,k, 2) = p3a(i,k,1,22)   ! dmi_1
                    p3o(i,j,k, 3) = p3a(i,k,1,23)   ! rhoi_1
                    p3o(i,j,k, 4) = p3a(i,k,1,24)   ! dhmax_1

                    p3o(i,j,k, 5) = p3a(i,k,1,34)   ! vmi_2
                    p3o(i,j,k, 6) = p3a(i,k,1,35)   ! dmi_2
                    p3o(i,j,k, 7) = p3a(i,k,1,36)   ! rhoi_2
                    p3o(i,j,k, 8) = p3a(i,k,1,37)   ! dhmax_2
                 enddo
                 enddo
              endif

          ! user-specified diagnostic output
              do i=1,ni
                 out2d(i,j,1) = p3a(i,1,1,38)
                 out2d(i,j,2) = p3a(i,1,1,39)
              enddo
              do k=1,nk
              do i=1,ni
                 out3d(i,j,k,1) = p3a(i,k,1,40)
                 out3d(i,j,k,2) = p3a(i,k,1,41)
                 out3d(i,j,k,3) = p3a(i,k,1,42)
              enddo
              enddo

           elseif ( ptype.eq.53 ) then
!---------------------------------------------------------------------------!
!   P3 config:  2-mom-cld, 2-mom-rain, 3-mom-ice, 3-ice-category, liq-frac  !
!---------------------------------------------------------------------------!

              do k=1,nk
              do i=1,ni
                 p3a(i,k,1, 1) = th0(i,j,k)+th3d(i,j,k)
                 p3a(i,k,1, 2) = q3d(i,j,k,nqv)
                 thten(i,j,k)    = p3a(i,k,1,1)
                 qten(i,j,k,nqv) = p3a(i,k,1,2)          ! qv (new)
                 p3a(i,k,1, 3) = q3d(i,j,k,nqc)          ! qc
                 p3a(i,k,1, 4) = q3d(i,j,k,nqr)          ! qr
                 p3a(i,k,1, 5) = q3d(i,j,k,nnr)          ! nr
                 p3a(i,k,1, 6) = th0(i,j,k)+tha(i,j,k)   ! th (old)
                 p3a(i,k,1, 7) = qa(i,j,k,nqv)           ! qv (old)
                 p3a(i,k,1, 8) = pi0(i,j,k)+pp3d(i,j,k)
                 p3a(i,k,1, 9) = prs(i,j,k)
                 p3a(i,k,1,10) = dz*rmh(i,j,k)
                 p3a(i,k,1,11) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
!p3a( ..12) for rainnc
!p3a( ..13) for snownc
!14,15,16 for diag_zdbz, _effc, _effi_ave

                 p3a(i,k,1,17) = q3d(i,j,k,nqi)          ! qit1
                 p3a(i,k,1,18) = q3d(i,j,k,nni1)         ! nit1
                 p3a(i,k,1,19) = q3d(i,j,k,nqir1)        ! qir1
                 p3a(i,k,1,20) = q3d(i,j,k,nbi1)         ! qib1
!21,22,23,24 for four diag_ (cat 1)

                 p3a(i,k,1,25) = q3d(i,j,k,nnc)          ! nc
                 p3a(i,k,1,26) = q3d(i,j,k,nzi1)         ! zit1
                 p3a(i,k,1,27) = q3d(i,j,k,nqil1)        ! qil1

                 p3a(i,k,1,28) = q3d(i,j,k,nqi2)         ! qit2
                 p3a(i,k,1,29) = q3d(i,j,k,nni2)         ! nii2
                 p3a(i,k,1,30) = q3d(i,j,k,nqir2)        ! qir2
                 p3a(i,k,1,31) = q3d(i,j,k,nbi2)         ! qib2
                 p3a(i,k,1,32) = q3d(i,j,k,nzi2)         ! zit2
                 p3a(i,k,1,33) = q3d(i,j,k,nqil2)        ! qil2
!34,35,36,37 for four diag_ (cat 2)

                 p3a(i,k,1,38) = q3d(i,j,k,nqi3)         ! qit3
                 p3a(i,k,1,39) = q3d(i,j,k,nni3)         ! nii3
                 p3a(i,k,1,40) = q3d(i,j,k,nqir3)        ! qir3
                 p3a(i,k,1,41) = q3d(i,j,k,nbi3)         ! qib3
                 p3a(i,k,1,42) = q3d(i,j,k,nzi3)         ! zit3
                 p3a(i,k,1,43) = q3d(i,j,k,nqil3)        ! qil3
!44,45,46,47 for four diag_ (cat 3)
!48,... for out2d and out3d diagnostic arrays
              enddo
              enddo

              call mp_p3_wrapper_wrf( th     = p3a(:,:,1,1),                   &
                               qv            = p3a(:,:,1,2),                   &
                               qc            = p3a(:,:,1,3),                   &
                               qr            = p3a(:,:,1,4),                   &
                               qnr           = p3a(:,:,1,5),                   &
                               th_old        = p3a(:,:,1,6),                   &
                               qv_old        = p3a(:,:,1,7),                   &
                               pii           = p3a(:,:,1,8),                   &
                               p             = p3a(:,:,1,9),                   &
                               dz            = p3a(:,:,1,10),                  &
                               w             = p3a(:,:,1,11),                  &
                               dt            = dt,                             &
                               itimestep     = nstep,                          &
                               rainnc        = p3a(:,1,1,12),                  &
                               snownc        = p3a(:,1,1,13),                  &
                              !snownc        = dum1(:,:,4),                    &
                               rainncv       = dum1(:,:,2),                    &
                               sr            = dum1(:,:,3),                    &
                               snowncv       = dum1(:,:,5),                    &
                               ids=1, ide=ni, jds=1, jde=1, kds=1, kde=nk ,    &
                               ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nk ,    &
                               its=1, ite=ni, jts=1, jte=1, kts=1, kte=nk ,    &
                               diag_zdbz     = p3a(:,:,1,14),                  &
                               diag_effc     = p3a(:,:,1,15),                  &
                               diag_effi_ave = p3a(:,:,1,16),                  &
                  n_iceCat      = 3,                              &

                               qit_1         = p3a(:,:,1,17),                  &
                               qni_1         = p3a(:,:,1,18),                  &
                               qir_1         = p3a(:,:,1,19),                  &
                               qib_1         = p3a(:,:,1,20),                  &
                               diag_vmi_1    = p3a(:,:,1,21),                  &
                               diag_dmi_1    = p3a(:,:,1,22),                  &
                               diag_rhoi_1   = p3a(:,:,1,23),                  &
                               diag_dhmax_1  = p3a(:,:,1,24),                  &

                               nc            = p3a(:,:,1,25),                  &
                               qzi_1         = p3a(:,:,1,26),                  &
                               qli_1         = p3a(:,:,1,27),                  &

                               qit_2         = p3a(:,:,1,28),                  &
                               qni_2         = p3a(:,:,1,29),                  &
                               qir_2         = p3a(:,:,1,30),                  &
                               qib_2         = p3a(:,:,1,31),                  &
                               qzi_2         = p3a(:,:,1,32),                  &
                               qli_2         = p3a(:,:,1,33),                  &
                               diag_vmi_2    = p3a(:,:,1,34),                  &
                               diag_dmi_2    = p3a(:,:,1,35),                  &
                               diag_rhoi_2   = p3a(:,:,1,36),                  &
                               diag_dhmax_2  = p3a(:,:,1,37),                  &

                               qit_3         = p3a(:,:,1,38),                  &
                               qni_3         = p3a(:,:,1,39),                  &
                               qir_3         = p3a(:,:,1,40),                  &
                               qib_3         = p3a(:,:,1,41),                  &
                               qzi_3         = p3a(:,:,1,42),                  &
                               qli_3         = p3a(:,:,1,43),                  &
                               diag_vmi_3    = p3a(:,:,1,44),                  &
                               diag_dmi_3    = p3a(:,:,1,45),                  &
                               diag_rhoi_3   = p3a(:,:,1,46),                  &
                               diag_dhmax_3  = p3a(:,:,1,47),                  &

                               diag2d_01     = p3a(:,1,1,48),                  &
                               diag2d_02     = p3a(:,1,1,49),                  &
                               diag3d_01     = p3a(:,:,1,50),                  &
                               diag3d_02     = p3a(:,:,1,51),                  &
                               diag3d_03     = p3a(:,:,1,52) )

              do k=1,nk
              do i=1,ni
                 if ( abs(p3a(i,k,1,1)-thten(i,j,k)).ge.tsmall .or.            &
                      abs(p3a(i,k,1,2)-qten(i,j,k,nqv)).ge.qsmall ) then
                    pp3d(i,j,k) = (rho(i,j,k)*(rd+rv*p3a(i,k,1,2))*p3a(i,k,1,1)*rp00)**rddcv
                    prs(i,j,k)  = p00*(pp3d(i,j,k)**cpdrd)
                    pp3d(i,j,k) = pp3d(i,j,k)-pi0(i,j,k)
                    th3d(i,j,k) = th3d(i,j,k)+(p3a(i,k,1,1)-thten(i,j,k))
                 endif
                 q3d(i,j,k,nqv)   = p3a(i,k,1, 2)    ! qv
                 q3d(i,j,k,nqc)   = p3a(i,k,1, 3)    ! qc
                 q3d(i,j,k,nqr)   = p3a(i,k,1, 4)    ! qr
                 q3d(i,j,k,nnr)   = p3a(i,k,1, 5)    ! nr
                 q3d(i,j,k,nqi)   = p3a(i,k,1,17)    ! qit_1
                 q3d(i,j,k,nni1)  = p3a(i,k,1,18)    ! nit_1
                 q3d(i,j,k,nqir1) = p3a(i,k,1,19)    ! qir_1
                 q3d(i,j,k,nbi1)  = p3a(i,k,1,20)    ! qib_1

                 q3d(i,j,k,nnc)   = p3a(i,k,1,25)    ! nc
                 q3d(i,j,k,nzi1)  = p3a(i,k,1,26)    ! zit_1
                 q3d(i,j,k,nqil1) = p3a(i,k,1,27)    ! qil_1

                 q3d(i,j,k,nqi2)  = p3a(i,k,1,28)    ! qit_2
                 q3d(i,j,k,nni2)  = p3a(i,k,1,29)    ! nit_2
                 q3d(i,j,k,nqir2) = p3a(i,k,1,30)    ! qir_2
                 q3d(i,j,k,nbi2)  = p3a(i,k,1,31)    ! bir_2
                 q3d(i,j,k,nzi2)  = p3a(i,k,1,32)    ! zit_2
                 q3d(i,j,k,nqil2) = p3a(i,k,1,33)    ! qil_2

                 q3d(i,j,k,nqi3)  = p3a(i,k,1,38)    ! qit_3
                 q3d(i,j,k,nni3)  = p3a(i,k,1,39)    ! nit_3
                 q3d(i,j,k,nqir3) = p3a(i,k,1,40)    ! qir_3
                 q3d(i,j,k,nbi3)  = p3a(i,k,1,41)    ! bir_3
                 q3d(i,j,k,nzi3)  = p3a(i,k,1,42)    ! zit_3
                 q3d(i,j,k,nqil3) = p3a(i,k,1,43)    ! qil_3
              enddo
              enddo

              if ( dowriteout .or. getdbz ) then
                 do k=1,nk
                 do i=1,ni
                    p3o(i,j,k, 1) = p3a(i,k,1,21)   ! vmi_1
                    p3o(i,j,k, 2) = p3a(i,k,1,22)   ! dmi_1
                    p3o(i,j,k, 3) = p3a(i,k,1,23)   ! rhoi_1
                    p3o(i,j,k, 4) = p3a(i,k,1,24)   ! dhmax_1

                    p3o(i,j,k, 5) = p3a(i,k,1,34)   ! vmi_2
                    p3o(i,j,k, 6) = p3a(i,k,1,35)   ! dmi_2
                    p3o(i,j,k, 7) = p3a(i,k,1,36)   ! rhoi_2
                    p3o(i,j,k, 8) = p3a(i,k,1,37)   ! dhmax_2

                    p3o(i,j,k, 9) = p3a(i,k,1,44)   ! vmi_3
                    p3o(i,j,k,10) = p3a(i,k,1,45)   ! dmi_3
                    p3o(i,j,k,11) = p3a(i,k,1,46)   ! rhoi_3
                    p3o(i,j,k,12) = p3a(i,k,1,47)   ! dhmax_3
                 enddo
                 enddo
              endif

          ! user-specified diagnostic output
              do i=1,ni
                 out2d(i,j,1) = p3a(i,1,1,48)
                 out2d(i,j,2) = p3a(i,1,1,49)
              enddo
              do k=1,nk
              do i=1,ni
                 out3d(i,j,k,1) = p3a(i,k,1,50)
                 out3d(i,j,k,2) = p3a(i,k,1,51)
                 out3d(i,j,k,3) = p3a(i,k,1,52)
              enddo
              enddo

           elseif (ptype.eq.54) then
!---------------------------------------------------------------------------!
!   P3 config:  2-mom-cld, 2-mom-rain, 3-mom-ice, 4-ice-category, liq-frac  !
!---------------------------------------------------------------------------!

              do k=1,nk
              do i=1,ni
                 p3a(i,k,1, 1) = th0(i,j,k)+th3d(i,j,k)
                 p3a(i,k,1, 2) = q3d(i,j,k,nqv)
                 thten(i,j,k)    = p3a(i,k,1,1)
                 qten(i,j,k,nqv) = p3a(i,k,1,2)          ! qv (new)
                 p3a(i,k,1, 3) = q3d(i,j,k,nqc)          ! qc
                 p3a(i,k,1, 4) = q3d(i,j,k,nqr)          ! qr
                 p3a(i,k,1, 5) = q3d(i,j,k,nnr)          ! nr
                 p3a(i,k,1, 6) = th0(i,j,k)+tha(i,j,k)   ! th (old)
                 p3a(i,k,1, 7) = qa(i,j,k,nqv)           ! qv (old)
                 p3a(i,k,1, 8) = pi0(i,j,k)+pp3d(i,j,k)
                 p3a(i,k,1, 9) = prs(i,j,k)
                 p3a(i,k,1,10) = dz*rmh(i,j,k)
                 p3a(i,k,1,11) = 0.5*(w3d(i,j,k)+w3d(i,j,k+1))
!p3a( ..12) for rainnc
!p3a( ..13) for snownc
!14,15,16 for diag_zdbz, _effc, _effi_ave

                 p3a(i,k,1,17) = q3d(i,j,k,nqi)          ! qit1
                 p3a(i,k,1,18) = q3d(i,j,k,nni1)         ! nit1
                 p3a(i,k,1,19) = q3d(i,j,k,nqir1)        ! qir1
                 p3a(i,k,1,20) = q3d(i,j,k,nbi1)         ! qib1
!21,22,23,24 for four diag_ (cat 1)

                 p3a(i,k,1,25) = q3d(i,j,k,nnc)          ! nc
                 p3a(i,k,1,26) = q3d(i,j,k,nzi1)         ! zit1
                 p3a(i,k,1,27) = q3d(i,j,k,nqil1)        ! qil1

                 p3a(i,k,1,28) = q3d(i,j,k,nqi2)         ! qit2
                 p3a(i,k,1,29) = q3d(i,j,k,nni2)         ! nii2
                 p3a(i,k,1,30) = q3d(i,j,k,nqir2)        ! qir2
                 p3a(i,k,1,31) = q3d(i,j,k,nbi2)         ! qib2
                 p3a(i,k,1,32) = q3d(i,j,k,nzi2)         ! zit2
                 p3a(i,k,1,33) = q3d(i,j,k,nqil2)        ! qil2
!34,35,36,37 for four diag_ (cat 2)

                 p3a(i,k,1,38) = q3d(i,j,k,nqi3)         ! qit3
                 p3a(i,k,1,39) = q3d(i,j,k,nni3)         ! nii3
                 p3a(i,k,1,40) = q3d(i,j,k,nqir3)        ! qir3
                 p3a(i,k,1,41) = q3d(i,j,k,nbi3)         ! qib3
                 p3a(i,k,1,42) = q3d(i,j,k,nzi3)         ! zit3
                 p3a(i,k,1,43) = q3d(i,j,k,nqil3)        ! qil3
!44,45,46,47 for four diag_ (cat 3)

                 p3a(i,k,1,48) = q3d(i,j,k,nqi4)         ! qit4
                 p3a(i,k,1,49) = q3d(i,j,k,nni4)         ! nii4
                 p3a(i,k,1,50) = q3d(i,j,k,nqir4)        ! qir4
                 p3a(i,k,1,51) = q3d(i,j,k,nbi4)         ! qib4
                 p3a(i,k,1,52) = q3d(i,j,k,nzi4)         ! zit4
                 p3a(i,k,1,53) = q3d(i,j,k,nqil4)        ! qil4
!54,55,56,57 for four diag_ (cat 4)
!58,... for out2d and out3d diagnostic arrays
              enddo
              enddo

              call mp_p3_wrapper_wrf( th     = p3a(:,:,1,1),                   &
                               qv            = p3a(:,:,1,2),                   &
                               qc            = p3a(:,:,1,3),                   &
                               qr            = p3a(:,:,1,4),                   &
                               qnr           = p3a(:,:,1,5),                   &
                               th_old        = p3a(:,:,1,6),                   &
                               qv_old        = p3a(:,:,1,7),                   &
                               pii           = p3a(:,:,1,8),                   &
                               p             = p3a(:,:,1,9),                   &
                               dz            = p3a(:,:,1,10),                  &
                               w             = p3a(:,:,1,11),                  &
                               dt            = dt,                             &
                               itimestep     = nstep,                          &
                               rainnc        = p3a(:,1,1,12),                  &
                               snownc        = p3a(:,1,1,13),                  &
                              !snownc        = dum1(:,:,4),                    &
                               rainncv       = dum1(:,:,2),                    &
                               sr            = dum1(:,:,3),                    &
                               snowncv       = dum1(:,:,5),                    &
                               ids=1, ide=ni, jds=1, jde=1, kds=1, kde=nk ,    &
                               ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nk ,    &
                               its=1, ite=ni, jts=1, jte=1, kts=1, kte=nk ,    &
                               diag_zdbz     = p3a(:,:,1,14),                  &
                               diag_effc     = p3a(:,:,1,15),                  &
                               diag_effi_ave = p3a(:,:,1,16),                  &
                  n_iceCat      = 4,                              &

                               qit_1         = p3a(:,:,1,17),                  &
                               qni_1         = p3a(:,:,1,18),                  &
                               qir_1         = p3a(:,:,1,19),                  &
                               qib_1         = p3a(:,:,1,20),                  &
                               diag_vmi_1    = p3a(:,:,1,21),                  &
                               diag_dmi_1    = p3a(:,:,1,22),                  &
                               diag_rhoi_1   = p3a(:,:,1,23),                  &
                               diag_dhmax_1  = p3a(:,:,1,24),                  &

                               nc            = p3a(:,:,1,25),                  &
                               qzi_1         = p3a(:,:,1,26),                  &
                               qli_1         = p3a(:,:,1,27),                  &

                               qit_2         = p3a(:,:,1,28),                  &
                               qni_2         = p3a(:,:,1,29),                  &
                               qir_2         = p3a(:,:,1,30),                  &
                               qib_2         = p3a(:,:,1,31),                  &
                               qzi_2         = p3a(:,:,1,32),                  &
                               qli_2         = p3a(:,:,1,33),                  &
                               diag_vmi_2    = p3a(:,:,1,34),                  &
                               diag_dmi_2    = p3a(:,:,1,35),                  &
                               diag_rhoi_2   = p3a(:,:,1,36),                  &
                               diag_dhmax_2  = p3a(:,:,1,37),                  &

                               qit_3         = p3a(:,:,1,38),                  &
                               qni_3         = p3a(:,:,1,39),                  &
                               qir_3         = p3a(:,:,1,40),                  &
                               qib_3         = p3a(:,:,1,41),                  &
                               qzi_3         = p3a(:,:,1,42),                  &
                               qli_3         = p3a(:,:,1,43),                  &
                               diag_vmi_3    = p3a(:,:,1,44),                  &
                               diag_dmi_3    = p3a(:,:,1,45),                  &
                               diag_rhoi_3   = p3a(:,:,1,46),                  &
                               diag_dhmax_3  = p3a(:,:,1,47),                  &

                               qit_4         = p3a(:,:,1,48),                  &
                               qni_4         = p3a(:,:,1,49),                  &
                               qir_4         = p3a(:,:,1,50),                  &
                               qib_4         = p3a(:,:,1,51),                  &
                               qzi_4         = p3a(:,:,1,52),                  &
                               qli_4         = p3a(:,:,1,53),                  &
                               diag_vmi_4    = p3a(:,:,1,54),                  &
                               diag_dmi_4    = p3a(:,:,1,55),                  &
                               diag_rhoi_4   = p3a(:,:,1,56),                  &
                               diag_dhmax_4  = p3a(:,:,1,57),                  &

                               diag2d_01     = p3a(:,1,1,58),                  &
                               diag2d_02     = p3a(:,1,1,59),                  &
                               diag3d_01     = p3a(:,:,1,60),                  &
                               diag3d_02     = p3a(:,:,1,61),                  &
                               diag3d_03     = p3a(:,:,1,62) )

              do k=1,nk
              do i=1,ni
                 if ( abs(p3a(i,k,1,1)-thten(i,j,k)).ge.tsmall .or.            &
                      abs(p3a(i,k,1,2)-qten(i,j,k,nqv)).ge.qsmall ) then
                    pp3d(i,j,k) = (rho(i,j,k)*(rd+rv*p3a(i,k,1,2))*p3a(i,k,1,1)*rp00)**rddcv
                    prs(i,j,k)  = p00*(pp3d(i,j,k)**cpdrd)
                    pp3d(i,j,k) = pp3d(i,j,k)-pi0(i,j,k)
                    th3d(i,j,k) = th3d(i,j,k)+(p3a(i,k,1,1)-thten(i,j,k))
                 endif
                 q3d(i,j,k,nqv)   = p3a(i,k,1, 2)    ! qv
                 q3d(i,j,k,nqc)   = p3a(i,k,1, 3)    ! qc
                 q3d(i,j,k,nqr)   = p3a(i,k,1, 4)    ! qr
                 q3d(i,j,k,nnr)   = p3a(i,k,1, 5)    ! nr
                 q3d(i,j,k,nqi)   = p3a(i,k,1,17)    ! qit_1
                 q3d(i,j,k,nni1)  = p3a(i,k,1,18)    ! nit_1
                 q3d(i,j,k,nqir1) = p3a(i,k,1,19)    ! qir_1
                 q3d(i,j,k,nbi1)  = p3a(i,k,1,20)    ! qib_1

                 q3d(i,j,k,nnc)   = p3a(i,k,1,25)    ! nc
                 q3d(i,j,k,nzi1)  = p3a(i,k,1,26)    ! zit_1
                 q3d(i,j,k,nqil1) = p3a(i,k,1,27)    ! qil_1

                 q3d(i,j,k,nqi2)  = p3a(i,k,1,28)    ! qit_2
                 q3d(i,j,k,nni2)  = p3a(i,k,1,29)    ! nit_2
                 q3d(i,j,k,nqir2) = p3a(i,k,1,30)    ! qir_2
                 q3d(i,j,k,nbi2)  = p3a(i,k,1,31)    ! bir_2
                 q3d(i,j,k,nzi2)  = p3a(i,k,1,32)    ! zit_2
                 q3d(i,j,k,nqil2) = p3a(i,k,1,33)    ! qil_2

                 q3d(i,j,k,nqi3)  = p3a(i,k,1,38)    ! qit_3
                 q3d(i,j,k,nni3)  = p3a(i,k,1,39)    ! nit_3
                 q3d(i,j,k,nqir3) = p3a(i,k,1,40)    ! qir_3
                 q3d(i,j,k,nbi3)  = p3a(i,k,1,41)    ! bir_3
                 q3d(i,j,k,nzi3)  = p3a(i,k,1,42)    ! zit_3
                 q3d(i,j,k,nqil3) = p3a(i,k,1,43)    ! qil_3

                 q3d(i,j,k,nqi4)  = p3a(i,k,1,48)    ! qit_4
                 q3d(i,j,k,nni4)  = p3a(i,k,1,49)    ! nit_4
                 q3d(i,j,k,nqir4) = p3a(i,k,1,50)    ! qir_4
                 q3d(i,j,k,nbi4)  = p3a(i,k,1,51)    ! bir_4
                 q3d(i,j,k,nzi4)  = p3a(i,k,1,52)    ! zit_4
                 q3d(i,j,k,nqil4) = p3a(i,k,1,53)    ! qil_4
              enddo
              enddo

              if ( dowriteout .or. getdbz ) then
                 do k=1,nk
                 do i=1,ni
                    p3o(i,j,k, 1) = p3a(i,k,1,21)   ! vmi_1
                    p3o(i,j,k, 2) = p3a(i,k,1,22)   ! dmi_1
                    p3o(i,j,k, 3) = p3a(i,k,1,23)   ! rhoi_1
                    p3o(i,j,k, 4) = p3a(i,k,1,24)   ! dhmax_1

                    p3o(i,j,k, 5) = p3a(i,k,1,34)   ! vmi_2
                    p3o(i,j,k, 6) = p3a(i,k,1,35)   ! dmi_2
                    p3o(i,j,k, 7) = p3a(i,k,1,36)   ! rhoi_2
                    p3o(i,j,k, 8) = p3a(i,k,1,37)   ! dhmax_2

                    p3o(i,j,k, 9) = p3a(i,k,1,44)   ! vmi_3
                    p3o(i,j,k,10) = p3a(i,k,1,45)   ! dmi_3
                    p3o(i,j,k,11) = p3a(i,k,1,46)   ! rhoi_3
                    p3o(i,j,k,12) = p3a(i,k,1,47)   ! dhmax_3

                    p3o(i,j,k,13) = p3a(i,k,1,54)   ! vmi_4
                    p3o(i,j,k,14) = p3a(i,k,1,55)   ! dmi_4
                    p3o(i,j,k,15) = p3a(i,k,1,56)   ! rhoi_4
                    p3o(i,j,k,16) = p3a(i,k,1,57)   ! dhmax_4
                 enddo
                 enddo
              endif

          ! user-specified diagnostic output
              do i=1,ni
                 out2d(i,j,1) = p3a(i,1,1,58)
                 out2d(i,j,2) = p3a(i,1,1,59)
              enddo
              do k=1,nk
              do i=1,ni
                 out3d(i,j,k,1) = p3a(i,k,1,60)
                 out3d(i,j,k,2) = p3a(i,k,1,61)
                 out3d(i,j,k,3) = p3a(i,k,1,62)
              enddo
              enddo

           endif p3_config

         !------   code common to all configs of P3:  ------

           if (axisymm.eq.1) then
              do i=1,ni
                 prate(i,j) = p3a(i,1,1,12)*rdt
                 qbudget(6) = qbudget(6) + p3a(i,1,1,12)*ruh(i)*rvh(j)*dx*dy*dum3(i,j,1)/rho(i,j,1)
              enddo
           else
              do i=1,ni
                 prate(i,j) = p3a(i,1,1,12)*rdt
                 qbudget(6) = qbudget(6) + p3a(i,1,1,12)*ruh(i)*rvh(j)*dx*dy
              enddo
           endif

           do n=1,nrain
           do i=1,ni
             ! convert from mm to cm:
              rain(i,j,n) = rain(i,j,n) + 0.1*p3a(i,1,1,12)  !actually, total precip
             !snow(i,j,n) = snow(i,j,n) + 0.1*p3a(i,1,1,13)
           enddo
           enddo

           if ( getdbz .and. qd_dbz.ge.1 ) then
              do k=1,nk
              do i=1,ni
                 qdiag(i,j,k,qd_dbz) = p3a(i,k,1,14)
              enddo
              enddo
           endif

           if (dorad) then
              do k=1,nk
              do i=1,ni
                 effc(i,j,k) = p3a(i,k,1,15)
                 effi(i,j,k) = p3a(i,k,1,16)
              enddo
              enddo
           endif


        enddo p3_jloop

        if(timestats.ge.1) time_microphy=time_microphy+mytime()

!                  --- END OF CODE FOR P3 SCHEME ---

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Jensen scheme   ccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ELSEIF(ptype.eq.55)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))

          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! cm1r17:  to make things easier to understand, use same arrays
            !          that are used for morrison code:
            ! dum1 = T  (this should have been calculated already)
            ! dum2 = pi (nondimensional pressure)
            ! dum4 = dz
            ! thten = copy of T  (this should have been calculated already)
            dum2(i,j,k)=pi0(i,j,k)+pp3d(i,j,k)
            dum4(i,j,k)=dz*rmh(i,j,k)
          enddo
          enddo
          enddo

!          write(6,*) 'got to call'
!          write(6,*) nstep, dt, prs(1,1,1), dum4(1,1,1), dum1(1,1,1)
!          write(6,*) q3d(ib,jb,kb,1),q3d(ib,jb,kb,2), q3d(ib,jb,kb,3), q3d(ib,jb,kb,4), &
!               q3d(ib,jb,kb,5), q3d(ib,jb,kb,6), q3d(ib,jb,kb,7), q3d(ib,jb,kb,8), &
!               q3d(ib,jb,kb,9), q3d(ib,jb,kb,10), q3d(ib,jb,kb,11), q3d(ib,jb,kb,12), &
!               q3d(ib,jb,kb,13), q3d(ib,jb,kb,14), q3d(ib,jb,kb,15), q3d(ib,jb,kb,16)

!          write(6,*) ni, nk, nj, ib, ie, jb, je, kb, ke, rain(ib,jb,1),dum5(ib,jb,1)
!          write(6,*) nrain, dx, dy, dz, qbudget(1), qbudget(2), qbudget(5), qbudget(6)
!          write(6,*)  ruh(1), rvh(1), rmh(1,1,1), dum3(1,1,1), rain(1,1,1), prate(1,1)

!               ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d, &
!                    nout3d=nout3d,out3d=out3d,eqtset=eqtset)

          call   mp_jensen_ishmael( &
          itimestep=nstep, dt_in=dt, p=prs, dz=dum4, t3d=dum1, &
               qv=q3d(ib,jb,kb,1), &
               qc=q3d(ib,jb,kb,2), &
               qr=q3d(ib,jb,kb,3), &
               nr=q3d(ib,jb,kb,7), &
               qi1=q3d(ib,jb,kb,4), &
               ni1=q3d(ib,jb,kb,8), &
               ai1=q3d(ib,jb,kb,11), &
               ci1=q3d(ib,jb,kb,14), &
               qi2=q3d(ib,jb,kb,5), &
               ni2=q3d(ib,jb,kb,9), &
               ai2=q3d(ib,jb,kb,12), &
               ci2=q3d(ib,jb,kb,15), &
               qi3=q3d(ib,jb,kb,6), &
               ni3=q3d(ib,jb,kb,10), &
               ai3=q3d(ib,jb,kb,13), &
               ci3=q3d(ib,jb,kb,16), &
               ids=1  ,ide=ni+1 , jds= 1 ,jde=nj+1 , kds=1  ,kde=nk+1 ,   &
               ims=ib ,ime=ie   , jms=jb ,jme=je   , kms=kb ,kme=ke ,     &
               its=1  ,ite=ni   , jts=1  ,jte=nj   , kts=1  ,kte=nk ,     &
               RAINNC=rain(ib,jb,1), RAINNCV=dum5(ib,jb,1), &
               diag_effc3d=effc(ibr,jbr,kbr), diag_effi3d=effi(ibr,jbr,kbr), diag_dbz3d=dum8,  &
       diag_vmi3d_1=p3o(1,1,1, 1), diag_di3d_1=p3o(1,1,1, 2), diag_rhopo3d_1=p3o(1,1,1, 3), diag_phii3d_1=p3o(1,1,1, 4), &
       diag_vmi3d_2=p3o(1,1,1, 5), diag_di3d_2=p3o(1,1,1, 6), diag_rhopo3d_2=p3o(1,1,1, 7), diag_phii3d_2=p3o(1,1,1, 8), &
       diag_vmi3d_3=p3o(1,1,1, 9), diag_di3d_3=p3o(1,1,1,10), diag_rhopo3d_3=p3o(1,1,1,11), diag_phii3d_3=p3o(1,1,1,12), &
               nrain=nrain,dx=dx, dy=dy, cm1dz=dz,                        &
               tcond=qbudget(1),tevac=qbudget(2),                         &
               tevar=qbudget(5),train=qbudget(6),                         &
               ruh=ruh,rvh=rvh,rmh=rmh,rr=dum3,rain=rain,prate=prate,     &
               ib3d=ib3d,ie3d=ie3d,jb3d=jb3d,je3d=je3d,kb3d=kb3d,ke3d=ke3d,eqtset=eqtset, &
               nout3d=nout3d,out3d=out3d,dowr=dowr,outfile=outfile,getdbz=getdbz,dorad=dorad)

            do k=1,nk
            do j=1,nj
            do i=1,ni
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                 th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
                 rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
              endif
            enddo
            enddo
            enddo

            if( getdbz .and. qd_dbz.ge.1 )then
              do k=1,nk
              do j=1,nj
              do i=1,ni
                qdiag(i,j,k,qd_dbz) = dum8(i,j,k)
              enddo
              enddo
              enddo
            endif


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc     WSM6          ccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.7)THEN

          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq( qsmall,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq( qsmall,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq( qsmall,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq( qsmall,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq( qsmall,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))

      ! NOTE: need to OpenMPize !

          wsmjloop:  &
          DO j=1,nj

            do k=1,nk
            do i=1,ni
              tmp2d(i,k, 1) = thten(i,j,k)      ! t
              tmp2d(i,k, 2) = q3d(i,j,k,nqv)    ! q
              tmp2d(i,k, 3) = q3d(i,j,k,nqc)    ! qc
              tmp2d(i,k, 4) = q3d(i,j,k,nqi)    ! qi
              tmp2d(i,k, 5) = q3d(i,j,k,nqr)    ! qr
              tmp2d(i,k, 6) = q3d(i,j,k,nqs)    ! qs
              tmp2d(i,k, 7) = q3d(i,j,k,nqg)    ! qg
              tmp2d(i,k, 8) = rho(i,j,k)        ! den
              tmp2d(i,k, 9) = prs(i,j,k)        ! p
              tmp2d(i,k,10) = dz*rmh(i,j,k)     ! delz
              tmp2d(i,k,11) = 0.0
              tmp2d(i,k,12) = 0.0
              tmp2d(i,k,13) = 0.0
              tmp2d(i,k,14) = 0.0
              tmp2d(i,k,15) = 0.0
            enddo
            enddo

            do i=1,ni
              tmp1d(i,1) = 0.0
              tmp1d(i,2) = 0.0
              tmp1d(i,3) = 0.0
              tmp1d(i,4) = 0.0
              tmp1d(i,5) = 0.0
              tmp1d(i,6) = 0.0
              tmp1d(i,7) = 0.0
            enddo

       call mp_wsm6_wrapper(t=tmp2d(ibtmp,kbtmp,1),        &
                        q=tmp2d(ibtmp,kbtmp,2),            &
                        qc=tmp2d(ibtmp,kbtmp,3),           &
                        qi=tmp2d(ibtmp,kbtmp,4),           &
                        qr=tmp2d(ibtmp,kbtmp,5),           &
                        qs=tmp2d(ibtmp,kbtmp,6),           &
                        qg=tmp2d(ibtmp,kbtmp,7),           &
                        den=tmp2d(ibtmp,kbtmp,8),          &
                        p=tmp2d(ibtmp,kbtmp,9),            &
                        delz=tmp2d(ibtmp,kbtmp,10),        &
                        delt=dt,                           &
                        g=g,                               &
                        cpd=cp,                            &
                        cpv=cpv,                           &
                        rd=rd,                             &
                        rv=rv,                             &
                        t0c=to,                            &
                        ep1=repsm1,                        &
                        ep2=eps,                           &
                        qmin=qsmall,                       &
                        xls=xls,                           &
                        xlv0=xlv,                          &
                        xlf0=xlf,                          &
                        den0=1.28,                         &
                        denr=1000.0,                       &
                        cliq=cpl,                          &
                        cice=cpi,                          &
                        psat=611.2,                        &
                        rain=tmp1d(ibtmp,1),               &
                        rainncv=tmp1d(ibtmp,2),            &
                        sr=tmp1d(ibtmp,3),                 &
                        snow=tmp1d(ibtmp,4),               &
                        snowncv=tmp1d(ibtmp,5),            &
                        graupel=tmp1d(ibtmp,6),            &
                        graupelncv=tmp1d(ibtmp,7),         &
                        rainprod2d=tmp2d(ibtmp,kbtmp,11),  &
                        evapprod2d=tmp2d(ibtmp,kbtmp,12),  &
                        its=1,                             &
                        ite=ni,                            &
                        kts=1,                             &
                        kte=nk,                            &
                        errmsg=errmsg,                     &
                        errflg=errflg                      &
                       )

            do k=1,nk
            do i=1,ni
              dum1(i,j,k)    = tmp2d(i,k,1)
              q3d(i,j,k,nqv) = tmp2d(i,k,2)
              q3d(i,j,k,nqc) = tmp2d(i,k,3)
              q3d(i,j,k,nqi) = tmp2d(i,k,4)
              q3d(i,j,k,nqr) = tmp2d(i,k,5)
              q3d(i,j,k,nqs) = tmp2d(i,k,6)
              q3d(i,j,k,nqg) = tmp2d(i,k,7)
              if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall .or.  &
                  abs(q3d(i,j,k,nqv)-qten(i,j,k,nqv)).ge.qsmall )then
                th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
                rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))  ! or, modify pressure (?)
              endif
            enddo
            enddo

            do n = 1,nrain
            do i=1,ni
              rain(i,j,n) = rain(i,j,n) + 0.1*tmp1d(i,1)  ! convert from mm to cm
            enddo
            enddo

            IF( dorad .and. radopt.ge.1 )THEN
              if(myid.eq.0) print *,'  getting wsm6 effective radii for radiation '
              call mp_wsm6_effectRad_wrapper(          &
                         do_microp_re=.true.,          &
                         t=tmp2d(ibtmp,kbtmp,1),       &
                         qc=tmp2d(ibtmp,kbtmp,3),      &
                         qi=tmp2d(ibtmp,kbtmp,4),      &
                         qs=tmp2d(ibtmp,kbtmp,6),      &
                         rho=tmp2d(ibtmp,kbtmp,8),     &
                         qmin=qsmall,                  &
                         t0c=to,                       &
                         re_qc_bg=2.49E-6,             &
                         re_qi_bg=4.99E-6,             &
                         re_qs_bg=9.99E-6,             &
                         re_qc_max=50.E-6,             &
                         re_qi_max=125.E-6,            &
                         re_qs_max=999.E-6,            &
                         re_qc=tmp2d(ibtmp,kbtmp,13),  &
                         re_qi=tmp2d(ibtmp,kbtmp,14),  &
                         re_qs=tmp2d(ibtmp,kbtmp,15),  &
                         its=1,                        &
                         ite=ni,                       &
                         kts=1,                        &
                         kte=nk,                       &
                         errmsg=errmsg,                &
                         errflg=errflg)
              do k=1,nk
              do i=1,ni
                effc(i,j,k) = tmp2d(i,k,13)
                effi(i,j,k) = tmp2d(i,k,14)
                effs(i,j,k) = tmp2d(i,k,15)
              enddo
              enddo
            ENDIF

            IF( getdbz .and. qd_dbz.ge.1 )THEN
              do i = 1,ni
                do k = 1,nk
                  p1d(k) = prs(i,j,k)
                  t1d(k) = (th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
                  qv1d(k)  = q3d(i,j,k,nqv)
                  qr1d(k)  = q3d(i,j,k,nqr)
                  qs1d(k)  = q3d(i,j,k,nqs)
                  qg1d(k)  = q3d(i,j,k,nqg)
                  dBZ1d(k) = -35._sp
                  zp(k) = zh(i,j,k)
                enddo
                call refl10cm_wsm6(qv1d=qv1d,  &
                                   qr1d=qr1d,  &
                                   qs1d=qs1d,  &
                                   qg1d=qg1d,  &
                                   t1d=t1d,    &
                                   p1d=p1d,    &
                                   dBZ=dbz1d,  &
                                   kts=1,      &
                                   kte=nk)
                do k = 1,nk
                  qdiag(i,j,k,qd_dbz) = max(-35._sp,dBZ1d(k))
                enddo
              enddo
            ENDIF

          ENDDO  wsmjloop

          if(timestats.ge.1) time_microphy=time_microphy+mytime()

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  insert new microphysics schemes here
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        ELSEIF(ptype.eq.8)THEN
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! otherwise, stop for undefined ptype

        ELSE

          print *,'  Undefined ptype!'
          call stopcm1

        ENDIF  ifptype

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   END microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


          if( dotbud .and. td_mp.ge.1 )then
            rdt = 1.0/dt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              tdiag(i,j,k,td_mp) = (th3d(i,j,k)-tdiag(i,j,k,td_mp))*rdt
            enddo
            enddo
            enddo
          endif

          if( doqbud .and. qd_mp.ge.1 )then
            rdt = 1.0/dt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              qdiag(i,j,k,qd_mp) = (q3d(i,j,k,nqv)-qdiag(i,j,k,qd_mp))*rdt
            enddo
            enddo
            enddo
          endif

        if(timestats.ge.1) time_microphy=time_microphy+mytime()


    end subroutine mp_driver

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  END MODULE mp_driver_module
