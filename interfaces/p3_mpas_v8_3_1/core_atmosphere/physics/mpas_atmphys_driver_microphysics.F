! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_microphysics
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop
 use mpas_dmpar, only : IO_NODE
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

 use mpas_atmphys_constants
 use mpas_atmphys_init_microphysics
 use mpas_atmphys_interface
 use mpas_atmphys_vars

!wrf physics:
 use module_mp_kessler
 use module_mp_thompson
 use module_mp_wsm6,only: wsm6
 use mp_wsm6,only: mp_wsm6_init,refl10cm_wsm6
! P3 aerosol-aware microphysics (MPAS wrapper + init)
 use module_mp_p3,        only : p3_init_for_mpas, mp_p3_aerosols
 implicit none
 private
 public:: allocate_microphysics,   &
          deallocate_microphysics, &
          driver_microphysics,     &
          init_microphysics


!MPAS driver for parameterization of cloud microphysics processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_microphysics:
! ------------------------------------------------
! allocate_microphysics     : allocate local arrays for parameterization of cloud microphysics.
! deallocate_microphysics   : deallocate local arrays for parameterization of cloud microphysics.
! microphysics_init         : initialization of individual cloud microphysics schemes.
! driver_microphysics       : main driver (called from mpas_atm_time_integration).
! precip_from_MPAS          : initialize timestep local arrays for precipitation.
! precip_to_MPAS            : copy local arrays to MPAS arrays.
! compute_radar_reflectivity: compute radar reflectivities.
! compute_relhum            : compute relative humidity.
!
! WRF physics called from microphysics_driver:
! --------------------------------------------
! * module_mp_kessler : Kessler cloud microphysics.
! * module_mp_thompson: Thompson cloud microphysics.
! * module_mp_wsm6    : WSM6 cloud microphysics.
!
! comments:
! unlike all the other physics parameterizations, parameterizations of cloud microphysics schemes
! are called at the bottom of subroutine atm_srk3 in module atm_time_integration.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed call to the Thompson cloud microphysics scheme until scheme is updated to that in WRF revision 3.5.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added subroutine compute_relhum to calculate the relative humidity using the functions rslf and rsif from
!   the Thompson cloud microphysics scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-12.
! * removed the argument tend from the call to microphysics_from_MPAS (not needed).
!   Laura D. Fowler (laura@ucar.edu) / 2013-11-07.
! * in call to subroutine wsm6, replaced the variable g (that originally pointed to gravity) with gravity,
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * moved the variable relhum from the diag_physics to the diag pool. Changed the argument
!   list for the subroutine compute_relhum accordingly.
!   Laura D. Fowler (laura@ucar.edu) / 2015-04-22.
! * added parameterization of the Thompson cloud microphysics from WRF version 3.8.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-28.
! * in subroutine compute_relhum, multiply relhum by 100. so that it has the same unit as in the initial
!   conditions.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-20.
! * added parameterization of the WSM6 cloud microphysics from WRF version 3.8.1. To initialize WSM6 as in its
!   original version, set the hail_option to 0.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.
! * since we removed the local variable microp_scheme from mpas_atmphys_vars.F, now defines microp_scheme as a
!   pointer to config_microp_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.

!--- initialization option for WSM6 from WRF version 3.8.1. this option could also be set as a namelist parameter.
 integer,parameter:: hail_opt = 0


 contains


!=================================================================================================================
 subroutine allocate_microphysics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme

! --- P3 config pointers ---
 integer, pointer :: p3_nCat
 logical, pointer :: p3_trplMomI
 logical, pointer :: p3_liqfrac

  logical, pointer :: p3_debug_on_ptr
  logical, save :: p3dbg_alloc_banner_done = .false.
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_p3_debug_on',p3_debug_on_ptr)
  if (.not. p3dbg_alloc_banner_done .and. trim(microp_scheme) == 'mp_p3_aerosols' &
      .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) then
     call mpas_log_write('P3DBG_ALLOC: scheme mp_p3_aerosols')
     p3dbg_alloc_banner_done = .true.
  end if


!sounding variables:
 if(.not.allocated(rho_p) ) allocate(rho_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(th_p)  ) allocate(th_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pi_p)  ) allocate(pi_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_p)) allocate(pres_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(z_p)   ) allocate(z_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(dz_p)  ) allocate(dz_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(w_p)   ) allocate(w_p(ims:ime,kms:kme,jms:jme)   )

!mass mixing ratios (always present):
 if(.not.allocated(qv_p)) allocate(qv_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qc_p)) allocate(qc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qr_p)) allocate(qr_p(ims:ime,kms:kme,jms:jme))

!surface precipitation (always present):
 if(.not.allocated(rainnc_p) ) allocate(rainnc_p(ims:ime,jms:jme) )
 if(.not.allocated(rainncv_p)) allocate(rainncv_p(ims:ime,jms:jme))

 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6")
       !mass mixing ratios:
       if(.not.allocated(qi_p)) allocate(qi_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qs_p)) allocate(qs_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qg_p)) allocate(qg_p(ims:ime,kms:kme,jms:jme))

       !surface precipitation:
       if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
       if(.not.allocated(snownc_p)    ) allocate(snownc_p(ims:ime,jms:jme)    )
       if(.not.allocated(snowncv_p)   ) allocate(snowncv_p(ims:ime,jms:jme)   )
       if(.not.allocated(graupelnc_p) ) allocate(graupelnc_p(ims:ime,jms:jme) )
       if(.not.allocated(graupelncv_p)) allocate(graupelncv_p(ims:ime,jms:jme))

       !cloud water,cloud ice,and snow effective radii:
       if(.not.allocated(recloud_p)) allocate(recloud_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(reice_p)  ) allocate(reice_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(resnow_p) ) allocate(resnow_p(ims:ime,kms:kme,jms:jme) )

       !precipitation flux:
       if(.not.allocated(rainprod_p)) allocate(rainprod_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(evapprod_p)) allocate(evapprod_p(ims:ime,kms:kme,jms:jme))

    microp2_select: select case(trim(microp_scheme))
       case("mp_thompson","mp_thompson_aerosols")
          if(.not.allocated(ntc_p)) allocate(ntc_p(ims:ime,jms:jme))
          if(.not.allocated(muc_p)) allocate(muc_p(ims:ime,jms:jme))
          if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme))

          microp3_select: select case(trim(microp_scheme))
             case("mp_thompson_aerosols")
                if(.not.allocated(nifa2d_p)) allocate(nifa2d_p(ims:ime,jms:jme))
                if(.not.allocated(nwfa2d_p)) allocate(nwfa2d_p(ims:ime,jms:jme))
                if(.not.allocated(nc_p)    ) allocate(nc_p(ims:ime,kms:kme,jms:jme)  )
                if(.not.allocated(nifa_p)  ) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
                if(.not.allocated(nwfa_p)  ) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))

             case default
          end select microp3_select

       case default
    end select microp2_select

    case ("mp_p3_aerosols")

       !--- P3 needs nCat and (optionally) 3-moment and liquid-fraction ice arrays,
       call mpas_pool_get_config(configs,'config_p3_nCat'    ,p3_nCat    )
       if (p3_nCat < 1 .or. p3_nCat > 3) then
          write(mpas_err_message,'(A,I0,A)') 'config_p3_nCat=', p3_nCat, &
               ' unsupported (must be 1..3 with current Registry scalars)'
          call mpas_log_write(mpas_err_message)
       end if
       call mpas_pool_get_config(configs,'config_p3_trplMomI',p3_trplMomI)
       call mpas_pool_get_config(configs,'config_p3_liqfrac' ,p3_liqfrac )
    
              !--- P3 category prognostics live only in state%scalars (MPAS-compatible).
       !    DO NOT allocate full-domain 4-D qit/qni/qir/qib/qzi/qli work arrays here.
       !    (Memory fix: avoid duplicating category tracers.)
       ! Rain number concentration work array kept (needed by MPAS/WRF interface packing)
       if(.not.allocated(nr_p)) allocate(nr_p(ims:ime,kms:kme,jms:jme))

! Bulk hydrometeor work arrays (needed for MPAS bulk-compatible scalar layout)
if(.not.allocated(qi_p)) allocate(qi_p(ims:ime,kms:kme,jms:jme))
if(.not.allocated(qs_p)) allocate(qs_p(ims:ime,kms:kme,jms:jme))
if(.not.allocated(qg_p)) allocate(qg_p(ims:ime,kms:kme,jms:jme))
if(.not.allocated(ni_p)) allocate(ni_p(ims:ime,kms:kme,jms:jme))

! Column diagnostics (match MPAS bulk diag_physics fields)
if(.not.allocated(rainprod_p)) allocate(rainprod_p(ims:ime,kms:kme,jms:jme))
if(.not.allocated(evapprod_p)) allocate(evapprod_p(ims:ime,kms:kme,jms:jme))
if (allocated(rainprod_p)) rainprod_p = 0._RKIND
if (allocated(evapprod_p)) evapprod_p = 0._RKIND
       
        ! P3 NOTE: Do not allocate full-domain 4D vmi/dmi/rhoi/dhmax arrays in MPAS.
        ! They are not required for MPAS bulk-compatible outputs and can be very memory-expensive.
        ! Column-local vmi/dmi/rhoi/dhmax buffers are still passed to mp_p3_aerosols.
       !--------------------------------------------------------------------------------
       ! P3 category prognostics / diagnostics are NOT present in standard MPAS IC files.
       ! For cold-start runs we must guarantee they start from a defined state; otherwise
       ! uninitialized values can propagate NaNs/negatives into P3 and later cause SIGFPE.
       !--------------------------------------------------------------------------------
       if (allocated(nr_p))    nr_p    = 0._RKIND
    
       !--- aerosol-aware: prognostic droplet number + aerosol number mixing ratios
       if(.not.allocated(nc_p)    ) allocate(nc_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nwfa_p)  ) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nifa_p)  ) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(nwfa2d_p)) allocate(nwfa2d_p(ims:ime,jms:jme))
       if(.not.allocated(nifa2d_p)) allocate(nifa2d_p(ims:ime,jms:jme))
    
       !--- surface precipitation bookkeeping (MPAS expects these arrays to exist for "snow"/"sr")
       if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
       if(.not.allocated(snownc_p)    ) allocate(snownc_p(ims:ime,jms:jme)    )
       if(.not.allocated(snowncv_p)   ) allocate(snowncv_p(ims:ime,jms:jme)   )
       if(.not.allocated(graupelnc_p) ) allocate(graupelnc_p(ims:ime,jms:jme) )
       if(.not.allocated(graupelncv_p)) allocate(graupelncv_p(ims:ime,jms:jme))
    
       !--- If your radiation/cloud optics coupling expects these “effective radius” arrays,
       !    allocate them here and fill them later from P3 diagnostics (effc/effi/etc).
       if(.not.allocated(recloud_p)) allocate(recloud_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(reice_p)  ) allocate(reice_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(resnow_p) ) allocate(resnow_p(ims:ime,kms:kme,jms:jme) )

    case default
 end select microp_select

 end subroutine allocate_microphysics

!=================================================================================================================
 subroutine deallocate_microphysics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

!sounding variables:
 if(allocated(rho_p) ) deallocate(rho_p )
 if(allocated(th_p)  ) deallocate(th_p  )
 if(allocated(pi_p)  ) deallocate(pi_p  )
 if(allocated(pres_p)) deallocate(pres_p)
 if(allocated(z_p)   ) deallocate(z_p   )
 if(allocated(dz_p)  ) deallocate(dz_p  )
 if(allocated(w_p)   ) deallocate(w_p   )

!mass mixing ratios:
 if(allocated(qv_p)) deallocate(qv_p)
 if(allocated(qc_p)) deallocate(qc_p)
 if(allocated(qr_p)) deallocate(qr_p)

 !surface precipitation:
 if(allocated(rainnc_p) ) deallocate(rainnc_p )
 if(allocated(rainncv_p)) deallocate(rainncv_p)

 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6")
       !mass mixing ratios:
       if(allocated(qi_p)) deallocate(qi_p)
       if(allocated(qs_p)) deallocate(qs_p)
       if(allocated(qg_p)) deallocate(qg_p)

       !surface precipitation:
       if(allocated(sr_p)        ) deallocate(sr_p        )
       if(allocated(snownc_p)    ) deallocate(snownc_p    )
       if(allocated(snowncv_p)   ) deallocate(snowncv_p   )
       if(allocated(graupelnc_p) ) deallocate(graupelnc_p )
       if(allocated(graupelncv_p)) deallocate(graupelncv_p)

       !cloud water,cloud ice,and snow effective radii:
       if(allocated(recloud_p)) deallocate(recloud_p)
       if(allocated(reice_p)  ) deallocate(reice_p  )
       if(allocated(resnow_p) ) deallocate(resnow_p )

       !precipitation flux:
       if(allocated(rainprod_p)) deallocate(rainprod_p)
       if(allocated(evapprod_p)) deallocate(evapprod_p)

    microp2_select: select case(trim(microp_scheme))
       case("mp_thompson","mp_thompson_aerosols")
          if(allocated(ntc_p)) deallocate(ntc_p)
          if(allocated(muc_p)) deallocate(muc_p)
          if(allocated(ni_p) ) deallocate(ni_p )
          if(allocated(nr_p) ) deallocate(nr_p )

          microp3_select: select case(trim(microp_scheme))
             case("mp_thompson_aerosols")
                if(allocated(nifa2d_p)) deallocate(nifa2d_p)
                if(allocated(nwfa2d_p)) deallocate(nwfa2d_p)
                if(allocated(nc_p)    ) deallocate(nc_p    )
                if(allocated(nifa_p)  ) deallocate(nifa_p  )
                if(allocated(nwfa_p)  ) deallocate(nwfa_p  )

             case default
          end select microp3_select

       case default
    end select microp2_select

    case ("mp_p3_aerosols")

       !--- P3 prognostic ice-category fields (4-D: i,k,j,nCat)

       !--- Optional P3 prognostics (depending on config)

       !--- P3 diagnostics / work arrays

       !--- Aerosol-aware fields (3-D + 2-D)
       if(allocated(nc_p)    ) deallocate(nc_p    )
       if(allocated(nwfa_p)  ) deallocate(nwfa_p  )
       if(allocated(nifa_p)  ) deallocate(nifa_p  )
       if(allocated(nwfa2d_p)) deallocate(nwfa2d_p)
       if(allocated(nifa2d_p)) deallocate(nifa2d_p)
       if(allocated(nr_p)) deallocate(nr_p)

       !--- P3 MPAS bulk-compatible bulk hydrometeor / number work arrays
       if(allocated(qi_p)) deallocate(qi_p)
       if(allocated(qs_p)) deallocate(qs_p)
       if(allocated(qg_p)) deallocate(qg_p)
       if(allocated(ni_p)) deallocate(ni_p)

       !--- P3 column diagnostic work arrays
       if(allocated(rainprod_p)) deallocate(rainprod_p)
       if(allocated(evapprod_p)) deallocate(evapprod_p)

       !--- Surface precipitation bookkeeping (MPAS expects these arrays for snow/sr/graupel channels)
       if(allocated(sr_p)        ) deallocate(sr_p        )
       if(allocated(snownc_p)    ) deallocate(snownc_p    )
       if(allocated(snowncv_p)   ) deallocate(snowncv_p   )
       if(allocated(graupelnc_p) ) deallocate(graupelnc_p )
       if(allocated(graupelncv_p)) deallocate(graupelncv_p)

       !--- If you allocated effective-radius arrays for radiation coupling, deallocate them here
       if(allocated(recloud_p)) deallocate(recloud_p)
       if(allocated(reice_p)  ) deallocate(reice_p  )
       if(allocated(resnow_p) ) deallocate(resnow_p )

    case default
 end select microp_select

end subroutine deallocate_microphysics

!=================================================================================================================
 subroutine init_microphysics(dminfo,configs,mesh,state,time_lev,sfc_input,diag_physics)

!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state

!local pointer:
 logical,pointer:: do_restart
 character(len=StrKIND),pointer:: microp_scheme

! --- P3 config pointers (must exist in MPAS config pool) ---
 character(len=StrKIND),pointer:: p3_lookup_file_dir
 integer,pointer:: p3_nCat
 logical,pointer:: p3_trplMomI
 logical,pointer:: p3_liqfrac
 logical,pointer:: p3_abort_on_err
 logical,pointer:: p3_dowr
 integer,pointer:: p3_autoAccr_param


 ! --- Additional WRF-P3 runtime option pointers (to mirror WRF options) ---
  logical,pointer:: p3_debug_on
  logical,pointer:: p3_scpf_on
  logical,pointer:: p3_aerosols
  logical,pointer:: p3_predictNc
  real,pointer::    p3_clbfact_dep
  real,pointer::    p3_clbfact_sub
  real,pointer::    p3_scpf_pfrac
  real,pointer::    p3_scpf_resfact
  integer,pointer:: p3_freq3Ddiag
!CCPP-compliant flags:
  character(len=StrKIND):: errmsg
 integer:: errflg

!local:
 integer :: p3_stat
 ! --- P3 scalar-index sanity check locals (used once in init; no timestep cost) ---
 integer, pointer :: idx_qit1, idx_qit2, idx_qit3
 integer, pointer :: idx_qni1, idx_qni2, idx_qni3
 integer, pointer :: idx_qir1, idx_qib1
 integer :: nCat_layout
!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine init_microphysics:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0
 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_do_restart'   ,do_restart   )

 microp_select: select case(trim(microp_scheme))
    case("mp_thompson","mp_thompson_aerosols")
       call thompson_init(l_mp_tables)
       call init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag_physics)

       microp2_select: select case(trim(microp_scheme))
          case("mp_thompson_aerosols")
             call init_thompson_aerosols_forMPAS(do_restart,dminfo,mesh,state,time_lev,diag_physics)

          case default
       end select microp2_select

    case("mp_wsm6")
       call mp_wsm6_init(den0=rho_a,denr=rho_r,dens=rho_s,cl=cliq,cpv=cpv, &
                         hail_opt=hail_opt,errmsg=errmsg,errflg=errflg)

    case("mp_p3_aerosols")

       ! Initialize aerosol tracer profiles (nwfa/nifa/nc) for P3 when not restarting
       call init_p3_aerosols_forMPAS(do_restart,dminfo,mesh,state,time_lev,diag_physics)
       ! P3 lookup tables + parameter initialization (must run before first p3_main call)
       call mpas_pool_get_config(configs,'config_p3_lookup_file_dir',p3_lookup_file_dir)
       call mpas_pool_get_config(configs,'config_p3_nCat'           ,p3_nCat           )
       call mpas_pool_get_config(configs,'config_p3_trplMomI'       ,p3_trplMomI       )
       call mpas_pool_get_config(configs,'config_p3_liqfrac'        ,p3_liqfrac        )
       call mpas_pool_get_config(configs,'config_p3_abort_on_err'   ,p3_abort_on_err   )
       call mpas_pool_get_config(configs,'config_p3_dowr'           ,p3_dowr           )
       call mpas_pool_get_config(configs,'config_p3_autoAccr_param' ,p3_autoAccr_param )

        call mpas_pool_get_config(configs,'config_p3_debug_on'     ,p3_debug_on     )
        call mpas_pool_get_config(configs,'config_p3_scpf_on'      ,p3_scpf_on      )
        call mpas_pool_get_config(configs,'config_p3_aerosols'     ,p3_aerosols     )
        call mpas_pool_get_config(configs,'config_p3_predictNc'    ,p3_predictNc    )
        call mpas_pool_get_config(configs,'config_p3_clbfact_dep'  ,p3_clbfact_dep  )
        call mpas_pool_get_config(configs,'config_p3_clbfact_sub'  ,p3_clbfact_sub  )
        call mpas_pool_get_config(configs,'config_p3_scpf_pfrac'   ,p3_scpf_pfrac   )
        call mpas_pool_get_config(configs,'config_p3_scpf_resfact' ,p3_scpf_resfact )
        call mpas_pool_get_config(configs,'config_p3_freq3Ddiag'   ,p3_freq3Ddiag   )
       p3_stat = 0
       ! --- SANITY CHECKS (run once at init; negligible cost) ---
       ! Derive the true nCat from Registry scalar index spacing (MPAS-compatible robustness)
       call mpas_pool_get_dimension(state,'index_p3_qit', idx_qit1)
       call mpas_pool_get_dimension(state,'index_p3_qni', idx_qni1)
       call mpas_pool_get_dimension(state,'index_p3_qir', idx_qir1)
       call mpas_pool_get_dimension(state,'index_p3_qib', idx_qib1)

       nCat_layout = idx_qni1 - idx_qit1

       if (nCat_layout < 1 .or. nCat_layout > 3) then
          write(mpas_err_message,'(A,I0,A,4(A,I0))') 'ERROR: Derived P3 nCat_layout=', nCat_layout, &
               ' from scalar index spacing; indices: qit=', idx_qit1, ' qni=', idx_qni1, ' qir=', idx_qir1, ' qib=', idx_qib1
          call physics_error_fatal(trim(mpas_err_message))
       end if

       if ( (idx_qir1 - idx_qni1) /= nCat_layout .or. (idx_qib1 - idx_qir1) /= nCat_layout ) then
          write(mpas_err_message,'(A,4(A,I0),A,I0)') 'ERROR: P3 scalar indices not evenly spaced by nCat_layout. ', &
               'qit=', idx_qit1, ' qni=', idx_qni1, ' qir=', idx_qir1, ' qib=', idx_qib1, ' nCat_layout=', nCat_layout
          call physics_error_fatal(trim(mpas_err_message))
       end if

       if (nCat_layout /= p3_nCat) then
          write(mpas_err_message,'(A,I0,A,I0)') 'P3DBG_FATAL: nCat_layout=', nCat_layout, ' config_p3_nCat=', p3_nCat
          call mpas_log_write(trim(mpas_err_message))
          call mpas_log_write('P3DBG_FATAL: P3 scalar layout must match config_p3_nCat; check Registry package activation.')
          call physics_error_fatal(trim(mpas_err_message))
       end if

if (len_trim(p3_lookup_file_dir) == 0) then
          call mpas_log_write('WARNING: config_p3_lookup_file_dir is empty; P3 lookup init may fail.')
       end if
               call p3_init_for_mpas(trim(p3_lookup_file_dir), p3_nCat, p3_trplMomI, p3_liqfrac, &
                     stat=p3_stat, abort_on_err=p3_abort_on_err, model='WRF', dowr=p3_dowr,    &
                     autoAccr_param_in=p3_autoAccr_param, &
                     debug_on=p3_debug_on, scpf_on=p3_scpf_on, p3_aerosols=p3_aerosols, &
                     clbfact_dep=p3_clbfact_dep, clbfact_sub=p3_clbfact_sub, &
                     scpf_pfrac=p3_scpf_pfrac, scpf_resfact=p3_scpf_resfact, &
                     freq3Ddiag_in=p3_freq3Ddiag, predictNc=p3_predictNc)


       if (p3_stat /= 0) then
          write(mpas_err_message,'(A,I0)') 'ERROR: p3_init_for_mpas failed, stat=', p3_stat
          call physics_error_fatal(trim(mpas_err_message))
       end if

    case default
 end select microp_select

!call mpas_log_write('--- end subroutine init_microphysics:')

 end subroutine init_microphysics

!=================================================================================================================
 subroutine driver_microphysics(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,tend,itimestep,its,ite)
!=================================================================================================================

 use mpas_constants, only : rvord

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: time_lev
 integer,intent(in):: itimestep
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(inout):: tend

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme

!local variables and arrays:
 integer :: istep
 integer :: i,j,k,c
 integer :: ii_loc, jj_loc

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg


 integer :: p3_badcount
 character(len=StrKIND) :: p3_dbgmsg
 real(kind=RKIND) :: pratio



  ! P3 debug throttling (to prevent multi-GB logs)
  logical, save :: p3dbg_drv_banner_printed = .false.
  logical, save :: p3dbg_drv_diag_printed   = .false.
  integer, parameter :: P3DBG_MAX_BAD_PI = 10
  integer, save :: p3dbg_bad_pi_n = 0
  logical, save :: p3dbg_supp_bad_pi_msg = .false.
  logical, pointer :: p3_debug_on_ptr
  logical, save :: p3dbg_alloc_banner_done = .false.
 integer, pointer :: idx_qit, idx_qni, idx_qir, idx_qib, idx_qzi, idx_qli
 real(kind=RKIND), dimension(:,:,:), pointer :: scalars
 integer :: p3_nCat_use
 integer :: nCat_layout_drv
! --- P3 config ---
 integer, pointer :: p3_nCat_ptr
 logical, pointer :: p3_trplMomI_ptr, p3_liqfrac_ptr
 integer :: p3_nCat
 logical :: p3_trplMomI, p3_liqfrac
  logical :: rolling_old
! --- P3 column buffers (packed 1D/2D) ---
 real(kind=RKIND), allocatable :: th1d(:), qv1d(:), qc1d(:), qr1d(:), nr1d(:)
real(kind=RKIND), allocatable :: qv_in(:), qc_in(:), qr_in(:), qi_in(:), qs_in(:), qg_in(:)
 real(kind=RKIND), allocatable :: p1d(:),  w1d(:),  dz1d(:)
 real(kind=RKIND), allocatable :: pi1d(:), t1d(:)
 real(kind=RKIND), allocatable :: qi1d(:), ni1d(:), qs1d(:), qg1d(:)
 real(kind=RKIND), allocatable :: thold1d(:), qvold1d(:)
 real(kind=RKIND), allocatable :: nc1d(:), nwfa1d(:), nifa1d(:)
 real(kind=RKIND), allocatable :: rainprod1d(:), evapprod1d(:)
 real(kind=RKIND), allocatable :: zdbz1d(:), effc1d(:), effiave1d(:)
 real(kind=RKIND), allocatable :: qit1d(:,:), qni1d(:,:), qir1d(:,:), qib1d(:,:)
 real(kind=RKIND), allocatable :: qzi1d(:,:), qli1d(:,:)
 real(kind=RKIND), allocatable :: vmi1d(:,:), dmi1d(:,:), rhoi1d(:,:), dhmax1d(:,:)
 real(kind=RKIND) :: qv_min, qv_max, qc_min, qc_max, qr_min, qr_max
real(kind=RKIND) :: p3_pmin, p3_pmax, p3_pimin, p3_pimax
real(kind=RKIND) :: rain_step_mm, snow_step_mm, graup_step_mm, ice_step_mm, sr_step
real(kind=RKIND) :: aerosol_source_w, aerosol_source_i

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('---enter subroutine driver_microphysics:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0



 p3_badcount = 0
 p3_badcount = 0
 nullify(p3_debug_on_ptr)
 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_p3_debug_on',p3_debug_on_ptr)

!... allocation of microphysics arrays:
!$OMP MASTER
 call allocate_microphysics(configs)
!$OMP END MASTER
!$OMP BARRIER

!... initialization of precipitation related arrays:
 call precip_from_MPAS(configs,diag_physics,its,ite)

!... initialization of soundings for non-hydrostatic dynamical cores.
 call microphysics_from_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,its,ite)


! --- P3 sanity: verify reconstructed p/pi are finite (first column only)
 if (trim(microp_scheme) == 'mp_p3_aerosols') then
    p3_pmin  = minval(pres_p(its,kts:kte,jts))
    p3_pmax  = maxval(pres_p(its,kts:kte,jts))
    p3_pimin = minval(pi_p(its,kts:kte,jts))
    p3_pimax = maxval(pi_p(its,kts:kte,jts))
    write(errmsg,'(A,4(1X,ES12.4))') 'P3DBG_AFTER_MPAS2PHY pmin,pmax,pimin,pimax=', p3_pmin,p3_pmax,p3_pimin,p3_pimax
    call mpas_log_write(trim(errmsg))
 end if
!... call to different cloud microphysics schemes:
 microp_select: select case(trim(microp_scheme))

    case ("mp_kessler")
       call mpas_timer_start('mp_kessler')
       call kessler( &
                  t        = th_p      , qv    = qv_p  , qc     = qc_p     ,                &
                  qr       = qr_p      , rho   = rho_p , pii    = pi_p     ,                &
                  dt_in    = dt_microp , z     = z_p   , xlv    = xlv      ,                &
                  cp       = cp        , ep2   = ep_2  , svp1   = svp1     ,                &
                  svp2     = svp2      , svp3  = svp3  , svpt0  = svpt0    ,                &
                  rhowater = rho_w     , dz8w  = dz_p  , rainnc = rainnc_p ,                &
                  rainncv  = rainncv_p ,                                                    &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde   , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme   , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                   )
       call mpas_timer_stop('mp_kessler')

    case ("mp_thompson")
       call mpas_timer_start('mp_thompson')
       istep = 1
       do while (istep .le. n_microp)
       call mp_gt_driver( &
                  th        = th_p        , qv         = qv_p         , qc         = qc_p         , &
                  qr        = qr_p        , qi         = qi_p         , qs         = qs_p         , &
                  qg        = qg_p        , ni         = ni_p         , nr         = nr_p         , &
                  pii       = pi_p        , p          = pres_p       , dz         = dz_p         , &
                  w         = w_p         , dt_in      = dt_microp    , itimestep  = itimestep    , &
                  rainnc    = rainnc_p    , rainncv    = rainncv_p    , snownc     = snownc_p     , &
                  snowncv   = snowncv_p   , graupelnc  = graupelnc_p  , graupelncv = graupelncv_p , &
                  sr        = sr_p        , rainprod   = rainprod_p   , evapprod   = evapprod_p   , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow    = resnow_p     , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs   = has_reqs     , &
                  ntc       = ntc_p       , muc        = muc_p                                    , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde           , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme           , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                           )
       istep = istep + 1
       enddo
       call mpas_timer_stop('mp_thompson')

    case ("mp_thompson_aerosols")
       call mpas_timer_start('mp_thompson_aerosols')
       istep = 1
       do while (istep .le. n_microp)
       call mp_gt_driver( &
                  th        = th_p        , qv         = qv_p         , qc         = qc_p         , &
                  qr        = qr_p        , qi         = qi_p         , qs         = qs_p         , &
                  qg        = qg_p        , ni         = ni_p         , nr         = nr_p         , &
                  pii       = pi_p        , p          = pres_p       , dz         = dz_p         , &
                  w         = w_p         , dt_in      = dt_microp    , itimestep  = itimestep    , &
                  rainnc    = rainnc_p    , rainncv    = rainncv_p    , snownc     = snownc_p     , &
                  snowncv   = snowncv_p   , graupelnc  = graupelnc_p  , graupelncv = graupelncv_p , &
                  sr        = sr_p        , rainprod   = rainprod_p   , evapprod   = evapprod_p   , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow    = resnow_p     , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs   = has_reqs     , &
                  nc        = nc_p        , nifa       = nifa_p       , nwfa       = nwfa_p       , &
                  nifa2d    = nifa2d_p    , nwfa2d     = nwfa2d_p     , ntc        = ntc_p        , &
                  muc       = muc_p       ,                                                         &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde           , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme           , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                           )
          istep = istep + 1
          enddo
          call mpas_timer_stop('mp_thompson_aerosols')

    case ("mp_wsm6")
       call mpas_timer_start('mp_wsm6')
       call wsm6( &
                  th        = th_p        , q          = qv_p         , qc        = qc_p      , &
                  qr        = qr_p        , qi         = qi_p         , qs        = qs_p      , &
                  qg        = qg_p        , den        = rho_p        , pii       = pi_p      , &
                  p         = pres_p      , delz       = dz_p         , delt      = dt_microp , &
                  g         = gravity     , cpd        = cp           , cpv       = cpv       , &
                  rd        = R_d         , rv         = R_v          , t0c       = svpt0     , &
                  ep1       = ep_1        , ep2        = ep_2         , qmin      = epsilon   , &
                  xls       = xls         , xlv0       = xlv          , xlf0      = xlf       , &
                  den0      = rho_a       , denr       = rho_w        , cliq      = cliq      , &
                  cice      = cice        , psat       = psat         , rain      = rainnc_p  , &
                  rainncv   = rainncv_p   , snow       = snownc_p     , snowncv   = snowncv_p , &
                  graupel   = graupelnc_p , graupelncv = graupelncv_p , sr        = sr_p      , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow   = resnow_p  , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs  = has_reqs  , &
                  re_qc_bg  = re_qc_bg    , re_qi_bg   = re_qi_bg     , re_qs_bg  = re_qs_bg  , &
                  re_qc_max = re_qc_max   , re_qi_max  = re_qi_max    , re_qs_max = re_qs_max , &
                  errmsg    = errmsg      , errflg     = errflg                               , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde       , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme       , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte         &
               )
       call mpas_timer_stop('mp_wsm6')

    case ("mp_p3_aerosols")
       call mpas_timer_start('mp_p3_aerosols')
        if (itimestep == 1 .and. its == 1 .and. jts == 1) then
           write(errmsg,'(A,6(I0,1X))') 'P3DBG_STAGE: P3 pack dims ims ime jms jme kts kte =', ims,ime,jms,jme,kts,kte
           call mpas_log_write(errmsg)
        end if

       call mpas_pool_get_config(configs,'config_p3_nCat',p3_nCat_ptr)
        p3_nCat = p3_nCat_ptr
call mpas_pool_get_dimension(state,'index_p3_qit', idx_qit)
call mpas_pool_get_dimension(state,'index_p3_qni', idx_qni)
call mpas_pool_get_dimension(state,'index_p3_qir', idx_qir)
call mpas_pool_get_dimension(state,'index_p3_qib', idx_qib)
nCat_layout_drv = idx_qni - idx_qit
        ! MPAS bulk-compatible + memory fix: read/write P3 category prognostics directly from state%scalars
        call mpas_pool_get_array(state,'scalars',scalars,time_lev)
        if (.not. associated(scalars)) then
           call physics_error_fatal('P3DBG_FATAL: state%scalars not associated in driver_microphysics')
        end if
        ! Optional P3 indices (only used when enabled)
        call mpas_pool_get_dimension(state,'index_p3_qzi', idx_qzi)
        call mpas_pool_get_dimension(state,'index_p3_qli', idx_qli)

        ! WRF-P3 parity: honor config_p3_nCat (dummy nCat passed into the P3 wrapper).
        ! The Registry package setup should now allocate exactly the active category count.
        p3_nCat_use = p3_nCat

        if (p3_nCat_use < 1 .or. p3_nCat_use > 3) then
           write(errmsg,'(A,I0,A)') 'P3DBG_FATAL: config_p3_nCat=', p3_nCat_use, ' is outside the supported range 1..3'
           call mpas_log_write(errmsg)
           call physics_error_fatal(trim(errmsg))
        end if

        if (nCat_layout_drv < 1 .or. nCat_layout_drv > 3) then
           write(errmsg,'(A,I0,A,4(A,I0))') 'P3DBG_FATAL: bad nCat_layout=', nCat_layout_drv, &
                ' indices qit/qni/qir/qib=', idx_qit, idx_qni, idx_qir, idx_qib
           call mpas_log_write(errmsg)
           call physics_error_fatal(trim(errmsg))
        end if

        if (nCat_layout_drv /= p3_nCat_use) then
           write(errmsg,'(A,I0,A,I0)') 'P3DBG_FATAL: scalar layout nCat_layout=', nCat_layout_drv, &
                ' != config_p3_nCat=', p3_nCat_use
           call mpas_log_write(errmsg)
           call mpas_log_write('P3DBG_FATAL: inactive P3 categories were allocated; check package setup before scalar allocation.')
           call physics_error_fatal(trim(errmsg))
        end if

        p3_nCat = p3_nCat_use

call mpas_pool_get_config(configs,'config_p3_trplMomI',p3_trplMomI_ptr)
        p3_trplMomI = p3_trplMomI_ptr
        call mpas_pool_get_config(configs,'config_p3_liqfrac',p3_liqfrac_ptr)
        p3_liqfrac = p3_liqfrac_ptr
        
        if (associated(p3_debug_on_ptr) .and. p3_debug_on_ptr .and. .not. p3dbg_drv_banner_printed) then
           p3dbg_drv_banner_printed = .true.
        end if
                if (.not. allocated(nc_p)) then
           stop 'nc_p not allocated in driver_microphysics'
        end if
        if (.not. allocated(nifa_p) .or. .not. allocated(nwfa_p)) then
           stop 'P3 aerosol number arrays not allocated in driver_microphysics'
        end if
       ! Column buffers
       allocate(th1d(kts:kte), qv1d(kts:kte), qc1d(kts:kte), qr1d(kts:kte), nr1d(kts:kte))
       allocate(qv_in(kts:kte), qc_in(kts:kte), qr_in(kts:kte), qi_in(kts:kte), qs_in(kts:kte), qg_in(kts:kte))
       allocate(pi1d(kts:kte), t1d(kts:kte))
       allocate(qi1d(kts:kte), ni1d(kts:kte), qs1d(kts:kte), qg1d(kts:kte))
       allocate(thold1d(kts:kte), qvold1d(kts:kte))    ! <<< NEW: "old" state buffers
       allocate(p1d (kts:kte), w1d (kts:kte), dz1d(kts:kte))
       allocate(nc1d(kts:kte), nwfa1d(kts:kte), nifa1d(kts:kte))
       allocate(rainprod1d(kts:kte), evapprod1d(kts:kte))
       allocate(zdbz1d(kts:kte), effc1d(kts:kte), effiave1d(kts:kte))

       ! Never carry gross process rates from a prior dynamics step.
       rainprod_p(its:ite,kts:kte,jts:jte) = 0._RKIND
       evapprod_p(its:ite,kts:kte,jts:jte) = 0._RKIND

       allocate(qit1d(kts:kte,p3_nCat), qni1d(kts:kte,p3_nCat), qir1d(kts:kte,p3_nCat), qib1d(kts:kte,p3_nCat))
       allocate(qzi1d(kts:kte,p3_nCat), qli1d(kts:kte,p3_nCat))
       allocate(vmi1d(kts:kte,p3_nCat), dmi1d(kts:kte,p3_nCat), rhoi1d(kts:kte,p3_nCat), dhmax1d(kts:kte,p3_nCat))
       ! Choose meaning of "old":
       !   - If .false.: old = state at start of the entire microphysics subcycling (baseline)
       rolling_old = .false.

       do j = jts, jte
          do i = its, ite

             !------------------------------------------
             ! Pack MPAS column ONCE for this (i,j)
             !------------------------------------------
             do k = kts, kte
                th1d(k) = th_p(i,k,j)
                qv1d(k) = qv_p(i,k,j)
                qc1d(k) = qc_p(i,k,j)
                qr1d(k) = qr_p(i,k,j)
                nr1d(k) = nr_p(i,k,j)

                p1d(k)  = pres_p(i,k,j)
                w1d(k)  = w_p(i,k,j)
                if (.not. ieee_is_finite(w1d(k))) w1d(k) = 0.0_RKIND
                dz1d(k) = dz_p(i,k,j)

                pi1d(k) = pi_p(i,k,j)
                ! Temperature (K) from theta * Exner
                t1d(k)  = th1d(k) * pi1d(k)

                !--------------------------------------------------------------------------------
                ! MPAS bulk-compatible cold-start:
                ! If MPAS bulk ice (qi_p/ni_p) is present but P3 category tracers in scalars are
                ! missing/zero (typical for IC files), seed category-1 and zero the rest.
                ! Also sanitize any non-finite/negative category fields to prevent NaN propagation.
                !--------------------------------------------------------------------------------
                block
                   real(kind=RKIND) :: sum_qit_cat, sum_qni_cat, vtmp, qi_seed
                   integer :: cc
                   sum_qit_cat = 0.0_RKIND
                   sum_qni_cat = 0.0_RKIND
                   do cc = 1, nCat_layout_drv
                      vtmp = scalars(idx_qit + (cc-1), k, i)
                      if (.not. ieee_is_finite(vtmp) .or. vtmp < 0._RKIND) vtmp = 0._RKIND
                      scalars(idx_qit + (cc-1), k, i) = vtmp
                      sum_qit_cat = sum_qit_cat + vtmp

                      vtmp = scalars(idx_qni + (cc-1), k, i)
                      if (.not. ieee_is_finite(vtmp) .or. vtmp < 0._RKIND) vtmp = 0._RKIND
                      scalars(idx_qni + (cc-1), k, i) = vtmp
                      sum_qni_cat = sum_qni_cat + vtmp

                      vtmp = scalars(idx_qir + (cc-1), k, i)
                      if (.not. ieee_is_finite(vtmp) .or. vtmp < 0._RKIND) vtmp = 0._RKIND
                      scalars(idx_qir + (cc-1), k, i) = vtmp

                      vtmp = scalars(idx_qib + (cc-1), k, i)
                      if (.not. ieee_is_finite(vtmp) .or. vtmp < 0._RKIND) vtmp = 0._RKIND
                      scalars(idx_qib + (cc-1), k, i) = vtmp

                      if (p3_trplMomI) then
                         vtmp = scalars(idx_qzi + (cc-1), k, i)
                         if (.not. ieee_is_finite(vtmp) .or. vtmp < 0._RKIND) vtmp = 0._RKIND
                         scalars(idx_qzi + (cc-1), k, i) = vtmp
                      else
                         scalars(idx_qzi + (cc-1), k, i) = 0._RKIND
                      end if

                      if (p3_liqfrac) then
                         vtmp = scalars(idx_qli + (cc-1), k, i)
                         if (.not. ieee_is_finite(vtmp) .or. vtmp < 0._RKIND) vtmp = 0._RKIND
                         scalars(idx_qli + (cc-1), k, i) = vtmp
                      else
                         scalars(idx_qli + (cc-1), k, i) = 0._RKIND
                      end if
                   end do

                   qi_seed = 0.0_RKIND
                   if (ieee_is_finite(qi_p(i,k,j)) .and. qi_p(i,k,j) > 0._RKIND) qi_seed = qi_seed + qi_p(i,k,j)
                   if (ieee_is_finite(qs_p(i,k,j)) .and. qs_p(i,k,j) > 0._RKIND) qi_seed = qi_seed + qs_p(i,k,j)
                   if (ieee_is_finite(qg_p(i,k,j)) .and. qg_p(i,k,j) > 0._RKIND) qi_seed = qi_seed + qg_p(i,k,j)
                   if (sum_qit_cat <= 1.0e-30_RKIND .and. qi_seed > 0._RKIND) then
                      scalars(idx_qit, k, i) = qi_seed
                      do cc = 2, nCat_layout_drv
                         scalars(idx_qit + (cc-1), k, i) = 0._RKIND
                      end do
                   end if
                   if (sum_qni_cat <= 1.0e-30_RKIND .and. ieee_is_finite(ni_p(i,k,j)) .and. ni_p(i,k,j) > 0._RKIND) then
                      scalars(idx_qni, k, i) = ni_p(i,k,j)
                      do cc = 2, nCat_layout_drv
                         scalars(idx_qni + (cc-1), k, i) = 0._RKIND
                      end do
                   end if

                   ! Force unused categories (beyond config nCat) to 0 every call (MPAS-compatible)
                   do cc = p3_nCat+1, nCat_layout_drv
                      scalars(idx_qit + (cc-1), k, i) = 0._RKIND
                      scalars(idx_qni + (cc-1), k, i) = 0._RKIND
                      scalars(idx_qir + (cc-1), k, i) = 0._RKIND
                      scalars(idx_qib + (cc-1), k, i) = 0._RKIND
                      scalars(idx_qzi + (cc-1), k, i) = 0._RKIND
                      scalars(idx_qli + (cc-1), k, i) = 0._RKIND
                   end do
                end block

                ! Default: single-category ice for MPAS wrapper
                                qi1d(k) = 0.0_RKIND
                ni1d(k) = 0.0_RKIND
                do c = 1, p3_nCat
                   qi1d(k) = qi1d(k) + max(scalars(idx_qit + (c-1), k, i), 0.0_RKIND)
                   ni1d(k) = ni1d(k) + max(scalars(idx_qni + (c-1), k, i), 0.0_RKIND)
                end do
                qs1d(k) = 0.0_RKIND
                qg1d(k) = 0.0_RKIND
                nc1d(k)   = nc_p(i,k,j)
                nwfa1d(k) = nwfa_p(i,k,j)
                nifa1d(k) = nifa_p(i,k,j)

             end do

                          do c = 1, p3_nCat
                do k = kts, kte
                   qit1d(k,c) = scalars(idx_qit + (c-1), k, i)
                   qni1d(k,c) = scalars(idx_qni + (c-1), k, i)
                   qir1d(k,c) = scalars(idx_qir + (c-1), k, i)
                   qib1d(k,c) = scalars(idx_qib + (c-1), k, i)
                   if (p3_trplMomI) then
                      qzi1d(k,c) = scalars(idx_qzi + (c-1), k, i)
                   else
                      qzi1d(k,c) = 0.0_RKIND
                   end if
                   if (p3_liqfrac) then
                      qli1d(k,c) = scalars(idx_qli + (c-1), k, i)
                   else
                      qli1d(k,c) = 0.0_RKIND
                   end if
vmi1d(k,c)   = 0.0_RKIND
                   dmi1d(k,c)   = 0.0_RKIND
                   rhoi1d(k,c)  = 0.0_RKIND
                   dhmax1d(k,c) = 0.0_RKIND
                end do
             end do

             !------------------------------------------
             ! Define "old" for this column
             !------------------------------------------
             thold1d(:) = th1d(:)
             qvold1d(:) = qv1d(:)

             ! Safety init for reflectivity diagnostic passed into P3
             do k = kts, kte
                zdbz1d(k) = -35.0_RKIND
             end do

             !------------------------------------------
             ! Subcycling loop for this single (i,j) column
             !------------------------------------------
             do istep = 1, n_microp

                rain_step_mm  = 0.0_RKIND
                snow_step_mm  = 0.0_RKIND
                graup_step_mm = 0.0_RKIND
                ice_step_mm   = 0.0_RKIND
                sr_step       = 0.0_RKIND
                rainprod1d    = 0.0_RKIND
                evapprod1d    = 0.0_RKIND

                !--------------------------------------------------------------------------------
                ! P3 safety: enforce finite / nonnegative inputs before calling P3.
                ! This prevents uninitialized / negative category fields from triggering SIGFPE
                ! (e.g., pow() invalid) inside P3 or downstream dynamics.
                !--------------------------------------------------------------------------------
                do k = kts, kte
                   if (.not. ieee_is_finite(qv1d(k))) qv1d(k) = qv_in(k)
                   if (.not. ieee_is_finite(qv1d(k))) qv1d(k) = 0._RKIND
                   if (qv1d(k) < 0._RKIND) qv1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qc1d(k))) qc1d(k) = qc_in(k)
                   if (.not. ieee_is_finite(qc1d(k))) qc1d(k) = 0._RKIND
                   if (qc1d(k) < 0._RKIND) qc1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qr1d(k))) qr1d(k) = qr_in(k)
                   if (.not. ieee_is_finite(qr1d(k))) qr1d(k) = 0._RKIND
                   if (qr1d(k) < 0._RKIND) qr1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qi1d(k))) qi1d(k) = qi_in(k)
                   if (.not. ieee_is_finite(qi1d(k))) qi1d(k) = 0._RKIND
                   if (qi1d(k) < 0._RKIND) qi1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qs1d(k))) qs1d(k) = qs_in(k)
                   if (.not. ieee_is_finite(qs1d(k))) qs1d(k) = 0._RKIND
                   if (qs1d(k) < 0._RKIND) qs1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qg1d(k))) qg1d(k) = qg_in(k)
                   if (.not. ieee_is_finite(qg1d(k))) qg1d(k) = 0._RKIND
                   if (qg1d(k) < 0._RKIND) qg1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nc1d(k))) nc1d(k) = 0._RKIND
                   if (nc1d(k) < 0._RKIND) nc1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nr1d(k))) nr1d(k) = 0._RKIND
                   if (nr1d(k) < 0._RKIND) nr1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(ni1d(k))) ni1d(k) = 0._RKIND
                   if (ni1d(k) < 0._RKIND) ni1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nwfa1d(k))) then
                      nwfa1d(k) = 0._RKIND
                   end if
                   if (nwfa1d(k) < 0._RKIND) then
                      nwfa1d(k) = 0._RKIND
                   end if
                   if (.not. ieee_is_finite(nifa1d(k))) then
                      nifa1d(k) = 0._RKIND
                   end if
                   if (nifa1d(k) < 0._RKIND) then
                      nifa1d(k) = 0._RKIND
                   end if
                   if (.not. ieee_is_finite(p1d(k))) p1d(k) = 1._RKIND
                   if (p1d(k) < 1._RKIND) p1d(k) = 1._RKIND
                   if (.not. ieee_is_finite(dz1d(k))) dz1d(k) = 1._RKIND
                   if (dz1d(k) <= 0._RKIND) dz1d(k) = 1._RKIND
                   if (.not. ieee_is_finite(t1d(k))) t1d(k) = 300._RKIND
                   if (t1d(k) < 100._RKIND) t1d(k) = 100._RKIND
                    if (t1d(k) > 450._RKIND) t1d(k) = 450._RKIND
                   do c = 1, p3_nCat
                      if (.not. ieee_is_finite(qit1d(k,c))) qit1d(k,c) = 0._RKIND
                      if (qit1d(k,c) < 0._RKIND) qit1d(k,c) = 0._RKIND
                      if (.not. ieee_is_finite(qni1d(k,c))) qni1d(k,c) = 0._RKIND
                      if (qni1d(k,c) < 0._RKIND) qni1d(k,c) = 0._RKIND
                      if (.not. ieee_is_finite(qir1d(k,c))) qir1d(k,c) = 0._RKIND
                      if (qir1d(k,c) < 0._RKIND) qir1d(k,c) = 0._RKIND
                      if (.not. ieee_is_finite(qib1d(k,c))) qib1d(k,c) = 0._RKIND
                      if (qib1d(k,c) < 0._RKIND) qib1d(k,c) = 0._RKIND
                       if (allocated(qzi1d)) then
                         if (.not. ieee_is_finite(qzi1d(k,c))) qzi1d(k,c) = 0._RKIND
                         if (qzi1d(k,c) < 0._RKIND) qzi1d(k,c) = 0._RKIND
                      end if
                      if (allocated(qli1d)) then
                         if (.not. ieee_is_finite(qli1d(k,c))) qli1d(k,c) = 0._RKIND
                         if (qli1d(k,c) < 0._RKIND) qli1d(k,c) = 0._RKIND
                      end if
                   end do
                end do
                 ! (memfix/physfix) Removed misplaced pre-call T->theta conversion block.
        ! Snapshot inputs for NaN recovery (MPAS-compatible: do not zero the column on non-finite returns)
        do k = kts, kte
           qv_in(k) = qv1d(k)
           qc_in(k) = qc1d(k)
           qr_in(k) = qr1d(k)
           qi_in(k) = qi1d(k)
           qs_in(k) = qs1d(k)
           qg_in(k) = qg1d(k)
        end do

                ii_loc = i
                jj_loc = j

                call mp_p3_aerosols( &
                     qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, nr1d, nc1d, &
                     nwfa1d, nifa1d, &
                     thold1d, qvold1d, &
                     t1d, p1d, w1d, dz1d, &
                     zdbz1d, effc1d, effiave1d, &
                     qit1d, qni1d, qir1d, qib1d, qzi1d, qli1d, &
                     vmi1d, dmi1d, rhoi1d, dhmax1d, &
                     rain_step_mm, snow_step_mm, graup_step_mm, ice_step_mm, &
                     itimestep=itimestep, kts=kts, kte=kte, nCat=p3_nCat, dt=dt_microp, &
                     ii=ii_loc, jj=jj_loc, rainprod1d=rainprod1d, evapprod1d=evapprod1d )


                 ! P3 pressure is treated as input-only; reset in case the kernel modifies it inadvertently.
                 do k = kts, kte
                    p1d(k) = pres_p(i,k,j)

                    ! Match the existing Thompson/MPAS contract: when microphysics is
                    ! subcycled, rainprod and evapprod contain the final substep's
                    ! post-limiter instantaneous rates.
                    if (.not. ieee_is_finite(rainprod1d(k)) .or. rainprod1d(k) < 0._RKIND) then
                       rainprod1d(k) = 0._RKIND
                    endif
                    if (.not. ieee_is_finite(evapprod1d(k)) .or. evapprod1d(k) < 0._RKIND) then
                       evapprod1d(k) = 0._RKIND
                    endif
                    rainprod_p(i,k,j) = rainprod1d(k)
                    evapprod_p(i,k,j) = evapprod1d(k)
                 end do
                !--------------------------------------------------------------------------------
                ! P3 safety: sanitize outputs (finite / nonnegative) before using them further.
                !--------------------------------------------------------------------------------
                do k = kts, kte
                   if (.not. ieee_is_finite(qv1d(k))) qv1d(k) = 0._RKIND
                   if (qv1d(k) < 0._RKIND) qv1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qc1d(k))) qc1d(k) = 0._RKIND
                   if (qc1d(k) < 0._RKIND) qc1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qr1d(k))) qr1d(k) = 0._RKIND
                   if (qr1d(k) < 0._RKIND) qr1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qi1d(k))) qi1d(k) = 0._RKIND
                   if (qi1d(k) < 0._RKIND) qi1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qs1d(k))) qs1d(k) = 0._RKIND
                   if (qs1d(k) < 0._RKIND) qs1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(qg1d(k))) qg1d(k) = 0._RKIND
                   if (qg1d(k) < 0._RKIND) qg1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nc1d(k))) nc1d(k) = 0._RKIND
                   if (nc1d(k) < 0._RKIND) nc1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nr1d(k))) nr1d(k) = 0._RKIND
                   if (nr1d(k) < 0._RKIND) nr1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(ni1d(k))) ni1d(k) = 0._RKIND
                   if (ni1d(k) < 0._RKIND) ni1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nwfa1d(k))) nwfa1d(k) = 0._RKIND
                   if (nwfa1d(k) < 0._RKIND) nwfa1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(nifa1d(k))) nifa1d(k) = 0._RKIND
                   if (nifa1d(k) < 0._RKIND) nifa1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(zdbz1d(k))) zdbz1d(k) = -999._RKIND
                   if (.not. ieee_is_finite(effc1d(k))) effc1d(k) = 0._RKIND
                   if (effc1d(k) < 0._RKIND) effc1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(effiave1d(k))) effiave1d(k) = 0._RKIND
                   if (effiave1d(k) < 0._RKIND) effiave1d(k) = 0._RKIND
                   if (.not. ieee_is_finite(pi1d(k))) pi1d(k) = 1._RKIND
                   if (pi1d(k) < 1.0e-6_RKIND) pi1d(k) = 1.0e-6_RKIND
                   if (.not. ieee_is_finite(t1d(k))) t1d(k) = 300._RKIND
                   if (t1d(k) < 100._RKIND) t1d(k) = 100._RKIND
                    if (t1d(k) > 450._RKIND) t1d(k) = 450._RKIND
                   do c = 1, p3_nCat
                      if (.not. ieee_is_finite(qit1d(k,c))) qit1d(k,c) = 0._RKIND
                      if (qit1d(k,c) < 0._RKIND) qit1d(k,c) = 0._RKIND
                      if (.not. ieee_is_finite(qni1d(k,c))) qni1d(k,c) = 0._RKIND
                      if (qni1d(k,c) < 0._RKIND) qni1d(k,c) = 0._RKIND
                      if (.not. ieee_is_finite(qir1d(k,c))) qir1d(k,c) = 0._RKIND
                      if (qir1d(k,c) < 0._RKIND) qir1d(k,c) = 0._RKIND
                      if (.not. ieee_is_finite(qib1d(k,c))) qib1d(k,c) = 0._RKIND
                      if (qib1d(k,c) < 0._RKIND) qib1d(k,c) = 0._RKIND
                       if (allocated(qzi1d)) then
                         if (.not. ieee_is_finite(qzi1d(k,c))) qzi1d(k,c) = 0._RKIND
                         if (qzi1d(k,c) < 0._RKIND) qzi1d(k,c) = 0._RKIND
                      end if
                      if (allocated(qli1d)) then
                         if (.not. ieee_is_finite(qli1d(k,c))) qli1d(k,c) = 0._RKIND
                         if (qli1d(k,c) < 0._RKIND) qli1d(k,c) = 0._RKIND
                      end if
                   end do
                end do

                ! MPAS bulk-compatible aerosol-aware source update:
                ! apply the diagnosed lowest-level aerosol source tendency after each P3 call.
                ! This mirrors MPAS Thompson aerosol-aware behavior and closes a missing part of
                ! the MPAS P3 aerosol-aware implementation.
                aerosol_source_w = 0._RKIND
                aerosol_source_i = 0._RKIND
                if (allocated(nwfa2d_p)) aerosol_source_w = nwfa2d_p(i,j)
                if (allocated(nifa2d_p)) aerosol_source_i = nifa2d_p(i,j)
                if (.not. ieee_is_finite(aerosol_source_w)) aerosol_source_w = 0._RKIND
                if (.not. ieee_is_finite(aerosol_source_i)) aerosol_source_i = 0._RKIND
                if (kts <= kte) then
                   nwfa1d(kts) = max(0._RKIND, nwfa1d(kts) + aerosol_source_w * dt_microp)
                   nifa1d(kts) = max(0._RKIND, nifa1d(kts) + aerosol_source_i * dt_microp)
                   if (i == its .and. j == jts .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) then
                      write(errmsg,'(A,4ES12.4)') 'P3DBG_AER_SRC: nwfa2d,nifa2d,nwfa(kts),nifa(kts)=', &
                           aerosol_source_w, aerosol_source_i, nwfa1d(kts), nifa1d(kts)
                      call mpas_log_write(trim(errmsg))
                   end if
                end if

                if (.not. ieee_is_finite(rain_step_mm)  .or. rain_step_mm  < 0._RKIND) rain_step_mm  = 0._RKIND
                if (.not. ieee_is_finite(snow_step_mm)  .or. snow_step_mm  < 0._RKIND) snow_step_mm  = 0._RKIND
                if (.not. ieee_is_finite(graup_step_mm) .or. graup_step_mm < 0._RKIND) graup_step_mm = 0._RKIND
                if (.not. ieee_is_finite(ice_step_mm)   .or. ice_step_mm   < 0._RKIND) ice_step_mm   = 0._RKIND

                ! --- Diagnose snow ratio (sr_step) for MPAS outputs ---
                ! WRF P3 defines total surface precipitation as (liquid + solid) and
                ! solid precipitation as prt_sol. Keep qs/qg diagnostics out of MPAS P3.
                block
                   real(kind=RKIND) :: solid_mm, total_mm
                   total_mm = rain_step_mm
                   solid_mm = snow_step_mm
                   if (total_mm > 1.0e-12_RKIND) then
                      sr_step = min(1.0_RKIND, max(0.0_RKIND, solid_mm / total_mm))
                   else
                      sr_step = 0.0_RKIND
                   end if
                end block

                ! Convert updated temperature back to theta (Exner assumed unchanged within dt_microp)
                 do k = kts, kte
                    if (.not. ieee_is_finite(pi1d(k)) .or. pi1d(k) <= 0._RKIND) then
                       if (p3_badcount < 10 .and. p3dbg_bad_pi_n < P3DBG_MAX_BAD_PI) then
                          write(p3_dbgmsg,'(A,3I10,1X,ES12.4,1X,ES12.4)') &
                               'P3DBG_BAD_PI1D_POST (j,i,k,pi,p)=', j,i,k, pi1d(k), p1d(k)
                          call mpas_log_write(trim(p3_dbgmsg))
                          p3_badcount = p3_badcount + 1
                          p3dbg_bad_pi_n = p3dbg_bad_pi_n + 1
                          if (p3dbg_bad_pi_n == P3DBG_MAX_BAD_PI .and. .not. p3dbg_supp_bad_pi_msg) then
                             call mpas_log_write('P3DBG_BAD_PI1D_POST: further occurrences suppressed')
                             p3dbg_supp_bad_pi_msg = .true.
                          end if
                       end if
                       if (ieee_is_finite(p1d(k)) .and. p1d(k) > 1._RKIND) then
                          pi1d(k) = (p1d(k)/P0)**rcp
                       else
                          pi1d(k) = 1._RKIND
                       end if
                    end if
                    th1d(k) = t1d(k) / max(1.0e-12_RKIND, pi1d(k))
                 end do


                ! --- Surface precipitation bookkeeping (mm) ---
                rainnc_p(i,j)  = rainnc_p(i,j)  + rain_step_mm
                rainncv_p(i,j) = rain_step_mm

                snownc_p(i,j)  = snownc_p(i,j)  + snow_step_mm
                snowncv_p(i,j) = snow_step_mm
                sr_p(i,j)      = sr_step

                graupelnc_p(i,j)  = graupelnc_p(i,j)  + graup_step_mm
                graupelncv_p(i,j) = graup_step_mm

                ! Optional: "rolling old" definition (old = previous substep state)
                if (rolling_old) then
                   thold1d(:) = th1d(:)
                   qvold1d(:) = qv1d(:)
                end if

             end do  ! istep

             !------------------------------------------

             !--------------------------------------------------------------------------------
             ! P3 carries a single total-ice category family for MPAS bulk outputs.
             ! Keep qi consistent with the sum of P3 category mass, set qs/qg to zero,
             ! and always recompute ni from the category number fields.
             !--------------------------------------------------------------------------------
             do k = kts, kte
                qi1d(k) = 0.0_RKIND
                ni1d(k) = 0.0_RKIND
                qs1d(k) = 0.0_RKIND
                qg1d(k) = 0.0_RKIND
                do c = 1, p3_nCat
                   if (.not. ieee_is_finite(qit1d(k,c)) .or. qit1d(k,c) < 0._RKIND) qit1d(k,c) = 0._RKIND
                   if (.not. ieee_is_finite(qni1d(k,c)) .or. qni1d(k,c) < 0._RKIND) qni1d(k,c) = 0._RKIND
                   qi1d(k) = qi1d(k) + qit1d(k,c)
                   ni1d(k) = ni1d(k) + qni1d(k,c)
                end do
                if (.not. ieee_is_finite(qi1d(k)) .or. qi1d(k) < 0._RKIND) qi1d(k) = 0._RKIND
                if (.not. ieee_is_finite(ni1d(k)) .or. ni1d(k) < 0._RKIND) ni1d(k) = 0._RKIND
             end do

             ! Unpack final column back to driver arrays ONCE
             !------------------------------------------
             do k = kts, kte
                ! Convert P3-updated temperature (t1d) back to dry potential temperature (theta)
                ! so MPAS state (theta_m/exner/pressure) remains thermodynamically consistent.
                if ((p1d(k) /= p1d(k)) .or. (p1d(k) <= 1._RKIND) .or. (abs(p1d(k)) > 1.0e9_RKIND)) p1d(k) = P0
                 ! (physfix) keep pi1d from dynamics; do not recompute from p1d here
                if ((pi1d(k) /= pi1d(k)) .or. (pi1d(k) < 1.0e-6_RKIND)) pi1d(k) = 1.0e-6_RKIND
                th1d(k) = t1d(k) / pi1d(k)

                th_p(i,k,j)   = th1d(k)
                 ! (physfix) do not overwrite pi_p from microphysics
                 ! (physfix) do not overwrite pres_p from microphysics
                qv_p(i,k,j) = qv1d(k)
                qc_p(i,k,j) = qc1d(k)
                qr_p(i,k,j) = qr1d(k)
                qi_p(i,k,j) = qi1d(k)
                qs_p(i,k,j) = qs1d(k)
                qg_p(i,k,j) = qg1d(k)
                ni_p(i,k,j) = ni1d(k)
                nr_p(i,k,j) = nr1d(k)

                nc_p(i,k,j)   = nc1d(k)
                nwfa_p(i,k,j) = nwfa1d(k)
                nifa_p(i,k,j) = nifa1d(k)

                ! Radiation coupling (WRF-P3 units: meters): map P3 effective radii to MPAS fields
                if (allocated(recloud_p)) recloud_p(i,k,j) = max(1.0e-9_RKIND, effc1d(k))
                if (allocated(reice_p))   reice_p(i,k,j)   = max(1.0e-9_RKIND, effiave1d(k))
                if (allocated(resnow_p))  resnow_p(i,k,j)  = max(1.0e-9_RKIND, effiave1d(k))


             end do

                          do c = 1, p3_nCat
                do k = kts, kte
                   scalars(idx_qit + (c-1), k, i) = qit1d(k,c)
                   scalars(idx_qni + (c-1), k, i) = qni1d(k,c)
                   scalars(idx_qir + (c-1), k, i) = qir1d(k,c)
                   scalars(idx_qib + (c-1), k, i) = qib1d(k,c)
                   if (p3_trplMomI) scalars(idx_qzi + (c-1), k, i) = qzi1d(k,c)
                   if (p3_liqfrac)  scalars(idx_qli + (c-1), k, i) = qli1d(k,c)
                   ! Per-category diags (vmi/dmi/rhoi/dhmax) are intentionally not stored full-domain (memory fix).
                end do
             end do

          end do
       end do

       deallocate(th1d, qv1d, qc1d, qr1d, nr1d)
       deallocate(qv_in, qc_in, qr_in, qi_in, qs_in, qg_in)
       deallocate(pi1d, t1d, qi1d, ni1d, qs1d, qg1d)
       deallocate(thold1d, qvold1d)  ! <<< NEW
       deallocate(p1d, w1d, dz1d, nc1d, nwfa1d, nifa1d, rainprod1d, evapprod1d, &
                  zdbz1d, effc1d, effiave1d)
       deallocate(qit1d, qni1d, qir1d, qib1d, qzi1d, qli1d, vmi1d, dmi1d, rhoi1d, dhmax1d)

       call mpas_timer_stop('mp_p3_aerosols')

    case default
 end select microp_select

!... calculate the 10cm radar reflectivity and relative humidity, if needed:
 if (l_diags) then
if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) call mpas_log_write('P3DBG_DRV: diagnostics block')
   !ensure that we only call compute_radar_reflectivity() if we are using an MPS that supports
   !the computation of simulated radar reflectivity:
   if(trim(microp_scheme) == "mp_wsm6"     .or. &
      trim(microp_scheme) == "mp_thompson" .or. &
      trim(microp_scheme) == "mp_thompson_aerosols" .or. &
      trim(microp_scheme) == "mp_p3_aerosols") then

if (trim(microp_scheme) == "mp_p3_aerosols") then
   if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) then
      call mpas_log_write('P3DBG_DRV: call compute_radar_reflectivity')
   end if
end if
call compute_radar_reflectivity(configs,diag_physics,its,ite)
if (trim(microp_scheme) == "mp_p3_aerosols") then
   if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) then
      call mpas_log_write('P3DBG_DRV: return compute_radar_reflectivity')
   end if
end if

   else
      call mpas_log_write('*** NOTICE: NOT computing simulated radar reflectivity')
      call mpas_log_write('            since WSM6/Thompson/P3 microphysics scheme was not selected')
   end if


if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) call mpas_log_write('P3DBG_DRV: call compute_relhum')
call compute_relhum(diag,its,ite)
if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) call mpas_log_write('P3DBG_DRV: return compute_relhum')
 end if

!... copy updated precipitation from the wrf-physics grid back to the geodesic-dynamics grid:

if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) call mpas_log_write('P3DBG_DRV: call precip_to_MPAS')
call precip_to_MPAS(configs,diag_physics,its,ite)
if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) call mpas_log_write('P3DBG_DRV: return precip_to_MPAS')
!... copy updated cloud microphysics variables from the wrf-physics grid back to the geodesic-
!    dynamics grid:
if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) then
   call mpas_log_write('P3DBG_DRV: call microphysics_to_MPAS')
end if
call microphysics_to_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,tend,its,ite)
if (itimestep == 1 .and. associated(p3_debug_on_ptr) .and. p3_debug_on_ptr) then
   call mpas_log_write('P3DBG_DRV: return microphysics_to_MPAS')
end if
!... deallocation of all microphysics arrays:
!$OMP BARRIER
!$OMP MASTER
 call deallocate_microphysics(configs)
!$OMP END MASTER

!call mpas_log_write('---enter subroutine driver_microphysics:')
!call mpas_log_write('')

 end subroutine driver_microphysics

!=================================================================================================================
 subroutine precip_from_MPAS(configs,diag_physics,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!output variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,pointer:: nCellsSolve
 real,dimension(:),pointer:: graupelncv,rainncv,snowncv,sr

!local variables and arrays:
 integer:: i,j

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )

!variables common to all cloud microphysics schemes:
 do j = jts, jte
 do i = its, ite
    rainncv_p(i,j) = 0._RKIND
    rainnc_p(i,j)  = 0._RKIND
 enddo
 enddo

 do i = its,ite
    rainncv(i) = 0._RKIND
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6","mp_p3_aerosols")
       do j = jts, jte
       do i = its, ite
          snowncv_p(i,j)    = 0._RKIND
          graupelncv_p(i,j) = 0._RKIND
          snownc_p(i,j)     = 0._RKIND
          graupelnc_p(i,j)  = 0._RKIND
          sr_p(i,j)         = 0._RKIND
       enddo
       enddo

       do i = its,ite
          snowncv(i)    = 0._RKIND
          graupelncv(i) = 0._RKIND
          sr(i)         = 0._RKIND
       enddo

    case default
 end select microp_select

 end subroutine precip_from_MPAS

!=================================================================================================================
 subroutine precip_to_MPAS(configs,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,dimension(:),pointer:: i_rainnc

 real(kind=RKIND),pointer:: config_bucket_rainnc
 real(kind=RKIND),dimension(:),pointer:: precipw
 real(kind=RKIND),dimension(:),pointer:: graupelnc,rainnc,snownc
 real(kind=RKIND),dimension(:),pointer:: graupelncv,rainncv,snowncv,sr

!local variables and arrays:
 integer:: i,j,k
 real(kind=RKIND):: rho_a

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme       )
 call mpas_pool_get_config(configs,'config_bucket_rainnc',config_bucket_rainnc)

 call mpas_pool_get_array(diag_physics,'i_rainnc'  ,i_rainnc  )
 call mpas_pool_get_array(diag_physics,'precipw'   ,precipw   )
 call mpas_pool_get_array(diag_physics,'graupelnc' ,graupelnc )
 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
 call mpas_pool_get_array(diag_physics,'rainnc'    ,rainnc    )
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'snownc'    ,snownc    )
 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )

 do i = its,ite
    precipw(i) = 0._RKIND
 enddo

!variables common to all cloud microphysics schemes:
 do j = jts,jte
 do i = its,ite

    !precipitable water:
    do k = kts,kte
       rho_a = rho_p(i,k,j) / (1._RKIND + qv_p(i,k,j))
       precipw(i) = precipw(i) + qv_p(i,k,j) * rho_a * dz_p(i,k,j)
    enddo

    !time-step precipitation:
    rainncv(i) = rainnc_p(i,j)

    !accumulated precipitation:
    rainnc(i) = rainnc(i) + rainncv(i)

   if(l_acrain .and. config_bucket_rainnc.gt.0._RKIND .and. &
      rainnc(i).gt.config_bucket_rainnc) then
      i_rainnc(i) = i_rainnc(i) + 1
      rainnc(i)   = rainnc(i) - config_bucket_rainnc
   endif

 enddo
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_wsm6")
       do j = jts,jte
       do i = its,ite
          !time-step precipitation:
          snowncv(i)    = snownc_p(i,j)
          graupelncv(i) = graupelnc_p(i,j)
          sr(i) = (snownc_p(i,j) + graupelnc_p(i,j)) / (rainnc_p(i,j)+1.e-12)

          !accumulated precipitation:
          snownc(i)    = snownc(i) + snowncv(i)
          graupelnc(i) = graupelnc(i) + graupelncv(i)
       enddo
       enddo

    case ("mp_p3_aerosols")
       do j = jts,jte
       do i = its,ite
          !time-step precipitation:
          snowncv(i)    = snownc_p(i,j)
          graupelncv(i) = graupelnc_p(i,j)   ! expected to be 0 for your current P3 wiring
          sr(i)         = sr_p(i,j)          ! P3 wrapper provides sr_step -> sr_p

          !accumulated precipitation:
          snownc(i)    = snownc(i) + snowncv(i)
          graupelnc(i) = graupelnc(i) + graupelncv(i)
       enddo
       enddo

    case default
 end select microp_select

 end subroutine precip_to_MPAS


!=================================================================================================================
 subroutine compute_radar_reflectivity(configs,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 real(kind=RKIND),dimension(:,:),pointer:: refl10cm
 real(kind=RKIND),dimension(:),pointer:: refl10cm_max,refl10cm_1km,refl10cm_1km_max

!local variables and arrays:
 integer:: i,j,k,kp
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qc1d,qr1d,qs1d,qg1d,t1d,p1d,nr1d,dBZ1d,zp
 real(kind=RKIND):: w1,w2

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call mpas_pool_get_array(diag_physics,'refl10cm',refl10cm)
 call mpas_pool_get_array(diag_physics,'refl10cm_max',refl10cm_max)
 call mpas_pool_get_array(diag_physics,'refl10cm_1km',refl10cm_1km)
 call mpas_pool_get_array(diag_physics,'refl10cm_1km_max',refl10cm_1km_max)

 microp_select: select case(trim(microp_scheme))
    case ("mp_kessler")
       call physics_error_fatal('--- calculation of radar reflectivity is not available' // &
                                 'with kessler cloud microphysics')

    case ("mp_wsm6")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j) + 0.5*dz_p(i,k,j) ! height AGL
          enddo

          call refl10cm_wsm6(qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             refl10cm(k,i) = dBZ1d(k)
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

    case ("mp_thompson","mp_thompson_aerosols")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qc1d) ) allocate(qc1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(nr1d) ) allocate(nr1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qc1d(k)  = qc_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             nr1d(k)  = nr_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j) + 0.5*dz_p(i,k,j) ! height AGL
          enddo

          call calc_refl10cm(qv1d,qc1d,qr1d,nr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             refl10cm(k,i) = dBZ1d(k)
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qc1d) ) deallocate(qc1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(nr1d) ) deallocate(nr1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

case ("mp_p3_aerosols")

   if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
   if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
   if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
   if(.not.allocated(qc1d) ) allocate(qc1d(kts:kte) )
   if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
   if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
   if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
   if(.not.allocated(nr1d) ) allocate(nr1d(kts:kte) )
   if(.not.allocated(dBZ1d)) allocate(dBZ1d(kts:kte))
   if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

   do j = jts,jte
   do i = its,ite

      do k = kts,kte
         p1d(k)  = pres_p(i,k,j)
         t1d(k)  = th_p(i,k,j) * pi_p(i,k,j)
         qv1d(k) = qv_p(i,k,j)
         qc1d(k) = qc_p(i,k,j)
         qr1d(k) = qr_p(i,k,j)
         qs1d(k) = qs_p(i,k,j)
         qg1d(k) = qg_p(i,k,j)
         nr1d(k) = nr_p(i,k,j)
         zp(k)   = z_p(i,k,j) - z_p(i,1,j) + 0.5_RKIND*dz_p(i,k,j)
      enddo

      call calc_refl10cm(qv1d,qc1d,qr1d,nr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

      kp = 1
      do k = kts,kte
         dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
         refl10cm(k,i) = dBZ1d(k)
         if(zp(k) .lt. 1000._RKIND) kp = k
      enddo

      refl10cm_max(i) = maxval(dBZ1d(:))
      kp = min(max(kp, kts), kte-1)
      w1 = (zp(kp+1)-1000._RKIND)/(zp(kp+1)-zp(kp))
      w2 = 1.0_RKIND - w1
      refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
      refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))

   enddo
   enddo

   if(allocated(p1d)  ) deallocate(p1d  )
   if(allocated(t1d)  ) deallocate(t1d  )
   if(allocated(qv1d) ) deallocate(qv1d )
   if(allocated(qc1d) ) deallocate(qc1d )
   if(allocated(qr1d) ) deallocate(qr1d )
   if(allocated(qs1d) ) deallocate(qs1d )
   if(allocated(qg1d) ) deallocate(qg1d )
   if(allocated(nr1d) ) deallocate(nr1d )
   if(allocated(dBZ1d)) deallocate(dBZ1d)
   if(allocated(zp)   ) deallocate(zp   )


    case default
 end select microp_select

 end subroutine compute_radar_reflectivity

!=================================================================================================================
 subroutine compute_relhum(diag,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag

!local pointers:
 real(kind=RKIND),dimension(:,:),pointer:: relhum

!local variables and arrays:
 integer:: i,j,k

 real(kind=RKIND):: tempc, rh
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qvs1d,t1d,p1d

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag,'relhum',relhum)
 do k = kts,kte
    do i = its,ite
       relhum(k,i) = 0._RKIND
    enddo
 enddo

 if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
 if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
 if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
 if(.not.allocated(qvs1d)) allocate(qvs1d(kts:kte))

 do j = jts,jte
 do i = its,ite
    do k = kts,kte
       p1d(k)  = pres_p(i,k,j)
       t1d(k)  = th_p(i,k,j) * pi_p(i,k,j)
       tempc   = t1d(k) - 273.16_RKIND
       qvs1d(k) = rslf(p1d(k),t1d(k))
       if(tempc .le. 0._RKIND) qvs1d(k) = rsif(p1d(k),t1d(k))
       qv1d(k) = qv_p(i,k,j)
       relhum(k,i) = qv1d(k) / qvs1d(k)
       relhum(k,i) = relhum(k,i) * 100._RKIND
    enddo
 enddo
 enddo

 if(allocated(p1d)  ) deallocate(p1d  )
 if(allocated(t1d)  ) deallocate(t1d  )
 if(allocated(qv1d) ) deallocate(qv1d )
 if(allocated(qvs1d)) deallocate(qvs1d)

 end subroutine compute_relhum

!=================================================================================================================
 end module mpas_atmphys_driver_microphysics
!=================================================================================================================
