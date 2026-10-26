   function atm_define_packages(packagePool) result(iErr)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      implicit none
      type (mpas_pool_type), intent(inout) :: packagePool !< Input: MPAS Pool for containing package logicals.

      integer :: iErr

      iErr = 0
      call mpas_pool_add_package(packagePool, 'mp_kessler_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_thompson_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_thompson_aers_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_p3_aers_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_p3_aers_cat2_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_p3_aers_cat3_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'mp_wsm6_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'cu_grell_freitas_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'cu_kain_fritsch_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'cu_ntiedtke_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'bl_ysu_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'bl_mynn_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'sf_noahmp_inActive', .false.)
      call mpas_pool_add_package(packagePool, 'iauActive', .false.)
      call mpas_pool_add_package(packagePool, 'limited_areaActive', .false.)
      call mpas_pool_add_package(packagePool, 'jedi_daActive', .false.)
      call mpas_pool_add_package(packagePool, 'no_invariant_streamActive', .false.)
      call mpas_pool_add_package(packagePool, 'ugwp_orog_streamActive', .false.)
      call mpas_pool_add_package(packagePool, 'ugwp_ngw_streamActive', .false.)
      call mpas_pool_add_package(packagePool, 'ugwp_diags_streamActive', .false.)
   end function atm_define_packages
