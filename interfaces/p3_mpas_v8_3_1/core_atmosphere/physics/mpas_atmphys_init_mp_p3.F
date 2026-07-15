! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!--------------------------------------------------------------------------------------------------
! MPAS microphysics initialization: P3 aerosol-aware (always-on) support.
!
! File name: mpas_atmphys_init_mp_p3.F
!
! This module provides:
!   1) init_p3_forMPAS          : initialize P3 lookup tables/constants (calls p3_init)
!   2) init_p3_aerosols_forMPAS : read aerosol climatology from netCDF and populate nwfa/nifa state
!
! Assumptions (by design, matching the request):
!   - P3 is ALWAYS aerosol-aware in MPAS.
!   - Aerosol climatology file path is hard-coded below.
!   - Climatology is expected on the MPAS mesh with dimensions matching nVertLevels and nCellsSolve.
!     Accepted variable shapes:
!       NWFA(nVertLevels,nCellsSolve) or NWFA(nCellsSolve,nVertLevels)
!       NIFA(nVertLevels,nCellsSolve) or NIFA(nCellsSolve,nVertLevels)
!     Optional 3rd (time/month) dimension is allowed; values are averaged over that dimension.
!--------------------------------------------------------------------------------------------------

#define DM_BCAST_R4(A) call mpas_dmpar_bcast_real4s(dminfo, size(A), A)

module mpas_atmphys_init_mp_p3

   use mpas_dmpar
   use mpas_kind_types
   use mpas_log
   use mpas_pool_routines
   use mpas_atmphys_utilities

   use netcdf
   use module_mp_p3, only : p3_init_for_mpas

   implicit none

   private
   public :: init_p3_forMPAS, init_p3_aerosols_forMPAS

   !------------------------------------------------------------------------------------------------
   ! Hard-coded input paths (MPAS-compatible).
   ! You can later replace these with config options if desired.
   !------------------------------------------------------------------------------------------------
   character(len=*), parameter :: P3_LOOKUP_TABLE_DIR = 'P3_LOOKUP_TABLES'
   character(len=*), parameter :: P3_AEROSOL_CLIMO_NC = 'P3_AEROSOL_CLIMO.nc'

contains

!==================================================================================================
subroutine init_p3_forMPAS(dminfo, mesh, state, time_lev, diag_physics)
!==================================================================================================
! Initialize P3 constants/lookup tables.
!
! Notes:
! - p3_init internally guards against multiple initialization (is_init SAVE).
! - We DO call this even in restart runs, because lookup tables/constants are not stored in MPAS
!   restart files and must be present in memory.
!
   type(dm_info),        intent(in)    :: dminfo
   type(mpas_pool_type), intent(in)    :: mesh
   type(mpas_pool_type), intent(inout) :: state
   type(mpas_pool_type), intent(inout) :: diag_physics
   integer,              intent(in)    :: time_lev

   integer :: stat
   logical :: owr

   ! ---- P3 configuration choices (hard-coded defaults; adjust as needed) ----
   integer, parameter :: nCat_default = 1
   logical, parameter :: trplMomI_default = .false.
   logical, parameter :: liqfrac_default  = .false.

   owr = (dminfo%my_proc_id == IO_NODE)

   call mpas_log_write('P3 init: calling p3_init_for_mpas (lookup dir="'//trim(P3_LOOKUP_TABLE_DIR)//'")')

   call p3_init_for_mpas( lookup_file_dir   = P3_LOOKUP_TABLE_DIR, &
                          nCat              = nCat_default,        &
                          trplMomI          = trplMomI_default,    &
                          liqfrac           = liqfrac_default,     &
                          stat              = stat,                &
                          abort_on_err      = .true.,              &
                          model             = 'WRF',               &
                          dowr              = owr )

   if (stat /= 0) then
      call physics_error_fatal('init_p3_forMPAS: p3_init returned error status.')
   end if

end subroutine init_p3_forMPAS


!==================================================================================================
subroutine init_p3_aerosols_forMPAS(do_restart, dminfo, mesh, state, time_lev)
!==================================================================================================
! Read aerosol climatology from netCDF and populate MPAS state scalars:
!   nwfa (hygroscopic / "water friendly")
!   nifa (ice nucleating / "ice friendly")
!
! Restart behavior:
!   - If do_restart=.true., return immediately (exact restartibility).
!
   logical,              intent(in)    :: do_restart
   type(dm_info),        intent(in)    :: dminfo
   type(mpas_pool_type), intent(in)    :: mesh
   type(mpas_pool_type), intent(inout) :: state
   integer,              intent(in)    :: time_lev

   integer, pointer :: nCellsSolve, nVertLevels
   integer, pointer :: index_nwfa, index_nifa

   real(kind=RKIND), dimension(:,:,:), pointer :: scalars
   real(kind=RKIND), dimension(:,:),   pointer :: nwfa, nifa

   real(kind=R4KIND), allocatable :: nwfa_in(:,:), nifa_in(:,:)   ! [k, iCell]

   if (do_restart) then
      call mpas_log_write('P3 aerosols: do_restart=.true.; skipping aerosol climatology init.')
      return
   end if

   ! Dimensions
   call mpas_pool_get_dimension(mesh,  'nCellsSolve',  nCellsSolve)
   call mpas_pool_get_dimension(mesh,  'nVertLevels',  nVertLevels)

   ! Scalar indices (must exist if P3 aerosol-aware is enabled)
   call mpas_pool_get_dimension(state, 'index_nwfa', index_nwfa)
   call mpas_pool_get_dimension(state, 'index_nifa', index_nifa)

   if (index_nwfa < 1 .or. index_nifa < 1) then
      call physics_error_fatal('init_p3_aerosols_forMPAS: index_nwfa/index_nifa not found in state%scalars.')
   end if

   ! Grab state scalar storage at this time_lev
   call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
   nwfa => scalars(index_nwfa, :, :)
   nifa => scalars(index_nifa, :, :)

   allocate(nwfa_in(nVertLevels, nCellsSolve))
   allocate(nifa_in(nVertLevels, nCellsSolve))
   nwfa_in = 0.0_R4KIND
   nifa_in = 0.0_R4KIND

   if (dminfo%my_proc_id == IO_NODE) then
      call mpas_log_write('P3 aerosols: reading netCDF climatology "'//trim(P3_AEROSOL_CLIMO_NC)//'"')
      call read_p3_aerosols_netcdf(P3_AEROSOL_CLIMO_NC, nVertLevels, nCellsSolve, nwfa_in, nifa_in)
   end if

   ! Broadcast to all tasks
   DM_BCAST_R4(nwfa_in)
   DM_BCAST_R4(nifa_in)

   ! Populate state scalars (convert real4 -> RKIND)
   nwfa(:,1:nCellsSolve) = real(nwfa_in(:,1:nCellsSolve), kind=RKIND)
   nifa(:,1:nCellsSolve) = real(nifa_in(:,1:nCellsSolve), kind=RKIND)

   call mpas_log_write('P3 aerosols: nwfa/nifa initialized from climatology.')

   deallocate(nwfa_in, nifa_in)

end subroutine init_p3_aerosols_forMPAS


!==================================================================================================
subroutine read_p3_aerosols_netcdf(fname, nVertLevels, nCellsSolve, nwfa_out, nifa_out)
!==================================================================================================
! IO-node helper: reads NWFA and NIFA from netCDF and writes into [k,iCell] arrays.
! Allows 2D or 3D (time/month) variables; 3D is averaged over the 3rd dimension.
!
   character(len=*),  intent(in)  :: fname
   integer,           intent(in)  :: nVertLevels, nCellsSolve
   real(kind=R4KIND), intent(out) :: nwfa_out(nVertLevels, nCellsSolve)
   real(kind=R4KIND), intent(out) :: nifa_out(nVertLevels, nCellsSolve)

   integer :: ncid
   integer :: varid_nwfa, varid_nifa
   integer :: ndims, dimids(NF90_MAX_VAR_DIMS)
   integer :: dimlen(3)
   integer :: i, ierr
   integer :: d1, d2, d3

   real(kind=R4KIND), allocatable :: tmp2(:,:)
   real(kind=R4KIND), allocatable :: tmp3(:,:,:)

   nwfa_out = 0.0_R4KIND
   nifa_out = 0.0_R4KIND

   ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
   call nc_check(ierr, 'open '//trim(fname))

   ! ---------- NWFA ----------
   ierr = nf90_inq_varid(ncid, 'NWFA', varid_nwfa)
   call nc_check(ierr, 'inq_varid NWFA')

   ierr = nf90_inquire_variable(ncid, varid_nwfa, ndims=ndims, dimids=dimids)
   call nc_check(ierr, 'inquire_variable NWFA')

   if (ndims < 2 .or. ndims > 3) then
      call physics_error_fatal('read_p3_aerosols_netcdf: NWFA must have 2 or 3 dimensions.')
   end if

   do i = 1, ndims
      ierr = nf90_inquire_dimension(ncid, dimids(i), len=dimlen(i))
      call nc_check(ierr, 'inquire_dimension NWFA dim')
   end do
   d1 = dimlen(1); d2 = dimlen(2)
   d3 = 1
   if (ndims == 3) d3 = dimlen(3)

   if (ndims == 2) then
      allocate(tmp2(d1,d2))
      ierr = nf90_get_var(ncid, varid_nwfa, tmp2)
      call nc_check(ierr, 'get_var NWFA (2D)')
      call pack_to_k_cell(tmp2, nVertLevels, nCellsSolve, nwfa_out)
      deallocate(tmp2)
   else
      allocate(tmp3(d1,d2,d3))
      ierr = nf90_get_var(ncid, varid_nwfa, tmp3)
      call nc_check(ierr, 'get_var NWFA (3D)')
      allocate(tmp2(d1,d2))
      tmp2 = sum(tmp3, dim=3) / real(d3, kind=R4KIND)
      call pack_to_k_cell(tmp2, nVertLevels, nCellsSolve, nwfa_out)
      deallocate(tmp2, tmp3)
   end if

   ! ---------- NIFA ----------
   ierr = nf90_inq_varid(ncid, 'NIFA', varid_nifa)
   call nc_check(ierr, 'inq_varid NIFA')

   ierr = nf90_inquire_variable(ncid, varid_nifa, ndims=ndims, dimids=dimids)
   call nc_check(ierr, 'inquire_variable NIFA')

   if (ndims < 2 .or. ndims > 3) then
      call physics_error_fatal('read_p3_aerosols_netcdf: NIFA must have 2 or 3 dimensions.')
   end if

   do i = 1, ndims
      ierr = nf90_inquire_dimension(ncid, dimids(i), len=dimlen(i))
      call nc_check(ierr, 'inquire_dimension NIFA dim')
   end do
   d1 = dimlen(1); d2 = dimlen(2)
   d3 = 1
   if (ndims == 3) d3 = dimlen(3)

   if (ndims == 2) then
      allocate(tmp2(d1,d2))
      ierr = nf90_get_var(ncid, varid_nifa, tmp2)
      call nc_check(ierr, 'get_var NIFA (2D)')
      call pack_to_k_cell(tmp2, nVertLevels, nCellsSolve, nifa_out)
      deallocate(tmp2)
   else
      allocate(tmp3(d1,d2,d3))
      ierr = nf90_get_var(ncid, varid_nifa, tmp3)
      call nc_check(ierr, 'get_var NIFA (3D)')
      allocate(tmp2(d1,d2))
      tmp2 = sum(tmp3, dim=3) / real(d3, kind=R4KIND)
      call pack_to_k_cell(tmp2, nVertLevels, nCellsSolve, nifa_out)
      deallocate(tmp2, tmp3)
   end if

   ierr = nf90_close(ncid)
   call nc_check(ierr, 'close '//trim(fname))

end subroutine read_p3_aerosols_netcdf


!==================================================================================================
subroutine pack_to_k_cell(src, nVertLevels, nCellsSolve, dst)
!==================================================================================================
! Map a 2D netCDF variable into MPAS [k,iCell] layout.
! Accepts src shaped either [nVertLevels,nCellsSolve] or [nCellsSolve,nVertLevels].
!
   integer,           intent(in)  :: nVertLevels, nCellsSolve
   real(kind=R4KIND), intent(in)  :: src(:,:)
   real(kind=R4KIND), intent(out) :: dst(nVertLevels, nCellsSolve)

   integer :: d1, d2

   d1 = size(src,1)
   d2 = size(src,2)

   if (d1 == nVertLevels .and. d2 == nCellsSolve) then
      dst(:,:) = src(:,:)
   else if (d1 == nCellsSolve .and. d2 == nVertLevels) then
      dst(:,:) = transpose(src(:,:))
   else
      call physics_error_fatal('pack_to_k_cell: netCDF variable dimensions do not match nVertLevels/nCellsSolve.')
   end if

end subroutine pack_to_k_cell


!==================================================================================================
subroutine nc_check(ierr, where)
!==================================================================================================
   integer,          intent(in) :: ierr
   character(len=*), intent(in) :: where
   character(len=StrKIND) :: msg

   if (ierr /= nf90_noerr) then
      write(msg,'(A,": ",A)') trim(where), trim(nf90_strerror(ierr))
      call physics_error_fatal(msg)
   end if
end subroutine nc_check


end module mpas_atmphys_init_mp_p3
