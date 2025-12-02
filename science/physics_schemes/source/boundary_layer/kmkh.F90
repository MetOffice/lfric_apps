! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To set the turbulent mixing coefficients KM and KH
!           (Note: should be used after any vertical interpolation
!                  but before any horizontal interpolation.)

!  Programming standard: UMDP 3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module kmkh_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'KMKH_MOD'
contains

subroutine kmkh (                                                              &
! in data
 bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                       &
 ntml,cumulus,ntdsc,sml_disc_inv,dsc_disc_inv,                                 &
 weight_1dbl,weight_1dbl_rho,                                                  &
! INOUT data
 rhokm,rhokh,rhokmz,rhokhz,rhokm_top,rhokh_top,tke_loc                         &
 )

use atm_fields_bounds_mod, only: pdims,pdims_s,tdims,ScmRowLen, ScmRow,        &
     tdims
use bl_option_mod, only:                                                       &
    Keep_Ri_FA, off, on, kprof_cu, except_disc_inv, zero, one
use bl_diags_mod, only: strnewbldiag
use cv_run_mod, only: l_param_conv
use model_domain_mod, only: model_type, mt_single_column
use s_scmop_mod,      only: default_streams,                                   &
                            t_avg, d_bl, scmdiag_bl

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! in arguments
integer, intent(in) ::                                                         &
 bl_levels

logical, intent(in) ::                                                         &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! in flag for Cu in the bl

! Additional variables for SCM diagnostics which are dummy in full UM
integer,intent(in) ::                                                          &
 nSCMDpkgs              ! No of SCM diagnostics packages

logical,intent(in) ::                                                          &
 L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

integer, intent(in) ::                                                         &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                    &
                              ! in Number of model levels in the
                              !    turbulently mixed layer.
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                   &
                              ! in Top level for turb mixing in
                              !    cloud layer
 sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! in Flags for whether discontinuous
 dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! in inversions are diagnosed

real(kind=r_bl), intent(in) ::                                                 &
 weight_1dbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
 weight_1dbl_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
                 bl_levels)
                              ! in Weighting applied to 1D BL scheme,
                              !    to blend with Smagorinsky scheme,
                              !    on theta and rho levels
! INOUT arguments
real(kind=r_bl), intent(in out) ::                                             &
 rhokmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                              ! INOUT Non-local turbulent mixing
                              !    coefficient for momentum.
 rhokhz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                              ! INOUT Non-local turbulent mixing
                              !    coefficient for heat and moisture
 rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           2:bl_levels),                                                       &
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for momentum.
 rhokh_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           2:bl_levels),                                                       &
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for heat
                              !    and moisture.
 rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,            &
       bl_levels),                                                             &
                              ! INOUT Layer K-1 - to - layer K
                              !       turbulent mixing coefficient
                              !       for momentum.
 rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                              ! INOUT Layer K-1 - to - layer K
                              !       turbulent mixing coefficient
                              !       for heat and moisture.
 tke_loc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,2:bl_levels)
                              ! INOUT Ri-based scheme diagnosed TKE
!----------------------------------------------------------------------
!  Define local storage.

character(len=*), parameter ::  RoutineName = 'KMKH'
! Scm diags
real(kind=r_bl) :: rhokh_diag(ScmRowLen,ScmRow,bl_levels)
real(kind=r_bl) :: rhokm_diag(ScmRowLen,ScmRow,bl_levels)
                           ! diffusivity of heat and momentum kg/(ms)

integer ::                                                                     &
 i,j,iScm,jScm,                                                                &
                     ! Loop counter (horizontal field index).
 k             ! Loop counter (vertical level index).

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (model_type == mt_single_column) then
  do k=1, bl_levels
    do j=pdims%j_start, pdims%j_end
      jScm = j - pdims%j_start + 1
      do i=pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1
        rhokh_diag(iScm,jScm,k) = rhokh(i,j,k)
        rhokm_diag(iScm,jScm,k) = rhokm(i,j,k)
      end do ! i
    end do ! j
  end do ! k
end if ! model_type

!$OMP PARALLEL DEFAULT(SHARED) private(k,j,i)
if (Keep_Ri_FA == on) then
  !-----------------------------------------------------------------------
  ! Set local K's to zero at the LCL in cumulus and at the
  ! top of a turbulent layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( ( ( cumulus(i,j) .or. sml_disc_inv(i,j) >= 1) .and.               &
             (k == ntml(i,j)+1 .or. k == ntml(i,j)+2) ) .or.                   &

             ( dsc_disc_inv(i,j)  >=  1 .and.                                  &
             (k == ntdsc(i,j)+1 .or. k == ntdsc(i,j)+2) ) ) then
          rhokh(i,j,k) = zero
          rhokm(i,j,k) = zero
          tke_loc(i,j,k) = zero
        end if

      end do ! P_POINTS,i
    end do ! P_POINTS,j
  end do ! BL_LEVELS
!$OMP end do NOWAIT

else if (Keep_Ri_FA == except_disc_inv) then
  !-----------------------------------------------------------------------
          ! Reduce local K's only at the top of a turbulent
          ! layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( ( sml_disc_inv(i,j)  >=  1 .and.                                  &
             (k == ntml(i,j)+1 .or. k == ntml(i,j)+2) ) .or.                   &

             ( dsc_disc_inv(i,j)  >=  1 .and.                                  &
             (k == ntdsc(i,j)+1 .or. k == ntdsc(i,j)+2) ) ) then
          rhokh(i,j,k) = (one-weight_1dbl(i,j,k))*rhokh(i,j,k)
          rhokm(i,j,k) = (one-weight_1dbl(i,j,k))*rhokm(i,j,k)
        end if

      end do ! P_POINTS,i
    end do ! P_POINTS,j
  end do ! BL_LEVELS
!$OMP end do NOWAIT

else
  !-----------------------------------------------------------------------
  ! Set local K's to zero from the LCL in cumulus and from the
  ! top of a turbulent layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        if ( (cumulus(i,j) .and. ( (l_param_conv .and. k >  ntml(i,j))         &
             .or. (.not. l_param_conv .and. k >= ntml(i,j)) )) .or.            &

             ( dsc_disc_inv(i,j)  >=  1 .and. k  >   ntdsc(i,j) ) .or.         &

             ( sml_disc_inv(i,j)  >=  1 .and. k  >   ntml(i,j) ) ) then
          !   This also means no local mixing within any DSC layer
          rhokh(i,j,k) = zero
          rhokm(i,j,k) = zero
          tke_loc(i,j,k) = zero
        end if

      end do ! P_POINTS,i
    end do ! P_POINTS,j
  end do ! BL_LEVELS
!$OMP end do NOWAIT

end if ! test on Keep_Ri_FA

if (kprof_cu == off) then
  !-----------------------------------------------------------------------
  ! Set non-local K's to zero at the LCL in cumulus layers,
  ! including level NTML if not l_param_conv convection scheme
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end

        if (cumulus(i,j) .and. ( (l_param_conv .and. k == ntml(i,j)+1)         &
             .or. (.not. l_param_conv .and.                                    &
                         k >= ntml(i,j) .and. k <  ntml(i,j)+2) )) then
          rhokhz(i,j,k)=zero
          rhokmz(i,j,k)=zero
          rhokh_top(i,j,k)=zero
          rhokm_top(i,j,k)=zero
        end if
      end do ! P_POINTS,i
    end do ! P_POINTS,j
  end do ! BL_LEVELS
!$OMP end do NOWAIT

end if  ! test on kprof_cu

! Need a barrier to ensure all previous possible loops have completed
!$OMP BARRIER
!-----------------------------------------------------------------------
! Save diffusion coefficients for diagnostics
!-----------------------------------------------------------------------
if (BL_diag%l_rhokmloc) then
!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokmloc(i,j,k)=rhokm(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokhloc) then
!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokhloc(i,j,k)=rhokh(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokmsurf) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokmsurf(i,j,k)=weight_1dbl(i,j,k)*rhokmz(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokhsurf) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokhsurf(i,j,k)=weight_1dbl_rho(i,j,k)*rhokhz(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokmsc) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokmsc(i,j,k)=weight_1dbl(i,j,k)*rhokm_top(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokhsc) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        BL_diag%rhokhsc(i,j,k)=weight_1dbl_rho(i,j,k)*rhokh_top(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
end if
!-----------------------------------------------------------------------
! Set KM and KH to be the maximum of the local and non-local
! values andstore RHO_KM on P-grid for output.
!-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end

      rhokh(i,j,k) = max( rhokh(i,j,k) ,                                       &
             weight_1dbl_rho(i,j,k)*(rhokhz(i,j,k)+rhokh_top(i,j,k)) )
      rhokm(i,j,k) = max( rhokm(i,j,k) ,                                       &
             weight_1dbl(i,j,k)*(rhokmz(i,j,k)+rhokm_top(i,j,k)) )

    end do ! P_POINTS,i
  end do ! P_POINTS,j
end do ! BL_LEVELS
!$OMP end do

!$OMP end PARALLEL

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine kmkh
end module kmkh_mod
