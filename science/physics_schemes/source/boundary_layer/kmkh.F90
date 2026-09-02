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

use atm_fields_bounds_mod, only: pdims,pdims_s,tdims,ScmRowLen,                &
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
 cumulus(pdims%i_start:pdims%i_end)
                              ! in flag for Cu in the bl

! Additional variables for SCM diagnostics which are dummy in full UM
integer,intent(in) ::                                                          &
 nSCMDpkgs              ! No of SCM diagnostics packages

logical,intent(in) ::                                                          &
 L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

integer, intent(in) ::                                                         &
 ntml(pdims%i_start:pdims%i_end),                                              &
                              ! in Number of model levels in the
                              !    turbulently mixed layer.
 ntdsc(pdims%i_start:pdims%i_end),                                             &
                              ! in Top level for turb mixing in
                              !    cloud layer
 sml_disc_inv(pdims%i_start:pdims%i_end),                                      &
                              ! in Flags for whether discontinuous
 dsc_disc_inv(pdims%i_start:pdims%i_end)
                              ! in inversions are diagnosed

real(kind=r_bl), intent(in) ::                                                 &
 weight_1dbl(pdims%i_start:pdims%i_end,bl_levels),                             &
 weight_1dbl_rho(pdims%i_start:pdims%i_end,                                    &
                 bl_levels)
                              ! in Weighting applied to 1D BL scheme,
                              !    to blend with Smagorinsky scheme,
                              !    on theta and rho levels
! INOUT arguments
real(kind=r_bl), intent(in out) ::                                             &
 rhokmz(tdims%i_start:tdims%i_end,                                             &
        2:bl_levels),                                                          &
                              ! INOUT Non-local turbulent mixing
                              !    coefficient for momentum.
 rhokhz(tdims%i_start:tdims%i_end,                                             &
        2:bl_levels),                                                          &
                              ! INOUT Non-local turbulent mixing
                              !    coefficient for heat and moisture
 rhokm_top(tdims%i_start:tdims%i_end,                                          &
           2:bl_levels),                                                       &
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for momentum.
 rhokh_top(tdims%i_start:tdims%i_end,                                          &
           2:bl_levels),                                                       &
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for heat
                              !    and moisture.
 rhokm(pdims_s%i_start:pdims_s%i_end,                                          &
       bl_levels),                                                             &
                              ! INOUT Layer K-1 - to - layer K
                              !       turbulent mixing coefficient
                              !       for momentum.
 rhokh(tdims%i_start:tdims%i_end,bl_levels),                                   &
                              ! INOUT Layer K-1 - to - layer K
                              !       turbulent mixing coefficient
                              !       for heat and moisture.
 tke_loc(pdims%i_start:pdims%i_end,2:bl_levels)
                              ! INOUT Ri-based scheme diagnosed TKE
!----------------------------------------------------------------------
!  Define local storage.

character(len=*), parameter ::  RoutineName = 'KMKH'
! Scm diags
real(kind=r_bl) :: rhokh_diag(ScmRowLen,bl_levels)
real(kind=r_bl) :: rhokm_diag(ScmRowLen,bl_levels)
                           ! diffusivity of heat and momentum kg/(ms)

integer ::                                                                     &
 i,iScm,                                                                       &
                     ! Loop counter (horizontal field index).
 k             ! Loop counter (vertical level index).

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (model_type == mt_single_column) then
  do k=1, bl_levels
    do i=pdims%i_start, pdims%i_end
      iScm = i - pdims%i_start + 1
      rhokh_diag(iScm,k) = rhokh(i,k)
      rhokm_diag(iScm,k) = rhokm(i,k)
    end do ! i
  end do ! k
end if ! model_type

!$OMP PARALLEL DEFAULT(SHARED) private(k,i)
if (Keep_Ri_FA == on) then
  !-----------------------------------------------------------------------
  ! Set local K's to zero at the LCL in cumulus and at the
  ! top of a turbulent layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      if ( ( ( cumulus(i) .or. sml_disc_inv(i) >= 1) .and.                     &
           (k == ntml(i)+1 .or. k == ntml(i)+2) ) .or.                         &

           ( dsc_disc_inv(i)  >=  1 .and.                                      &
           (k == ntdsc(i)+1 .or. k == ntdsc(i)+2) ) ) then
        rhokh(i,k) = zero
        rhokm(i,k) = zero
        tke_loc(i,k) = zero
      end if

    end do ! P_POINTS,i
  end do ! BL_LEVELS
!$OMP end do NOWAIT

else if (Keep_Ri_FA == except_disc_inv) then
  !-----------------------------------------------------------------------
          ! Reduce local K's only at the top of a turbulent
          ! layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      if ( ( sml_disc_inv(i)  >=  1 .and.                                      &
           (k == ntml(i)+1 .or. k == ntml(i)+2) ) .or.                         &

           ( dsc_disc_inv(i)  >=  1 .and.                                      &
           (k == ntdsc(i)+1 .or. k == ntdsc(i)+2) ) ) then
        rhokh(i,k) = (one-weight_1dbl(i,k))*rhokh(i,k)
        rhokm(i,k) = (one-weight_1dbl(i,k))*rhokm(i,k)
      end if

    end do ! P_POINTS,i
  end do ! BL_LEVELS
!$OMP end do NOWAIT

else
  !-----------------------------------------------------------------------
  ! Set local K's to zero from the LCL in cumulus and from the
  ! top of a turbulent layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      if ( (cumulus(i) .and. ( (l_param_conv .and. k >  ntml(i))               &
           .or. (.not. l_param_conv .and. k >= ntml(i)) )) .or.                &

           ( dsc_disc_inv(i)  >=  1 .and. k  >   ntdsc(i) ) .or.               &

           ( sml_disc_inv(i)  >=  1 .and. k  >   ntml(i) ) ) then
        !   This also means no local mixing within any DSC layer
        rhokh(i,k) = zero
        rhokm(i,k) = zero
        tke_loc(i,k) = zero
      end if

    end do ! P_POINTS,i
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
    do i = pdims%i_start, pdims%i_end

      if (cumulus(i) .and. ( (l_param_conv .and. k == ntml(i)+1)               &
           .or. (.not. l_param_conv .and.                                      &
                       k >= ntml(i) .and. k <  ntml(i)+2) )) then
        rhokhz(i,k)=zero
        rhokmz(i,k)=zero
        rhokh_top(i,k)=zero
        rhokm_top(i,k)=zero
      end if
    end do ! P_POINTS,i
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
    do i = pdims%i_start, pdims%i_end
      BL_diag%rhokmloc(i,k)=rhokm(i,k)
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokhloc) then
!$OMP do SCHEDULE(STATIC)
  do k = 1, bl_levels
    do i = pdims%i_start, pdims%i_end
      BL_diag%rhokhloc(i,k)=rhokh(i,k)
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokmsurf) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      BL_diag%rhokmsurf(i,k)=weight_1dbl(i,k)*rhokmz(i,k)
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokhsurf) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      BL_diag%rhokhsurf(i,k)=weight_1dbl_rho(i,k)*rhokhz(i,k)
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokmsc) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      BL_diag%rhokmsc(i,k)=weight_1dbl(i,k)*rhokm_top(i,k)
    end do
  end do
!$OMP end do NOWAIT
end if

if (BL_diag%l_rhokhsc) then
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do i = pdims%i_start, pdims%i_end
      BL_diag%rhokhsc(i,k)=weight_1dbl_rho(i,k)*rhokh_top(i,k)
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
  do i = pdims%i_start, pdims%i_end

    rhokh(i,k) = max( rhokh(i,k) ,                                             &
           weight_1dbl_rho(i,k)*(rhokhz(i,k)+rhokh_top(i,k)) )
    rhokm(i,k) = max( rhokm(i,k) ,                                             &
           weight_1dbl(i,k)*(rhokmz(i,k)+rhokm_top(i,k)) )

  end do ! P_POINTS,i
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
