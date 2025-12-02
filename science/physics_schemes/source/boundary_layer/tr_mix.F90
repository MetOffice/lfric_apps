! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    subroutine TR_MIX ------------------------------------------------

!    Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

!    Purpose: Calculate tracer flux and pass through to IMP_MIX to solve

!    Programming standard: UMDP 3

!    Documentation: UM Documentation Paper No 24.

! ----------------------------------------------------------------------
module tr_mix_mod
use um_types, only: real_umphys

implicit none
character(len=*), parameter, private :: ModuleName='TR_MIX_MOD'

contains
subroutine tr_mix (                                                            &
! in fields
  r_theta_levels, r_rho_levels, r_dims,                                        &
  bl_levels,gamma_in,rhokh_rdz,rhokh_1,dtrdz,surf_em,res_factor,               &
  kent, we_lim, t_frac, zrzi,                                                  &
  kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,zhnl ,zhsc, z_uv,                 &
! INOUT / out fields
  field,f_field,surf_dep_flux                                                  &
  )

use atm_fields_bounds_mod, only: pdims, array_dims
use imp_mix_mod, only: imp_mix
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! arguments passed in
integer, intent(in) ::                                                         &
  bl_levels   ! in No. of atmospheric levels for boundary layer

type(array_dims), intent(in) :: r_dims ! in size of r_X_levels fields

real(kind=real_umphys), intent(in) ::                                          &
 gamma_in(bl_levels),                                                          &
 r_theta_levels(r_dims%i_start:r_dims%i_end,r_dims%j_start:r_dims%j_end,       &
                 0:bl_levels),                                                 &
 r_rho_levels(r_dims%i_start:r_dims%i_end,r_dims%j_start:r_dims%j_end,         &
               bl_levels),                                                     &
                               ! in height of model rho and theta levels
 rhokh_rdz(pdims%i_start:,pdims%j_start:, 2:),                                 &
                               ! in Mixing coeff. above surface
                               !    = RHOKH(,K)*RDZ(K)
                               !    for K>=2 (from IMP_SOLVER).
 rhokh_1(pdims%i_start:,pdims%j_start:),                                       &
                               ! in  Surface exchange coeff.
                               !     from P243 (SF_EXCH)
 dtrdz(pdims%i_start:,pdims%j_start:, :),                                      &
                               ! in  dt/(rho*r*r*dz) for scalar
                               !     flux divergence
 surf_em(pdims%i_start:,pdims%j_start:),                                       &
                               ! in, Surface emissions in kg/m2/s
 res_factor(pdims%i_start:,pdims%j_start:),                                    &
                               ! in, dry dep coeff=Ra/(Ra+Rb+Rc)
  we_lim(pdims%i_start:,pdims%j_start:,:),                                     &
                               ! in rho*entrainment rate implied by
                               !     placing of subsidence
  zrzi(pdims%i_start:,pdims%j_start:,:),                                       &
                               ! in (z-z_base)/(z_i-z_base)
  t_frac(pdims%i_start:,pdims%j_start:,:),                                     &
                               ! in a fraction of the timestep
  we_lim_dsc(pdims%i_start:,pdims%j_start:,:),                                 &
                               ! in rho*entrainment rate implied by
                               !     placing of subsidence
  zrzi_dsc(pdims%i_start:,pdims%j_start:,:),                                   &
                               ! in (z-z_base)/(z_i-z_base)
  t_frac_dsc(pdims%i_start:,pdims%j_start:,:),                                 &
                               ! in a fraction of the timestep
  z_uv(pdims%i_start:,pdims%j_start:,:),                                       &
                               ! in Z_uv(*,K) is height of half
                               !    level k-1/2.
  zhsc(pdims%i_start:,pdims%j_start:),                                         &
                               ! in Top of decoupled layer
  zhnl(pdims%i_start:,pdims%j_start:)
                               ! in Top of surface mixed layer

integer, intent(in) ::                                                         &
  kent(pdims%i_start:,pdims%j_start:),                                         &
                               ! in grid-level of SML inversion
  kent_dsc(pdims%i_start:,pdims%j_start:)
                               ! in grid-level of DSC inversion

! INOUT arguments
real(kind=real_umphys), intent(in out) ::                                      &
 field(pdims%i_start:,pdims%j_start:,:)
                               ! INOUT Tracer amount in kg/kg.

! out arguments
real(kind=real_umphys), intent(out) ::                                         &
 f_field(pdims%i_start:,pdims%j_start:,:),                                     &
                               ! out Flux of tracer in kg/m2/s.
 surf_dep_flux(pdims%i_start:,pdims%j_start:)
                               ! out, surface deposn flux (kg/m2/s)

!    Local and other symbolic constants :-

real(kind=real_umphys),parameter:: one=1.0
real(kind=real_umphys),parameter:: smallp=tiny(one)

!   Workspace :-

real(kind=real_umphys) ::                                                      &
 rhok_dep(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                                ! Surface deposition coeficient
 gamma_rhokh_rdz(pdims%i_start:pdims%i_end,                                    &
                 pdims%j_start:pdims%j_end,2:bl_levels),                       &
                         ! gamma*RHOKH_RDZ
 dz_disc,                                                                      &
                         ! Temporary in subgrid zi calculation
 dzlkp1,                                                                       &
                         ! Thickness of theta-level K+1
 f_field_ent,                                                                  &
                         ! Time-level n entrainment flux
 dfield_inv            ! inversion jump

! Arrays for vectorisation
real(kind=real_umphys) ::                                                      &
        dfield_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                         ! Jump in field across SML inversion
        dfield_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                         ! Jump in field across DSC inversion

!  Local scalars :-
integer ::                                                                     &
 i,j,                                                                          &
            ! Loop counter (horizontal field index).
 k,                                                                            &
            ! Loop counter (vertical index).
 km1,                                                                          &
            ! Max(K-1,2)
 ient     ! Loop counter for entrainment

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='TR_MIX'

!-----------------------------------------------------------------------
!   0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!$OMP PARALLEL DEFAULT(SHARED) private(km1,k,j,i,dzlkp1,dz_disc,               &
!$OMP                          f_field_ent, dfield_inv, ient)

!-----------------------------------------------------------------------
!  1.  Calculate flux of tracer:
!-----------------------------------------------------------------------
!  1.1 Above the surface
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  if ( k>1 ) then
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        gamma_rhokh_rdz(i,j,k) = gamma_in(k) * rhokh_rdz(i,j,k)
        f_field(i,j,k) = - rhokh_rdz(i,j,k) *                                  &
                           (field(i,j,k) - field(i,j,k-1))
      end do
    end do
  else

    !-----------------------------------------------------------------------
    !  1.2 At the surface: (i) set surface flux equal to input emissions
    !                    (should be passed in as ancillary file, else ZERO)
    !                      (ii) Use input resistance factors to calculate
    !                    surface deposition (if ZERO then no dry deposition)
    !-----------------------------------------------------------------------

    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        f_field(i,j,1) = surf_em(i,j)  ! Inject surface emissions
        rhok_dep(i,j) = res_factor(i,j) * rhokh_1(i,j)
        surf_dep_flux(i,j) = -rhok_dep(i,j) * field(i,j,1)
        rhok_dep(i,j) = gamma_in(1) * rhok_dep(i,j)
      end do
    end do
  end if
end do
!$OMP end do NOWAIT

!-----------------------------------------------------------------------
! Add on explicit entrainment fluxes
! These are calculated from the time-level n profile and parametrized
! entrainment rate to give a time-level n flux (F_FIELD_ENT).  This is
! then implemented using an equivalent RHOKH as the diagnosed inversion
! jump can change significantly across the timestep and so, therefore,
! should the parametrized flux.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    k = kent(i,j)-1     ! equal to originally diagnosed NTML
    !----------------------------
    ! diagnose SML inversion jump
    !----------------------------
    if (k  ==  bl_levels-1) then
      dfield_sml(i,j) = field(i,j,k+1) - field(i,j,k)
    else
      dzlkp1  = z_uv(i,j,k+2) - z_uv(i,j,k+1)
      dz_disc = z_uv(i,j,k+2) - zhnl(i,j)
      if (dz_disc/dzlkp1  >   0.1) then
        dfield_sml(i,j) = (field(i,j,k+1)-field(i,j,k))                        &
                            * dzlkp1 /dz_disc

        if ( field(i,j,k+2)  >   field(i,j,k+1) .and.                          &
                   field(i,j,k+1)  >   field(i,j,k) ) then
          dfield_sml(i,j) = min( field(i,j,k+2)-field(i,j,k),                  &
                                 dfield_sml(i,j) )
        else if ( field(i,j,k+2)  <   field(i,j,k+1) .and.                     &
                  field(i,j,k+1)  <   field(i,j,k) ) then
          dfield_sml(i,j) = max( field(i,j,k+2)-field(i,j,k),                  &
                                 dfield_sml(i,j) )
        else  ! FIELD non-monotonic
          dfield_sml(i,j) = field(i,j,k+1)-field(i,j,k)
        end if
      else
        dfield_sml(i,j) = field(i,j,k+2) - field(i,j,k)
      end if
    end if
  end do

  do i = pdims%i_start, pdims%i_end

    k = kent_dsc(i,j)-1     ! equal to originally diagnosed NTDSC

    !----------------------------
    ! diagnose DSC inversion jump
    !----------------------------
    if (k  ==  bl_levels-1) then
      dfield_dsc(i,j) = field(i,j,k+1) - field(i,j,k)
    else if ( kent_dsc(i,j)  >   2 ) then

      dzlkp1  = z_uv(i,j,k+2) - z_uv(i,j,k+1)
      dz_disc = z_uv(i,j,k+2) - zhsc(i,j)
      if (dz_disc/dzlkp1  >   0.1) then
        dfield_dsc(i,j) = (field(i,j,k+1)-field(i,j,k))                        &
                            * dzlkp1 /dz_disc

        if ( field(i,j,k+2)  >   field(i,j,k+1) .and.                          &
             field(i,j,k+1)  >   field(i,j,k) ) then
          dfield_dsc(i,j) = min( field(i,j,k+2)-field(i,j,k),                  &
                                 dfield_dsc(i,j) )
        else if ( field(i,j,k+2)  <   field(i,j,k+1) .and.                     &
                  field(i,j,k+1)  <   field(i,j,k) ) then
          dfield_dsc(i,j) = max( field(i,j,k+2)-field(i,j,k),                  &
                                 dfield_dsc(i,j) )
        else  ! FIELD non-monotonic
          dfield_dsc(i,j) = field(i,j,k+1)-field(i,j,k)
        end if
      else
        dfield_dsc(i,j) = field(i,j,k+2) - field(i,j,k)
      end if
    end if
  end do
end do
!$OMP end do
! OpenMP needs to (implicitly) barrier here (by not having NOWAIT)
! as dfield_dsc and dfield_sml used below

      !------------------------------------
      ! calculate entrainment fluxes and KH
      !------------------------------------
do ient = 1, 3

!$OMP do SCHEDULE(STATIC)
  do j = pdims%j_start, pdims%j_end
    !CDIR nodep
    do i = pdims%i_start, pdims%i_end

      k = kent(i,j)-2+ient
      if ( k > 1 .and. k <= bl_levels .and.                                    &
                       t_frac(i,j,ient) > 0.0) then
        if ( abs(field(i,j,k)-field(i,j,k-1)) >= smallp ) then
          dfield_inv = dfield_sml(i,j)
          ! DFIELD_INV must have same sign as
          ! local gradient to get right sign of flux
          if ( dfield_sml(i,j) / (field(i,j,k)-field(i,j,k-1))                 &
                <   0.0 ) dfield_inv = field(i,j,k)-field(i,j,k-1)
          f_field_ent = t_frac(i,j,ient) * ( f_field(i,j,1)                    &
            - ( we_lim(i,j,ient) * dfield_inv + f_field(i,j,1) )               &
              * zrzi(i,j,ient) )
          ! interpolation to surface flux must not change the sign
          ! of the entrainment flux otherwise KH will be <0!
          if ( f_field_ent / (field(i,j,k)-field(i,j,k-1))  > 0.0 )            &
            f_field_ent = - t_frac(i,j,ient) *                                 &
                      we_lim(i,j,ient) * dfield_inv * zrzi(i,j,ient)
          ! Restrict size of RHOKH for numerical safety
          km1 = max( k-1, 2 )
          gamma_rhokh_rdz(i,j,k) = gamma_rhokh_rdz(i,j,k) +                    &
                  min( - gamma_in(k)*f_field_ent/                              &
                         ( field(i,j,k)-field(i,j,k-1) ),                      &
                       10.0*gamma_rhokh_rdz(i,j,km1) )
          ! Recalculate explicit flux using entrainment KH
          f_field(i,j,k) = - (gamma_rhokh_rdz(i,j,k) / gamma_in(k)) *          &
                             (field(i,j,k) - field(i,j,k-1))
        end if
      end if
    end do

    !CDIR nodep
    do i = pdims%i_start, pdims%i_end

      k = kent_dsc(i,j)-2+ient
      if ( kent_dsc(i,j) >= 3 .and. k <= bl_levels ) then
        if ( t_frac_dsc(i,j,ient) > 0.0                                        &
             .and. abs(field(i,j,k)-field(i,j,k-1)) >= smallp ) then
          dfield_inv = dfield_dsc(i,j)
          ! DFIELD_INV must have same sign as
          ! local gradient to get right sign of flux
          if ( dfield_dsc(i,j) / (field(i,j,k)-field(i,j,k-1))                 &
                <   0.0 ) dfield_inv = field(i,j,k)-field(i,j,k-1)
          f_field_ent = - t_frac_dsc(i,j,ient) *                               &
             we_lim_dsc(i,j,ient) * dfield_inv * zrzi_dsc(i,j,ient)
          ! Restrict size of RHOKH for numerical safety
          km1 = max( k-1, 2 )
          gamma_rhokh_rdz(i,j,k) = gamma_rhokh_rdz(i,j,k) +                    &
                   min( - gamma_in(k)*f_field_ent/                             &
                                    ( field(i,j,k)-field(i,j,k-1) ),           &
                        10.0*gamma_rhokh_rdz(i,j,km1) )
          ! Recalculate explicit flux using entrainment KH
          f_field(i,j,k) = - (gamma_rhokh_rdz(i,j,k) / gamma_in(k)) *          &
                           (field(i,j,k) - field(i,j,k-1))
        end if
      end if

    end do
  end do
!$OMP end do NOWAIT
end do ! IENT

!$OMP end PARALLEL

!-----------------------------------------------------------------------
!  2.  Call routine IMPL_CAL to calculate incrememnts to tracer field
!      and suface deposition flux for output
!-----------------------------------------------------------------------

call imp_mix (                                                                 &
 bl_levels,dtrdz,r_theta_levels,r_rho_levels,r_dims,                           &
 gamma_rhokh_rdz,rhok_dep,f_field,surf_dep_flux,field                          &
 )

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine tr_mix
end module tr_mix_mod
