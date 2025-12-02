! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Various tests to determine whether ascent is convective
!
module cumulus_test_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName='CUMULUS_TEST_MOD'

contains

! Subroutine Interface:
subroutine cumulus_test(npnts,mbl,                                             &
                        ntml, k_plume,                                         &
                        zhpar, land_frac, qw, cloud_fraction, z_theta,         &
                        z_lcl, nlcl, cumulus )

! ------------------------------------------------------------------------------
! Description:
!   Works out which unstable points are convective using top of parcel.
!   If parcel top is < 3000m performs a series of tests to check whether
!   the point is convective or stratocumulus.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
use cv_run_mod, only: l_cvdiag_ctop_qmax
use cv_diag_param_mod, only:                                                   &
    sc_cftol
use bl_option_mod, only: zero
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
use nlsizes_namelist_mod, only: model_levels

implicit none

! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Number of points
 ,mbl                    ! maximum number of boundary levels bl_levels-1

integer, intent(in) ::                                                         &
  ntml(npnts)          & ! Top of parcel ascent
 ,k_plume(npnts)         ! Starting model level for plume ascent

real(kind=r_bl), intent(in) ::                                                 &
  zhpar(npnts)                       & ! height of parcel top (m)
 ,land_frac(npnts)                   & ! fraction of land
 ,qw(npnts,model_levels)             & ! total water (kg/kg)
 ,cloud_fraction(npnts,model_levels) & ! layer cloud fraction
 ,z_theta(npnts,model_levels)        & ! height of model theta levels(m)
 ,z_lcl(npnts)                         ! height of LCL (m)

integer, intent(in out) ::                                                     &
  nlcl(npnts)                 ! Lifting condensation level

logical, intent(in out) ::                                                     &
  cumulus(npnts)              ! true if convective point


!-----------------------------------------------------------------------
! Local variables

integer ::                                                                     &
  ii,k,kk           ! loop counter

integer ::                                                                     &
  kcucheck(npnts)   ! Level for calculating cloud-layer qw gradient up to

real(kind=r_bl) ::                                                             &
  grad_cld        & ! SVL gradient in layer above LCL.
 ,grad_sub          ! SVL gradient in layer below LCL.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CUMULUS_TEST'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) private(ii, grad_cld, grad_sub, k, kk)

! Find level just below 2.5km, for use in gradient check
!$OMP do SCHEDULE(STATIC)
do ii=1, npnts
  kcucheck(ii) = 1
end do
!$OMP end do
do k = 2, model_levels
!$OMP do SCHEDULE(STATIC)
  do ii=1, npnts
    if ( z_theta(ii,k) > 2500.0 .and. kcucheck(ii) == 1 )  kcucheck(ii) = k-1
  end do
!$OMP end do
end do

if ( l_cvdiag_ctop_qmax ) then
  ! Optional extra; move kcucheck down if a lower value of qw occurs
  ! below the cloud-top.  This is to avoid a problem where strong
  ! convective detrainment creates a local maximum in the qw profile at
  ! the cloud-top, and we spuriously diagnose it to be well-mixed in qw
  ! based on the gradient between cloud-top and the SML, even though
  ! a much stronger qw gradient occurs between the SML-top and
  ! a local minimum below cloud-top.
!$OMP do SCHEDULE(STATIC)
  do ii=1, npnts
    ! First move kcucheck down to the parcel-top height if above
    ! (on its own this would be redundant because separate if tests
    !  use ntml instead of kcucheck if ntml is lower, but need to do this
    !  to restrict the search for low qw to below cloud-top)
    kcucheck(ii) = min( kcucheck(ii), ntml(ii) )
    if ( kcucheck(ii) > nlcl(ii) + 2 ) then
      ! Don't allow to move kcucheck below nlcl+2
      kk = kcucheck(ii)
      do k = nlcl(ii) + 2, kk
        if ( qw(ii,k) < qw(ii,kcucheck(ii)) )  kcucheck(ii) = k
      end do
    end if
  end do
!$OMP end do
end if  ! ( l_cvdiag_ctop_qmax )

!$OMP do SCHEDULE(STATIC)
do ii=1, npnts

  !-----------------------------------------------------------------------
  !     Check lifting condensation level against height of parcel ascent.
  !     If lifting condensation level lower than parcel ascent, and is
  !     within bl_levels, then decide on type of cloudy layer.
  !     Gradient tests and cumulus parametriztion require a minimum number
  !     of grid-levels between LCL and top of parcel ascent, otherwise
  !     define as stratocumulus.  To avoid resolution sensitivity, also
  !     require this depth to be more than 400m, i.e. cumulus clouds are
  !     deeper than 400m.  Note this depth is physically plausible and
  !     close to the two grid-level requirement in the original L38
  !     implementation.
  !-----------------------------------------------------------------------
  ! Note the following test still uses MBL making it possibly resolution
  ! dependent. This test is only likely to change results if the number of
  ! boundary layer levels is reduced. This is unlikely to be done.
  !-----------------------------------------------------------------------

  if ( ntml(ii)-nlcl(ii) >= 3 .and. nlcl(ii)  <   mbl-1                        &
                              .and. zhpar(ii)-z_lcl(ii) >= 400.0_r_bl ) then

    !-----------------------------------------------------------------------
    !     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
    !     For stratocumulus top of mixed layer = zh
    !     For cumulus top of mixed layer = Z_LCL
    !     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top
    !     above 3000m indicates convection.
    !     Diagnosis is done by comparing gradients
    !-----------------------------------------------------------------------

    if (zhpar(ii) >= 3000.0_r_bl) then

      cumulus(ii) = .true.

      ! Added for low LCL levels
      ! nlcl is not permitted to be less than 2 if Cu is diagnosed
      nlcl(ii) = max(2, nlcl(ii))


      ! Tops < 3000m and nlcl > kplume
    else if ( nlcl(ii)  >   k_plume(ii)) then

      ! Current test is against a height of ~<2.5km
      ! This could be replaced by a scale height if a suitable method
      ! for determining a sensible height was possible from profile/cumulus
      ! depth information available in this routine

      if (ntml(ii)  >   kcucheck(ii)                                           &
                  .and. nlcl(ii)  <=  kcucheck(ii)-2) then

        grad_cld = abs( qw(ii,kcucheck(ii)) - qw(ii,nlcl(ii)) )                &
                    /( z_theta(ii,kcucheck(ii)) - z_theta(ii,nlcl(ii)) )
      else
        grad_cld = abs( qw(ii,ntml(ii))     - qw(ii,nlcl(ii)) )                &
                    /( z_theta(ii,ntml(ii)) - z_theta(ii,nlcl(ii)) )
      end if

      grad_sub   =  abs( qw(ii,nlcl(ii)) - qw(ii,k_plume(ii)) )                &
                    /( z_theta(ii,nlcl(ii)) - z_theta(ii,k_plume(ii)) )

      if (grad_cld  >   1.10_r_bl*grad_sub) then
        !-----------------------------------------------------------------------
        !     Not well mixed, however it is possible that the depth of a well
        !     mixed boundary layer has increased but not yet been mixed yet so
        !     test gradient from next level DOwn.
        !     Note typical cumulus profiles are expected to have a fairly
        !     uniform q profile from the surface to the cloud base and then a
        !     decreasing profile of q above this in the cloud. Typical the
        !     decreasing gradient from the cloud base to 2.5km will be the
        !     order of > 1.10 the below cloud value.
        !-----------------------------------------------------------------------

        ! test against a height ~ 2.5km

        if (ntml(ii)  <=  kcucheck(ii)) then
          grad_cld = abs( qw(ii,ntml(ii)-1) - qw(ii,nlcl(ii)) )                &
                /( z_theta(ii,ntml(ii)-1) - z_theta(ii,nlcl(ii)) )
        end if

        if ( grad_cld  >   1.10_r_bl*grad_sub) then
          !-----------------------------------------------------------------------
          !      Diagnose a cumulus layer
          !-----------------------------------------------------------------------
          cumulus(ii) = .true.
        end if

      else

        ! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
        ! and not yet been mixed (so could have been erroneously identIfied as
        ! well-mixed)

        ! First check using level below (recalculate grad_sub)

        if (nlcl(ii) - k_plume(ii)  >=  2) then

          grad_sub = abs( qw(ii,nlcl(ii)-1) - qw(ii,k_plume(ii)) )             &
                     /( z_theta(ii,nlcl(ii)-1) - z_theta(ii,k_plume(ii)) )

          if ( grad_cld  >   1.10_r_bl*grad_sub) then
            cumulus(ii) =.true.
          end if

        end if

        ! If still diagnosing well-mixed, check using level above
        ! (recalculate grad_cld)

        if (.not. cumulus(ii) ) then

          if (ntml(ii)  >   kcucheck(ii)                                       &
                   .and. nlcl(ii)  <=  kcucheck(ii)-2) then

            grad_cld = abs( qw(ii,kcucheck(ii)) - qw(ii,nlcl(ii)+1) )          &
                       /( z_theta(ii,kcucheck(ii)) - z_theta(ii,nlcl(ii)+1) )
          else
            grad_cld = abs( qw(ii,ntml(ii)) - qw(ii,nlcl(ii)+1) )              &
                       /( z_theta(ii,ntml(ii)) - z_theta(ii,nlcl(ii)+1) )
          end if

          if ( grad_cld  >   1.10_r_bl*grad_sub) then
            cumulus(ii) =.true.
          end if

        end if   ! not cumulus
      end if     ! test on cloud gradient
    end if       ! test on cloud top height
  end if         ! tests on nlcl

end do           ! ii loop
!$OMP end do

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are DOne on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land. The wording of this comment implies
!      the problem occurred for situation where the gradient tests were
!      done i.e. zhpar < 3000m.
!-----------------------------------------------------------------------
! The following code still uses a MBL test. I think this test was designed
! with MBL ~ 3000m in mind and so should really be replaced with
!   zhpar(ii) < 3000.0 instead of ntml(ii)  <   MBL
! Changing this may alter results so I am not doing it yet.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do ii=1, npnts
  k=nlcl(ii)

  if ( (Land_frac(ii) > zero) .and. cumulus(ii) .and.                          &
                                       ntml(ii)  <   mbl ) then
    do while( k  <=  ntml(ii) .and. cloud_fraction(ii,k)                       &
                                                  >=  sc_cftol )
      k = k + 1
    end do
    if (k  ==  ntml(ii)+1) cumulus(ii) = .false.
  end if
end do       ! ii loop
!$OMP end do
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
end subroutine cumulus_test

end module cumulus_test_mod
