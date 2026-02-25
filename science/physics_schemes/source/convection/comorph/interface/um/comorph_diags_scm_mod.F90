! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module comorph_diags_scm_mod

use diag_type_mod, only: diag_list_type

use um_types, only: real_umphys

implicit none

save

! Super-arrays for SCM diagnostics output by CoMorph:

! 4-D diagostics on x,y,layers,updraft_types domain
real(kind=real_umphys), target, allocatable :: scm_diags_xylu(:,:,:,:,:)
! 4-D diagostics on x,y,layers,dndraft_types domain
real(kind=real_umphys), target, allocatable :: scm_diags_xyld(:,:,:,:,:)
! 3-D diagnostics on x,y,z domain
real(kind=real_umphys), target, allocatable :: scm_diags_xyz(:,:,:,:)
! 4-D diagnostics on x,y,z,updraft_types domain
real(kind=real_umphys), target, allocatable :: scm_diags_xyzu(:,:,:,:,:)
! 4-D diagnostics on x,y,z,downdraft_types domain
real(kind=real_umphys), target, allocatable :: scm_diags_xyzd(:,:,:,:,:)
! 5-D diagnostics on x,y,z,layers,updraft_types domain
real(kind=real_umphys), target, allocatable :: scm_diags_xyzlu(:,:,:,:,:,:)
! 5-D diagnostics on x,y,z,layers,downdraft_types domain
real(kind=real_umphys), target, allocatable :: scm_diags_xyzld(:,:,:,:,:,:)

! Number of diagnostics in each category
integer :: n_diags_xylt
integer :: n_diags_xyz
integer :: n_diags_xyzt
integer :: n_diags_xyzlt
! Assuming same number of diags from updrafts and downdrafts

! Lists of pointers to the active diagnostic structures within
! the comorph_diags structure
type(diag_list_type), allocatable ::  scm_list_xylu(:)
type(diag_list_type), allocatable ::  scm_list_xyld(:)
type(diag_list_type), allocatable ::  scm_list_xyz(:)
type(diag_list_type), allocatable ::  scm_list_xyzu(:)
type(diag_list_type), allocatable ::  scm_list_xyzd(:)
type(diag_list_type), allocatable ::  scm_list_xyzlu(:)
type(diag_list_type), allocatable ::  scm_list_xyzld(:)


contains


!----------------------------------------------------------------
! Subroutine to set diagnostic requests and assign comorph's
! diagnostic pointers to memory.
!----------------------------------------------------------------
subroutine comorph_diags_scm_reqs( comorph_diags, l_tracer )

use comorph_diags_type_mod, only: comorph_diags_type
use draft_diags_type_mod, only: draft_diags_type
use comorph_constants_mod, only: n_updraft_types, n_dndraft_types,             &
                                 n_conv_layers_diag, n_tracers,                &
                                 l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,    &
                                 l_cv_cloudfrac, l_par_core,                   &
                                 l_turb_par_gen
use subregion_mod, only: n_regions
use nlsizes_namelist_mod, only: row_length, rows, model_levels
use raise_error_mod, only: raise_fatal

implicit none

! Structure storing meta-data and pointers to fields for all
! diagnostics calculated by CoMorph
type(comorph_diags_type), target, intent(in out) :: comorph_diags

! Flag for whether tracers are available in this call
logical, intent(in) :: l_tracer

! Number of primary fields for which the SCM outputs diagnostics
integer :: n_fields
! Number of condensed water species
integer :: n_cond
! Number of other parcel properties
integer :: n_par

! Pointer to either updraft or downdraft diags structure and array
type(draft_diags_type), pointer :: draft_diags => null()
! Number of convection types for updrafts or downdrafts
integer :: n_conv_types
! Number of convective drafts for xyzlt domain diagnostics
integer :: n_drafts

! Counters
integer :: i, j, k, i_type, i_layr, i_diag, i_tracer, i_draft, i_region

character(len=*), parameter :: routinename                                     &
                               = "COMORPH_DIAGS_SCM_REQS"


! Work out how many primary fields there are

n_cond = 1  ! qcl non-optional
if ( l_cv_rain )   n_cond = n_cond + 1
if ( l_cv_cf )     n_cond = n_cond + 1
if ( l_cv_snow )   n_cond = n_cond + 1
if ( l_cv_graup )  n_cond = n_cond + 1

n_fields = 5 + n_cond  ! u,v,w,T,q + q_cond
if ( l_cv_cloudfrac )  n_fields = n_fields + 3
if ( l_tracer )        n_fields = n_fields + n_tracers

! Set number of non-primary-field parcel properties
n_par = 2  ! mass-flux, radius



! Diagnostics on x,y,layer,updraft/downdraft types domain...

! Work out how many diagnostics there'll be
               ! CAPE + mass-flux-weighted CAPE
               ! term height, ccb height, ccb area, ccb massflux,
               ! ccb mean Tv excess, ccb core Tv excess
n_diags_xylt = 8

if ( n_diags_xylt > 0 ) then

  if ( n_updraft_types > 0 ) then
    ! Diagnostics for updraft types

    ! Allocate super-array and list
    allocate( scm_diags_xylu( row_length, rows,                                &
                              n_conv_layers_diag, n_updraft_types,             &
                              n_diags_xylt ) )
    allocate( scm_list_xylu( n_diags_xylt ) )
    ! Initialise to zero
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, i_layr, i_type, i_diag )        &
!$OMP SHARED( row_length, rows, n_diags_xylt, scm_diags_xylu)
    do i_diag = 1, n_diags_xylt
      do i_type = 1, n_updraft_types
        do i_layr = 1, n_conv_layers_diag
          do j = 1, rows
            do i = 1, row_length
              scm_diags_xylu(i,j,i_layr,i_type,i_diag) = 0.0
            end do
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

    ! Initialise diagnostic counter to zero
    i_diag = 0

    ! CAPE and mass-flux-weighted CAPE
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % cape )
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % mfw_cape )
    ! Convective cloud top and base properties
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % term_height )
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % ccb_height )
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % ccb_area )
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % ccb_massflux_d )
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % ccb_mean_tv_ex )
    call assign_diag_xylt( i_diag, 1,                                          &
                           comorph_diags % updraft_diags_2d % ccb_core_tv_ex )

    ! Check the number of diagnostics counted equals the number
    ! specified and used as the list length
    if ( i_diag < n_diags_xylt ) then
      call raise_fatal( routinename,                                           &
       "Number of CoMorph SCM diags is less than the "      //                 &
       "length of the list.  Try decreasing n_diags_xylt." )
    end if

  end if  ! ( n_updraft_types > 0 )

  if ( n_dndraft_types > 0 ) then
    ! Diagnostics for dndraft types

    ! Allocate super-array and list
    allocate( scm_diags_xyld( row_length, rows,                                &
                              n_conv_layers_diag, n_dndraft_types,             &
                              n_diags_xylt ) )
    allocate( scm_list_xyld( n_diags_xylt ) )
    ! Initialise to zero
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, i_layr, i_type, i_diag )        &
!$OMP SHARED( row_length, rows, n_diags_xylt, scm_diags_xyld, n_dndraft_types )
    do i_diag = 1, n_diags_xylt
      do i_type = 1, n_dndraft_types
        do i_layr = 1, n_conv_layers_diag
          do j = 1, rows
            do i = 1, row_length
              scm_diags_xyld(i,j,i_layr,i_type,i_diag) = 0.0
            end do
          end do
        end do
      end do
    end do
!$OMP END PARALLEL DO

    ! Initialise diagnostic counter to zero
    i_diag = 0

    ! CAPE and mass-flux-weighted CAPE
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % cape )
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % mfw_cape )
    ! Convective cloud top and base properties
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % term_height )
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % ccb_height )
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % ccb_area )
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % ccb_massflux_d )
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % ccb_mean_tv_ex )
    call assign_diag_xylt( i_diag, 2,                                          &
                           comorph_diags % dndraft_diags_2d % ccb_core_tv_ex )

    ! Check the number of diagnostics counted equals the number
    ! specified and used as the list length
    if ( i_diag < n_diags_xylt ) then
      call raise_fatal( routinename,                                           &
       "Number of CoMorph SCM diags is less than the "      //                 &
       "length of the list.  Try decreasing n_diags_xylt." )
    end if

  end if  ! ( n_dndraft_types > 0 )

end if  ! ( n_diags_xylt > 0 )



! Diagnostics on x,y,z domain...

! Work out how many diagnostics there'll be
              ! Primary field increments and input values
n_diags_xyz = n_fields*2                                                       &
              ! layer_mass
            + 1
              ! sub-grid clear, cloudy, and rainy/icy regions
n_diags_xyz = n_diags_xyz + n_regions*5


if ( n_diags_xyz > 0 ) then

  ! Allocate super-array and list
  allocate( scm_diags_xyz( row_length, rows, model_levels,                     &
                           n_diags_xyz ) )
  allocate( scm_list_xyz( n_diags_xyz ) )

  ! Initialise super-array to zero
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, k, i_diag )                     &
!$OMP SHARED( row_length, rows, model_levels, n_diags_xyz, scm_diags_xyz )
  do i_diag = 1, n_diags_xyz
    do k = 1, model_levels
      do j = 1, rows
        do i = 1, row_length
          scm_diags_xyz(i,j,k,i_diag) = 0.0
        end do
      end do
    end do
  end do
!$OMP END PARALLEL DO

  ! Initialise diagnostic counter to zero
  i_diag = 0

  ! Increments for each primary field
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_incr % wind_u )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_incr % wind_v )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_incr % wind_w )
  call assign_diag_xyz( i_diag,                                                &
                   comorph_diags % fields_incr % temperature )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_incr % q_vap )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_incr % q_cl )
  if ( l_cv_rain )  call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_incr % q_rain )
  if ( l_cv_cf )    call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_incr % q_cf )
  if ( l_cv_snow )  call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_incr % q_snow )
  if ( l_cv_graup ) call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_incr % q_graup )
  if ( l_cv_cloudfrac ) then
    call assign_diag_xyz( i_diag,                                              &
                          comorph_diags % fields_incr % cf_liq )
    call assign_diag_xyz( i_diag,                                              &
                          comorph_diags % fields_incr % cf_ice )
    call assign_diag_xyz( i_diag,                                              &
                          comorph_diags % fields_incr % cf_bulk )
  end if
  if ( l_tracer .and. n_tracers > 0 ) then
    allocate( comorph_diags % fields_incr % tracers(n_tracers) )
    do i_tracer = 1, n_tracers
      call assign_diag_xyz( i_diag,                                            &
             comorph_diags % fields_incr % tracers(i_tracer) )
    end do
  end if

  ! Input values for each primary field
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_inp % wind_u )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_inp % wind_v )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_inp % wind_w )
  call assign_diag_xyz( i_diag,                                                &
                   comorph_diags % fields_inp % temperature )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_inp % q_vap )
  call assign_diag_xyz( i_diag,                                                &
                        comorph_diags % fields_inp % q_cl )
  if ( l_cv_rain )  call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_inp % q_rain )
  if ( l_cv_cf )    call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_inp % q_cf )
  if ( l_cv_snow )  call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_inp % q_snow )
  if ( l_cv_graup ) call assign_diag_xyz( i_diag,                              &
                        comorph_diags % fields_inp % q_graup )
  if ( l_cv_cloudfrac ) then
    call assign_diag_xyz( i_diag,                                              &
                          comorph_diags % fields_inp % cf_liq )
    call assign_diag_xyz( i_diag,                                              &
                          comorph_diags % fields_inp % cf_ice )
    call assign_diag_xyz( i_diag,                                              &
                          comorph_diags % fields_inp % cf_bulk )
  end if
  if ( l_tracer .and. n_tracers > 0 ) then
    allocate( comorph_diags % fields_inp % tracers(n_tracers) )
    do i_tracer = 1, n_tracers
      call assign_diag_xyz( i_diag,                                            &
             comorph_diags % fields_inp % tracers(i_tracer) )
    end do
  end if

  ! Miscellaneous
  call assign_diag_xyz( i_diag,                                                &
                    comorph_diags % layer_mass )

  ! Sub-grid clear, cloudy and icy/rainy properties...
  do i_region = 1, n_regions
    ! Fractional area of each region
    call assign_diag_xyz( i_diag,                                              &
           comorph_diags % genesis_diags % subregion_diags(i_region) % frac )
    ! Temperature of each region
    call assign_diag_xyz( i_diag,                                              &
           comorph_diags % genesis_diags % subregion_diags(i_region) % fields  &
                         % temperature )
    ! Water-vapour in each region
    call assign_diag_xyz( i_diag,                                              &
           comorph_diags % genesis_diags % subregion_diags(i_region) % fields  &
                         % q_vap )
    ! Relative humidity in each region
    call assign_diag_xyz( i_diag,                                              &
           comorph_diags % genesis_diags % subregion_diags(i_region) % fields  &
                         % rel_hum_liq )
    ! Moist static stability for lifted parcel
    call assign_diag_xyz( i_diag,                                              &
           comorph_diags % genesis_diags % subregion_diags(i_region) % Nsq_up )
  end do

  ! Check the number of diagnostics counted equals the number
  ! specified and used as the list length
  if ( i_diag < n_diags_xyz ) then
    call raise_fatal( routinename,                                             &
         "Number of CoMorph SCM diags is less than the "      //               &
         "length of the list.  Try decreasing n_diags_xyz." )
  end if

end if



! Diagnostics on x,y,z,updraft/downdraft types domain...

n_fields = 5 + n_cond + 1  ! u,v,w,T,q + q_cond + Tv'

! Work out how many diagnostics there'll be
               ! In-plume mean primary-field values,
               ! detrained air properties,
               ! initiating air properties
n_diags_xyzt = 3*n_fields                                                      &
               ! Other parcel properties
             + n_par                                                           &
               ! Initiating mass, entrainment, detrainment
             + 3
               ! In-plume and initiating core primary field values
if ( l_par_core )  n_diags_xyzt = n_diags_xyzt + 2*n_fields

if ( n_diags_xyzt > 0 ) then

  ! Loop over draft types (updraft, downdraft)
  n_drafts = 2
  do i_draft = 1, n_drafts

    ! Assign pointer and allocate arrays for either
    ! updraft or downdraft
    select case (i_draft)
    case (1)
      n_conv_types = n_updraft_types
      if ( n_conv_types > 0 ) then
        draft_diags => comorph_diags % updraft
        allocate( scm_diags_xyzu( row_length, rows, model_levels,              &
                                  n_conv_types, n_diags_xyzt ) )
        allocate( scm_list_xyzu( n_diags_xyzt ) )
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, k, i_type, i_diag )             &
!$OMP SHARED( row_length, rows, model_levels, n_conv_types, n_diags_xyzt,      &
!$OMP         scm_diags_xyzu )
        do i_diag = 1, n_diags_xyzt
          do i_type = 1, n_conv_types
            do k = 1, model_levels
              do j = 1, rows
                do i = 1, row_length
                  scm_diags_xyzu(i,j,k,i_type,i_diag) = 0.0
                end do
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if
    case (2)
      n_conv_types = n_dndraft_types
      if ( n_conv_types > 0 ) then
        draft_diags => comorph_diags % dndraft
        allocate( scm_diags_xyzd( row_length, rows, model_levels,              &
                                  n_conv_types, n_diags_xyzt ) )
        allocate( scm_list_xyzd( n_diags_xyzt ) )
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, k, i_type, i_diag )             &
!$OMP SHARED( row_length, rows, model_levels, n_conv_types, n_diags_xyzt,      &
!$OMP         scm_diags_xyzd )
        do i_diag = 1, n_diags_xyzt
          do i_type = 1, n_conv_types
            do k = 1, model_levels
              do j = 1, rows
                do i = 1, row_length
                  scm_diags_xyzd(i,j,k,i_type,i_diag) = 0.0
                end do
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if
    end select

    if ( n_conv_types > 0 ) then

      ! Initialise diagnostic counter to zero
      i_diag = 0

      ! In-plume mean primary-field values
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % mean % wind_u )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % mean % wind_v )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % mean % wind_w )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                        draft_diags % par % mean % temperature )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % mean % q_vap )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % mean % q_cl )
      if ( l_cv_rain )  call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % par % mean % q_rain )
      if ( l_cv_cf )    call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % par % mean % q_cf )
      if ( l_cv_snow )  call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % par % mean % q_snow )
      if ( l_cv_graup ) call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % par % mean % q_graup )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                   draft_diags % par % mean % virt_temp_excess )

      ! In-plume core primary-field values, if used
      if ( l_par_core ) then
        call assign_diag_xyzt( i_diag, i_draft,                                &
                             draft_diags % par % core % wind_u )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                             draft_diags % par % core % wind_v )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                             draft_diags % par % core % wind_w )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                        draft_diags % par % core % temperature )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                               draft_diags % par % core % q_vap )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                               draft_diags % par % core % q_cl )
        if ( l_cv_rain )  call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % par % core % q_rain )
        if ( l_cv_cf )    call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % par % core % q_cf )
        if ( l_cv_snow )  call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % par % core % q_snow )
        if ( l_cv_graup ) call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % par % core % q_graup )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                   draft_diags % par % core % virt_temp_excess )
      end if  ! ( l_par_core )

      ! Detrained primary-field values
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % plume_model % det_fields % wind_u )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % plume_model % det_fields % wind_v )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % plume_model % det_fields % wind_w )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                        draft_diags % plume_model % det_fields % temperature )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % plume_model % det_fields % q_vap )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % plume_model % det_fields % q_cl )
      if ( l_cv_rain )  call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % plume_model % det_fields % q_rain )
      if ( l_cv_cf )    call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % plume_model % det_fields % q_cf )
      if ( l_cv_snow )  call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % plume_model % det_fields % q_snow )
      if ( l_cv_graup ) call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % plume_model % det_fields % q_graup )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                   draft_diags % plume_model % det_fields % virt_temp_excess )

      ! Initiating parcel primary-field values
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % wind_u )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % wind_v )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % wind_w )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % temperature )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % q_vap )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % q_cl )
      if ( l_cv_rain )  call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % gen % mean % q_rain )
      if ( l_cv_cf )    call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % gen % mean % q_cf )
      if ( l_cv_snow )  call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % gen % mean % q_snow )
      if ( l_cv_graup ) call assign_diag_xyzt( i_diag, i_draft,                &
                             draft_diags % gen % mean % q_graup )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % mean % virt_temp_excess )

      ! Initiating core primary-field values, if used
      if ( l_par_core ) then
        call assign_diag_xyzt( i_diag, i_draft,                                &
                             draft_diags % gen % core % wind_u )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                             draft_diags % gen % core % wind_v )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                             draft_diags % gen % core % wind_w )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                        draft_diags % gen % core % temperature )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                               draft_diags % gen % core % q_vap )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                               draft_diags % gen % core % q_cl )
        if ( l_cv_rain )  call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % gen % core % q_rain )
        if ( l_cv_cf )    call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % gen % core % q_cf )
        if ( l_cv_snow )  call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % gen % core % q_snow )
        if ( l_cv_graup ) call assign_diag_xyzt( i_diag, i_draft,              &
                             draft_diags % gen % core % q_graup )
        call assign_diag_xyzt( i_diag, i_draft,                                &
                   draft_diags % gen % core % virt_temp_excess )
      end if  ! ( l_par_core )

      ! Other parcel properties
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % massflux_d )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % par % radius )


      ! Initiating mass
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                             draft_diags % gen % massflux_d )

      ! Entrained and detrained mass
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                        draft_diags % plume_model % ent_mass_d )
      call assign_diag_xyzt( i_diag, i_draft,                                  &
                        draft_diags % plume_model % det_mass_d )

      ! Check the number of diagnostics counted equals the number
      ! specified and used as the list length
      if ( i_diag < n_diags_xyzt ) then
        call raise_fatal( routinename,                                         &
         "Number of CoMorph SCM diags is less than the "      //               &
         "length of the list.  Try decreasing n_diags_xyzt." )
      end if

      draft_diags => null()

    end if  ! ( n_conv_types > 0 )

  end do  ! i_draft = 1, n_drafts

end if  ! ( n_diags_xyzt > 0 )


! Diagnostics on x,y,z,layer,updraft/downdraft types domain...

! Work out how many diagnostics there'll be
                ! mass-flux, radius, mean Tv excess, core/mean ratio
n_diags_xyzlt = 4
                ! Add core Tv excess if using core
if ( l_par_core )  n_diags_xyzlt = n_diags_xyzlt + 1

if ( n_diags_xyzlt > 0 ) then

  n_drafts = 2

  ! Loop over draft types (updraft, downdraft)
  do i_draft = 1, n_drafts

    ! Assign pointer and allocate arrays for either
    ! updraft or downdraft
    select case (i_draft)
    case (1)
      n_conv_types = n_updraft_types
      if ( n_conv_types > 0 ) then
        draft_diags => comorph_diags % updraft
        allocate( scm_diags_xyzlu( row_length, rows, model_levels,             &
                                   n_conv_layers_diag, n_conv_types,           &
                                   n_diags_xyzlt ) )
        allocate( scm_list_xyzlu( n_diags_xyzlt ) )
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, k, i_layr, i_type, i_diag )     &
!$OMP SHARED( row_length, rows, model_levels, n_conv_types, n_diags_xyzlt,     &
!$OMP         scm_diags_xyzlu )
        do i_diag = 1, n_diags_xyzlt
          do i_type = 1, n_conv_types
            do i_layr = 1, n_conv_layers_diag
              do k = 1, model_levels
                do j = 1, rows
                  do i = 1, row_length
                    scm_diags_xyzlu(i,j,k,i_layr,i_type,i_diag) = 0.0
                  end do
                end do
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if
    case (2)
      n_conv_types = n_dndraft_types
      if ( n_conv_types > 0 ) then
        draft_diags => comorph_diags % dndraft
        allocate( scm_diags_xyzld( row_length, rows, model_levels,             &
                                   n_conv_layers_diag, n_conv_types,           &
                                   n_diags_xyzlt ) )
        allocate( scm_list_xyzld( n_diags_xyzlt ) )
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE( i, j, k, i_layr, i_type, i_diag )     &
!$OMP SHARED( row_length, rows, model_levels, n_conv_types, n_diags_xyzlt,     &
!$OMP         scm_diags_xyzld )
        do i_diag = 1, n_diags_xyzlt
          do i_type = 1, n_conv_types
            do i_layr = 1, n_conv_layers_diag
              do k = 1, model_levels
                do j = 1, rows
                  do i = 1, row_length
                    scm_diags_xyzld(i,j,k,i_layr,i_type,i_diag) = 0.0
                  end do
                end do
              end do
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end if
    end select

    if ( n_conv_types > 0 ) then

      ! Initialise diagnostic counter to zero
      i_diag = 0

      ! Mass-flux and parcel radius
      call assign_diag_xyzlt( i_diag, i_draft,                                 &
                              draft_diags % par % massflux_d )
      call assign_diag_xyzlt( i_diag, i_draft,                                 &
                              draft_diags % par % radius )

      ! Parcel mean Tv excess
      call assign_diag_xyzlt( i_diag, i_draft,                                 &
                   draft_diags % par % mean % virt_temp_excess )

      ! Core/mean ratio
      call assign_diag_xyzlt( i_diag, i_draft,                                 &
                   draft_diags % plume_model % core_mean_ratio )

      ! Core Tv excess if used
      if ( l_par_core ) then
        call assign_diag_xyzlt( i_diag, i_draft,                               &
                   draft_diags % par % core % virt_temp_excess )
      end if

      ! Check the number of diagnostics counted equals the number
      ! specified and used as the list length
      if ( i_diag < n_diags_xyzlt ) then
        call raise_fatal( routinename,                                         &
         "Number of CoMorph SCM diags is less than the "      //               &
         "length of the list.  Try decreasing n_diags_xyzlt." )
      end if

      draft_diags => null()

    end if  ! ( n_conv_types > 0 )

  end do  ! i_draft = 1, n_drafts

end if  ! ( n_diags_xyzlt > 0 )


return
end subroutine comorph_diags_scm_reqs


!----------------------------------------------------------------
! Subroutine to process SCM diags from CoMorph; calls scmoutput
! for each diag, then deallocates the work arrays
!----------------------------------------------------------------
subroutine comorph_diags_scm_proc( )

use comorph_constants_mod, only: n_updraft_types, n_dndraft_types,             &
                                 n_conv_layers_diag
use s_scmop_mod, only: default_streams,                                        &
                       t_inst, d_all, d_point
use scmoutput_mod, only: scmoutput

implicit none

! Strings to store convection type / layer
character(len=3) :: ctype
character(len=3) :: clayr
! String to store current diag name with convection type appended
character(len=100) :: diag_name

! Counters
integer :: i_diag, i_type, i_layr

character(len=*), parameter :: RoutineName                                     &
                               = "COMORPH_DIAGS_SCM_PROC"


! If any diagnostics on n_diags_xylt
! (updraft and downdraft type and layer 2D) domains
if ( n_diags_xylt > 0 ) then

  ! If any for updrafts
  if ( n_updraft_types > 0 ) then

    ! Loop over updraft diagnostics
    do i_diag = 1, n_diags_xylt
      ! Loop over updraft types
      do i_type = 1, n_updraft_types
        ! Loop over updraft layers
        do i_layr = 1, n_conv_layers_diag

          ! Create diag name with current type and layer indices appended
          write(ctype,"(I3)") i_type
          write(clayr,"(I3)") i_layr
          diag_name = "updraft" // "_" // trim(adjustl(ctype)) //              &
                                   "_" // trim(adjustl(clayr)) //              &
            trim(adjustl(scm_list_xylu(i_diag)%pt%diag_name(8:)))

          ! Output data for current type
          call scmoutput( scm_diags_xylu(:,:,i_layr,i_type,i_diag),            &
                          diag_name, diag_name,                                &
                          " ", t_inst, d_point, default_streams, '',           &
                          routinename )

        end do
      end do
    end do

    deallocate( scm_list_xylu )
    deallocate( scm_diags_xylu )

  end if

  ! If any for downdrafts
  if ( n_dndraft_types > 0 ) then

    ! Loop over dndraft diagnostics
    do i_diag = 1, n_diags_xylt
      ! Loop over dndraft types
      do i_type = 1, n_dndraft_types
        ! Loop over dndraft layers
        do i_layr = 1, n_conv_layers_diag

          ! Create diag name with current type and layer indices appended
          write(ctype,"(I3)") i_type
          write(clayr,"(I3)") i_layr
          diag_name = "dndraft" // "_" // trim(adjustl(ctype)) //              &
                                   "_" // trim(adjustl(clayr)) //              &
            trim(adjustl(scm_list_xyld(i_diag)%pt%diag_name(8:)))

          ! Output data for current type
          call scmoutput( scm_diags_xyld(:,:,i_layr,i_type,i_diag),            &
                          diag_name, diag_name,                                &
                          " ", t_inst, d_point, default_streams, '',           &
                          routinename )

        end do
      end do
    end do

    deallocate( scm_list_xyld )
    deallocate( scm_diags_xyld )

  end if

end if  ! ( n_diags_xylt > 0 )


! If any diagnostics on x,y,z domain
if ( n_diags_xyz > 0 ) then

  ! Loop over the requested diags in the list and call scmoutput
  ! for each one.  Using diagnostic names stored in the list.
  do i_diag = 1, n_diags_xyz
    call scmoutput( scm_diags_xyz(:,:,:,i_diag),                               &
                    scm_list_xyz(i_diag)%pt % diag_name,                       &
                    scm_list_xyz(i_diag)%pt % diag_name,                       &
                    " ", t_inst, d_all, default_streams, '',                   &
                    routinename )
  end do

  deallocate( scm_list_xyz )
  deallocate( scm_diags_xyz )

end if


! If any diagnostics on n_diags_xyzt
! (updraft and downdraft type) domains
if ( n_diags_xyzt > 0 ) then

  ! Note: in the calls below, scmoutput is called separately for
  ! each convection type.  i.e. a separate diagnostic field is
  ! created for each type.  This is really annoying; it would be
  ! much more sensible to create one super-diagnostic with a
  ! 4th dimension for convection type.  However, the SCM
  ! diagnostic system is hardwired to allow only 3 dimensions.
  ! Outputting 4-D diags from the SCM would take a lot of
  ! work modifying the scm diag system.

  ! If any for updrafts
  if ( n_updraft_types > 0 ) then

    ! Loop over updraft diagnostics
    do i_diag = 1, n_diags_xyzt
      ! Loop over updraft types
      do i_type = 1, n_updraft_types

        ! Create diag name with current type index appended
        write(ctype,"(I3)") i_type
        diag_name = "updraft" // trim(adjustl(ctype)) //                       &
          trim(adjustl(scm_list_xyzu(i_diag)%pt%diag_name(8:)))

        ! Output data for current type
        call scmoutput( scm_diags_xyzu(:,:,:,i_type,i_diag),                   &
                        diag_name, diag_name,                                  &
                        " ", t_inst, d_all, default_streams, '',               &
                        routinename )

      end do
    end do

    deallocate( scm_list_xyzu )
    deallocate( scm_diags_xyzu )

  end if

  ! If any for downdrafts
  if ( n_dndraft_types > 0 ) then

    ! Loop over updraft diagnostics
    do i_diag = 1, n_diags_xyzt
      ! Loop over downdraft types
      do i_type = 1, n_dndraft_types

        ! Create diag name with current type index appended
        write(ctype,"(I3)") i_type
        diag_name = "dndraft" // trim(adjustl(ctype)) //                       &
          trim(adjustl(scm_list_xyzd(i_diag)%pt%diag_name(8:)))

        ! Output data for current type
        call scmoutput( scm_diags_xyzd(:,:,:,i_type,i_diag),                   &
                        diag_name, diag_name,                                  &
                        " ", t_inst, d_all, default_streams, '',               &
                        routinename )

      end do
    end do

    deallocate( scm_list_xyzd )
    deallocate( scm_diags_xyzd )

  end if

end if  ! ( n_diags_xyzt > 0 )


! If any diagnostics on n_diags_xyzlt
! (updraft and downdraft type and layer) domains
if ( n_diags_xyzlt > 0 ) then

  ! If any for updrafts
  if ( n_updraft_types > 0 ) then

    ! Loop over updraft diagnostics
    do i_diag = 1, n_diags_xyzlt
      ! Loop over updraft types
      do i_type = 1, n_updraft_types
        ! Loop over updraft layers
        do i_layr = 1, n_conv_layers_diag

          ! Create diag name with current type and layer indices appended
          write(ctype,"(I3)") i_type
          write(clayr,"(I3)") i_layr
          diag_name = "updraft" // "_" // trim(adjustl(ctype)) //              &
                                   "_" // trim(adjustl(clayr)) //              &
            trim(adjustl(scm_list_xyzlu(i_diag)%pt%diag_name(8:)))

          ! Output data for current type
          call scmoutput( scm_diags_xyzlu(:,:,:,i_layr,i_type,i_diag),         &
                          diag_name, diag_name,                                &
                          " ", t_inst, d_all, default_streams, '',             &
                          routinename )

        end do
      end do
    end do

    deallocate( scm_list_xyzlu )
    deallocate( scm_diags_xyzlu )

  end if

  ! If any for downdrafts
  if ( n_dndraft_types > 0 ) then

    ! Loop over dndraft diagnostics
    do i_diag = 1, n_diags_xyzlt
      ! Loop over dndraft types
      do i_type = 1, n_dndraft_types
        ! Loop over dndraft layers
        do i_layr = 1, n_conv_layers_diag

          ! Create diag name with current type and layer indices appended
          write(ctype,"(I3)") i_type
          write(clayr,"(I3)") i_layr
          diag_name = "dndraft" // "_" // trim(adjustl(ctype)) //              &
                                   "_" // trim(adjustl(clayr)) //              &
            trim(adjustl(scm_list_xyzld(i_diag)%pt%diag_name(8:)))

          ! Output data for current type
          call scmoutput( scm_diags_xyzld(:,:,:,i_layr,i_type,i_diag),         &
                          diag_name, diag_name,                                &
                          " ", t_inst, d_all, default_streams, '',             &
                          routinename )

        end do
      end do
    end do

    deallocate( scm_list_xyzld )
    deallocate( scm_diags_xyzld )

  end if

end if  ! ( n_diags_xyzlt > 0 )


return
end subroutine comorph_diags_scm_proc


!----------------------------------------------------------------
! Routine to switch on a diagnostic and assign its pointer to
! the next address in the diags super-array, for x,y,layer,type
! domain
!----------------------------------------------------------------
subroutine assign_diag_xylt( i_diag, i_draft, diag )

use diag_type_mod, only: diag_type
use raise_error_mod, only: raise_fatal

implicit none

! Current value of diagnostic counter
integer, intent(in out) :: i_diag
! Indicator for updraft or downdraft
integer, intent(in) :: i_draft
! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

character(len=*), parameter :: routinename = "ASSIGN_DIAG_XYLT"

! Set flag to request diag in 4-D in x,y,l,t
diag % request % x_y_lay_typ = .true.
! Increment diagnostic counter
i_diag = i_diag + 1

! Check not gone beyond the end of the list
if ( i_diag > n_diags_xylt ) then
  call raise_fatal( routinename,                                               &
         "Number of CoMorph SCM diags exceeds the length "    //               &
         "of the list.  Try increasing n_diags_xylt." )
end if

! Assign pointers for either updraft or downdraft accordingly
select case (i_draft)
case (1)
  scm_list_xylu(i_diag)%pt => diag
  diag % field_4d => scm_diags_xylu(:,:,:,:,i_diag)
case (2)
  scm_list_xyld(i_diag)%pt => diag
  diag % field_4d => scm_diags_xyld(:,:,:,:,i_diag)
end select

return
end subroutine assign_diag_xylt


!----------------------------------------------------------------
! Routine to switch on a diagnostic and assign its pointer to
! the next address in the diags super-array, for x,y,z domain
!----------------------------------------------------------------
subroutine assign_diag_xyz( i_diag, diag )

use diag_type_mod, only: diag_type
use raise_error_mod, only: raise_fatal

implicit none

! Current value of diagnostic counter
integer, intent(in out) :: i_diag
! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

character(len=*), parameter :: routinename = "ASSIGN_DIAG_XYZ"

! Set flag to request diag in 3-D in x,y,z
diag % request % x_y_z = .true.
! Increment diagnostic counter
i_diag = i_diag + 1

! Check not gone beyond the end of the list
if ( i_diag > n_diags_xyz ) then
  call raise_fatal( routinename,                                               &
         "Number of CoMorph SCM diags exceeds the length "    //               &
         "of the list.  Try increasing n_diags_xyz." )
end if

! Point element of diag meta-data pointer list at structure
scm_list_xyz(i_diag)%pt => diag
! Assign pointer to the next address in the diags super-array
diag % field_3d => scm_diags_xyz(:,:,:,i_diag)

return
end subroutine assign_diag_xyz


!----------------------------------------------------------------
! Version of the above, for diags that are for each conv type
!----------------------------------------------------------------
subroutine assign_diag_xyzt( i_diag, i_draft, diag )

use diag_type_mod, only: diag_type
use raise_error_mod, only: raise_fatal

implicit none

! Current value of diagnostic counter
integer, intent(in out) :: i_diag
! Indicator for updraft or downdraft
integer, intent(in) :: i_draft
! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

character(len=*), parameter :: routinename = "ASSIGN_DIAG_XYZT"

! Set flag to request diag in 4-D in x,y,z,t
diag % request % x_y_z_typ = .true.
! Increment diagnostic counter
i_diag = i_diag + 1

! Check not gone beyond the end of the list
if ( i_diag > n_diags_xyzt ) then
  call raise_fatal( routinename,                                               &
         "Number of CoMorph SCM diags exceeds the length "    //               &
         "of the list.  Try increasing n_diags_xyzt." )
end if

! Assign pointers for either updraft or downdraft accordingly
select case (i_draft)
case (1)
  scm_list_xyzu(i_diag)%pt => diag
  diag % field_4d => scm_diags_xyzu(:,:,:,:,i_diag)
case (2)
  scm_list_xyzd(i_diag)%pt => diag
  diag % field_4d => scm_diags_xyzd(:,:,:,:,i_diag)
end select

return
end subroutine assign_diag_xyzt


!----------------------------------------------------------------
! Version of the above, for diags that are for each conv type and layer
!----------------------------------------------------------------
subroutine assign_diag_xyzlt( i_diag, i_draft, diag )

use diag_type_mod, only: diag_type
use raise_error_mod, only: raise_fatal

implicit none

! Current value of diagnostic counter
integer, intent(in out) :: i_diag
! Indicator for updraft or downdraft
integer, intent(in) :: i_draft
! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

character(len=*), parameter :: routinename = "ASSIGN_DIAG_XYZLT"

! Set flag to request diag in 5-D in x,y,z,l,t
diag % request % x_y_z_lay_typ = .true.
! Increment diagnostic counter
i_diag = i_diag + 1

! Check not gone beyond the end of the list
if ( i_diag > n_diags_xyzlt ) then
  call raise_fatal( routinename,                                               &
         "Number of CoMorph SCM diags exceeds the length "    //               &
         "of the list.  Try increasing n_diags_xyzlt." )
end if

! Assign pointers for either updraft or downdraft accordingly
select case (i_draft)
case (1)
  scm_list_xyzlu(i_diag)%pt => diag
  diag % field_5d => scm_diags_xyzlu(:,:,:,:,:,i_diag)
case (2)
  scm_list_xyzld(i_diag)%pt => diag
  diag % field_5d => scm_diags_xyzld(:,:,:,:,:,i_diag)
end select

return
end subroutine assign_diag_xyzlt


end module comorph_diags_scm_mod
#endif
