! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          This module contains additional variables and key routines
!          for the water tracer scheme.
!
! Water tracers are held in arrays in the following order:
!   1. Normal water tracer
!   2. H218O water isotope if present
!   3. HDO water isotope if present
!   4. Non-isotopic water tracers if present
!
! Which type of water tracer passes through JULES?:
!      Isotopes     - Always
!      Non-isotopic - Either they all do or they all don't (set by user)
!      Normal       - Yes if there are other types of water tracers passing
!                     through JULES, otherwise no.
!
! This file belongs in section: Water_Tracers
!
! Code description:
! Code Owner: Please refer to the UM file CodeOwners.txt
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

module water_tracers_mod

use um_types,                only: real_umphys
use atm_fields_bounds_mod,   only: tdims
use free_tracers_inputs_mod, only: max_wtrac

! DrHook parameters
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Water tracer constant info derived type
type :: wtrac_info_type
  character(len=7) :: wt_class             ! Type of water tracer
  integer :: num                           ! Water tracer number
  real(kind=real_umphys) :: standard_ratio ! Standard ratio
  real(kind=real_umphys) :: qlimit         ! Water tracer vapour limit
  real(kind=real_umphys) :: qcf_limit      ! Water tracer ice limit
  real(kind=real_umphys) :: ch4            ! Water tracer methane value
  logical :: is_iso                        ! Is the water tracer an isotope?
end type

! Store water tracer information here
type(wtrac_info_type) :: wtrac_info(max_wtrac)

!Water tracer derived type
type :: wtrac_type
  real, pointer :: q(:,:,:) => null()
  real, pointer :: qcl(:,:,:) => null()
  real, pointer :: qcf(:,:,:) => null()
  real, pointer :: qcf2(:,:,:) => null()
  real, pointer :: qr(:,:,:) => null()
  real, pointer :: qgr(:,:,:) => null()
  real, pointer :: mv(:,:,:) => null()
  real, pointer :: mcl(:,:,:) => null()
  real, pointer :: mcf(:,:,:) => null()
  real, pointer :: mcf2(:,:,:) => null()
  real, pointer :: mr(:,:,:) => null()
  real, pointer :: mgr(:,:,:) => null()
  real, pointer :: conv_prog_dq(:,:,:) => null()
  real, allocatable :: tag_mask(:,:)
  ! Remove the next two items in the future and just use values in wtrac_info
  ! (But leave for now for use in JULES)
  character(len=7) :: wt_class
  real :: standard_ratio
end type wtrac_type

! Water tracer derived type for surface prognostic fields (which are used
! in JULES)
type :: sfc_wtrac_type
  real, pointer :: snow_surft(:,:)      => null() ! Snow on tiles
                                                  ! (land_pts,ntiles)
  real, pointer :: snow_grnd_surft(:,:) => null() ! Snow under canopy
                                                  ! on tiles
                                                  ! (land_pts,ntiles)
  real, pointer :: sice_surft(:,:,:)    => null() ! Snow layer ice mass
                                                  ! on tiles
                                                  ! (land_pts,ntiles,nsmax)
  real, pointer :: sliq_surft(:,:,:)    => null() ! Snow layer liquid mass
                                                  ! on tiles
                                                  ! (land_pts,ntiles,nsmax)
  real, pointer :: canopy_surft(:,:)    => null() ! Canopy water on tiles
                                                  ! (land_pts,ntiles)
  real, pointer :: sthu_soilt(:,:)      => null() ! Unfrozen soil moisture
                                                  ! as frac of saturation
                                                  ! (land_pts,sm_levels)
  real, pointer :: sthf_soilt(:,:)      => null() ! Frozen soil moisture
                                                  ! as frac of saturation
                                                  ! (land_pts,sm_levels)
  real, pointer :: smcl_soilt(:,:)      => null() ! Soil layer moisture
                                                  ! (land_pts,sm_levels)
  real, pointer :: sthzw_soilt(:)       => null() ! Soil layer moisture in
                                                  ! deep-zw layer
                                                  ! (land_pts)
  real, pointer :: tot_surf_runoff_gb(:)=> null() ! Surface runoff accumulated
                                                  ! over river routing timestep
                                                  ! (land_pts)
  real, pointer :: tot_sub_runoff_gb(:) => null() ! Sub-surf runoff accumulated
                                                  ! over river routing timestep
                                                  ! (land_pts)
  real, pointer :: acc_lake_evap_gb(:,:)=> null() ! Lake evaporation accumulated
                                                  ! over river routing timestep
                                                  ! (row_length,rows)
  real, pointer :: twatstor(:,:)        => null() ! Store in rivers on RIVER
                                                  ! grid
                                                  !(river_row_length,river_rows)
  real, pointer :: inlandout_atm_gb(:)  => null() ! Inland base outflow
                                                  ! (land_pts)

end type sfc_wtrac_type

! Central control on whether phase changes are limited by the amount of
! available source water tracer.

logical, parameter :: l_wtrac_limit_chg=.false.

! Central control on whether to ensure there are no negative water tracer
! ratios calculated in wtrac_calc_ratio

logical, parameter :: l_wtrac_no_neg_ratio = .false.

! Minimum water value used in water tracer ratio calculation
! (Standard ratio used for values below this minimum.)

real(kind=real_umphys), parameter :: min_q_ratio = 1.0e-18

character(len=*), private, parameter :: ModuleName = 'WATER_TRACERS_MOD'

contains

subroutine wtrac_init(wtrac)


use free_tracers_inputs_mod, only: n_wtrac, norm_class_wtrac,                  &
    noniso_class_wtrac, h218o_class_wtrac, hdo_class_wtrac, class_wtrac,       &
    qlimit_h218o_wtrac, qlimit_hdo_wtrac

use turb_diff_mod, only: qlimit

use pc2_constants_mod, only: condensate_limit

use UM_ParCore,                       only: mype

use umPrintMgr,                       only:                                    &
    PrintStatus,                                                               &
    PrStatus_Normal,                                                           &
    umPrint,                                                                   &
    umMessage


implicit none

type(wtrac_type) :: wtrac(n_wtrac)

integer :: i_wt

! Error handling
character(len=*), parameter  :: routinename='WTRAC_INIT'

!DrHook variables
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!Fill in metadata in wtrac_info structure
do i_wt=1,n_wtrac
  select case(class_wtrac(i_wt))

  case (norm_class_wtrac)
    ! Normal water tracers
    wtrac_info(i_wt)%wt_class=norm_class_wtrac
    wtrac_info(i_wt)%num=i_wt
    wtrac_info(i_wt)%standard_ratio=1.0
    wtrac_info(i_wt)%qlimit=qlimit
    wtrac_info(i_wt)%qcf_limit=condensate_limit
    wtrac_info(i_wt)%ch4=1.0
    wtrac_info(i_wt)%is_iso=.false.

  case (h218o_class_wtrac)
    ! Isotope H218O
    wtrac_info(i_wt)%wt_class=h218o_class_wtrac
    wtrac_info(i_wt)%num=i_wt
    wtrac_info(i_wt)%standard_ratio=1.0
    if (qlimit_h218o_wtrac < 0.0) then
      wtrac_info(i_wt)%qlimit=qlimit
    else
      wtrac_info(i_wt)%qlimit=qlimit_h218o_wtrac
    end if
    wtrac_info(i_wt)%qcf_limit=condensate_limit
    wtrac_info(i_wt)%ch4=1.0
    wtrac_info(i_wt)%is_iso=.true.

  case (hdo_class_wtrac)
    ! Isotope HDO
    wtrac_info(i_wt)%wt_class=hdo_class_wtrac
    wtrac_info(i_wt)%num=i_wt
    wtrac_info(i_wt)%standard_ratio=1.0
    if (qlimit_hdo_wtrac < 0.0) then
      wtrac_info(i_wt)%qlimit=qlimit
    else
      wtrac_info(i_wt)%qlimit=qlimit_hdo_wtrac
    end if
    wtrac_info(i_wt)%qcf_limit=condensate_limit
    wtrac_info(i_wt)%ch4=1.0
    wtrac_info(i_wt)%is_iso=.true.

  case (noniso_class_wtrac)
    ! Non-isotopic water tracers
    wtrac_info(i_wt)%wt_class=noniso_class_wtrac
    wtrac_info(i_wt)%num=i_wt
    wtrac_info(i_wt)%standard_ratio=0.0
    wtrac_info(i_wt)%qlimit=0.0
    wtrac_info(i_wt)%qcf_limit=1.0e-3*condensate_limit
    wtrac_info(i_wt)%ch4=0.0
    wtrac_info(i_wt)%is_iso=.false.
  end select
  ! Required by all classes
  allocate (wtrac(i_wt)%tag_mask(tdims%i_start:tdims%i_end,                    &
                                 tdims%j_start:tdims%j_end))
  wtrac(i_wt)%tag_mask(:,:)=1.0

  ! Remove this in future when no longer used in JULES
  wtrac(i_wt)%wt_class = wtrac_info(i_wt)%wt_class
  wtrac(i_wt)%standard_ratio = wtrac_info(i_wt)%standard_ratio

end do

!print out info if required.

if (PrintStatus >= PrStatus_Normal .and. mype == 0) then
  do i_wt=1,n_wtrac
    write(umMessage,'(A,I0,A,A7,A,I0,A,E13.5,A,E13.5,A,E13.5,A,E13.5)')        &
       'Water tracer #',i_wt,                                                  &
       ' wt_class:',wtrac_info(i_wt)%wt_class,                                 &
       ' num:',wtrac_info(i_wt)%num,                                           &
       ' standard_ratio:',wtrac_info(i_wt)%standard_ratio,                     &
       ' qlimit:',wtrac_info(i_wt)%qlimit,                                     &
       ' qcf_limit:',wtrac_info(i_wt)%qcf_limit,                               &
       ' ch4: ',wtrac_info(i_wt)%ch4
    call umPrint(umMessage,src=routinename)
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine wtrac_init


! ----------------------------------------------


end module water_tracers_mod
