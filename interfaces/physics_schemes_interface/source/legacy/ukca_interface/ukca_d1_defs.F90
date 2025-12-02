! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to define fields from D1 needed for UKCA
!  these are set in ukca_setd1defs
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code Description:
!  Language:  Fortran 2003
!  This code is written to UMDP3 programming standards.
!
! ######################################################################
!
module ukca_d1_defs

!!!! Temporary use of ukca_fieldname_mod. Should use ukca_api_mod when
!!!! ukca_api_mod is no longer dependent on ukca_d1_defs.
use ukca_fieldname_mod, only: ukca_maxlen_fieldname => maxlen_fieldname

implicit none

!     No of Prognostics and Diagnostics required
integer, save :: Nukca_D1items      ! Size of UkcaD1codes array
integer, save :: n_ntp              ! No of non-transported prognostics
integer, save :: intp_first         ! item for 1st non-transported prognostic
integer, save :: istrat_first       ! item for 1st strat flux diag
integer, save :: imode_first        ! item for 1st MODE diag
integer, parameter :: ukca_sect=34       ! stash section for UKCA
integer, parameter :: ukca_glomap_diag_sect=38
                                         ! stash section for UKCA GLOMAP diags
integer, parameter :: ukca_diag_sect=50  ! stash section for UKCA diags
integer, parameter :: ukca_s34_plev=51   ! stash section 34 on pressure levels
integer, parameter :: ukca_s50_plev=52   ! stash section 50 on pressure levels
integer, parameter :: item1_stratflux  = 1    ! 1st item No. Strat fluxes (s38)
integer, parameter :: item1_mode_diags = 201  ! 1st item No. for MODE
                                              ! diags (s38)
integer, parameter :: item1_plume_diags = 900 ! 1st item No. for plume
                                              ! scavenging diags (s38)

!  Diagnostic range for nitrate scheme
!  includes the CMIP diagnostics but not the PM or plume scav diags
integer, parameter :: item1_nitrate_diags = 578
integer, parameter :: item1_nitrate_noems = 584
integer, parameter :: itemN_nitrate_diags = 668

!  Diagnostic range for dust 3 insoluble mode scheme and ccn diags
!  includes the CMIP diagnostics but not the PM or plume scav diags
integer, parameter :: item1_dust3mode_diags = 675
integer, parameter :: itemN_dust3mode_diags = 701

!  Diagnostic range for microplastic scheme
!  includes the CMIP diagnostics but not the PM or plume scav diags
integer, parameter :: item1_microplastic_diags = 706
integer, parameter :: itemN_microplastic_diags = 761

type :: code
  integer :: section        ! section code
  integer :: item           ! item code
  integer :: n_levels       ! number of levels
  integer :: address        ! address in D1
  integer :: length         ! length of field
  integer :: halo_type      ! field halo type
  integer :: grid_type      ! grid type
  integer :: field_type     ! field grid type
  integer :: len_dim1       ! length of array dim1
  integer :: len_dim2       ! length of array dim2
  integer :: len_dim3       ! length of array dim3
  logical :: prognostic     ! prognostic t/f
  logical :: required       ! t/f
  character(len=ukca_maxlen_fieldname) :: fieldname ! UKCA field name
end type code

!     Number of tracers, emissions and diagnostics, set according
!     to GUI choices in ukca_setd1defs

integer, save   :: nmax_plume_diags = 43 ! Max no plume scavenging diags
integer, save   :: n_strat_fluxdiags  ! No. strat flux diags
integer, save   :: n_MODE_diags       ! No. diagnostics for MODE

!     list of prognostic/diagnostics from/to D1

type(code), allocatable, save :: UkcaD1Codes(:)

!     Some of the fields extracted from D1 are potentially processed on the UM
!     side of the UKCA interface instead of/ as well as by UKCA.
!     The logicals below are used to track the actual requirements.

!     The D1 field frac_types may be required by UKCA and/or may be required
!     by the UM for use in determining an active tile indicator input to UKCA.
logical, save :: l_ukca_use_frac_types = .false.
logical, save :: l_um_use_frac_types = .false.

!     The D1 fields net_surf_sw and tot_surf_sw may be required by the UM
!     for calculating UKCA's input surface albedo.
logical, save :: l_um_use_surf_sw = .false.

!     The D1 field rho_r2 may be required by UKCA and/or may be required by
!     the UM for plume scavenging diagnostics calculations.
logical, save :: l_ukca_use_rho_r2 = .false.
logical, save :: l_um_use_rho_r2 = .false.

!     The D1 field exner_theta_levels may be required by UKCA and/or
!     may be required by the UM for post-processing of UKCA diagnostics.
logical, save :: l_ukca_use_exner_theta_levels = .false.
logical, save :: l_um_use_exner_theta_levels = .false.

!      The flux logicals are now turned on via STASH panel requests.
logical, save :: L_ukca_stratflux = .false.  ! T for stratospheric fluxes
logical, save :: L_ukca_mode_diags= .false.  ! T for MODE diags.
logical, save :: l_ukca_plume_diags= .false. ! T for plume scavenging diags.

! STASH Codes of species used for  coupling of UKCA with radiation scheme
! and sulphur cycle of UM

!! Note: to BE UPDATED FOR any SEC 34 ITEM CHANGES     !!

! species used for feedback to radiation scheme
integer, parameter :: i_ukca_grg_o3  =  1, i_ukca_grg_ch4 = 9 ,                &
                      i_ukca_grg_n2o = 49, i_ukca_grg_f11 = 55,                &
                      i_ukca_grg_f12 = 56, i_ukca_grg_f113= 64,                &
                      i_ukca_grg_h22 = 65, i_ukca_grg_cf3chf2 = -1,            &
                      i_ukca_grg_chf2chf2 = -1

! Item numbers for UKCA oxidants used in the CLASSIC sulphur cycle (i.e. O3,
! HONO2, H2O2, OH, HO2)
! Set in atmos_ukca_setup - address of OH and HO2 varies with solver since
! these are non-transported prognostics when used with the B-E solver.
! Item numbers for OH and HO2 and OH when non-transported prognostics are
! specified as parameters to ensure consistency between atmos_ukca_setup and
! the ntp_name2stashitem function.
integer, parameter :: i_oh_be  = 995
integer, parameter :: i_ho2_be = 993
integer, save :: ukca_item_sulpc(5)

character(len=*), parameter, private :: ModuleName='UKCA_D1_DEFS'

contains

! ----------------------------------------------------------------------
integer function stash2ntpindex(item)
! ----------------------------------------------------------------------
! Description:
! Given a STASH item number look up where it is in the ntp array.
! All NTPs are assumed to be in the UKCA section.
! Abort here if failed to look up.
! ----------------------------------------------------------------------

use um_stashcode_mod,        only: stashcode_ukca_sec
use ereport_mod,             only: ereport
use errormessagelength_mod,  only: errormessagelength
use parkind1,                only: jprb, jpim
use yomhook,                 only: lhook, dr_hook

implicit none

! Function argument
integer, intent(in)        :: item   ! Item number in UKCA section to look up

! Local variables

integer :: i                                         ! Loop counter
integer :: j                                         ! UkcaD1Codes index
integer :: errcode                                   ! error code
character(len=errormessagelength)   :: cmessage      ! Error return message

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName = 'STASH2NTPINDEX'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set a default value. If it is still this on exit, the search has failed.
stash2ntpindex = -999

! Search all entries in ntp array to find the one we want
ntp_loop: do i = 1, n_ntp
  j = intp_first + i - 1
  if (UkcaD1Codes(j)%section == stashcode_ukca_sec .and.                       &
      UkcaD1Codes(j)%item == item) then
    stash2ntpindex = i
    exit ntp_loop
  end if
end do ntp_loop

! If stash2ntpindex is -999 then call ereport
if (stash2ntpindex == -999) then
  write(cmessage,'(A,2I6,A)') 'Failed to find: ', stashcode_ukca_sec, item,    &
    ' in UkcaD1Codes structure.'
  errcode = 1
  call ereport(RoutineName,errcode,cmessage)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return

end function stash2ntpindex
! ----------------------------------------------------------------------

end module ukca_d1_defs
