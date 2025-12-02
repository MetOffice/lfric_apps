! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_conv_mod

use um_types, only: real_umphys

implicit none

! Description:
! Water tracer (working array) structures used in convection code plus
! related routines
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

logical :: l_wtrac_conv = .false.    ! Controls water tracer calcs in convection

! Structure to hold water tracer fields that are used in ni_conv_ctl
! NOTE: this structure differs from the others as the arrays have an
! extra dimension of the number of water tracers.  This makes it easier
! to pass the fields in glue_conv.

type :: conv_ni_wtrac_type

  real(kind=real_umphys), allocatable :: q_inc(:,:,:,:)   ! Conv increments
  real(kind=real_umphys), allocatable :: qcl_inc(:,:,:,:) ! Conv increments
  real(kind=real_umphys), allocatable :: qcf_inc(:,:,:,:) ! Conv increments
  real(kind=real_umphys), allocatable :: q_conv(:,:,:,:)   ! Local copy of q
  real(kind=real_umphys), allocatable :: qcl_conv(:,:,:,:) ! Local copy of qcl
  real(kind=real_umphys), allocatable :: qcf_conv(:,:,:,:) ! Local copy of qcf
  real(kind=real_umphys), allocatable :: qrain_conv(:,:,:,:)  ! Copy of qrain
  real(kind=real_umphys), allocatable :: qgraup_conv(:,:,:,:) ! Copy of qgraup
  real(kind=real_umphys), allocatable :: qcf2_conv(:,:,:,:)   ! Copy of qcf2
  real(kind=real_umphys), allocatable :: dqbydt(:,:,:,:)
  real(kind=real_umphys), allocatable :: dqclbydt(:,:,:,:)
  real(kind=real_umphys), allocatable :: dqcfbydt(:,:,:,:)
  real(kind=real_umphys), allocatable :: dq_add(:,:,:,:)   ! Added q to ensure
                                                           !    q > min value
  real(kind=real_umphys), allocatable :: it_rain(:,:,:)    ! Copy of conv_rain
  real(kind=real_umphys), allocatable :: it_snow(:,:,:)    ! Copy of conv_snow

end type

! Structure to hold water tracer fields for the different convection types
! for the environment calculations.
! These are set in glue_conv and have size n_dp, n_sh or n_md.

type :: conv_e_wtrac_type

  real(kind=real_umphys), allocatable :: q(:,:)
                                  ! Water tracer specific humidity (kg/kg)
  real(kind=real_umphys), allocatable :: qcl(:,:)
                                  ! Water tracer liquid condensate (kg/kg)
  real(kind=real_umphys), allocatable :: qcf(:,:)
                                  ! Water tracer ice condensate (kg/kg)
  real(kind=real_umphys), allocatable :: dqbydt(:,:)
                                  ! Rate of change of water tracer q (kg/kg/s)
  real(kind=real_umphys), allocatable :: dqclbydt(:,:)
                                  ! Rate of change of water tracer qcl (kg/kg/s)
  real(kind=real_umphys), allocatable :: dqcfbydt(:,:)
                                  ! Rate of change of water traer qcf (kg/kg/s)
  real(kind=real_umphys), allocatable :: rain(:)
                                  ! Water tracer surface rainfall (kg/m2/s)
  real(kind=real_umphys), allocatable :: snow(:)
                                  ! Water tracer surface snowfall (kg/m2/s)

end type

! Structure to hold water tracer fields for the parcel calculations
! These are set in deep_conv, shallow_conv or mid_conv and have size n_dp,
! n_sh or n_md.

type :: conv_p_wtrac_type

  real(kind=real_umphys), allocatable :: q(:,:)
                               ! Water tracer parcel specific humidity (kg/kg)
  real(kind=real_umphys), allocatable :: qcl(:,:)
                               ! Water tracer parcel liquid condensate (kg/kg)
  real(kind=real_umphys), allocatable :: qcf(:,:)
                                ! Water tracer parcel ice condensate (kg/kg)
  real(kind=real_umphys), allocatable :: qpert(:)
                                ! Water tracer content of q parcel perturbation
  real(kind=real_umphys), allocatable :: precip(:,:)
                                ! Water tracer precip from each layer (kg/m2/s)
  real(kind=real_umphys), allocatable :: Qlkp1(:)
                                ! Amount of water tracer in condensation
                                ! to liquid water in the parcel (kg/kg)
  real(kind=real_umphys), allocatable :: Qfkp1(:)
                                ! Amount of water tracer in deposition to
                                ! in ice water in the parcel (kg/kg)
  real(kind=real_umphys), allocatable :: Frezkp1(:)
                                ! Amount of water tracer in freezing from
                                ! liq to ice water in the parcel (kg/kg)
end type

! Structure to hold water tracer fields for the downdraught calculations
! These are set in dd_all_call and dd_call and have size ndd.

type :: conv_dd_wtrac_type

  real(kind=real_umphys), allocatable :: precip_k(:)
                                ! Water tracer precip added when
                                ! descending from layer K to K-1 (kg/m**2/s)
  real(kind=real_umphys), allocatable :: q_k(:)
                                ! Parcel water tracer mix ratio of layer
                                ! k (kg/kg) (i.e. updraught value)
  real(kind=real_umphys), allocatable :: dqbydt_k(:)
                                ! Increment to water tracer mixing ratio
                                ! of layer k (kg/kg/s)
  real(kind=real_umphys), allocatable :: qdd_k(:)
                                ! Water tracer mixing ratio of downdraught
                                ! in layer K (kg/kg)
  real(kind=real_umphys), allocatable :: delqd(:)
                                ! Water tracer content of moistening
                                ! necessary to achieve saturation (kg/kg)
  real(kind=real_umphys), allocatable :: qe_k(:)
                                ! Environment water tracer mixing ratio
                                ! of layer k (kg/kg)
  real(kind=real_umphys), allocatable :: qe_km1(:)
                                ! Environment water tracer mixing ratio
                                ! of layer k-1 (kg/kg)
  real(kind=real_umphys), allocatable :: rain(:)
                                ! Water tracer sfc rainfall (kg/m**2/s)
  real(kind=real_umphys), allocatable :: snow(:)
                                ! Water tracer sfc snowfall (kg/m**2/s)
  real(kind=real_umphys), allocatable :: rain_dd(:)
                                ! Amount of water tracer rainfall passing
                                ! through downdraught (kg/m**2/s)
  real(kind=real_umphys), allocatable :: snow_dd(:)
                                ! Amount of water tracer snowfall passing
                                ! through downdraught (kg/m**2/s)
end type

! Structure to hold water tracer fields used to calculate phase change and
! evaporation for precip falling through environment part of the
! downdraught calculations and also below the cloud base.
! These are set in:
! 1) dd_all_call and dd_call and have size ndd
! 2) evap_bcb_nodd with size n_nodd
! 3) evap_bcb_nodd_all with size nsofar.

type :: conv_ev_wtrac_type

  real(kind=real_umphys), allocatable :: dqbydt_km1(:)
                                     ! Increment to water tracer mixing
                                     ! ratio of layer k-1(kg/kg/s)
  real(kind=real_umphys), allocatable :: rain_env(:)
                                     ! Amount of water tracer rainfall passing
                                     ! through environment (kg/m**2/s)
  real(kind=real_umphys), allocatable :: snow_env(:)
                                     ! Amount of water tracer snowfall passing
                                     ! through environment (kg/m**2/s)
end type


! Structure to hold water tracer fields following the 2nd compression
! for the downdraught calculations
! These are set in downd and have size nddon.

type :: conv_dd2_wtrac_type

  real(kind=real_umphys), allocatable :: qdd_k(:)
                                     ! Water tracer content in downdraught
                                     ! at layer k compressed (kg/kg)
  real(kind=real_umphys), allocatable :: qe_k(:)
                                     ! Water tracer content of environment
                                     ! layer k compressed (kg/kg)
  real(kind=real_umphys), allocatable :: qe_km1(:)
                                     ! Water tracer content of environment
                                     ! layer k-1 compressed (kg/kg)
  real(kind=real_umphys), allocatable :: dqbydt_k(:)
                                     ! Increment to water tracer in
                                     ! layer k compressed (kg/kg)
  real(kind=real_umphys), allocatable :: dqbydt_km1(:)
                                     ! Increment to water tracer in
                                     ! layer k-1 compressed (kg/kg)
  real(kind=real_umphys), allocatable :: delqd(:)
                                     ! Water tracer content of moistening
                                     ! necessary to achieve saturation
                                     ! compressed (kg/kg)
  real(kind=real_umphys), allocatable :: rain(:)
                                     ! Water tracer downdraught rainfall
                                     ! compressed (kg/m**2/S)
  real(kind=real_umphys), allocatable :: snow(:)
                                     ! Water tracer downdraught snowfall
                                     ! compressed (kg/m**2/S)
end type

character(len=*), parameter, private :: ModuleName='WTRAC_CONV_MOD'

contains

! Subroutine Interface:
subroutine wtrac_alloc_conv_ni(row_length, rows, nlev,                         &
                               n_wtrac, n_wtrac_conv,                          &
                               wtrac_as, wtrac_ni)
!
! Description:
! Allocate arrays used by water tracers in the individual convection schemes
!

use mphys_inputs_mod,        only: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
use wtrac_atm_step_mod,      only: atm_step_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: row_length        ! Row length
integer, intent(in) :: rows              ! No. of rows
integer, intent(in) :: nlev              ! No. of vertical levels USED in
                                         ! CONVECTION
integer, intent(in) :: n_wtrac           ! No. of water tracers
integer, intent(in) :: n_wtrac_conv      ! No. of water tracers USED in
                                         ! CONVECTION

! Water tracer structure (containing 'star' fields)
type(atm_step_wtrac_type), intent(in) :: wtrac_as(n_wtrac)

type(conv_ni_wtrac_type), intent(in out) :: wtrac_ni

! Local variables
integer :: i,j,k,i_wt      ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_NI'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_wtrac_conv) then
  allocate(wtrac_ni%q_inc(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qcl_inc(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qcf_inc(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%q_conv(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qcl_conv(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qcf_conv(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qrain_conv(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qgraup_conv(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%qcf2_conv(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%dqbydt(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%dqclbydt(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%dqcfbydt(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%dq_add(row_length,rows,nlev,n_wtrac_conv))
  allocate(wtrac_ni%it_rain(row_length,rows,n_wtrac_conv))
  allocate(wtrac_ni%it_snow(row_length,rows,n_wtrac_conv))
else
  allocate(wtrac_ni%q_inc(1,1,1,1))
  allocate(wtrac_ni%qcl_inc(1,1,1,1))
  allocate(wtrac_ni%qcf_inc(1,1,1,1))
  allocate(wtrac_ni%q_conv(1,1,1,1))
  allocate(wtrac_ni%qcl_conv(1,1,1,1))
  allocate(wtrac_ni%qcf_conv(1,1,1,1))
  allocate(wtrac_ni%qrain_conv(1,1,1,1))
  allocate(wtrac_ni%qgraup_conv(1,1,1,1))
  allocate(wtrac_ni%qcf2_conv(1,1,1,1))
  allocate(wtrac_ni%dqbydt(1,1,1,1))
  allocate(wtrac_ni%dqclbydt(1,1,1,1))
  allocate(wtrac_ni%dqcfbydt(1,1,1,1))
  allocate(wtrac_ni%dq_add(1,1,1,1))
  allocate(wtrac_ni%it_rain(1,1,1))
  allocate(wtrac_ni%it_snow(1,1,1))
end if

! Initialise fields
if (l_wtrac_conv) then

  do i_wt = 1, n_wtrac_conv
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(i_wt, nlev, rows, row_length, wtrac_ni) private(k, i, j)
    do k= 1, nlev
      do j = 1, rows
        do i = 1, row_length
          wtrac_ni%dqbydt(i,j,k,i_wt)   = 0.0
          wtrac_ni%dqclbydt(i,j,k,i_wt) = 0.0
          wtrac_ni%dqcfbydt(i,j,k,i_wt) = 0.0
          wtrac_ni%q_inc(i,j,k,i_wt)    = 0.0
          wtrac_ni%qcl_inc(i,j,k,i_wt)  = 0.0
          wtrac_ni%qcf_inc(i,j,k,i_wt)  = 0.0
          wtrac_ni%dq_add(i,j,k,i_wt)   = 0.0
        end do
      end do
    end do
!$OMP end PARALLEL do
  end do

  do i_wt = 1, n_wtrac_conv
!$OMP PARALLEL DEFAULT(none) private(k, i, j)                                  &
!$OMP SHARED(i_wt, nlev, rows, row_length, wtrac_ni, wtrac_as, l_mcr_qrain,    &
!$OMP l_mcr_qgraup, l_mcr_qcf2)
!$OMP do SCHEDULE(STATIC)
    do k= 1, nlev
      do j = 1, rows
        do i = 1, row_length
          wtrac_ni%q_conv(i,j,k,i_wt)   = wtrac_as(i_wt)%q_star(i,j,k)
          wtrac_ni%qcl_conv(i,j,k,i_wt) = wtrac_as(i_wt)%qcl_star(i,j,k)
          wtrac_ni%qcf_conv(i,j,k,i_wt) = wtrac_as(i_wt)%qcf_star(i,j,k)
        end do ! i
      end do ! j
    end do ! k
!$OMP end do NOWAIT

    if (l_mcr_qrain) then
!$OMP do SCHEDULE(STATIC)
      do k = 1, nlev
        do j = 1, rows
          do i = 1, row_length
            wtrac_ni%qrain_conv(i,j,k,i_wt) = wtrac_as(i_wt)%qr_star(i,j,k)
          end do ! i
        end do ! j
      end do ! k
!$OMP end do NOWAIT
    end if

    if (l_mcr_qgraup) then
!$OMP do SCHEDULE(STATIC)
      do k = 1, nlev
        do j = 1, rows
          do i = 1, row_length
            wtrac_ni%qgraup_conv(i,j,k,i_wt) = wtrac_as(i_wt)%qgr_star(i,j,k)
          end do ! i
        end do ! j
      end do ! k
!$OMP end do NOWAIT
    end if

    if (l_mcr_qcf2) then
!$OMP do SCHEDULE(STATIC)
      do k = 1, nlev
        do j = 1, rows
          do i = 1, row_length
            wtrac_ni%qcf2_conv(i,j,k,i_wt) = wtrac_as(i_wt)%qcf2_star(i,j,k)
          end do ! i
        end do ! j
      end do ! k
!$OMP end do NOWAIT
    end if

!$OMP end PARALLEL

  end do  ! i_wt
end if ! l_wtrac_conv

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_ni

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_ni(wtrac_ni)

!
! Description:
! Deallocate arrays used by water tracers in ni_conv_ctl
!
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
type(conv_ni_wtrac_type), intent(in out) :: wtrac_ni

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_NI'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

deallocate(wtrac_ni%it_snow)
deallocate(wtrac_ni%it_rain)
deallocate(wtrac_ni%dq_add)
deallocate(wtrac_ni%dqcfbydt)
deallocate(wtrac_ni%dqclbydt)
deallocate(wtrac_ni%dqbydt)
deallocate(wtrac_ni%qcf2_conv)
deallocate(wtrac_ni%qgraup_conv)
deallocate(wtrac_ni%qrain_conv)
deallocate(wtrac_ni%qcf_conv)
deallocate(wtrac_ni%qcl_conv)
deallocate(wtrac_ni%q_conv)
deallocate(wtrac_ni%qcf_inc)
deallocate(wtrac_ni%qcl_inc)
deallocate(wtrac_ni%q_inc)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_ni

! --------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_conv_e(np_c, nlev, n_wtrac, conv_e_wtrac)
!
! Description:
!   Allocate arrays used by water tracers for the environment
!   calculations in the convection scheme.
!
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: np_c
integer, intent(in) :: nlev
integer, intent(in) :: n_wtrac

type(conv_e_wtrac_type), intent(in out) :: conv_e_wtrac(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_E'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    allocate(conv_e_wtrac(i_wt)%q(np_c,nlev))
    allocate(conv_e_wtrac(i_wt)%qcl(np_c,nlev))
    allocate(conv_e_wtrac(i_wt)%qcf(np_c,nlev))
    allocate(conv_e_wtrac(i_wt)%dqbydt(np_c,nlev))
    allocate(conv_e_wtrac(i_wt)%dqclbydt(np_c,nlev))
    allocate(conv_e_wtrac(i_wt)%dqcfbydt(np_c,nlev))
    allocate(conv_e_wtrac(i_wt)%rain(np_c))
    allocate(conv_e_wtrac(i_wt)%snow(np_c))
  end do
else
  allocate(conv_e_wtrac(1)%q(1,1))
  allocate(conv_e_wtrac(1)%qcl(1,1))
  allocate(conv_e_wtrac(1)%qcf(1,1))
  allocate(conv_e_wtrac(1)%dqbydt(1,1))
  allocate(conv_e_wtrac(1)%dqclbydt(1,1))
  allocate(conv_e_wtrac(1)%dqcfbydt(1,1))
  allocate(conv_e_wtrac(1)%rain(1))
  allocate(conv_e_wtrac(1)%snow(1))
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_e

! -----------------------------------------------------------------
! Subroutine Interface:
subroutine wtrac_dealloc_conv_e(n_wtrac, conv_e_wtrac)

!
! Description:
!   Deallocate arrays used by water tracers for the environment
!   calculations in the convection scheme.
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: n_wtrac

type(conv_e_wtrac_type), intent(in out) :: conv_e_wtrac(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_E'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  deallocate(conv_e_wtrac(i_wt)%snow)
  deallocate(conv_e_wtrac(i_wt)%rain)
  deallocate(conv_e_wtrac(i_wt)%dqcfbydt)
  deallocate(conv_e_wtrac(i_wt)%dqclbydt)
  deallocate(conv_e_wtrac(i_wt)%dqbydt)
  deallocate(conv_e_wtrac(i_wt)%qcf)
  deallocate(conv_e_wtrac(i_wt)%qcl)
  deallocate(conv_e_wtrac(i_wt)%q)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_e

! --------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_conv_p(np_c, nlev, n_wtrac, wtrac_p)

!
! Description:
!   Allocate arrays used by water tracers for the parcel calculations in
!   the convection scheme.
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: np_c
integer, intent(in) :: nlev
integer, intent(in) :: n_wtrac

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_P'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (l_wtrac_conv) then
  do i_wt = 1, n_wtrac
    allocate(wtrac_p(i_wt)%q(np_c,nlev))
    allocate(wtrac_p(i_wt)%qcl(np_c,nlev))
    allocate(wtrac_p(i_wt)%qcf(np_c,nlev))
    allocate(wtrac_p(i_wt)%qpert(np_c))
    allocate(wtrac_p(i_wt)%precip(np_c,nlev))
    allocate(wtrac_p(i_wt)%Qlkp1(np_c))
    allocate(wtrac_p(i_wt)%Qfkp1(np_c))
    allocate(wtrac_p(i_wt)%Frezkp1(np_c))
  end do
else
  do i_wt = 1, n_wtrac
    allocate(wtrac_p(i_wt)%q(1,1))
    allocate(wtrac_p(i_wt)%qcl(1,1))
    allocate(wtrac_p(i_wt)%qcf(1,1))
    allocate(wtrac_p(i_wt)%qpert(1))
    allocate(wtrac_p(i_wt)%precip(1,1))
    allocate(wtrac_p(i_wt)%Qlkp1(1))
    allocate(wtrac_p(i_wt)%Qfkp1(1))
    allocate(wtrac_p(i_wt)%Frezkp1(1))
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_p

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_p(n_wtrac, wtrac_p)

!
! Description:
!   Deallocate arrays used by water tracers for parcel calculations in the
!   convection scheme.
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: n_wtrac

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_P'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  deallocate(wtrac_p(i_wt)%Frezkp1)
  deallocate(wtrac_p(i_wt)%Qfkp1)
  deallocate(wtrac_p(i_wt)%Qlkp1)
  deallocate(wtrac_p(i_wt)%precip)
  deallocate(wtrac_p(i_wt)%qpert)
  deallocate(wtrac_p(i_wt)%qcf)
  deallocate(wtrac_p(i_wt)%qcl)
  deallocate(wtrac_p(i_wt)%q)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_p

! --------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_conv_dd(np_c, n_wtrac, wtrac_dd)

!
! Description:
!   Allocate arrays used by water tracers in the downdraught part of the
!   convection scheme
!   Note, this routine is only called if  l_wtrac_conv=T.  This is because it
!   was found that allocating these arrays (to size 1) caused an increase in
!   run time in model runs without water tracers. This approach is acceptable
!   as the whole structure, rather than individual arrays within the structure,
!   are passed into lower-level routines.
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: np_c
integer, intent(in) :: n_wtrac

type(conv_dd_wtrac_type), intent(in out) :: wtrac_dd(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_DD'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  allocate(wtrac_dd(i_wt)%precip_k(np_c))
  allocate(wtrac_dd(i_wt)%q_k(np_c))
  allocate(wtrac_dd(i_wt)%dqbydt_k(np_c))
  allocate(wtrac_dd(i_wt)%qdd_k(np_c))
  allocate(wtrac_dd(i_wt)%delqd(np_c))
  allocate(wtrac_dd(i_wt)%qe_k(np_c))
  allocate(wtrac_dd(i_wt)%qe_km1(np_c))
  allocate(wtrac_dd(i_wt)%rain(np_c))
  allocate(wtrac_dd(i_wt)%snow(np_c))
  allocate(wtrac_dd(i_wt)%rain_dd(np_c))
  allocate(wtrac_dd(i_wt)%snow_dd(np_c))
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_dd

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_dd(n_wtrac, wtrac_dd)

!
! Description:
!   Deallocate arrays used by water tracers in the downdraught calculations
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: n_wtrac

type(conv_dd_wtrac_type), intent(in out) :: wtrac_dd(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_DD'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  deallocate(wtrac_dd(i_wt)%snow_dd)
  deallocate(wtrac_dd(i_wt)%rain_dd)
  deallocate(wtrac_dd(i_wt)%snow)
  deallocate(wtrac_dd(i_wt)%rain)
  deallocate(wtrac_dd(i_wt)%qe_km1)
  deallocate(wtrac_dd(i_wt)%qe_k)
  deallocate(wtrac_dd(i_wt)%delqd)
  deallocate(wtrac_dd(i_wt)%qdd_k)
  deallocate(wtrac_dd(i_wt)%dqbydt_k)
  deallocate(wtrac_dd(i_wt)%q_k)
  deallocate(wtrac_dd(i_wt)%precip_k)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_dd
! --------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_conv_ev(np_c, n_wtrac, wtrac_ev)

!
! Description:
!   Allocate arrays used by water tracers in the phase change and
!   evaporation of precip in environment or below cloud base calculations
!   Note, this routine is only called if  l_wtrac_conv=T.  This is because it
!   was found that allocating these arrays (to size 1) caused an increase in
!   run time in model runs without water tracers. This approach is acceptable
!   as the whole structure, rather than individual arrays within the structure,
!   are passed into lower-level routines.
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: np_c
integer, intent(in) :: n_wtrac

type(conv_ev_wtrac_type), intent(in out) :: wtrac_ev(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_EV'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


do i_wt = 1, n_wtrac
  allocate(wtrac_ev(i_wt)%dqbydt_km1(np_c))
  allocate(wtrac_ev(i_wt)%rain_env(np_c))
  allocate(wtrac_ev(i_wt)%snow_env(np_c))
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_ev

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_ev(n_wtrac, wtrac_ev)

!
! Description:
!   Deallocate arrays used by water tracers in the phase change and
!   evaporation of precip in environment or below cloud base calculations
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: n_wtrac

type(conv_ev_wtrac_type), intent(in out) :: wtrac_ev(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_EV'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  deallocate(wtrac_ev(i_wt)%snow_env)
  deallocate(wtrac_ev(i_wt)%rain_env)
  deallocate(wtrac_ev(i_wt)%dqbydt_km1)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_ev

! --------------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_alloc_conv_dd2(np_c, n_wtrac, wtrac_dd2)

!
! Description:
!   Allocate 2nd compression arrays used by water tracers in the downdraught
!   part of the convection scheme
!   Note, this routine is only called if  l_wtrac_conv=T.  This is because it
!   was found that allocating these arrays (to size 1) caused an increase in
!   run time in model runs without water tracers. This approach is acceptable
!   as the whole structure, rather than individual arrays within the structure,
!   are passed into lower-level routines.
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: np_c
integer, intent(in) :: n_wtrac

type(conv_dd2_wtrac_type), intent(in out) :: wtrac_dd2(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_ALLOC_CONV_DD2'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  allocate(wtrac_dd2(i_wt)%qdd_k(np_c))
  allocate(wtrac_dd2(i_wt)%qe_k(np_c))
  allocate(wtrac_dd2(i_wt)%qe_km1(np_c))
  allocate(wtrac_dd2(i_wt)%dqbydt_k(np_c))
  allocate(wtrac_dd2(i_wt)%dqbydt_km1(np_c))
  allocate(wtrac_dd2(i_wt)%delqd(np_c))
  allocate(wtrac_dd2(i_wt)%rain(np_c))
  allocate(wtrac_dd2(i_wt)%snow(np_c))
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_alloc_conv_dd2

! -----------------------------------------------------------------

! Subroutine Interface:
subroutine wtrac_dealloc_conv_dd2(n_wtrac, wtrac_dd2)

!
! Description:
!   Deallocate 2nd compression arrays used by water tracers in the
!   downdraught calculations
!

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: n_wtrac

type(conv_dd2_wtrac_type), intent(in out) :: wtrac_dd2(n_wtrac)

! Local variables
integer :: i_wt      ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_DEALLOC_CONV_DD2'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  deallocate(wtrac_dd2(i_wt)%snow)
  deallocate(wtrac_dd2(i_wt)%rain)
  deallocate(wtrac_dd2(i_wt)%delqd)
  deallocate(wtrac_dd2(i_wt)%dqbydt_km1)
  deallocate(wtrac_dd2(i_wt)%dqbydt_k)
  deallocate(wtrac_dd2(i_wt)%qe_km1)
  deallocate(wtrac_dd2(i_wt)%qe_k)
  deallocate(wtrac_dd2(i_wt)%qdd_k)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_dealloc_conv_dd2

end module wtrac_conv_mod
