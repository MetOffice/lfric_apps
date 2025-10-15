!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Random seed generator

module random_seed_gen_alg_mod

  use constants_mod,       only: r_def, i_def, str_def

  implicit none

  private

  public random_seed_gen_alg

 contains
  !>@brief Generation of a reproducible random seed for use by ensemble runs.
  !>@details Use full model date in randomising seed to reduce chance of
  !> recycling from one run to the next. The formula calculates the days since
  !> ~2000AD and adds in time suitably inflated to fully change the seed.
  !> Only use last two digits of year to prevent numerical overflow at some date
  !> in the future. A random number generated from this seed is used to
  !> multiply the seed again.
  !>@param[out] iranseed      : Random seed controlled by ensemble member
  !>                            number for stochastic physics
  !>@param[in] ensemble_number: Allows different ensemble members to have
  !>                            different random number seeds

  subroutine random_seed_gen_alg(iranseed, ensemble_number)
    use xios, only: xios_date, xios_get_current_date, &
                    xios_date_get_day_of_year, xios_date_get_second_of_day

    implicit none

    integer(i_def), intent(out) :: iranseed(:)
    integer(i_def), intent(in) :: ensemble_number

    type(xios_date) :: datetime
    integer(i_def) :: i
    integer(i_def) :: year, month, day, utc_shift, hour, minute, &
                      milli_ensemble_number
    integer(i_def) :: random_seed_size, iarg, max_iarg
    integer(i_def), allocatable :: prevseed(:)
    real(r_def), allocatable :: rnum(:)

    ! Use the datetime from the XIOS clock, should eventually be replaced by
    ! calls to the model clock.
    call xios_get_current_date(datetime)

    ! Values required for the random seed function.
    year = int(datetime%year, i_def)
    month = int(datetime%month, i_def)
    day = int(datetime%day, i_def)
    hour = int(datetime%hour, i_def) + 1
    minute = int(datetime%minute, i_def)
    milli_ensemble_number = ensemble_number + 100

    ! Fetch random seed array size from intrinsic and allocate arrays.
    random_seed_size = size(iranseed)
    if(.not. allocated(prevseed)) allocate(prevseed(random_seed_size))
    if(.not. allocated(rnum)) allocate(rnum(random_seed_size))

    ! Take only the last 2 digits of the year.
    year = year - 100 * nint(0.01_r_def * year)

    ! Random seed function that is both reproducible and can be recycled between
    ! runs. The formulation used here is ported directly from the UM code and
    ! calculates the date in Julian calendar form
    ! (Crossref DOI link: https://doi.org/10.1145/364096.364097) then adds some
    ! additional variation based on the ensemble number, hours, and minutes.
    iarg = int((day - 32075 + 1461*(year + 4800 + (month - 14)/12)/4             &
           + 367 * (month - 2 - (month - 14)/12*12)/12 - 3*((year + 4900         &
           + (month - 14)/12)/100)/4)*1000 + milli_ensemble_number**2.86_r_def   &
           + ensemble_number**3.79_r_def + hour**5.12_r_def + minute**3.24_r_def)

    ! Constrain iarg in a range to prevent numerical overflow.
    max_iarg = floor(sqrt(real(huge(iarg), r_def)))
    iarg = mod(iarg, max_iarg)
    iarg = max(iarg, 256)
    prevseed = iarg
    ! Generate initial random number set to use.
    call random_seed(put = prevseed(1:random_seed_size))
    call random_number(rnum)

    ! Range of seed from 0 to 2**31 (32-bit Int).
    iranseed = int(iarg*rnum)

  end subroutine random_seed_gen_alg

end module random_seed_gen_alg_mod
