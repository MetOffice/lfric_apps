!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Routines for calculating semi-Lagrangian interpolation coefficients.
!!
!> @details This module contains functions and subroutines which allow the
!!          semi-Lagrangian transport scheme to pre-compute interpolation
!!          coefficients, and to apply monotonicity to the transported field.
!------------------------------------------------------------------------------

module sl_support_mod

  use constants_mod,                      only: i_def, r_tran
  use transport_enumerated_types_mod,     only: vertical_monotone_none,           &
                                                vertical_monotone_strict,         &
                                                vertical_monotone_relaxed,        &
                                                vertical_monotone_order_constant, &
                                                vertical_monotone_order_linear,   &
                                                vertical_monotone_order_high
  use vertical_mass_remapping_kernel_mod, only: local_point_1d_array

  implicit none

  private

  public :: monotone_cubic_sl
  public :: monotone_quintic_sl
  public :: compute_cubic_coeffs
  public :: compute_quintic_coeffs
  public :: compute_cubic_hermite_coeffs

  contains

  !-------------------------------------------------------------------------------
  !> @brief   Limits the cubic semi-Lagrangian advective transport.
  !> @details This subroutine identifies the non-monotone cubic-interpolated values
  !!          and modifies them using the specified order/option of the modification.
  !> @param[inout] fi        The interpolated field
  !> @param[in]    f         The grid-data field
  !> @param[in]    sc        The cubic-stencil used for interpolation for each point
  !> @param[in]    cl        Linear weights
  !> @param[in]    vertical_monotone
  !!                         Option to identify non-monotone points
  !> @param[in]    vertical_monotone_order
  !!                         Option of the modification for non-monotone points
  !> @param[in]    ns        First index of data-in
  !> @param[in]    nf        End index of data-in
  !-------------------------------------------------------------------------------
  subroutine monotone_cubic_sl(fi,f,sc,cl,vertical_monotone,vertical_monotone_order,ns,nf)

    implicit none

    ! Arguments
    integer(kind=i_def),                     intent(in)    :: vertical_monotone
    integer(kind=i_def),                     intent(in)    :: vertical_monotone_order
    integer(kind=i_def),                     intent(in)    :: ns
    integer(kind=i_def),                     intent(in)    :: nf
    real(kind=r_tran),   dimension(ns:nf),   intent(inout) :: fi(ns:nf)
    real(kind=r_tran),   dimension(ns:nf),   intent(in)    :: f(ns:nf)
    integer(kind=i_def), dimension(4,ns:nf), intent(in)    :: sc
    real(kind=r_tran),   dimension(2,ns:nf), intent(in)    :: cl

    ! Local variables
    integer(kind=i_def), dimension(ns:nf) :: no_mono_id
    integer(kind=i_def) :: k
    real(kind=r_tran)   :: test1, test2, test3, test4
    real(kind=r_tran)   :: minv, maxv, sigma, xi

    ! Identify points/values that are non-monotone using
    ! either the strict or relaxed conditions
    ! if no_mono_id(:) /= 0, fi is non-monotone
    !    no_mono_id(:) == 1, flat-left  for mod_high
    !    no_mono_id(:) == 2, flat-right for mod_high

    no_mono_id(:) = 0_i_def

    select case( vertical_monotone )

    case ( vertical_monotone_strict )
      do k = ns, nf
        test1 = (fi(k)-f(sc(2,k)))*(f(sc(3,k))-fi(k))
        test2 = (fi(k)-f(sc(2,k)))*(f(sc(3,k))-f(sc(2,k)))
        if ( test1 < 0.0_r_tran ) then
          if (test2  < 0.0_r_tran ) then
            no_mono_id(k) = 1_i_def
          else
            no_mono_id(k) = 2_i_def
          end if
        end if
      end do
    case ( vertical_monotone_relaxed )
      do k = ns, nf
        test1 = (fi(k)-f(sc(2,k)))*(f(sc(3,k))-fi(k))
        test2 = (f(sc(2,k))-f(sc(1,k)))*(f(sc(3,k))-f(sc(4,k)))
        test3 = (fi(k)-f(sc(2,k)))*(f(sc(2,k))-f(sc(1,k)))
        test4 = (fi(k)-f(sc(2,k)))*(f(sc(3,k))-f(sc(2,k)))
        if ( test1 < 0.0_r_tran .and. (test2 < 0.0_r_tran .or. test3 < 0.0_r_tran) ) then
          if (test4  < 0.0_r_tran ) then
            no_mono_id(k) = 1_i_def
          else
            no_mono_id(k) = 2_i_def
          end if
        end if
      end do
    end select

    ! Perform modification for identified no-mono points/values

    select case( vertical_monotone_order )

    case ( vertical_monotone_order_constant )
      do k = ns, nf
        if ( no_mono_id(k) /= 0_i_def ) then
          minv = min(f(sc(2,k)), f(sc(3,k)))
          maxv = max(f(sc(2,k)), f(sc(3,k)))
          fi(k) = max( minv, min(fi(k), maxv) )
        end if
      end do
    case ( vertical_monotone_order_linear )
      do k = ns, nf
        if ( no_mono_id(k) /= 0_i_def ) then
          fi(k) = cl(1,k)*f(sc(2,k)) + cl(2,k)*f(sc(3,k))
        end if
      end do
    case ( vertical_monotone_order_high )
      do k = ns, nf
        if (no_mono_id(k) == 1_i_def ) then
          sigma = f(sc(3,k)) -  f(sc(2,k))
          xi = cl(2,k)
          fi(k) =  sigma*xi**3 + f(sc(2,k))
        elseif (no_mono_id(k) == 2_i_def ) then
          sigma = f(sc(3,k)) -  f(sc(2,k))
          xi = cl(1,k)
          fi(k) = -sigma*xi**3 + f(sc(3,k))
        end if
      end do
    end select

  end subroutine monotone_cubic_sl

  !-------------------------------------------------------------------------------
  !> @brief   Limits the quintic semi-Lagrangian advective transport.
  !> @details This subroutine identifies the non-monotone quintic-interpolated values
  !!          and modifies them using the specified order/option of the modification.
  !> @param[inout] fi        The interpolated field
  !> @param[in]    f         The grid-data field
  !> @param[in]    sq        The quintic-stencil used for interpolation for each point
  !> @param[in]    cl        Linear weights
  !> @param[in]    vertical_monotone
  !!                         Option to identify non-monotone points
  !> @param[in]    vertical_monotone_order
  !!                         Option of the modification for non-monotone points
  !> @param[in]    ns        Start index of data-in
  !> @param[in]    nf        End index of data-in
  !-------------------------------------------------------------------------------
  subroutine monotone_quintic_sl(fi,f,sq,cl,vertical_monotone,vertical_monotone_order,ns,nf)

    implicit none

    ! Arguments
    integer(kind=i_def),                     intent(in)    :: vertical_monotone
    integer(kind=i_def),                     intent(in)    :: vertical_monotone_order
    integer(kind=i_def),                     intent(in)    :: ns
    integer(kind=i_def),                     intent(in)    :: nf
    real(kind=r_tran),   dimension(ns:nf),   intent(inout) :: fi(ns:nf)
    real(kind=r_tran),   dimension(ns:nf),   intent(in)    :: f(ns:nf)
    integer(kind=i_def), dimension(6,ns:nf), intent(in)    :: sq
    real(kind=r_tran),   dimension(2,ns:nf), intent(in)    :: cl

    ! Local variables
    integer(kind=i_def), dimension(ns:nf) :: no_mono_id
    integer(kind=i_def) :: k
    real(kind=r_tran)   :: test1, test2, test3, test4
    real(kind=r_tran)   :: minv, maxv, sigma, xi

    ! Identify points/values that are non-monotone using
    ! either the strict or relaxed conditions
    ! if no_mono_id(:) # 0, fi is non-monotone
    !    no_mono_id(:) = 1, flat-left  for mod_high
    !    no_mono_id(:) = 2, flat-right for mod_high

    no_mono_id(:) = 0_i_def

    select case( vertical_monotone )

    case ( vertical_monotone_strict )
      do k = ns, nf
        test1 = (fi(k)-f(sq(3,k)))*(f(sq(4,k))-fi(k))
        test2 = (fi(k)-f(sq(3,k)))*(f(sq(4,k))-f(sq(3,k)))
        if ( test1 < 0.0_r_tran ) then
          if (test2  < 0.0_r_tran ) then
            no_mono_id(k) = 1_i_def
          else
            no_mono_id(k) = 2_i_def
          end if
        end if
      end do
    case ( vertical_monotone_relaxed )
      do k = ns, nf
        test1 = (fi(k)-f(sq(3,k)))*(f(sq(4,k))-fi(k))
        test2 = (f(sq(3,k))-f(sq(2,k)))*(f(sq(4,k))-f(sq(5,k)))
        test3 = (fi(k)-f(sq(3,k)) )*(f(sq(3,k))-f(sq(2,k)))
        test4 = (fi(k)-f(sq(3,k)) )*(f(sq(4,k))-f(sq(3,k)))
        if ( test1 < 0.0_r_tran .and. (test2 < 0.0_r_tran .or. test3 < 0.0_r_tran) ) then
          if ( test4  < 0.0_r_tran ) then
            no_mono_id(k) = 1_i_def
          else
            no_mono_id(k) = 2_i_def
          end if
        end if
      end do
    end select

    ! Perform modification for identified no-mono points/values

    select case( vertical_monotone_order )

    case ( vertical_monotone_order_constant )
      do k = ns, nf
        if ( no_mono_id(k) /= 0_i_def ) then
          minv = min(f(sq(3,k)), f(sq(4,k)))
          maxv = max(f(sq(3,k)), f(sq(4,k)))
          fi(k) = max( minv, min(fi(k), maxv) )
        end if
      end do
    case ( vertical_monotone_order_linear )
      do k = ns, nf
        if ( no_mono_id(k) /= 0_i_def ) then
          fi(k) = cl(1,k)*f(sq(3,k)) + cl(2,k)*f(sq(4,k))
        end if
      end do
    case ( vertical_monotone_order_high )
      do k = ns, nf
        if ( no_mono_id(k) == 1_i_def ) then
          sigma = f(sq(4,k)) -  f(sq(3,k))
          xi = cl(2,k)
          fi(k) =  sigma*xi**5 + f(sq(3,k))
        elseif ( no_mono_id(k) == 2_i_def ) then
          sigma = f(sq(4,k)) -  f(sq(3,k))
          xi = cl(1,k)
          fi(k) = -sigma*xi**5 + f(sq(4,k))
        end if
      end do
    end select

  end subroutine monotone_quintic_sl

  !-------------------------------------------------------------------------------
  !> @brief   This subroutine computes the cubic-Lagrange weights.
  !> @details Compute the cubic-Lagrange interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  zi  The interpolation points
  !> @param[in]  zg  The grid points where the data is located
  !> @param[out] sc  The cubic-stencil used for interpolation for each point
  !> @param[out] cc  The cubic interpolation weights
  !> @param[out] cl  Linear weights
  !> @param[in]  nzi The size of interpolation points
  !> @param[in]  nzg The size of the grid-data
  !-------------------------------------------------------------------------------
  subroutine compute_cubic_coeffs(zi,zg,dz,sc,cc,cl,nzi,nzg)

    implicit none

    ! Arguments
    integer(kind=i_def),                   intent(in)  :: nzi, nzg
    real(kind=r_tran),   dimension(nzi),   intent(in)  :: zi
    real(kind=r_tran),   dimension(nzg),   intent(in)  :: zg
    real(kind=r_tran),   dimension(nzg),   intent(in)  :: dz
    real(kind=r_tran),   dimension(4,nzi), intent(out) :: cc
    integer(kind=i_def), dimension(4,nzi), intent(out) :: sc
    real(kind=r_tran),   dimension(2,nzi), intent(out) :: cl

    ! Local variables
    real(kind=r_tran)   :: z1, z2, z3, z4, xi
    real(kind=r_tran)   :: d1, d2, d3, d4
    real(kind=r_tran)   :: n1, n2, n3, n4
    integer(kind=i_def) :: k, km

    do k = 1, nzi
      km = local_point_1d_array(zi(k), zg, 1, nzg)
      xi = (zi(k) - zg(km))/dz(km)

      sc(1,k) = max(1 , km - 1 )
      sc(2,k) = km
      sc(3,k) = min(nzg, km + 1 )
      sc(4,k) = min(nzg, km + 2 )

      z1 = 0.0_r_tran
      z2 = z1 + dz(sc(1,k))
      z3 = z2 + dz(sc(2,k))
      z4 = z3 + dz(sc(3,k))
      xi = z2 + xi*(z3-z2)

      d1 = (z1-z2)*(z1-z3)*(z1-z4)
      d2 = (z2-z1)*(z2-z3)*(z2-z4)
      d3 = (z3-z1)*(z3-z2)*(z3-z4)
      d4 = (z4-z1)*(z4-z2)*(z4-z3)

      n1 = (xi-z2)*(xi-z3)*(xi-z4)
      n2 = (xi-z1)*(xi-z3)*(xi-z4)
      n3 = (xi-z1)*(xi-z2)*(xi-z4)
      n4 = (xi-z1)*(xi-z2)*(xi-z3)

      ! cubic weights
      cc(1,k) = n1/d1
      cc(2,k) = n2/d2
      cc(3,k) = n3/d3
      cc(4,k) = n4/d4

      ! linear weights
      cl(1,k) = (z3-xi)/(z3-z2)
      cl(2,k) = 1.0_r_tran - cl(1,k)

      ! Next to boundaries there are not enough points for cubic
      ! so revert to linear

      if( sc(1,k) == sc(2,k) .or. sc(3,k) == sc(4,k) ) then
        cc(1,k) = 0.0_r_tran
        cc(2,k) = cl(1,k)
        cc(3,k) = cl(2,k)
        cc(4,k) = 0.0_r_tran
      end if
    end do

  end subroutine compute_cubic_coeffs

  !-------------------------------------------------------------------------------
  !> @breif   This subroutine computes cubic-Hermite weights.
  !> @details Compute the cubic-Hermite interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  zi  The interpolation points
  !> @param[in]  zg  The grid points where the data is located
  !> @param[out] sc  The cubic-stencil used for interpolation for each point
  !> @param[out] cc  The cubic interpolation weights
  !> @param[out] cl  Linear weights
  !> @param[in]  nzi The size of interpolation points
  !> @param[in]  nzg The size of the grid-data
  !-------------------------------------------------------------------------------
  subroutine compute_cubic_hermite_coeffs(zi,zg,dz,sc,cc,cl,nzi,nzg)

    implicit none

    ! Arguments
    integer(kind=i_def),                   intent(in)  :: nzi, nzg
    real(kind=r_tran),   dimension(nzi),   intent(in)  :: zi
    real(kind=r_tran),   dimension(nzg),   intent(in)  :: zg
    real(kind=r_tran),   dimension(nzg),   intent(in)  :: dz
    real(kind=r_tran),   dimension(4,nzi), intent(out) :: cc
    integer(kind=i_def), dimension(4,nzi), intent(out) :: sc
    real(kind=r_tran),   dimension(2,nzi), intent(out) :: cl

    ! Local variables
    real(kind=r_tran)   :: xi, alfa, beta, inv_1p_alfa, inv_1p_beta
    real(kind=r_tran)   :: xip2, xip3, c1, c2, c3, c4
    integer(kind=i_def) :: k, km

    do k = 1, nzi
      km = local_point_1d_array(zi(k), zg, 1, nzg)
      xi = (zi(k) - zg(km))/dz(km)

      sc(1,k) = max(1 , km - 1 )
      sc(2,k) = km
      sc(3,k) = min(nzg, km + 1 )
      sc(4,k) = min(nzg, km + 2 )

      alfa = dz(sc(1,k))/dz(sc(2,k))
      beta = dz(sc(3,k))/dz(sc(2,k))
      inv_1p_alfa=1.0_r_tran/(1.0_r_tran + alfa)
      inv_1p_beta=1.0_r_tran/(1.0_r_tran + beta)
      xip2 = xi**2
      xip3 = xi**3
      c1 = 2.0_r_tran*xip3 - 3.0_r_tran*xip2 + 1.0_r_tran
      c2 = 1.0_r_tran - c1
      c3 = xip3 - 2.0_r_tran*xip2 + xi
      c4 = xip3 - xip2

      ! Cubic-hermite weights
      if( sc(1,k) == sc(2,k) ) then
        cc(1,k) = 0.0
        cc(2,k) = c1 - c3 - c4*inv_1p_beta
        cc(3,k) = c2 + c3
        cc(4,k) = c4*inv_1p_beta
      elseif ( sc(3,k) == sc(4,k) ) then
        cc(1,k) = -c3*inv_1p_alfa
        cc(2,k) = c1 - c4
        cc(3,k) = c2 + c4 + c3*inv_1p_alfa
        cc(4,k) = 0.0
      else
        cc(1,k) = -c3*inv_1p_alfa
        cc(2,k) = c1 - c4*inv_1p_beta
        cc(3,k) = c2 + c3*inv_1p_alfa
        cc(4,k) = c4*inv_1p_beta
      end if

      ! linear weights
      cl(2,k) = xi
      cl(1,k) = 1.0_r_tran - cl(2,k)
    end do

  end subroutine compute_cubic_hermite_coeffs

  !-------------------------------------------------------------------------------
  !> @details This subroutine computes the quintic-Lagrange weights.
  !> @details Compute the quintic-Lagrange interpolation weights for vertical
  !!          semi-Lagrangian advective transport.
  !> @param[in]  zi        The interpolation points
  !> @param[in]  zg        The grid points where the data is located
  !> @param[out] sq        The quintic-stencil used for interpolation for each point
  !> @param[out] cq        The quintic interpolation weights
  !> @param[out] cl        Linear weights
  !> @param[in]  nzi       The size of interpolation points
  !> @param[in]  nzg       The size of the grid-data
  !-------------------------------------------------------------------------------
  subroutine compute_quintic_coeffs(zi,zg,dz,sq,cq,cl,nzi,nzg)

    implicit none

    ! Arguments
    integer(kind=i_def),                   intent(in)  :: nzi, nzg
    real(kind=r_tran),   dimension(nzi),   intent(in)  :: zi
    real(kind=r_tran),   dimension(nzg),   intent(in)  :: zg
    real(kind=r_tran),   dimension(nzg),   intent(in)  :: dz
    real(kind=r_tran),   dimension(6,nzi), intent(out) :: cq
    integer(kind=i_def), dimension(6,nzi), intent(out) :: sq
    real(kind=r_tran),   dimension(2,nzi), intent(out) :: cl

    ! Local variables
    real(kind=r_tran)   :: z1, z2, z3, z4, z5, z6, xi
    real(kind=r_tran)   :: d1, d2, d3, d4, d5, d6
    real(kind=r_tran)   :: n1, n2, n3, n4, n5, n6
    integer(kind=i_def) :: k, km

    do k = 1, nzi
      km = local_point_1d_array(zi(k), zg, 1, nzg)
      xi = (zi(k) - zg(km))/dz(km)

      sq(1,k) = max(1 , km - 2 )
      sq(2,k) = max(1 , km - 1 )
      sq(3,k) = km
      sq(4,k) = min(nzg, km + 1 )
      sq(5,k) = min(nzg, km + 2 )
      sq(6,k) = min(nzg, km + 3 )

      z1 = 0.0_r_tran
      z2 = z1 + dz(sq(1,k))
      z3 = z2 + dz(sq(2,k))
      z4 = z3 + dz(sq(3,k))
      z5 = z4 + dz(sq(4,k))
      z6 = z5 + dz(sq(5,k))
      xi = z3 + xi*(z4-z3)

      d1 = (z1-z2)*(z1-z3)*(z1-z4)*(z1-z5)*(z1-z6)
      d2 = (z2-z1)*(z2-z3)*(z2-z4)*(z2-z5)*(z2-z6)
      d3 = (z3-z1)*(z3-z2)*(z3-z4)*(z3-z5)*(z3-z6)
      d4 = (z4-z1)*(z4-z2)*(z4-z3)*(z4-z5)*(z4-z6)
      d5 = (z5-z1)*(z5-z2)*(z5-z3)*(z5-z4)*(z5-z6)
      d6 = (z6-z1)*(z6-z2)*(z6-z3)*(z6-z4)*(z6-z5)

      n1 = (xi-z2)*(xi-z3)*(xi-z4)*(xi-z5)*(xi-z6)
      n2 = (xi-z1)*(xi-z3)*(xi-z4)*(xi-z5)*(xi-z6)
      n3 = (xi-z1)*(xi-z2)*(xi-z4)*(xi-z5)*(xi-z6)
      n4 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z5)*(xi-z6)
      n5 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z4)*(xi-z6)
      n6 = (xi-z1)*(xi-z2)*(xi-z3)*(xi-z4)*(xi-z5)

      ! quintic weights
      cq(1,k) = n1/d1
      cq(2,k) = n2/d2
      cq(3,k) = n3/d3
      cq(4,k) = n4/d4
      cq(5,k) = n5/d5
      cq(6,k) = n6/d6

      ! linear weights
      cl(1,k) = (z4-xi)/(z4-z3)
      cl(2,k) = 1.0_r_tran - cl(1,k)

      if( sq(1,k) == sq(2,k) .or. sq(5,k) == sq(6,k) ) then
        ! Revert to cubic weights
        d1 = (z2-z3)*(z2-z4)*(z2-z5)
        d2 = (z3-z2)*(z3-z4)*(z3-z5)
        d3 = (z4-z2)*(z4-z3)*(z4-z5)
        d4 = (z5-z2)*(z5-z3)*(z5-z4)

        n1 = (xi-z3)*(xi-z4)*(xi-z5)
        n2 = (xi-z2)*(xi-z4)*(xi-z5)
        n3 = (xi-z2)*(xi-z3)*(xi-z5)
        n4 = (xi-z2)*(xi-z3)*(xi-z4)

        cq(1,k) = 0.0_r_tran
        cq(2,k) = n1/d1
        cq(3,k) = n2/d2
        cq(4,k) = n3/d3
        cq(5,k) = n4/d4
        cq(6,k) = 0.0_r_tran
      end if

      if( sq(2,k) == sq(3,k) .or. sq(4,k) == sq(5,k) ) then
        ! Revert to linear weights
        cq(1,k) = 0.0_r_tran
        cq(2,k) = 0.0_r_tran
        cq(3,k) = cl(1,k)
        cq(4,k) = cl(2,k)
        cq(5,k) = 0.0_r_tran
        cq(6,k) = 0.0_r_tran
      end if
    end do

  end subroutine compute_quintic_coeffs

end module sl_support_mod