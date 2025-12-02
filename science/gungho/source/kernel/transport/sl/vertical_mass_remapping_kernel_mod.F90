!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a vertical mass remapping of tracer-mass at W3.
!! @details The SLICE scheme (using different order reconstructions) is
!!          used to perform a vertical mass remapping of tracer mass.

module vertical_mass_remapping_kernel_mod

use argument_mod,         only : arg_type,              &
                                 GH_FIELD, GH_SCALAR,   &
                                 GH_REAL, GH_INTEGER,   &
                                 GH_READWRITE, GH_READ, &
                                 GH_WRITE, CELL_COLUMN, GH_LOGICAL
use fs_continuity_mod,    only : W3, W2v
use constants_mod,        only : r_tran, i_def, EPS_R_TRAN, l_def
use kernel_mod,           only : kernel_type
! TODO #3011: these config options should be passed through as arguments
use transport_config_mod, only : slice_order_constant, slice_order_linear, &
                                 slice_order_parabola, slice_order_cubic
use transport_enumerated_types_mod, only : monotone_none,                    &
                                           monotone_strict,                  &
                                           monotone_relaxed,                 &
                                           vertical_monotone_order_constant, &
                                           vertical_monotone_order_linear,   &
                                           vertical_monotone_order_high
implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed
!>                                      by the PSy layer.
type, public, extends(kernel_type) :: vertical_mass_remapping_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                      &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v), & ! departure points
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),  & ! mass
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     W2v), & ! mass flux
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! order of construction
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! monotone scheme
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),           & ! monotone order
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)            & ! enforce min-val
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vertical_mass_remapping_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
 public :: vertical_mass_remapping_code
 public :: local_point_1d_array
private :: piecewise_reconstruction
private :: remapp_mass
private :: monotone_quadratic_coeffs
private :: monotone_cubic_coeffs
private :: monotone_edge_values

contains

  !-------------------------------------------------------------------------------
  !> @details This kernel remaps mass from one grid to another, for
  !!          vertical motion only, then interpolates theta at the
  !!          departure point using 1D cubic-Lagrange interpolation.
  !> @param[in]     nlayers      The number of layers
  !> @param[in]     dep_dist_z   The vertical departure distance used for SL advection
  !> @param[in,out] mass         The mass field that needs to be remapped
  !> @param[in,out] flux         The mass flux corresponding to the transport
  !> @param[in]     order        Order of the reconstruction of underlying function
  !> @param[in] vertical_monotone        The monotone scheme to be used
  !> @param[in] vertical_monotone_order  The order of monotone reconstruction option
  !> @param[in] enforce_min_value        Option to enforce a minimum value
  !> @param[in]     ndf_w2v      The number of degrees of freedom per cell
  !!                             on W2v space
  !> @param[in]     undf_w2v     The number of unique degrees of freedom
  !!                             on W2v space for this partition
  !> @param[in]     map_w2v      The dofmap for the cell at the base of the column
  !!                             on W2v space
  !> @param[in]     ndf_w3       The number of degrees of freedom per cell
  !!                             on w3 space
  !> @param[in]     undf_w3      The number of unique degrees of freedom
  !!                             on w3 space
  !> @param[in]     map_w3       The dofmap for the cell at the base of the column
  !!                             on w3 space
  !-------------------------------------------------------------------------------

  subroutine vertical_mass_remapping_code( nlayers,                     &
                                           dep_dist_z,                  &
                                           mass,                        &
                                           flux,                        &
                                           order,                       &
                                           vertical_monotone,           &
                                           vertical_monotone_order,     &
                                           enforce_min_value,           &
                                           ndf_w2v, undf_w2v, map_w2v,  &
                                           ndf_w3, undf_w3, map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                         :: nlayers
  integer(kind=i_def), intent(in)                         :: ndf_w3
  integer(kind=i_def), intent(in)                         :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)      :: map_w3
  integer(kind=i_def), intent(in)                         :: ndf_w2v
  integer(kind=i_def), intent(in)                         :: undf_w2v
  integer(kind=i_def), dimension(ndf_w2v), intent(in)     :: map_w2v

  real(kind=r_tran), dimension(undf_w2v), intent(in)      :: dep_dist_z
  real(kind=r_tran), dimension(undf_w3),  intent(inout)   :: mass
  real(kind=r_tran), dimension(undf_w2v), intent(inout)   :: flux
  integer(kind=i_def), intent(in)                         :: order
  logical(kind=l_def), intent(in)  :: enforce_min_value
  integer(kind=i_def), intent(in)  :: vertical_monotone,       &
                                      vertical_monotone_order
  ! Local variables
  integer(kind=i_def)                     :: k, nz, nzl
  real(kind=r_tran), dimension(nlayers+1) :: dist, zl, zd, m_flux
  real(kind=r_tran), dimension(nlayers)   :: dz, m0, mn
  real(kind=r_tran), dimension(nlayers)   :: a0, a1, a2, a3

  nz  = nlayers
  nzl = nlayers + 1_i_def

  ! Extract and fill local column from global data
  do k=0,nlayers
    dist(k+1) = dep_dist_z(map_w2v(1)+k)
  end do
  do k=0,nlayers - 1
    m0(k+1) = mass(map_w3(1)+k)
  end do

  do k=1, nzl
    zl(k) = real(k-1,r_tran)
  end do
  do k=1, nzl
    zd(k) = zl(k) - dist(k)
    zd(k) = min( max(zl(1),zd(k)), zl(nzl) )
  end do
  do k=1,nz
    dz(k)  = 1.0_r_tran
  end do

  call piecewise_reconstruction(zl,dz,m0,a0,a1,a2,a3,order,1,nz,            &
                                vertical_monotone, vertical_monotone_order, &
                                enforce_min_value                           )
  call remapp_mass(zd,zl,dz,m0,a0,a1,a2,a3,1,nz,mn,m_flux )

  ! Export the local new mass back to the global mass

  do k=0,nlayers - 1
    mass(map_w3(1)+k) = mn(k+1)
    flux(map_w2v(1)+k) = m_flux(k+1)
  end do

  flux(map_w2v(1)+nlayers) = m_flux(nlayers+1)

  end subroutine vertical_mass_remapping_code

  !-------------------------------------------------------------------------------
  !> @brief   This routine returns the cubic coefficents (a0,a1,a2,a3) for each cell.
  !> @details These coefficents define a piece-wise function for each cell. This
  !>          function is defined such that the integral of the function over the
  !>          cell gives the mass of that cell.
  !> @param[in]  xl     The grid-edges of the given mass (mass(i) is the mass of
  !!                    cell [xl(i),xl(i+1)])
  !> @param[in]  dx     Size of the cells, dx(i)=xl(i+1)-xl(i)
  !> @param[in]  mass   The mass of grid cells
  !> @param[out] a0     The first coefficent of the cell reconstruction.
  !>                    Each cell has a cubic function f(x) = a0 + a1*x + a2*x^2 + a3*x^3
  !> @param[out] a1     The second coefficent of the cell reconstruction.
  !> @param[out] a2     The third  coefficent of the cell reconstruction.
  !> @param[out] a3     The fourth coefficent of the cell reconstruction.
  !> @param[in] order   Order the piecewise function (constant, linear, quadratic, or cubic)
  !> @param[in] ns      Index of the starting cell
  !> @param[in] nf      Index of the end cell
  !> @param[in] vertical_monotone        The monotone scheme option to be used
  !> @param[in] vertical_monotone_order  Order of sub-cell monotone reconstruction
  !> @param[in] enforce_min_value        Option to enforcce a minimum value
  !-------------------------------------------------------------------------------
  subroutine piecewise_reconstruction(xl,dx,mass,a0,a1,a2,a3,order,ns,nf,          &
                                      vertical_monotone, vertical_monotone_order,  &
                                      enforce_min_value)
   implicit none

   integer(kind=i_def), intent(in)                  :: ns, nf, order
   integer(kind=i_def), intent(in)                  :: vertical_monotone,  &
                                                       vertical_monotone_order
   logical(kind=l_def), intent(in)                  :: enforce_min_value
   real(kind=r_tran), dimension(ns:nf+1), intent(in) :: xl
   real(kind=r_tran), dimension(ns:nf), intent(in)   :: dx, mass
   real(kind=r_tran), dimension(ns:nf), intent(out)  :: a0, a1, a2, a3

   ! Local variables
   real(kind=r_tran), dimension(ns:nf) :: rho_left_cv, rho_right_cv, slope_rho
   real(kind=r_tran), dimension(ns:nf) :: y_centre_cv, rho_bar
   integer(kind=i_def), parameter :: min_cvs_required = 4_i_def
   real(kind=r_tran), dimension(ns:nf+1) :: rho_left_cv_all,slope_im,slope_ip
   integer(kind=i_def) :: j, j0, j1, cv_start
   real(kind=r_tran)    :: y0, y1, y2, y3, y4, y, m1, m2, m3, m4
   real(kind=r_tran), parameter :: c2=2.0_r_tran, c3=3.0_r_tran, c4=4.0_r_tran,  &
                                  c6=6.0_r_tran, c9=9.0_r_tran, c12=12.0_r_tran

    j0 = int(min_cvs_required/2,i_def)
    j1 = min_cvs_required - 1_i_def

    do j = ns, nf
      y_centre_cv(j) = 0.5_r_tran * ( xl(j) + xl(j+1) )
      a1(j) = 0.0_r_tran
      a2(j) = 0.0_r_tran
      a3(j) = 0.0_r_tran
      rho_bar(j) = mass(j) / dx(j)
    end do

   ! Use piecewise-constant for less than 4 cells
    if ( ((nf - ns + 1 ) < min_cvs_required) .or.    &
                   (order == slice_order_constant) ) then
      do j = ns, nf
        a0(j) = rho_bar(j)
      end do

    else

      do j = ns, nf + 1
        cv_start = min( max( j - j0, ns ) + j1, nf ) - j1
        y0 = xl(cv_start)
        y1 = xl(cv_start+1) - y0
        y2 = xl(cv_start+2) - y0
        y3 = xl(cv_start+3) - y0
        y4 = xl(cv_start+4) - y0
        m1 =   mass(cv_start)
        m2 =   mass(cv_start+1) + m1
        m3 =   mass(cv_start+2) + m2
        m4 =   mass(cv_start+3) + m3
        m1 = m1 / ( y1*(y1-y2)*(y1-y3)*(y1-y4) )
        m2 = m2 / ( y2*(y2-y1)*(y2-y3)*(y2-y4) )
        m3 = m3 / ( y3*(y3-y1)*(y3-y2)*(y3-y4) )
        m4 = m4 / ( y4*(y4-y1)*(y4-y2)*(y4-y3) )

        y  = xl(j) - y0
        rho_left_cv_all(j) =                                                             &
             ( y*( y*(c4*y - c3*(y2+y3+y4)) + c2*(y3*y2+y2*y4+y3*y4) ) - y3*y2*y4 )*m1 + &
             ( y*( y*(c4*y - c3*(y1+y3+y4)) + c2*(y1*y3+y1*y4+y3*y4) ) - y1*y3*y4 )*m2 + &
             ( y*( y*(c4*y - c3*(y2+y4+y1)) + c2*(y2*y4+y1*y2+y1*y4) ) - y1*y2*y4 )*m3 + &
                  (y*(y*(c4*y-c3*(y2+y3+y1))+c2*(y3*y2+y1*y2+y1*y3))-y1*y2*y3)*m4
        y  = y_centre_cv( max( ns, j - 1 ) ) - y0
        slope_im(j) =                                                           &
                  ( y*(c12*y - c6*(y2+y3+y4)) + c2*(y3*y2+y2*y4+y3*y4) )*m1 +   &
                  ( y*(c12*y - c6*(y1+y3+y4)) + c2*(y1*y3+y1*y4+y3*y4) )*m2 +   &
                  ( y*(c12*y - c6*(y2+y4+y1)) + c2*(y2*y4+y1*y2+y1*y4) )*m3 +   &
                  ( y*(c12*y - c6*(y2+y3+y1)) + c2*(y3*y2+y1*y2+y1*y3) )*m4
        y  = y_centre_cv( min( j, nf    ) ) - y0
        slope_ip(j) =                                                           &
                  ( y*(c12*y - c6*(y2+y3+y4)) + c2*(y3*y2+y2*y4+y3*y4) )*m1 +   &
                  ( y*(c12*y - c6*(y1+y3+y4)) + c2*(y1*y3+y1*y4+y3*y4) )*m2 +   &
                  ( y*(c12*y - c6*(y2+y4+y1)) + c2*(y2*y4+y1*y2+y1*y4) )*m3 +   &
                  ( y*(c12*y - c6*(y2+y3+y1)) + c2*(y3*y2+y1*y2+y1*y3) )*m4
      end do


      if ( (vertical_monotone /= monotone_none) .or. enforce_min_value ) then
         call monotone_edge_values(rho_left_cv_all,rho_bar,dx, &
                   vertical_monotone,enforce_min_value, ns, nf )
      end if

      do j = ns, nf
        slope_rho(j)   = 0.5_r_tran * ( slope_im(j+1) + slope_ip(j) )
        slope_rho(j)   = slope_rho(j) * dx(j)
        rho_left_cv(j) = rho_left_cv_all(j)
        rho_right_cv(j)= rho_left_cv_all(j+1)
        !
        ! If rho_bar is outside [rho_left_cv,rho_right_cv] use constant
        ! In this case a constant is the only valid representation
        !
        m1 = (rho_bar(j)-rho_left_cv(j))*(rho_right_cv(j)-rho_bar(j))
        if (m1 < 0.0_r_tran) then
           rho_left_cv(j)  = rho_bar(j)
           rho_right_cv(j) = rho_bar(j)
           slope_rho(j)    = 0.0_r_tran
        end if
      end do

      select case( order )
      case ( slice_order_linear )
         !
         ! For linear we satisfy the mass/integral and one value left or right
         ! depending on whether the mass is close to the left or right values
         !
         do j = ns, nf
           if ( abs(rho_bar(j)-rho_left_cv(j)) < abs(rho_right_cv(j)-rho_bar(j)) ) then
             a0(j) = rho_left_cv(j)
             a1(j) = c2*(rho_bar(j)-rho_left_cv(j))
           else
             a0(j) = c2*rho_bar(j)-rho_right_cv(j)
             a1(j) = c2*(rho_right_cv(j)-rho_bar(j))
           end if
           m1 = (rho_bar(j)-rho_left_cv(j))*(rho_right_cv(j)-rho_bar(j))
           if (m1 < 0.0_r_tran) then
               a0(j) = rho_bar(j)
               a1(j) = 0.0_r_tran
           end if
         end do
      case ( slice_order_parabola )
         do j = ns, nf
            a0(j) = rho_left_cv(j)
            a1(j) = -c4*rho_left_cv(j) - c2*rho_right_cv(j) + c6*rho_bar(j)
            a2(j) =  c3*rho_left_cv(j) + c3*rho_right_cv(j) - c6*rho_bar(j)
         end do
         if ( vertical_monotone /= monotone_none )  then
            call monotone_quadratic_coeffs(rho_left_cv, rho_right_cv, rho_bar,    &
                                           a0,a1,a2,ns,nf,vertical_monotone_order )
         end if
      case ( slice_order_cubic )
         do j = ns, nf
            a0(j) = rho_left_cv(j)
            a1(j) = -c6*rho_left_cv(j) + c6*rho_bar(j)   - c2*slope_rho(j)
            a2(j) =  c9*rho_left_cv(j) - c3*rho_right_cv(j)  &
                   - c6*rho_bar(j)     + c6*slope_rho(j)
            a3(j) = -c4*rho_left_cv(j) + c4*rho_right_cv(j) - c4*slope_rho(j)
         end do
         if ( vertical_monotone /= monotone_none )  then
            call monotone_cubic_coeffs(rho_left_cv, rho_right_cv, rho_bar,       &
                                       a0,a1,a2,a3,ns,nf,vertical_monotone_order )
         end if
      end select

    end if

  end subroutine piecewise_reconstruction

  !-------------------------------------------------------------------------------
  !> @brief   This routine remaps conservatively the mass from a grid xl to a
  !>          superposed grid xld.
  !> @details The mass(i) is the mass of cell [xl(i),xl(i+1)] and
  !>           mass_d(i) is the mass of the cell [xld(i),xld(i+1)]
  !>           if xl and xld are bounded by the same boundaries then
  !>           sum{mass} = sum{mass_d}
  !> @param[in] xld  Departure points of grid xl
  !> @param[in] xl   The grid-edges of the given mass (mass(i) is the mass of
  !!                 cell [xl(i),xl(i+1)])
  !> @param[in] dx   Size of the cells, dx(i)=xl(i+1)-xl(i)
  !> @param[in] mass The mass of grid cells
  !> @param[in] a0   The first coefficent of the cell reconstruction.
  !>                 Each cell has a cubic function f(x) = a0 + a1*x + a2*x^2 + a3*x^3
  !> @param[in] a1   The second coefficent of the cell reconstruction.
  !> @param[in] a2   The third  coefficent of the cell reconstruction.
  !> @param[in] a3   The fourth coefficent of the cell reconstruction.
  !> @param[in] ns   Index of the starting cell
  !> @param[in] nf   Index of the end cell
  !> @param[out] mass_d     The updated mass / or mass of the Lagrangian cell
  !> @param[out] mass_flux  The mass flux (or mass swept from xld(i) to xl(i))
  !-------------------------------------------------------------------------------
  subroutine remapp_mass( xld, xl, dx, mass, a0, a1, a2, a3, ns, nf, mass_d, mass_flux )
    implicit none
    integer(kind=i_def), intent(in)                     :: ns,nf
    real(kind=r_tran), dimension(ns:nf+1), intent(in)   :: xl, xld
    real(kind=r_tran), dimension(ns:nf),   intent(in)   :: dx, mass, a0,a1,a2,a3
    real(kind=r_tran), dimension(ns:nf),   intent(out)  :: mass_d
    real(kind=r_tran), dimension(ns:nf+1), intent(out)  :: mass_flux

    ! Local variables
    real(kind=r_tran)   :: m1,m2,hd,s,s1,s2
    integer(kind=i_def) :: j,i,jd,sig1,sig2
    real(kind=r_tran)   :: c1,c2,c3,sig,sigm1

    c1 = 0.5_r_tran
    c2 = 1.0_r_tran/3.0_r_tran
    c3 = 0.25_r_tran

    do j = ns + 1, nf
       m2 = 0.0_r_tran
       sig  = sign(1.0_r_tran, xl(j) - xld(j) )
       s = (1.0_r_tran - sig)/2.0_r_tran
       sig1 = int(sig, i_def)
       sig2 = int(s , i_def)

       jd = local_point_1d_array(xld(j), xl, ns, nf)
       hd = 1.0_r_tran/dx(jd)

       s = ( xld(j) - xl(jd) ) * hd
       s1 = max(sig*s,0.0_r_tran)
       s2 = max(sig,s)

       m1 = a0(jd)*(s2-s1) + c1*a1(jd)*(s2**2-s1**2)  &
                           + c2*a2(jd)*(s2**3-s1**3)  &
                           + c3*a3(jd)*(s2**4-s1**4)
       m1 = m1 * dx(jd)

       ! Make sure that m1 is of the right-signe
       !      and discard the wrong mass.
       ! This to maintain the right mass bounds
       !
       sigm1 = sign(1.0_r_tran,mass(jd))
       m1 = sigm1*max(sigm1*m1,0.0_r_tran)

       do i = jd+sig1, j-sig1-sig2, sig1
          m2 = m2 + mass(i)
       end do

       mass_flux(j) = (m1 + m2)*sig

    end do
    mass_flux(ns  ) = 0.0_r_tran
    mass_flux(nf+1) = 0.0_r_tran

  ! Compute the new/updated arrival mass (or mass of the lagrangian cell)
  do j = ns, nf
      mass_d(j)= mass(j) + mass_flux(j) - mass_flux(j+1)
   end do

  end subroutine remapp_mass

  !-------------------------------------------------------------------------------
  !> @brief Compute location index of the point p in a 1d array x.
  !> @param[in] p Coordinate of the point
  !> @param[in] x The 1d array of grid-points
  !> @param[in] ns Start of index of grid-points
  !> @param[in] nf End of index of grid-points
  !> @return location index "l" such as x(l) < p < x(l+1)
  !-------------------------------------------------------------------------------
  function local_point_1d_array(p, x, ns, nf)

  implicit none

  integer(kind=i_def), intent(in) :: ns, nf
  real(kind=r_tran), dimension(ns:nf), intent(in) :: x
  real(kind=r_tran),   intent(in) :: p
  integer(kind=i_def)             :: jl,jm,ju
  integer(kind=i_def)             :: local_point_1d_array

  jl=ns-1_i_def
  ju=nf+1_i_def
  do while ((ju-jl) > 1_i_def )
    jm=(ju+jl)/2_i_def
    if( (x(nf) > x(ns)) .eqv. (p > x(jm)) ) then
      jl=jm
    else
      ju=jm
    end if
  end do
  local_point_1d_array = max(ns,min(jl,nf))

  end function local_point_1d_array

  !-------------------------------------------------------------------------------
  !> @brief   This routine to modify the quadratic-coefficient to achieve
  !>          monotone reconstruction f(x) = a0 + a1*x + a2*x^2.
  !> @param[in]    vl   The left-value of the cell reconstruction
  !> @param[in]    vr   The right-value of the cell reconstruction
  !> @param[in]    vb   The average value mass/dx  of the cell
  !> @param[inout] a0   The first coefficent of the cell reconstruction
  !> @param[inout] a1   The second coefficent of the cell reconstruction
  !> @param[inout] a2   The third  coefficent of the cell reconstruction
  !> @param[in]    ns   Index of the starting cell
  !> @param[in]    nf   Index of the end cell
  !> @param[in] vertical_monotone_order   The order/option of the modified coefficients
  !-------------------------------------------------------------------------------
   subroutine  monotone_quadratic_coeffs(vl,vr,vb,a0,a1,a2,ns,nf,vertical_monotone_order)
   implicit none
   integer(kind=i_def),                 intent(in)    :: ns, nf, vertical_monotone_order
   real(kind=r_tran), dimension(ns:nf), intent(in)    :: vl, vr, vb
   real(kind=r_tran), dimension(ns:nf), intent(inout) :: a0, a1, a2

   integer(kind=i_def) :: i
   real(kind=r_tran)   :: xm, vh, test1, test2
   integer(kind=i_def), dimension(ns:nf) :: mod_id

  ! Identify the cells that needs modified and
  ! if the modified reconstruction preserves:
  !    the left value     "mod_id=1"
  ! or the right-value    "mod_id=2"

   mod_id(:) = 0_i_def
   do i = ns, nf
      xm = -a1(i)/(2.0_r_tran*a2(i)+EPS_R_TRAN)
      vh = 0.5_r_tran*(vr(i) + vl(i))
      test1 = xm*(1.0_r_tran - xm)
      test2 = (vb(i) - vl(i))*(vh - vb(i))
      !
      ! test1 = test if the turning point is inside the cell
      ! test2 = test if integral (vb) is close to left cell value (vl)
      !
      if ( test1 > 0.0_r_tran .and. test2 >= 0.0_r_tran ) then     ! flat left edge
          mod_id(i) = 1_i_def
      elseif ( test1 > 0.0_r_tran .and. test2 < 0.0_r_tran ) then  ! flat right edge
          mod_id(i) = 2_i_def
      end if
   end do

   select case( vertical_monotone_order )
   case ( vertical_monotone_order_constant )
        do i = ns, nf
          if ( mod_id(i) /= 0_i_def ) then
            a0(i) = vb(i)
            a1(i) = 0.0_r_tran
            a2(i) = 0.0_r_tran
          end if
        end do

   case ( vertical_monotone_order_linear )
      do i = ns, nf
        if ( mod_id(i) == 1_i_def ) then
           ! Linear left satisfying vl and vb
             a0(i) = vl(i)
             a1(i) = 2.0_r_tran*(vb(i)-vl(i))
             a2(i) = 0.0_r_tran
         elseif ( mod_id(i) == 2_i_def ) then
            ! Linear right satisfying vr and vb
             a0(i) = 2.0_r_tran*vb(i) - vr(i)
             a1(i) = 2.0_r_tran*(vr(i) - vb(i))
             a2(i) = 0.0_r_tran
        end if
      end do

   case ( vertical_monotone_order_high ) !high-order modification
     do i = ns, nf
        if ( mod_id(i) == 1_i_def ) then     ! flat left edge
           !-------------------------------------------------
           ! In this case the modified quadratic is flattened
           ! at x=0, with 3 conditions:
           ! (1) df(0)=0; (2) f(0)=vl,(3) Int(f,x=0..1) = m
           !-------------------------------------------------
            a1(i)=0.0_r_tran
            a2(i)=3.0_r_tran*(vb(i) - vl(i))
        elseif ( mod_id(i) == 2_i_def ) then  ! flat right edge
           !-------------------------------------------------
           ! In this case the modified quadratic is flatened
           ! at x=1, with 3 conditions:
           ! (1) df(1)=0; (2) f(1)=vr,(3) Int(f,x=0..1) = vb
           !-------------------------------------------------
           a0(i)=-2.0_r_tran*vr(i) + 3.0_r_tran*vb(i)
           a1(i)= 6.0_r_tran*vr(i) - 6.0_r_tran*vb(i)
           a2(i)=-3.0_r_tran*vr(i) + 3.0_r_tran*vb(i)
        end if
     end do
   end select

   end subroutine  monotone_quadratic_coeffs

  !-------------------------------------------------------------------------------
  !> @brief   This routine to modify the cubic-coefficient to achieve
  !>          monotone reconstruction f(x) = a0 + a1*x + a2*x^2 + a3*x^3.
  !> @param[in]    vl   The left-value of the cell reconstruction
  !> @param[in]    vr   The right-value of the cell reconstruction
  !> @param[in]    vb   The average value mass/dx  of the cell
  !> @param[inout] a0   The first coefficent of the cell reconstruction
  !> @param[inout] a1   The second coefficent of the cell reconstruction
  !> @param[inout] a2   The third  coefficent of the cell reconstruction
  !> @param[inout] a3   The fourth coefficent of the cell reconstruction
  !> @param[in]    ns   Index of the starting cell
  !> @param[in]    nf   Index of the end cell
  !> @param[in] vertical_monotone_order   The order/option of the modified coefficients
  !-------------------------------------------------------------------------------
  subroutine  monotone_cubic_coeffs(vl,vr,vb,a0,a1,a2,a3,ns,nf,vertical_monotone_order)
  implicit none
  integer(kind=i_def),                 intent(in)    :: ns, nf, vertical_monotone_order
  real(kind=r_tran), dimension(ns:nf), intent(in)    :: vl, vr, vb
  real(kind=r_tran), dimension(ns:nf), intent(inout) :: a0, a1, a2, a3
  integer(kind=i_def) :: i
  real(kind=r_tran)   :: a, b, c, d, x1, x2, v1, v2, vh
  real(kind=r_tran)   :: test1, test2, test3
  integer(kind=i_def), dimension(ns:nf) :: mod_id

  ! Identify the cells that needs modified and
  ! if the modified reconstruction preserves:
  !    the left value     "mod_id=1"
  ! or the right-value    "mod_id=2"

  mod_id(:) = 0_i_def
  do i = ns, nf
     a = 3.0_r_tran*a3(i)
     b = 2.0_r_tran*a2(i)
     c = a1(i)
     d = b**2 - 4.0_r_tran*a*c
     if ( d > 0.0_r_tran ) then  ! possible turning points within the cell
        x1 = (-b+sqrt(d))/(2.0_r_tran*a + EPS_R_TRAN)
        x2 = (-b-sqrt(d))/(2.0_r_tran*a + EPS_R_TRAN)
        v1 = a0(i) + a1(i)*x1 + a2(i)*(x1**2) + a3(i)*(x1**3)
        v2 = a0(i) + a1(i)*x2 + a2(i)*(x2**2) + a3(i)*(x2**3)
        test1 = x1*(1.0_r_tran - x1)
        test2 = x2*(1.0_r_tran - x2)
        if ( test1 > 0.0_r_tran .or. test2 > 0.0_r_tran ) then
            vh = 0.5_r_tran*(vr(i) + vl(i))
            test3 = (vb(i) - vl(i))*(vh - vb(i))
            if ( test3 >= 0.0_r_tran ) then  ! flat left
               mod_id(i) = 1_i_def
            else                            ! flat right
               mod_id(i) = 2_i_def
            end if
        end if
     end if
  end do


  select case( vertical_monotone_order )
  case ( vertical_monotone_order_constant )
      do i = ns, nf
        if ( mod_id(i) /= 0_i_def ) then
          a0(i) = vb(i)
          a1(i) = 0.0_r_tran
          a2(i) = 0.0_r_tran
          a3(i) = 0.0_r_tran
        end if
      end do
   case ( vertical_monotone_order_linear )
      do i = ns, nf
        if ( mod_id(i) == 1_i_def ) then
           ! Linear left satisfying vl and vb
             a0(i) = vl(i)
             a1(i) = 2.0_r_tran*(vb(i)-vl(i))
             a2(i) = 0.0_r_tran
             a3(i) = 0.0_r_tran
         elseif ( mod_id(i) == 2_i_def ) then
            ! Linear right satisfying vr and vb
             a0(i) = 2.0_r_tran*vb(i) - vr(i)
             a1(i) = 2.0_r_tran*(vr(i) - vb(i))
             a2(i) = 0.0_r_tran
             a3(i) = 0.0_r_tran
        end if
      end do
   case ( vertical_monotone_order_high ) !high-order modification
     do i = ns, nf
         if (  mod_id(i) == 1_i_def ) then  ! flatten left edge
            !---------------------------------------------
            ! In this case the modified cubic is flatened
            ! at x=0, with 4 conditions:
            ! (1) df(0)=0; (2) d2f(0)=0, (3) f(0)=vl, and
            ! (4) Int(f,x=0..1) = vb
            !----------------------------------------------
            a1(i) = 0.0_r_tran
            a2(i) = 0.0_r_tran
            a3(i) = -4.0_r_tran*vl(i) + 4.0_r_tran*vb(i)
        elseif ( mod_id(i) == 2_i_def ) then  ! flatten right edge
            !---------------------------------------------
            ! In this case the modified cubic is flatened
            ! at x=1, with 4 conditions:
            ! (1) df(1)=0; (2) d2f(1)=0, (3) f(1)=vr, and
            ! (4) Int(f,x=0..1) = vb
            !----------------------------------------------
            a0(i) =  -3.0_r_tran*vr(i) +  4.0_r_tran*vb(i)
            a1(i) =  12.0_r_tran*vr(i) - 12.0_r_tran*vb(i)
            a2(i) = - a1(i)
            a3(i) =   4.0_r_tran*vr(i) -  4.0_r_tran*vb(i)
        end if
     end do
   end select

  end subroutine  monotone_cubic_coeffs
  !-------------------------------------------------------------------------------
  !> @brief   This makes the edges values of cell monotone within the existing data
  !> @param[inout] fl           The edge function/reconstruction values of cells.
  !> @param[in]    fb           The intgeral/f_bar of the cell
  !> @param[in]    dx           cells size/lengh
  !> @param[in]    vertical_monotone   No-monotone identification option
  !> @param[in]    enforce_min_value   Enforce_min_value or not option
  !> @param[in]    ns           Index of the starting cell
  !> @param[in]    nf           Index of the end cell
  !-------------------------------------------------------------------------------
  subroutine monotone_edge_values(fl,fb,dx,vertical_monotone,enforce_min_value,ns,nf )
  implicit none
  integer(kind=i_def), intent(in) :: ns,nf,vertical_monotone
  logical(kind=l_def), intent(in) :: enforce_min_value
  real(kind=r_tran),   dimension(ns:nf+1), intent(inout) :: fl
  real(kind=r_tran),   dimension(ns:nf), intent(in)      :: fb, dx
  ! Local variables
  integer(kind=i_def) :: i, im1, im2, ip1
  real(kind=r_tran)   :: a, b, minv, maxv, test1, test2, test3
  real(kind=r_tran), dimension(ns:nf)     :: df
  real(kind=r_tran), dimension(ns-1:nf+1) :: fbe
  real(kind=r_tran) :: relative_min_value

  ! Add two extra points beyond boundaries with linear extrapolation
  ! so every edge value is within two average-value points
  fbe(ns:nf) = fb(ns:nf)
  ip1 = min(ns+1, nf)
  a = dx(ns )/(dx(ns)+dx(ip1))
  b = dx(ip1)/(dx(ns)+dx(ip1))
  fbe(ns-1) = a*(3.0_r_tran*fb(ns) - 2.0_r_tran*fb(ip1)) + b*fb(ns)
  im1 = max(nf-1, ns)
  a = dx(nf )/(dx(nf)+dx(im1))
  b = dx(im1)/(dx(nf)+dx(im1))
  fbe(nf+1) = a*(3.0_r_tran*fb(nf) - 2.0_r_tran*fb(im1)) + b*fb(nf)

  select case( vertical_monotone )

   case ( monotone_strict )
     do i = ns, nf + 1_i_def
        im1 = i - 1_i_def
        minv = min( fbe(i), fbe(im1)  )
        maxv = max( fbe(i), fbe(im1) )
        fl(i) = max( minv, min(fl(i), maxv) )
     end do

   case ( monotone_relaxed )
     do i = ns - 1, nf
        df(i) = fbe(i+1) - fbe(i)
     end do

     do i = ns, nf + 1_i_def
        im1 = i - 1_i_def
        im2 = max(ns - 1_i_def, i - 2_i_def)

        test1 = (fl(i)-fbe(im1))*(fbe(i)-fl(i))
        test2 = df(im2)*(fl(i)-fbe(im1))
        test3 = df(im2)*df(i)

        if ( ( test1 < 0.0_r_tran ) .and.  &
             ( test2 <= 0.0_r_tran .or. test3 <= 0.0_r_tran ) ) then
           minv = min( fbe(i), fbe(im1) )
           maxv = max( fbe(i), fbe(im1) )
           fl(i) = max( minv, min(fl(i), maxv) )
        end if
     enddo

  end select
  !
  ! Enforce min-value if required
  ! Here we are enforcing min_value= zero+epsilon (a relative zero)
  !
  if ( enforce_min_value )  then
       relative_min_value = maxval(abs(fl(:)))
       relative_min_value = max(EPS_R_TRAN, relative_min_value*EPS_R_TRAN)
       fl(:) = max( fl(:), relative_min_value )
  end if
  end subroutine monotone_edge_values

end module vertical_mass_remapping_kernel_mod
