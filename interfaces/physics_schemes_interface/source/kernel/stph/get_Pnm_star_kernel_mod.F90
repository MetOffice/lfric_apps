!-------------------------------------------------------------------------------
!(c) Crown copyright 2021 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Create matrix of Legendre Polynomial coefficients
!> @details This kernel populates the Pnm_star variable with the associated
!!          Legendre Polynomial (including the Spherical harmonic amplitude,
!!          hence the "star" tag) with wavenumber n and m numbers using
!!          recurrence relationship. See documentation at
!!          SH_decomp_cubedsphere.pdf in TicketDetails of LFRIC ticket #2682
module get_Pnm_star_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, &
                               GH_WRITE, GH_INTEGER, GH_REAL, &
                               ANY_DISCONTINUOUS_SPACE_1,     &
                               ANY_DISCONTINUOUS_SPACE_2,     &
                               GH_READ, CELL_COLUMN

  use constants_mod,     only: r_def, i_def, pi
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: get_Pnm_star_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                   &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! Pnm_star
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! latitude
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          & ! stph_n_max
        /)
        integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass ::  get_Pnm_star_code
  end type get_Pnm_star_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public get_Pnm_star_code
contains

  !> @brief Populate Pnm_star matrix with values of the Assoc. Legendre Polynomials
  !> @param[in]      nlayers     The number of layers
  !> @param[in,out]  Pnm_star    Legendre Polynomials x SH amplitude
  !> @param[in]      latitude    2D field of latitudes
  !> @param[in]      stph_n_max  Stochastic Physics maximum wavenumber
  !> @param[in]      ndf_sp      Number of degrees of freedom per cell for spectral space
  !> @param[in]      undf_sp     Number of unique degrees of freedom for spectral space
  !> @param[in]      map_sp      Dofmap for the cell at the base of the column for spectral space
  !> @param[in]      ndf_2d      Number of degrees of freedom per cell for 2D space
  !> @param[in]      undf_2d     Number of unique degrees of freedom for 2D space
  !> @param[in]      map_2d      Dofmap for the cell at the base of the column for 2D space

  subroutine  get_Pnm_star_code(nlayers,   &
                                Pnm_star,  &
                                latitude,  &
                                stph_n_max,&
                                ndf_sp,    &
                                undf_sp,   &
                                map_sp,    &
                                ndf_2d,    &
                                undf_2d,   &
                                map_2d     &
                                )

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_sp, ndf_2d
    integer(kind=i_def), intent(in) :: undf_sp, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_sp)  ::  map_sp
    integer(kind=i_def), intent(in), dimension(ndf_2d)  ::  map_2d

    !fields
    real(kind=r_def), intent(inout), dimension(undf_sp) :: Pnm_star
    real(kind=r_def), intent(in),    dimension(undf_2d) :: latitude

    !SPT scalar
    integer(kind=i_def), intent(in) :: stph_n_max


    ! Equivalent Spherical values
    real(kind=r_def) :: sin_esp_lat_rad, cos_esp_lat_rad

    ! Factors for the recurrence relationships
    real(kind=r_def) :: fact_1st, fact_2nd, coeff

    ! Integers for iteration
    integer(kind=i_def) :: n, m, n_row, np1_row

    !!!! Compute the Pnm_star coefficients:
    ! Legendre Polynomials from -pi/2  to pi/2 defined as P[sin phi]
    ! in the range of [-1,1]

    ! Compute the Legendre Polynomials and the
    ! amplitude of the spherical harmonics together using adapted
    ! recurrence relationships for the Legendre Polynomials.
    !
    ! As Ynm = sqrt ((2n + 1)/4Pi * (n - m)!/(n + m)! )  Pnm [sin phi] exp(i m lambda)
    ! We create
    ! Pnm_star = sqrt ((2n + 1)/4Pi * (n - m)!/(n + m)! )  Pnm [sin phi]
    !
    ! Pnm_star includes the amplitude to simplify the operations
    ! done in the recurrence relationships, where factorials are
    ! removed. This avoids the division of very large numbers.

    ! Indexing of the Pnm_star matrix:
    ! It is a triangular matrix, so the indexing can be
    !      | m=0 | m=1 | m=2 | m=3 | m=4 |
    !  --------------------------------------
    ! n=0 |  0  |     |     |     |     |
    ! n=1 |  1  |  2  |     |     |     |
    ! n=2 |  3  |  4  |  5  |     |     |
    ! n=3 |  6  |  7  |  8  |  9  |     |
    ! n=4 | 10  | 11  | 12  | 13  | 14  |
    !
    ! Therefore the (n,m) location is given by
    ! m + sum (k, k-1, ...,n )

    ! Set values of the cosine and sine of latitude
    sin_esp_lat_rad=SIN(latitude(map_2d(1)))
    cos_esp_lat_rad=COS(latitude(map_2d(1)))

    ! Set Pnm_star (0, 0)
    Pnm_star(map_sp(1)) = 1.0_r_def/sqrt(4.0_r_def*pi)

    ! --- Set Pnm_star(n+1, n+1)

    ! Use recurrence relationship  Pnm(n+1, n+1) = - (2*n + 1) * cos phi * Pnm(n, n)
    ! transformed into
    ! Pnm_star(n+1, n+1) = - sqrt( (2*n + 3) / (2*n + 2) * cos phi * Pnm_star(n, n)
    !
    ! Indices Pnm_star(n+1, n+1) -> n_row + 2*(n + 1); where n_row = Sum (n)
    !         Pnm_star(n, n) -> n_row + n

    n_row=0
    do n = 0,stph_n_max-1
      n_row = n_row + n
      np1_row= n_row + n + 1
      coeff=(-1.0_r_def)*sqrt((2.0_r_def*n + 3.0_r_def)/(2.0_r_def*n + 2.0_r_def))
      Pnm_star(map_sp(1) + np1_row + n+1) =  coeff * cos_esp_lat_rad * Pnm_star(map_sp(1) + (n_row + n) )
    end do

    ! --- Set Pnm_star(n+1,n)

    ! Use recurrence relationship  Pnm(n+1, n) =  ( 2*n+1 ) * sin phi * Pnm(n, m)
    ! transformed into
    ! Pnm_star(n+1, n) = sqrt( 2*n + 3 ) sin phi * Pnm_star(n, n)
    !
    ! Indices Pnm_star(n+1, n) -> n_row + (n+1) + n; where n_row = Sum (n)
    !         Pnm_star(n,n) -> n_row + n
    n_row=0
    do n = 0,stph_n_max-1
      n_row = n_row + n
      np1_row= n_row + n + 1
      coeff=sqrt(2.0_r_def*n + 3.0_r_def)
      Pnm_star(map_sp(1) + np1_row + n) = coeff * sin_esp_lat_rad * Pnm_star(map_sp(1) + (n_row + n) )
    end do

    ! --- Set Pnm_star(n,m)

    ! Use recurrence relationship
    ! (n - m) * Pnm(n,m) = (2*n-1) * sin phi * Pnm(n-1, m)  - &
    !                           (n + m - 1) * Pnm(n-2, m)
    !
    ! Going from 2 (0 and 1 done in previous loop)
    ! Tansformed into
    ! P_star(n,m)  = fact_1st * sin phi * P_star(n-1,m)  -
    !                fact_2nd * P_star(n-2,m)
    !
    ! Where
    !
    ! fact_1st = sqrt( (4*n*n - 1) / (n*n m *m) )
    ! fact_2nd = sqrt( (2*n + 1) / (2*n - 3)  * ( (n-1)**2 - m*m ) / ( n*n - m*m ) )
    !
    ! Indices Pnm_star(n, m) -> n_row + m; where n_row = Sum (n)
    !         Pnm_star(n-1, m) -> (n_row-n) + m
    !         Pnm_star(n-2, m) -> (n_row-n-(n-1)) + m

    n_row=1 ! starting at one as 0 n loop start at 2
    do n=2,stph_n_max
      n_row = n_row + n
      do m=0,n-2
        fact_1st= sqrt((4.0_r_def*n*n - 1)/(n*n - m*m))
        fact_2nd= (-1.0_r_def)*sqrt((2.0_r_def*n + 1)/(2.0_r_def*n - 3.0_r_def)* &
                  ((n-1.0_r_def)**2.0_r_def-m*m)/(n*n - m*m))
        ! Compute SH weigth factor
        Pnm_star(map_sp(1) + n_row + m) = fact_1st * sin_esp_lat_rad* Pnm_star(map_sp(1) + n_row - n + m) + &
                                          fact_2nd * Pnm_star(map_sp(1) + n_row - n - (n-1) + m)
      end do
    end do

    ! Normalization factor for Spherical Harmonics for n != 0
    ! Indices Pnm_star(n,m) -> n_row + m; where n_row = Sum (n)
    n_row=0
    do n=1, stph_n_max
      n_row = n_row + n
      do m=1, n
        Pnm_star(map_sp(1)+ n_row + m) = Pnm_star(map_sp(1)+ n_row +m)*sqrt(2.0_r_def)
      end do
    end do

  end subroutine  get_Pnm_star_code

end module get_Pnm_star_kernel_mod
