!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface with tnuc code in UM to pass tnuc into LFRic
module lsp_prognostic_tnuc_kernel_mod

  use argument_mod,      only : arg_type,          &
                                GH_FIELD, GH_REAL, &
                                GH_READ, GH_WRITE, &
                                DOMAIN,            &
                                GH_INTEGER,        &
                                ANY_DISCONTINUOUS_SPACE_1

  use fs_continuity_mod, only:  Wtheta
  use constants_mod,     only : r_def, i_def, r_um
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: lsp_prognostic_tnuc_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                  &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                  & ! n_acc_ins (name in LFRic)
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                  & ! n_cor_ins
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                  & ! theta_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                  & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA)                   & ! tnuc (out)
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: lsp_prognostic_tnuc_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: lsp_prognostic_tnuc_code

contains

  !> @param[in]     nlayers             Number of layers
  !> @param[in]     n_acc_ins           Accumulation mode
  !> @param[in]     n_cor_ins           Coarse mode
  !> @param[in]     theta_in_wth        Potential temperature field
  !> @param[in]     exner_in_wth        Exner pressure in potential temperature space
  !> @param[out]    tnuc                Temperature of nucleation (K)
  !> @param[in]     ndf_wth             Number of degrees of freedom per cell for theta space
  !> @param[in]     undf_wth            Number of unique degrees of freedom in segment for theta space
  !> @param[in]     map_wth             Dofmap for a segment of the theta space

  subroutine lsp_prognostic_tnuc_code(nlayers,             &
                                      seg_len,             &
                                      n_acc_ins,           &
                                      n_cor_ins,           &
                                      theta_in_wth,        &
                                      exner_in_wth,        &
                                      tnuc,                &
                                      ndf_wth,             &
                                      undf_wth,            &
                                      map_wth)

  !---------------------------------------
  ! UM modules
  !---------------------------------------

  use lsp_prognostic_tnuc_mod, only: lsp_prognostic_tnuc
  use planet_constants_mod,    only: p_zero, kappa

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, seg_len
  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in), dimension(ndf_wth, seg_len) :: map_wth
  integer(kind=i_def) :: i, j, k

  real(kind=r_def), intent(in),  dimension(undf_wth)  ::  n_acc_ins
  real(kind=r_def), intent(in),  dimension(undf_wth)  ::  n_cor_ins
  real(kind=r_def), intent(in),  dimension(undf_wth)  ::  theta_in_wth
  real(kind=r_def), intent(in),  dimension(undf_wth)  ::  exner_in_wth

  real(kind=r_def), intent(inout),  dimension(undf_wth) ::  tnuc

  real(r_um), allocatable :: dust_tot_nd(:,:,:)
  real(r_um), allocatable :: tnuc_new(:,:,:)

  real(r_um), pointer :: arcldust_b1(:,:,:)
  real(r_um), pointer :: arcldust_b2(:,:,:)
  real(r_um), pointer :: arcldust_b3(:,:,:)
  real(r_um), pointer :: arcldust_b4(:,:,:)
  real(r_um), pointer :: arcldust_b5(:,:,:)
  real(r_um), pointer :: arcldust_b6(:,:,:)

  real(r_um), target :: arcldust(seg_len,1,0:nlayers)

  !------------------------------------------
  ! Local variables for the kernel
  !------------------------------------------
  real(r_um), dimension(seg_len,1,0:nlayers) ::      &
              tracer_ukca_iacc, tracer_ukca_icor,    &
              theta, exner_theta_levels,             &
              dust_div1, dust_div2,                  &
              dust_div3, dust_div4,                  &
              dust_div5, dust_div6,                  &
              p_theta_levels

  ! Set dust variables to 0 as not used in LFRic
  dust_div1       = 0.0_r_um
  dust_div2       = 0.0_r_um
  dust_div3       = 0.0_r_um
  dust_div4       = 0.0_r_um
  dust_div5       = 0.0_r_um
  dust_div6       = 0.0_r_um

  ! Set arcldust variables to target size arcldust(seg_len,1,0:nlayers)
  arcldust_b1 => arcldust
  arcldust_b2 => arcldust
  arcldust_b3 => arcldust
  arcldust_b4 => arcldust
  arcldust_b5 => arcldust
  arcldust_b6 => arcldust

  !---------------------------------------------
  ! Mapping of LFRic fields into UM variables
  !--------------------------------------------

  ! allocate size of dust_tot_nd and tnuc_new
  allocate(dust_tot_nd(seg_len,1,0:nlayers))
  allocate(tnuc_new(seg_len,1,0:nlayers))

  ! Map LFRic iacc onto UM iacc to pass to UM
  j = 1
  do k = 1, nlayers
    do i = 1, seg_len
      tracer_ukca_iacc(i,j,k) = real(n_acc_ins(map_wth(1,i) + k),kind=r_um)
    end do ! i
  end do ! k

  ! Map LFRic iacc onto UM iacc to pass to UM
  j = 1
  do k = 1, nlayers
    do i = 1, seg_len
      tracer_ukca_icor(i,j,k) = real(n_cor_ins(map_wth(1,i) + k),kind=r_um)
    end do ! i
  end do ! k


  ! Map LFRic theta onto UM theta to pass to UM
  j = 1
  do k = 1, nlayers
    do i = 1, seg_len
      theta(i,j,k) = real(theta_in_wth(map_wth(1,i) + k),kind=r_um)
    end do ! i
  end do ! k

  ! Map LFRic exner onto UM exner to pass to UM
  j = 1
  do k = 1, nlayers
    do i = 1, seg_len
      exner_theta_levels(i,j,k) = real(exner_in_wth(map_wth(1,i) + k),kind=r_um)
    end do ! i
  end do ! k

 ! Calculate pressure on theta levels
  j = 1
  do k = 1, nlayers
    do i = 1, seg_len
      p_theta_levels(i,j,k)    = real((p_zero*(exner_in_wth(                 &
                                      map_wth(1,i) + k))                     &
                                      **(1.0_r_um/kappa)),kind=r_um)
    end do ! i
  end do ! k

  !---------------------------------------------------
  ! Initialisation of prognostic variables and arrays
  !---------------------------------------------------

  call lsp_prognostic_tnuc(dust_div1                                        &
                          ,dust_div2                                        &
                          ,dust_div3                                        &
                          ,dust_div4                                        &
                          ,dust_div5                                        &
                          ,dust_div6                                        &
                          ,tracer_ukca_iacc                                 &
                          ,tracer_ukca_icor                                 &
                          ,theta                                            &
                          ,exner_theta_levels                               &
                          ,p_theta_levels                                   &
                          ,dust_tot_nd                                      &
                          ,tnuc_new                                         &
                          ,arcldust_b1                                      &
                          ,arcldust_b2                                      &
                          ,arcldust_b3                                      &
                          ,arcldust_b4                                      &
                          ,arcldust_b5                                      &
                          ,arcldust_b6 )

  ! Map tnuc_new (UM) back to tnuc (LFRic) to be passed out of kernel
  j = 1
  do k = 1, nlayers
    do i = 1, seg_len
      tnuc(map_wth(1,i) + k) = real(tnuc_new(i,j,k),r_def)
     end do ! i
  end do ! k

  ! deallocate variables
  deallocate(dust_tot_nd)
  deallocate(tnuc_new)

  end subroutine lsp_prognostic_tnuc_code

end module lsp_prognostic_tnuc_kernel_mod
