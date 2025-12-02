!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2018.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief Routine to process and dump time series of maxima and minima of fields to file
module minmax_tseries_mod

  use constants_mod,                     only: r_def, str_max_filename, i_def
  use field_mod,                         only: field_type, field_proxy_type
  use diagnostic_alg_mod,                only: scalar_nodal_diagnostic_alg, &
                                               vector_nodal_diagnostic_alg
  use files_config_mod,                  only: diag_stem_name
  use lfric_mpi_mod,                     only: global_mpi
  use mesh_mod,                          only: mesh_type
  use fs_continuity_mod,                 only: W1, W2

  implicit none

  integer(i_def), private :: unitno

  private
  public :: minmax_tseries, minmax_tseries_init, minmax_tseries_final

contains


 subroutine minmax_tseries_init(field_name)

    use io_utility_mod, only: claim_io_unit

    implicit none

    character(len=*),    intent(in)    :: field_name
    character(len=str_max_filename)    :: fname

    if ( global_mpi%get_comm_rank() == 0 ) then

      unitno = claim_io_unit()
      fname = trim(diag_stem_name) // "_" // &
              trim("nodal_minmax_") // trim(field_name) // trim(".m")

      open(unit=unitno, file=fname, status='replace')

    end if

  end subroutine minmax_tseries_init


!> @brief Routine to process and dump time series of maxima and minima of fields to file
!> @details Writes field to a .m formatted file by dumping the values of maxima
!>          and minima on nodal points.
!> @param[in] field       Field to output
!> @param[in] field_name Name of field to output
!> @param[in] mesh       Mesh all fields are on
 subroutine minmax_tseries(field, field_name, mesh)
   use scalar_mod, only : scalar_type

   implicit none

   type(field_type), intent(in)    :: field
   character(len=*), intent(in)    :: field_name

   type(mesh_type),  intent(in), pointer :: mesh

   ! Internal variables
   type(field_type)                   :: nodal_output(3)
   type(field_type)                   :: nodal_coordinates(3)
   type(field_type)                   :: level
   integer(i_def)                     :: output_dim
   integer(i_def)                     :: fs

   type(field_proxy_type)             :: n_p(3)

   integer(i_def)                     :: i

   character(len=str_max_filename)    :: fname

   fname = trim(diag_stem_name) // "_" // &
           trim("nodal_minmax_") // trim(field_name) // trim(".m")

   fs = field%which_function_space()

   if ( fs == W1 .or. fs == W2) then
     ! Vector field
     call vector_nodal_diagnostic_alg(nodal_output, output_dim, &
                                      nodal_coordinates, level, &
                                      field_name, field)

   else
     ! Scalar field
     call scalar_nodal_diagnostic_alg( nodal_output, nodal_coordinates, level, &
                                       field_name, field, mesh, .false. )
   end if

   do i = 1,3
     n_p(i) = nodal_output(i)%get_proxy()
   end do

   if ( global_mpi%get_comm_rank() == 0 ) then
     do i=1,3
       write(unitno,'(2e16.8)') n_p(i)%get_max(), n_p(i)%get_min()
     enddo
   end if

 end subroutine minmax_tseries

 subroutine minmax_tseries_final()

    use io_utility_mod, only: release_io_unit

    implicit none

    if ( global_mpi%get_comm_rank() == 0 ) then
      close(unitno)
      call release_io_unit( unitno )
    end if

  end subroutine minmax_tseries_final

end module minmax_tseries_mod
