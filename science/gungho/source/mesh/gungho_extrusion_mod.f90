!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!>
!> @brief Sets up vertical extrusion of mesh
!> @details This code contains two functions, create_extrusion() and
!>          create_shifted_extrusion(). The function create_extrusion() generates
!>          the standard vertical extrusion with different options, including
!>          uniform, quadratic, geometric and dcmip spacing.
!>          create_shifted_extrusion() creates a vertical extrusion with the same
!>          options as create_extrusion() but the top and bottom layers are
!>          half the normal cell height.
!>          There are technical infrastructure limitations which mean two different
!>          functions have been used to create the normal vertical mesh and the shifted
!>          vertical mesh. Tickets #1645 and #1659 deal with this issue of multiple
!>          instances of the mesh with different vertical extrusion.

module gungho_extrusion_mod

  use base_mesh_config_mod, only : key_from_geometry, &
                                   geometry_planar,   &
                                   geometry_spherical
  use constants_mod,        only : r_def, i_def
  use extrusion_mod,        only : extrusion_type,             &
                                   PRIME_EXTRUSION,            &
                                   uniform_extrusion_type,     &
                                   quadratic_extrusion_type,   &
                                   geometric_extrusion_type,   &
                                   shifted_extrusion_type,     &
                                   double_level_extrusion_type
  use extrusion_config_mod, only : key_from_method,            &
                                   method_specified_values,    &
                                   method_uniform,             &
                                   method_quadratic,           &
                                   method_geometric,           &
                                   method_dcmip,               &
                                   method_um_L38_29t_9s_40km,  &
                                   method_um_L85_50t_35s_85km, &
                                   method_um_L70_61t_9s_40km,  &
                                   method_um_L120_99t_21s_40km,&
                                   method_um_L140_122t_18s_40km,&
                                   method_um_L70_50t_20s_80km, &
                                   eta_values
  use log_mod,              only : log_event,       &
                                   log_level_error, &
                                   log_scratch_space

  implicit none

  private
  public create_extrusion
  public create_shifted_extrusion
  public create_double_level_extrusion

  character(*), parameter :: module_name = 'gungho_extrusion_mod'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with user-specified values
  !>
  type, public, extends(extrusion_type) :: specified_values_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => specified_values_extrude
  end type specified_values_extrusion_type

  interface specified_values_extrusion_type
    module procedure specified_values_extrusion_constructor
  end interface specified_values_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with specific UM configuration L38_29t_9s_40km
  !>
  type, public, extends(extrusion_type) :: um_L38_29t_9s_40km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L38_29t_9s_40km_extrude
  end type um_L38_29t_9s_40km_extrusion_type

  interface um_L38_29t_9s_40km_extrusion_type
    module procedure um_L38_29t_9s_40km_extrusion_constructor
  end interface um_L38_29t_9s_40km_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with specific UM configuration L70_50t_20s_80km
  !>
  type, public, extends(extrusion_type) :: um_L70_50t_20s_80km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L70_50t_20s_80km_extrude
  end type um_L70_50t_20s_80km_extrusion_type

  interface um_L70_50t_20s_80km_extrusion_type
    module procedure um_L70_50t_20s_80km_extrusion_constructor
  end interface um_L70_50t_20s_80km_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with specific UM configuration L85_50t_35s_85km
  !>
  type, public, extends(extrusion_type) :: um_L85_50t_35s_85km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L85_50t_35s_85km_extrude
  end type um_L85_50t_35s_85km_extrusion_type

  interface um_L85_50t_35s_85km_extrusion_type
    module procedure um_L85_50t_35s_85km_extrusion_constructor
  end interface um_L85_50t_35s_85km_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with specific UM configuration L70_61t_9s_40km
  !>
  type, public, extends(extrusion_type) :: um_L70_61t_9s_40km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L70_61t_9s_40km_extrude
  end type um_L70_61t_9s_40km_extrusion_type

  interface um_L70_61t_9s_40km_extrusion_type
    module procedure um_L70_61t_9s_40km_extrusion_constructor
  end interface um_L70_61t_9s_40km_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with specific UM configuration L120_99t_21s_40km
  !>
  type, public, extends(extrusion_type) :: um_L120_99t_21s_40km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L120_99t_21s_40km_extrude
  end type um_L120_99t_21s_40km_extrusion_type

  interface um_L120_99t_21s_40km_extrusion_type
    module procedure um_L120_99t_21s_40km_extrusion_constructor
  end interface um_L120_99t_21s_40km_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes with specific UM configuration L140_122t_18s_40km
  !>
  type, public, extends(extrusion_type) :: um_L140_122t_18s_40km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L140_122t_18s_40km_extrude
  end type um_L140_122t_18s_40km_extrusion_type

  interface um_L140_122t_18s_40km_extrusion_type
    module procedure um_L140_122t_18s_40km_extrusion_constructor
  end interface um_L140_122t_18s_40km_extrusion_type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes using DCMIP scheme.
  !>
  type, public, extends(extrusion_type) :: dcmip_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => dcmip_extrude
  end type dcmip_extrusion_type

  interface dcmip_extrusion_type
    module procedure dcmip_extrusion_constructor
  end interface dcmip_extrusion_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a specified_values_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function specified_values_extrusion_constructor( atmosphere_bottom, &
                                                   atmosphere_top,    &
                                                   number_of_layers,  &
                                                   extrusion_id ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(specified_values_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function specified_values_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh with specified values
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine specified_values_extrude( this, eta )

    implicit none

    class(specified_values_extrusion_type), intent(in)  :: this
    real(r_def),                            intent(out) :: eta(0:)

    eta(0) = 0.0000000_r_def
    eta(1:this%get_number_of_layers()-1) = eta_values ! from the namelist
    eta(this%get_number_of_layers()) = 1.0000000_r_def

  end subroutine specified_values_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a um_L38_29t_9s_40km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function um_L38_29t_9s_40km_extrusion_constructor( atmosphere_bottom, &
                                                     atmosphere_top,    &
                                                     number_of_layers,  &
                                                     extrusion_id ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L38_29t_9s_40km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L38_29t_9s_40km_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh with specific UM configuration L38_29t_9s_40km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine um_L38_29t_9s_40km_extrude( this, eta )

    implicit none

    class(um_L38_29t_9s_40km_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:)

    if (this%get_number_of_layers() /= 38)then
      call log_event( "Extrusion L38_29t_9s_40km reqires 38 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/                                    &
        0.0000000_r_def,  0.0005095_r_def,  0.0020380_r_def,  0.0045854_r_def, &
        0.0081519_r_def,  0.0127373_r_def,  0.0183417_r_def,  0.0249651_r_def, &
        0.0326074_r_def,  0.0412688_r_def,  0.0509491_r_def,  0.0616485_r_def, &
        0.0733668_r_def,  0.0861040_r_def,  0.0998603_r_def,  0.1146356_r_def, &
        0.1304298_r_def,  0.1472430_r_def,  0.1650752_r_def,  0.1839264_r_def, &
        0.2037966_r_def,  0.2246857_r_def,  0.2465938_r_def,  0.2695209_r_def, &
        0.2934670_r_def,  0.3184321_r_def,  0.3444162_r_def,  0.3714396_r_def, &
        0.3998142_r_def,  0.4298913_r_def,  0.4620737_r_def,  0.4968308_r_def, &
        0.5347160_r_def,  0.5763897_r_def,  0.6230643_r_def,  0.6772068_r_def, &
        0.7443435_r_def,  0.8383348_r_def,  1.0000000_r_def /)

  end subroutine um_L38_29t_9s_40km_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a um_L85_50t_35s_85km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function um_L85_50t_35s_85km_extrusion_constructor( atmosphere_bottom, &
                                                      atmosphere_top,    &
                                                      number_of_layers,  &
                                                      extrusion_id ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L85_50t_35s_85km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L85_50t_35s_85km_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh with specific UM configuration L85_50t_35s_85km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine um_L85_50t_35s_85km_extrude( this, eta )

    implicit none

    class(um_L85_50t_35s_85km_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:)

    if (this%get_number_of_layers() /= 85)then
      call log_event( "Extrusion L85_50t_35s_85km reqires 85 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/ &
        0.0000000E+00_r_def,   0.2352941E-03_r_def,   0.6274510E-03_r_def, &
        0.1176471E-02_r_def,   0.1882353E-02_r_def,   0.2745098E-02_r_def, &
        0.3764706E-02_r_def,   0.4941176E-02_r_def,   0.6274510E-02_r_def, &
        0.7764705E-02_r_def,   0.9411764E-02_r_def,   0.1121569E-01_r_def, &
        0.1317647E-01_r_def,   0.1529412E-01_r_def,   0.1756863E-01_r_def, &
        0.2000000E-01_r_def,   0.2258823E-01_r_def,   0.2533333E-01_r_def, &
        0.2823529E-01_r_def,   0.3129411E-01_r_def,   0.3450980E-01_r_def, &
        0.3788235E-01_r_def,   0.4141176E-01_r_def,   0.4509804E-01_r_def, &
        0.4894118E-01_r_def,   0.5294117E-01_r_def,   0.5709804E-01_r_def, &
        0.6141176E-01_r_def,   0.6588235E-01_r_def,   0.7050980E-01_r_def, &
        0.7529411E-01_r_def,   0.8023529E-01_r_def,   0.8533333E-01_r_def, &
        0.9058823E-01_r_def,   0.9600001E-01_r_def,   0.1015687E+00_r_def, &
        0.1072942E+00_r_def,   0.1131767E+00_r_def,   0.1192161E+00_r_def, &
        0.1254127E+00_r_def,   0.1317666E+00_r_def,   0.1382781E+00_r_def, &
        0.1449476E+00_r_def,   0.1517757E+00_r_def,   0.1587633E+00_r_def, &
        0.1659115E+00_r_def,   0.1732221E+00_r_def,   0.1806969E+00_r_def, &
        0.1883390E+00_r_def,   0.1961518E+00_r_def,   0.2041400E+00_r_def, &
        0.2123093E+00_r_def,   0.2206671E+00_r_def,   0.2292222E+00_r_def, &
        0.2379856E+00_r_def,   0.2469709E+00_r_def,   0.2561942E+00_r_def, &
        0.2656752E+00_r_def,   0.2754372E+00_r_def,   0.2855080E+00_r_def, &
        0.2959203E+00_r_def,   0.3067128E+00_r_def,   0.3179307E+00_r_def, &
        0.3296266E+00_r_def,   0.3418615E+00_r_def,   0.3547061E+00_r_def, &
        0.3682416E+00_r_def,   0.3825613E+00_r_def,   0.3977717E+00_r_def, &
        0.4139944E+00_r_def,   0.4313675E+00_r_def,   0.4500474E+00_r_def, &
        0.4702109E+00_r_def,   0.4920571E+00_r_def,   0.5158098E+00_r_def, &
        0.5417201E+00_r_def,   0.5700686E+00_r_def,   0.6011688E+00_r_def, &
        0.6353697E+00_r_def,   0.6730590E+00_r_def,   0.7146671E+00_r_def, &
        0.7606701E+00_r_def,   0.8115944E+00_r_def,   0.8680208E+00_r_def, &
        0.9305884E+00_r_def,   0.1000000E+01_r_def /)

  end subroutine um_L85_50t_35s_85km_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a um_L70_50t_20s_80km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function um_L70_50t_20s_80km_extrusion_constructor( atmosphere_bottom, &
                                                      atmosphere_top,    &
                                                      number_of_layers,  &
                                                      extrusion_id ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L70_50t_20s_80km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L70_50t_20s_80km_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh with specific UM configuration L70_50t_20s_80km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine um_L70_50t_20s_80km_extrude( this, eta )

    implicit none

    class(um_L70_50t_20s_80km_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:)


    if (this%get_number_of_layers() /= 70)then
      call log_event( "Extrusion L70_50t_20s_80km reqires 70 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/  &
       0.0000000_r_def,  0.0002500_r_def,  0.0006667_r_def,  0.0012500_r_def, &
       0.0020000_r_def,  0.0029167_r_def,  0.0040000_r_def,  0.0052500_r_def, &
       0.0066667_r_def,  0.0082500_r_def,  0.0100000_r_def,  0.0119167_r_def, &
       0.0140000_r_def,  0.0162500_r_def,  0.0186667_r_def,  0.0212500_r_def, &
       0.0240000_r_def,  0.0269167_r_def,  0.0300000_r_def,  0.0332500_r_def, &
       0.0366667_r_def,  0.0402500_r_def,  0.0440000_r_def,  0.0479167_r_def, &
       0.0520000_r_def,  0.0562500_r_def,  0.0606667_r_def,  0.0652500_r_def, &
       0.0700000_r_def,  0.0749167_r_def,  0.0800000_r_def,  0.0852500_r_def, &
       0.0906668_r_def,  0.0962505_r_def,  0.1020017_r_def,  0.1079213_r_def, &
       0.1140113_r_def,  0.1202745_r_def,  0.1267154_r_def,  0.1333406_r_def, &
       0.1401592_r_def,  0.1471838_r_def,  0.1544313_r_def,  0.1619238_r_def, &
       0.1696895_r_def,  0.1777643_r_def,  0.1861929_r_def,  0.1950307_r_def, &
       0.2043451_r_def,  0.2142178_r_def,  0.2247466_r_def,  0.2360480_r_def, &
       0.2482597_r_def,  0.2615432_r_def,  0.2760868_r_def,  0.2921094_r_def, &
       0.3098631_r_def,  0.3296378_r_def,  0.3517651_r_def,  0.3766222_r_def, &
       0.4046373_r_def,  0.4362943_r_def,  0.4721379_r_def,  0.5127798_r_def, &
       0.5589045_r_def,  0.6112759_r_def,  0.6707432_r_def,  0.7382500_r_def, &
       0.8148403_r_def,  0.9016668_r_def,  1.0000000_r_def /)

  end subroutine um_L70_50t_20s_80km_extrude

  !> @brief Extrudes the mesh with specific UM configuration L70_61t_9s_40km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine um_L70_61t_9s_40km_extrude( this, eta )

    implicit none

    class(um_L70_61t_9s_40km_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:)

    if (this%get_number_of_layers() /= 70)then
      call log_event( "Extrusion L70_61t_9s_40km reqires 70 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/ &
       0.0000000E+00_r_def,   0.1250000E-03_r_def,   0.5416666E-03_r_def,   0.1125000E-02_r_def,   0.1875000E-02_r_def, &
       0.2791667E-02_r_def,   0.3875000E-02_r_def,   0.5125000E-02_r_def,   0.6541667E-02_r_def,   0.8125000E-02_r_def, &
       0.9875000E-02_r_def,   0.1179167E-01_r_def,   0.1387500E-01_r_def,   0.1612500E-01_r_def,   0.1854167E-01_r_def, &
       0.2112500E-01_r_def,   0.2387500E-01_r_def,   0.2679167E-01_r_def,   0.2987500E-01_r_def,   0.3312500E-01_r_def, &
       0.3654167E-01_r_def,   0.4012500E-01_r_def,   0.4387500E-01_r_def,   0.4779167E-01_r_def,   0.5187500E-01_r_def, &
       0.5612501E-01_r_def,   0.6054167E-01_r_def,   0.6512500E-01_r_def,   0.6987500E-01_r_def,   0.7479167E-01_r_def, &
       0.7987500E-01_r_def,   0.8512500E-01_r_def,   0.9054167E-01_r_def,   0.9612500E-01_r_def,   0.1018750E+00_r_def, &
       0.1077917E+00_r_def,   0.1138750E+00_r_def,   0.1201250E+00_r_def,   0.1265417E+00_r_def,   0.1331250E+00_r_def, &
       0.1398750E+00_r_def,   0.1467917E+00_r_def,   0.1538752E+00_r_def,   0.1611287E+00_r_def,   0.1685623E+00_r_def, &
       0.1761954E+00_r_def,   0.1840590E+00_r_def,   0.1921980E+00_r_def,   0.2006732E+00_r_def,   0.2095645E+00_r_def, &
       0.2189729E+00_r_def,   0.2290236E+00_r_def,   0.2398690E+00_r_def,   0.2516917E+00_r_def,   0.2647077E+00_r_def, &
       0.2791699E+00_r_def,   0.2953717E+00_r_def,   0.3136506E+00_r_def,   0.3343919E+00_r_def,   0.3580330E+00_r_def, &
       0.3850676E+00_r_def,   0.4160496E+00_r_def,   0.4515977E+00_r_def,   0.4924007E+00_r_def,   0.5392213E+00_r_def, &
       0.5929016E+00_r_def,   0.6543679E+00_r_def,   0.7246365E+00_r_def,   0.8048183E+00_r_def,   0.8961251E+00_r_def, &
       0.1000000E+01_r_def &
       /)

  end subroutine um_L70_61t_9s_40km_extrude

  !> @brief Extrudes the mesh with specific UM configuration L120_99t_21s_40km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine um_L120_99t_21s_40km_extrude( this, eta )

    implicit none

    class(um_L120_99t_21s_40km_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:)

    if (this%get_number_of_layers() /= 120)then
      call log_event( "Extrusion L120_99t_21s_40km reqires 120 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/ &
       0.0000000E+00_r_def,   0.1250000E-03_r_def,   0.3000000E-03_r_def,   0.5250000E-03_r_def,   0.8000000E-03_r_def, &
       0.1125000E-02_r_def,   0.1500000E-02_r_def,   0.1925000E-02_r_def,   0.2400000E-02_r_def,   0.2925000E-02_r_def, &
       0.3500000E-02_r_def,   0.4125000E-02_r_def,   0.4800000E-02_r_def,   0.5525000E-02_r_def,   0.6300000E-02_r_def, &
       0.7125000E-02_r_def,   0.8000000E-02_r_def,   0.8925000E-02_r_def,   0.9900000E-02_r_def,   0.1092500E-01_r_def, &
       0.1200000E-01_r_def,   0.1312500E-01_r_def,   0.1430000E-01_r_def,   0.1552500E-01_r_def,   0.1680000E-01_r_def, &
       0.1812500E-01_r_def,   0.1950000E-01_r_def,   0.2092500E-01_r_def,   0.2240000E-01_r_def,   0.2392500E-01_r_def, &
       0.2550000E-01_r_def,   0.2712500E-01_r_def,   0.2880000E-01_r_def,   0.3052500E-01_r_def,   0.3230000E-01_r_def, &
       0.3412500E-01_r_def,   0.3600000E-01_r_def,   0.3792500E-01_r_def,   0.3990000E-01_r_def,   0.4192500E-01_r_def, &
       0.4400000E-01_r_def,   0.4612500E-01_r_def,   0.4830000E-01_r_def,   0.5052500E-01_r_def,   0.5280000E-01_r_def, &
       0.5512500E-01_r_def,   0.5750000E-01_r_def,   0.5992500E-01_r_def,   0.6240000E-01_r_def,   0.6492500E-01_r_def, &
       0.6750000E-01_r_def,   0.7012500E-01_r_def,   0.7280000E-01_r_def,   0.7554308E-01_r_def,   0.7838037E-01_r_def, &
       0.8132010E-01_r_def,   0.8436962E-01_r_def,   0.8753612E-01_r_def,   0.9082688E-01_r_def,   0.9424930E-01_r_def, &
       0.9781096E-01_r_def,   0.1015196E+00_r_def,   0.1053834E+00_r_def,   0.1094103E+00_r_def,   0.1136090E+00_r_def, &
       0.1179882E+00_r_def,   0.1225567E+00_r_def,   0.1273237E+00_r_def,   0.1322988E+00_r_def,   0.1374915E+00_r_def, &
       0.1429119E+00_r_def,   0.1485701E+00_r_def,   0.1544765E+00_r_def,   0.1606420E+00_r_def,   0.1670773E+00_r_def, &
       0.1737938E+00_r_def,   0.1808028E+00_r_def,   0.1881162E+00_r_def,   0.1957460E+00_r_def,   0.2037043E+00_r_def, &
       0.2120037E+00_r_def,   0.2206570E+00_r_def,   0.2296771E+00_r_def,   0.2390774E+00_r_def,   0.2488715E+00_r_def, &
       0.2590731E+00_r_def,   0.2696964E+00_r_def,   0.2807557E+00_r_def,   0.2922657E+00_r_def,   0.3042411E+00_r_def, &
       0.3166971E+00_r_def,   0.3296493E+00_r_def,   0.3431131E+00_r_def,   0.3571047E+00_r_def,   0.3716401E+00_r_def, &
       0.3867360E+00_r_def,   0.4024089E+00_r_def,   0.4186760E+00_r_def,   0.4355546E+00_r_def,   0.4530622E+00_r_def, &
       0.4712166E+00_r_def,   0.4900359E+00_r_def,   0.5095385E+00_r_def,   0.5297430E+00_r_def,   0.5506685E+00_r_def, &
       0.5723339E+00_r_def,   0.5947588E+00_r_def,   0.6179630E+00_r_def,   0.6419664E+00_r_def,   0.6667892E+00_r_def, &
       0.6924521E+00_r_def,   0.7189759E+00_r_def,   0.7463816E+00_r_def,   0.7746907E+00_r_def,   0.8039246E+00_r_def, &
       0.8341056E+00_r_def,   0.8652556E+00_r_def,   0.8973972E+00_r_def,   0.9305530E+00_r_def,   0.9647462E+00_r_def, &
       0.1000000E+01_r_def &
       /)

  end subroutine um_L120_99t_21s_40km_extrude

  !> @brief Extrudes the mesh with specific UM configuration L140_122t_18s_40km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine um_L140_122t_18s_40km_extrude( this, eta )

    implicit none

    class(um_L140_122t_18s_40km_extrusion_type), intent(in)  :: this
    real(r_def),                   intent(out) :: eta(0:)

    if (this%get_number_of_layers() /= 140)then
      call log_event( "Extrusion L140_122t_18s_40km reqires 140 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/ &
       0.0000000E+00_r_def,   0.5000001E-04_r_def,   0.1333334E-03_r_def,   0.2500001E-03_r_def,   0.4000001E-03_r_def, &
       0.5833335E-03_r_def,   0.8000002E-03_r_def,   0.1050000E-02_r_def,   0.1333334E-02_r_def,   0.1650000E-02_r_def, &
       0.2000001E-02_r_def,   0.2383334E-02_r_def,   0.2800001E-02_r_def,   0.3250001E-02_r_def,   0.3733334E-02_r_def, &
       0.4250001E-02_r_def,   0.4800001E-02_r_def,   0.5383335E-02_r_def,   0.6000001E-02_r_def,   0.6650002E-02_r_def, &
       0.7333335E-02_r_def,   0.8050002E-02_r_def,   0.8800002E-02_r_def,   0.9583336E-02_r_def,   0.1040000E-01_r_def, &
       0.1125000E-01_r_def,   0.1213334E-01_r_def,   0.1305000E-01_r_def,   0.1400000E-01_r_def,   0.1498334E-01_r_def, &
       0.1600000E-01_r_def,   0.1705000E-01_r_def,   0.1813334E-01_r_def,   0.1925001E-01_r_def,   0.2040001E-01_r_def, &
       0.2158334E-01_r_def,   0.2280001E-01_r_def,   0.2405001E-01_r_def,   0.2533334E-01_r_def,   0.2665001E-01_r_def, &
       0.2800001E-01_r_def,   0.2938334E-01_r_def,   0.3080001E-01_r_def,   0.3225001E-01_r_def,   0.3373334E-01_r_def, &
       0.3525001E-01_r_def,   0.3680001E-01_r_def,   0.3838334E-01_r_def,   0.4000001E-01_r_def,   0.4165001E-01_r_def, &
       0.4333334E-01_r_def,   0.4505001E-01_r_def,   0.4680001E-01_r_def,   0.4858335E-01_r_def,   0.5040001E-01_r_def, &
       0.5225001E-01_r_def,   0.5413335E-01_r_def,   0.5605002E-01_r_def,   0.5800001E-01_r_def,   0.5998335E-01_r_def, &
       0.6200002E-01_r_def,   0.6405002E-01_r_def,   0.6613335E-01_r_def,   0.6825002E-01_r_def,   0.7040001E-01_r_def, &
       0.7258336E-01_r_def,   0.7480001E-01_r_def,   0.7705002E-01_r_def,   0.7933335E-01_r_def,   0.8165002E-01_r_def, &
       0.8400003E-01_r_def,   0.8638336E-01_r_def,   0.8880002E-01_r_def,   0.9125002E-01_r_def,   0.9373336E-01_r_def, &
       0.9625003E-01_r_def,   0.9880003E-01_r_def,   0.1013834E+00_r_def,   0.1040000E+00_r_def,   0.1066500E+00_r_def, &
       0.1093334E+00_r_def,   0.1120500E+00_r_def,   0.1148000E+00_r_def,   0.1175834E+00_r_def,   0.1204000E+00_r_def, &
       0.1232500E+00_r_def,   0.1261334E+00_r_def,   0.1290500E+00_r_def,   0.1320000E+00_r_def,   0.1349834E+00_r_def, &
       0.1380000E+00_r_def,   0.1410500E+00_r_def,   0.1441334E+00_r_def,   0.1472500E+00_r_def,   0.1504001E+00_r_def, &
       0.1535837E+00_r_def,   0.1568019E+00_r_def,   0.1600566E+00_r_def,   0.1633512E+00_r_def,   0.1666905E+00_r_def, &
       0.1700811E+00_r_def,   0.1735319E+00_r_def,   0.1770539E+00_r_def,   0.1806609E+00_r_def,   0.1843694E+00_r_def, &
       0.1881994E+00_r_def,   0.1921741E+00_r_def,   0.1963207E+00_r_def,   0.2006706E+00_r_def,   0.2052595E+00_r_def, &
       0.2101280E+00_r_def,   0.2153219E+00_r_def,   0.2208926E+00_r_def,   0.2268972E+00_r_def,   0.2333995E+00_r_def, &
       0.2404696E+00_r_def,   0.2481849E+00_r_def,   0.2566306E+00_r_def,   0.2658995E+00_r_def,   0.2760931E+00_r_def, &
       0.2873217E+00_r_def,   0.2997051E+00_r_def,   0.3133729E+00_r_def,   0.3284651E+00_r_def,   0.3451325E+00_r_def, &
       0.3635374E+00_r_def,   0.3838541E+00_r_def,   0.4062692E+00_r_def,   0.4309826E+00_r_def,   0.4582078E+00_r_def, &
       0.4881724E+00_r_def,   0.5211188E+00_r_def,   0.5573052E+00_r_def,   0.5970055E+00_r_def,   0.6405105E+00_r_def, &
       0.6881285E+00_r_def,   0.7401855E+00_r_def,   0.7970266E+00_r_def,   0.8590163E+00_r_def,   0.9265389E+00_r_def, &
       0.1000000E+01_r_def &
       /)

  end subroutine um_L140_122t_18s_40km_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a um_L70_61t_9s_40km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function um_L70_61t_9s_40km_extrusion_constructor( atmosphere_bottom, &
                                                     atmosphere_top,    &
                                                     number_of_layers,   &
                                                     extrusion_id ) result(new)
    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L70_61t_9s_40km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L70_61t_9s_40km_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a um_L120_99t_21s_40km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function um_L120_99t_21s_40km_extrusion_constructor( atmosphere_bottom, &
                                                       atmosphere_top,    &
                                                       number_of_layers,   &
                                                       extrusion_id ) result(new)
    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L120_99t_21s_40km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L120_99t_21s_40km_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a um_L140_122t_18s_40km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  !>
  function um_L140_122t_18s_40km_extrusion_constructor( atmosphere_bottom, &
                                                       atmosphere_top,    &
                                                       number_of_layers,   &
                                                       extrusion_id ) result(new)
    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L140_122t_18s_40km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L140_122t_18s_40km_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a dcmip_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New dcmip_extrusion_type object.
  !>
  function dcmip_extrusion_constructor( atmosphere_bottom, &
                                        atmosphere_top,    &
                                        number_of_layers,  &
                                        extrusion_id ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(dcmip_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function dcmip_extrusion_constructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Extrudes the mesh using the DCMIP scheme.
  !>
  !> For more information see DCMIP-TestCaseDocument_v1.7.pdf,
  !> Appendix F.2. - Eq. 229.
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  !>
  subroutine dcmip_extrude( this, eta )

    implicit none

    class(dcmip_extrusion_type), intent(in)  :: this
    real(r_def),                 intent(out) :: eta(0:)

    real(r_def), parameter :: phi_flatten = 15.0_r_def

    integer(i_def) :: k

    do k = 0, this%get_number_of_layers()
      eta(k) = dcmip_func(real(k,r_def)/real(this%get_number_of_layers(),r_def))
    end do

  end subroutine dcmip_extrude

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Helper function for generating DCMIP extrusion
  !>
  !> @param[in] eta_uni   Input value which increases incrementally with level number
  !> @return    eta       Vertical eta coordinate
  !>
  function dcmip_func(eta_uni) result(eta)
    implicit none

    real(r_def), intent(in) :: eta_uni
    real(r_def) :: eta

    real(r_def), parameter :: phi_flatten = 15.0_r_def

    eta = ( sqrt(phi_flatten*(eta_uni**2_i_def) + 1.0_r_def) &
                    - 1.0_r_def ) / &
                  ( sqrt(phi_flatten + 1.0_r_def) - 1.0_r_def )

  end function dcmip_func

  !> @brief Creates vertical mesh extrusion
  !> @details Creates vertical mesh with nlayers.
  !> @param  method            Extrusion method to apply
  !> @param  geometry          Enumeration integer for 2D mesh geometrical shape
  !> @param  number_of_layers  Layers in extrusion
  !> @param  domain_height     Height above surface level (m)
  !> @param  scaled_radius     Radius of surface level, spherical geometries only
  !> @return extrusion         Extrusion class
  function create_extrusion( method,           &
                             geometry,         &
                             number_of_layers, &
                             domain_height,    &
                             scaled_radius ) result(extrusion)

    implicit none

    integer,        intent(in) :: method
    integer,        intent(in) :: geometry
    integer(i_def), intent(in) :: number_of_layers
    real(r_def),    intent(in) :: domain_height
    real(r_def),    intent(in) :: scaled_radius

    class(extrusion_type), allocatable :: extrusion

    real(r_def) :: atmosphere_bottom

    if (allocated(extrusion)) deallocate(extrusion)

    select case (geometry)
      case (geometry_planar)
        atmosphere_bottom = 0.0_r_def
      case (geometry_spherical)
        atmosphere_bottom = scaled_radius
      case default
        write( log_scratch_space,                      &
               '(A, ": Unrecognised geometry: ", A)' ) &
             module_name, key_from_geometry( geometry )
        call log_event( log_scratch_space, log_level_error )
    end select

    select case (method)
      case (method_specified_values)
        allocate( extrusion, source=specified_values_extrusion_type(   &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )
      case (method_uniform)
        allocate( extrusion, source=uniform_extrusion_type(              &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )

      case (method_um_l38_29t_9s_40km)
        allocate( extrusion, source=um_l38_29t_9s_40km_extrusion_type(   &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )

      case (method_um_l85_50t_35s_85km)
        allocate( extrusion, source=um_l85_50t_35s_85km_extrusion_type(  &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )

      case (method_um_l70_50t_20s_80km)
        allocate( extrusion, source=um_l70_50t_20s_80km_extrusion_type(  &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )

      case (method_um_l70_61t_9s_40km)
        allocate( extrusion, source=um_l70_61t_9s_40km_extrusion_type(   &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )

      case (method_um_L120_99t_21s_40km)
        allocate( extrusion, source=um_L120_99t_21s_40km_extrusion_type( &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )
      case (method_um_L140_122t_18s_40km)
        allocate( extrusion, source=um_L140_122t_18s_40km_extrusion_type(&
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )
      case (method_quadratic)
        allocate( extrusion, source=quadratic_extrusion_type(            &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )
      case (method_geometric)
        allocate( extrusion, source=geometric_extrusion_type(            &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )
      case (method_dcmip)
        allocate( extrusion, source=dcmip_extrusion_type(                &
                                       atmosphere_bottom, domain_height, &
                                       number_of_layers, prime_extrusion ) )

      case default
        write( log_scratch_space, '(A)' )                              &
            trim(module_name) // ': Unrecognised extrusion method: ' //&
            trim(key_from_method( method ))
        call log_event( log_scratch_space, log_level_error )

    end select

  end function create_extrusion

  !> @brief Creates vertical mesh extrusion for vertically shifted mesh.
  !> @details Creates vertically shifted mesh with nlayers+1 with the top and
  !>          bottom levels having half the cell height of the normal extrusion.
  !> @param[in] old The original extrusion.
  !> @return new     Extrusion class
  function create_shifted_extrusion(old) result(new)

    implicit none

    class(extrusion_type),  intent(in) :: old
    class(extrusion_type), allocatable :: new

    if (allocated(new)) deallocate(new)

    allocate(new, source=shifted_extrusion_type(old))

  end function create_shifted_extrusion

  !> @brief Creates vertical mesh extrusion for double level mesh.
  !> @details Creates vertically double level mesh from normal and shifted meshes
  !> @return new     Extrusion class
  function create_double_level_extrusion(old) result(new)

    implicit none

    class(extrusion_type),  intent(in) :: old
    class(extrusion_type), allocatable :: new

    if (allocated(new)) deallocate(new)

    allocate(new, source=double_level_extrusion_type(old))

  end function create_double_level_extrusion

end module gungho_extrusion_mod
