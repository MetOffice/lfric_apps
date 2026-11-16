.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _coordinates_and_vectors_vector_component_transforms:

Vector Component Transforms
===========================

This page describes transformations implemented in LFRic between computational
wind fields and their components, in different function spaces.


For the mathematical background, see
:ref:`coordinates_and_vectors_vector_components`. In brief, the values stored
in a :math:`\mathbb{W}_2` (or :math:`\mathbb{W}_1`) field are *computational*
values, related to physical components through the Jacobian :math:`J` of the
coordinate mapping. The transforms below fall into four groups.

Paths are relative to ``lfric_core/components/science/source/`` (``CORE``),
``lfric_core/infrastructure/source/`` (``INFRA``) and
``lfric_apps/science/gungho/source/`` (``GUNGHO``).

Conventions
-----------

* The coordinate field :math:`\boldsymbol\chi` is either in the ``native``
  system (cubed-sphere :math:`(\alpha,\beta,r)` on global meshes,
  :math:`(\lambda,\phi,r)` on limited-area meshes, optionally rotated and
  stretched) or in global Cartesian :math:`(X,Y,Z)` (``coord_system_xyz``).
  On a plane :math:`\boldsymbol\chi` is Cartesian.
* The Jacobian routines (``CORE/kernel/geometry/sci_coordinate_jacobian_mod.F90``,
  and the column-optimised ``sci_native_jacobian_mod.F90`` for native
  coordinates) return the Jacobian from the reference element to **global
  Cartesian** space, including any rotation and stretching. Therefore the
  Piola transforms below always return global Cartesian components on the
  sphere; converting to :math:`(\lambda,\phi,r)` components is a separate step.
* Lowest-order :math:`\mathbb{W}_2` DoFs are ordered W, S, E, N, bottom, top.

Coordinate-system transforms
----------------------------

Pure maths, in ``INFRA/utilities/coord_transform_mod.F90``:

.. list-table::
    :header-rows: 1
    :widths: 35 65

    * - Routine
      - Purpose
    * - ``sphere2cart_vector``
      - :math:`(\lambda,\phi,r)` components to Cartesian components
    * - ``cart2sphere_vector``, ``cart2sphere_scalar``
      - Cartesian to :math:`(\lambda,\phi,r)` components, in m/s. The scalar
        version takes and returns the three components as separate scalars
    * - ``alphabetar2xyz_vector``, ``xyz2alphabetar_vector``
      - Cubed-sphere :math:`(\alpha,\beta,r)` components to and from Cartesian,
        via the panel rotation matrices (:math:`\alpha,\beta` basis vectors are
        non-orthogonal)

Position transforms (``chi2xyz``, ``chi2llr``, ``chi2abr``) are in
``CORE/kernel/geometry/sci_chi_transform_mod.F90``, and account for the
rotation and stretching of the mesh.

Field-level kernels in ``CORE/kernel/geometry/``:

.. list-table::
    :header-rows: 1
    :widths: 40 60

    * - Kernel
      - Purpose
    * - ``sci_nodal_xyz_coordinates_kernel_mod``
      - Global :math:`(X,Y,Z)` coordinates at the nodal points of a target space
    * - ``sci_convert_cart2sphere_vector_kernel_mod``
      - In place: three Cartesian component fields to
        :math:`(\lambda,\phi,r)` components, given the nodal :math:`(X,Y,Z)`
        fields. Singular at the poles
    * - ``sci_pointwise_convert_xyz2llr_kernel_mod``
      - In place: nodal :math:`(X,Y,Z)` fields to :math:`(\lambda,\phi,r)`

:math:`\mathbb{W}_2` and :math:`\mathbb{W}_1` to components
-----------------------------------------------------------

Piola transforms, in ``CORE/kernel/inter_function_space/``. Each returns three
fields holding Cartesian components at the nodal points of a target space.

.. list-table::
    :header-rows: 1
    :widths: 40 60

    * - Kernel
      - Transform
    * - ``sci_convert_hdiv_field_kernel_mod``
      - :math:`\mathbb{W}_2`, any target space, any order and coordinate
        system: :math:`\mathbf{u}=J\hat{\mathbf{u}}/\det J`
    * - ``sci_convert_hdiv_native_kernel_mod``
      - As above, with :math:`\mathbb{W}_2` as the target; lowest order and
        native coordinates only (faster, whole column at once)
    * - ``sci_dg_convert_hdiv_field_kernel_mod``,
        ``sci_dg_convert_hdiv_native_kernel_mod``
      - As above, with horizontally discontinuous targets
        (:math:`\mathbb{W}_3`, :math:`\mathbb{W}_\theta`; native version
        :math:`\mathbb{W}_3` only)
    * - ``sci_convert_hcurl_field_kernel_mod``
      - :math:`\mathbb{W}_1`: :math:`\boldsymbol\xi=J^{-T}\hat{\boldsymbol\xi}`

Callers choose between the general and ``_native`` kernels with the same test:
``coord_system == coord_system_native`` and lowest-order in both directions.

.. important::

    The continuous-target kernels (``convert_hdiv_*``, ``convert_hcurl``)
    increment the target, so values at shared DoFs are summed over cells. The
    result **must be multiplied by the reciprocal multiplicity**
    (``get_rmultiplicity_fv`` for lowest order, ``_fe`` otherwise). The ``dg_``
    kernels overwrite and must **not** be rescaled.

Components to :math:`\mathbb{W}_2`
----------------------------------

.. list-table::
    :header-rows: 1
    :widths: 40 60

    * - Routine
      - Purpose
    * - ``sci_convert_phys_to_hdiv_kernel_mod``
      - Pointwise: :math:`(u_\lambda,u_\phi,u_r)` held in three
        :math:`\mathbb{W}_2`-sized fields to a :math:`\mathbb{W}_2` field,
        :math:`\hat{\mathbf{u}}=\det J\,J^{-1}\mathbf{u}`
    * - ``sci_compute_map_u_operators_kernel_mod``
      - Operators for a Galerkin projection of lon/lat components in
        :math:`\mathbb{W}_3` and radial in :math:`\mathbb{W}_\theta` to
        :math:`\mathbb{W}_2`. Lowest order only
    * - ``sci_compute_sample_u_ops_kernel_mod``
      - Operators for *sampling* the same components into
        :math:`\mathbb{W}_{2}` (broken): :math:`\det J\,\mathbf{n}\cdot
        J^{-1}\mathbf{e}_{\lambda|\phi|r}` at each DoF. Lowest order only
    * - ``sci_project_ws_to_w2_operator_kernel_mod``
      - Operators projecting Cartesian components in :math:`\mathbb{W}_3`
        to :math:`\mathbb{W}_2`
    * - ``sci_project_ws_to_w1_operator_kernel_mod``
      - As above, to :math:`\mathbb{W}_1`

The operators are created and cached per mesh in
``CORE/algorithm/sci_mapping_constants_mod.x90``
(``get_project_{x,y,z}dot_to_w2``, ``get_u_{lon,lat,up}_map``,
``get_u_{lon,lat,up}_sample``). The Galerkin and sampling versions are
selected by the ``sample_physics_winds`` option. The sampling approach applies
the operators with ``dg_matrix_vector_kernel`` and then
``average_w2b_to_w2_kernel`` to restore continuity; the Galerkin approach
assembles a right-hand side and applies the :math:`\mathbb{W}_2` mass-matrix
solver.

Splitting by direction
----------------------

These do not change basis; they select the horizontal and vertical DoFs
(lowest order only, using the ``face_selector_ew/ns`` fields):
``split_w2_field_kernel_mod`` and ``combine_w2_field_kernel_mod``
(wrapped in ``CORE/algorithm/sci_split_combine_w2_alg_mod.x90``),
``sci_extract_uv_kernel_mod`` (to :math:`\mathbb{W}_{2H}`) and
``sci_extract_w_kernel_mod`` (to :math:`\mathbb{W}_\theta`).

Algorithms using the transforms
-------------------------------

.. list-table::
    :header-rows: 1
    :widths: 35 65

    * - Algorithm
      - Sequence
    * - ``map_physics_winds``
        (``GUNGHO/algorithm/physics/physics_mappings_alg_mod.x90``)
      - ``dg_convert_hdiv_*`` to :math:`\mathbb{W}_3`/:math:`\mathbb{W}_\theta`,
        then (spherical) ``nodal_xyz_coordinates`` and
        ``convert_cart2sphere_vector``
    * - ``split_wind_alg`` (same file)
      - ``convert_hdiv_*`` to :math:`\mathbb{W}_2` components, rescale by
        multiplicity, ``convert_cart2sphere_vector``, then ``extract_w`` and
        ``extract_uv``
    * - ``set_wind``
        (``GUNGHO/algorithm/physics/map_fd_to_prognostics_alg_mod.x90``)
      - Physics components to :math:`\mathbb{W}_2` by sampling, or by Galerkin
        projection. With ``sample_physics_winds_correction`` (spherical,
        fully periodic only) a correction is built with
        ``w3_to_w2_correction`` and ``convert_phys_to_hdiv``
    * - ``interp_w2_to_w3wth_alg``, ``interp_w3wth_to_w2_alg``
        (``GUNGHO/algorithm/interpolation/interpolation_alg_mod.x90``)
      - The two recipes above, in each direction
    * - Advective-form wind transport
        (``GUNGHO/algorithm/transport/control/wind_transport_alg_mod.x90``)
      - ``dg_convert_hdiv_*`` to Cartesian components in :math:`\mathbb{W}_3`
        (not converted to lon/lat), transport as scalars, then project the
        increment back with ``project_{x,y,z}dot_to_w2``
    * - ``diagnostic_alg_mod``
        (``GUNGHO/algorithm/diagnostics/``)
      - ``convert_hcurl_field``/``convert_hdiv_*``, rescale by multiplicity,
        ``convert_cart2sphere_vector`` and ``pointwise_convert_xyz2llr``

Analytic initial winds (e.g. ``initial_u_kernel``) are built in
:math:`(\lambda,\phi,r)` components and converted with ``sphere2cart_vector``
before being projected onto the model space.

Adjoint versions of the kernels are in ``lfric_apps/science/adjoint/source/``.
