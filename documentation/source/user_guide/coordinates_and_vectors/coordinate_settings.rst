.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
     Some of the content of this file has been produced with the assistance of
     Met Office GitHub Copilot Enterprise.
   ------------------------------------------------------------------------------

.. _coordinates_and_vectors_coordinate_settings:

Coordinate settings
===================

The coordinates of the mesh are stored in the :math:`\chi` coordinate fields,
which are used for computations throughout the model. Three options in the
``finite_element`` namelist control how these fields are represented:
``coord_system``, ``coord_space`` and ``coord_order``. A fourth,
``coord_order_nonprime``, sets the order on meshes other than the prime mesh.

For how vector components are converted between coordinate systems, see
:ref:`coordinates_and_vectors_vector_component_transforms`.

``coord_system``
----------------

Selects which coordinate system is stored in the :math:`\chi` fields.

.. list-table::
    :header-rows: 1
    :widths: 20 80

    * - Value
      - Meaning
    * - ``'xyz'``
      - Global Cartesian coordinates :math:`(X,Y,Z)`
    * - ``'native'``
      - The system native to the mesh. This is :math:`(\alpha,\beta,h)` on
        cubed-sphere meshes and :math:`(\lambda,\phi,h)` on spherical
        limited-area meshes, where :math:`h` is the radius minus the planet
        radius

On planar meshes the coordinates are Cartesian in either case.

The choice also determines the vertical order of the coordinate space (see
``coord_order`` below).

``coord_space``
---------------

Selects the finite element space that holds the :math:`\chi` fields. The
spaces differ in their continuity:

.. list-table::
    :header-rows: 1
    :widths: 15 40 45

    * - Value
      - Continuity
      - Restrictions
    * - ``'Wchi'``
      - Fully discontinuous
      - None
    * - ``'Wtheta'``
      - Discontinuous in the horizontal, continuous in the vertical
      - Only appropriate for ``'native'`` coordinates, which separate the
        horizontal and vertical
    * - ``'W0'``
      - Fully continuous
      - Only valid for non-periodic planar domains, or for ``'xyz'``
        coordinates on a spherical domain

The model stops with an error if ``'W0'`` is chosen for any other
combination of domain and ``coord_system``.

``coord_order``
---------------

The horizontal polynomial order of the coordinate space on the prime mesh.
It must be at least 1. Order 2 is recommended for accurate representation of
orography.

The vertical order is not set directly. It follows from ``coord_system`` and
``coord_space``:

.. list-table::
    :header-rows: 1
    :widths: 25 25 50

    * - ``coord_system``
      - ``coord_space``
      - Vertical order
    * - ``'xyz'``
      - any
      - Same as the horizontal order
    * - ``'native'``
      - ``'Wchi'``
      - 1
    * - ``'native'``
      - ``'Wtheta'``
      - 0

For ``coord_space='W0'`` the value of ``coord_order`` is reduced by one before
use, because a W0 polynomial is one order above the W3 polynomial of the same
nominal order.

``'native'`` with ``'W0'`` (only permitted on non-periodic planar domains)
has vertical order 0.

``coord_order_nonprime``
------------------------

Horizontal order for the coordinate spaces of all meshes other than the prime
mesh, for example coarser meshes where higher-order coordinates may be
unnecessary. The vertical order is set in the same way as for ``coord_order``.
