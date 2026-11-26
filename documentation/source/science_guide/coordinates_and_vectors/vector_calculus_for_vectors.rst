.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
     Some of the content of this file has been produced with the assistance of
     Met Office GitHub Copilot Enterprise.
   ------------------------------------------------------------------------------

.. _coordinates_and_vectors_vector_calculus:

Vector Calculus for Vectors
===========================

Transporting or diffusing a *vector* is more involved than doing so for a
scalar, because the basis vectors of a curvilinear coordinate system change
from place to place. This page derives the evolution equation for each
contravariant component of a vector :math:`\mathbf{u}` that is both advected
and diffused, and shows how the extra terms (the Christoffel symbols) are
obtained from the coordinate Jacobian and the metric tensor described in
:ref:`coordinates_and_vectors_coordinate_systems`.

The equation
------------

We consider the vector transport equation with a conservative diffusion term,

.. math::
    :label: vector_equation

    \frac{\partial\mathbf{u}}{\partial t}
    + \left(\mathbf{a}\cdot\nabla\right)\mathbf{u}
    = \frac{1}{\rho}\nabla\cdot\left(\rho\nu\nabla\mathbf{u}\right),

where :math:`\mathbf{a}` is the advecting wind (which may be :math:`\mathbf{u}`
itself), :math:`\rho` is the dry density and :math:`\nu` is a diffusion
coefficient. For a Cartesian coordinate system, each Cartesian component of
:math:`\mathbf{u}` simply obeys the scalar transport-diffusion equation. In
curvilinear coordinates this is not true for the components of
:math:`\mathbf{u}`, and the equation for each component contains extra terms.

Notation
--------

Coordinates are written :math:`\chi_i` (the coordinate system label is dropped
here), with :math:`\partial_i = \partial/\partial\chi_i`. As in
:ref:`coordinates_and_vectors_coordinate_systems`:

* the basis vectors :math:`\mathbf{e}_i` are the columns of the Jacobian
  :math:`J` from the coordinates to Cartesian coordinates, and the dual basis
  vectors :math:`\mathbf{e}^i` are the rows of :math:`J^{-1}`;
* the metric tensor is :math:`g_{ij}=\mathbf{e}_i\cdot\mathbf{e}_j`, i.e.
  :math:`g = J^TJ`, and its inverse is :math:`g^{ij}=\mathbf{e}^i\cdot\mathbf{e}^j
  = (J^{-1}J^{-T})_{ij}`;
* :math:`\sqrt{g}=\det J` is the volume of a cell with unit coordinate
  extents.

Vectors are written with contravariant components,
:math:`\mathbf{u}=\sum_i u^i\mathbf{e}_i`, and similarly
:math:`\mathbf{a}=\sum_i a^i\mathbf{e}_i`. Summations are written out
explicitly, rather than using the Einstein summation convention.

The contravariant component :math:`u^i` is obtained by taking the dot product
of a vector with :math:`\mathbf{e}^i`. The component form of
equation :eq:`vector_equation` is found by taking this dot product with every
term in the equation.

Christoffel symbols
-------------------

The basis vectors vary in space, so the derivative of a vector has
contributions from the change in the components *and* from the change in the
basis vectors:

.. math::

    \partial_j\mathbf{u} = \sum_k\left(\partial_j u^k\,\mathbf{e}_k
    + u^k\,\partial_j\mathbf{e}_k\right).

The derivative of a basis vector is itself a vector, so can be expanded in
the basis. The expansion coefficients are the *Christoffel symbols*
:math:`\Gamma^i_{jk}`,

.. math::
    :label: christoffel_definition

    \partial_j\mathbf{e}_k = \sum_i\Gamma^i_{jk}\,\mathbf{e}_i
    \qquad\Longleftrightarrow\qquad
    \Gamma^i_{jk} = \mathbf{e}^i\cdot\partial_j\mathbf{e}_k .

Substituting this in gives the :math:`i`-th component of
:math:`\partial_j\mathbf{u}`, which is the *covariant derivative* of the
component :math:`u^i`:

.. math::
    :label: covariant_derivative

    \nabla_j u^i = \partial_j u^i + \sum_k\Gamma^i_{jk}u^k .

The extra term is the correction for the changing basis vectors. The
Christoffel symbols are symmetric in their lower indices, since
:math:`\partial_j\mathbf{e}_k = \partial_j\partial_k\boldsymbol{\chi}
= \partial_k\mathbf{e}_j`.

Relation to the Jacobian and metric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since :math:`[\mathbf{e}_k]_l=J_{lk}`, and the dual basis vectors are the rows
of :math:`J^{-1}`, equation :eq:`christoffel_definition` gives the Christoffel
symbols directly in terms of the *second* derivatives of the mapping to
Cartesian coordinates:

.. math::
    :label: christoffel_jacobian

    \Gamma^i_{jk} = \sum_l\left(J^{-1}\right)_{il}\,\partial_j J_{lk}
    = \sum_l\left(J^{-1}\right)_{il}
    \frac{\partial^2\chi^{(XYZ)}_l}{\partial\chi_j\,\partial\chi_k}.

Equivalently, they can be written using only the metric tensor,

.. math::
    :label: christoffel_metric

    \Gamma^i_{jk} = \frac{1}{2}\sum_m g^{im}\left(
    \partial_j g_{km} + \partial_k g_{jm} - \partial_m g_{jk}\right).

The two forms agree because :math:`\partial_j g_{km} = \partial_j
\left(\mathbf{e}_k\cdot\mathbf{e}_m\right)` can be expanded using
equation :eq:`christoffel_definition`. Two useful consequences are:

* If the Jacobian to Cartesian coordinates is constant (e.g. Cartesian
  coordinates on a plane), all Christoffel symbols vanish.
* Contracting over two indices gives (Jacobi's formula for the derivative of a
  determinant)

  .. math::
      :label: christoffel_trace

      \sum_j\Gamma^j_{jk} = \frac{1}{2}\sum_{j,m}g^{jm}\partial_k g_{jm}
      = \partial_k\ln\sqrt{g} = \partial_k\ln\det J,

  which is what allows divergences to be written in terms of :math:`\sqrt{g}`
  below.

**Example: spherical polar coordinates.** For :math:`(\lambda,\varphi,r)`
(longitude, latitude, radius) the metric is
:math:`g=\mathrm{diag}\left(r^2\cos^2\varphi,\;r^2,\;1\right)`
and :math:`\sqrt{g}=r^2\cos\varphi`. The non-zero Christoffel symbols are

.. math::

    \Gamma^\lambda_{\lambda\varphi} = -\tan\varphi, \quad
    \Gamma^\lambda_{\lambda r} = \frac{1}{r}, \quad
    \Gamma^\varphi_{\lambda\lambda} = \sin\varphi\cos\varphi, \quad
    \Gamma^\varphi_{\varphi r} = \frac{1}{r}, \quad
    \Gamma^r_{\lambda\lambda} = -r\cos^2\varphi, \quad
    \Gamma^r_{\varphi\varphi} = -r,

together with those related by the symmetry of the lower indices. For the
equiangular cubed-sphere coordinates :math:`(\alpha,\beta,r)` the metric is
not diagonal, and the Christoffel symbols are best evaluated from
equation :eq:`christoffel_jacobian` or :eq:`christoffel_metric`, using the
Jacobian given in :ref:`coordinates_and_vectors_equiangular_cubed_sphere_coordinates`.

Transport term
--------------

The directional derivative along :math:`\mathbf{a}` is
:math:`(\mathbf{a}\cdot\nabla)=\sum_j a^j\partial_j` (see
:ref:`coordinates_and_vectors_vector_components`). Applying this to the vector
:math:`\mathbf{u}` and using equation :eq:`covariant_derivative`, the
:math:`i`-th component of the transport term is

.. math::
    :label: transport_term

    \left[\left(\mathbf{a}\cdot\nabla\right)\mathbf{u}\right]^i
    = \sum_j a^j\nabla_j u^i
    = \sum_j a^j\partial_j u^i + \sum_{j,k}\Gamma^i_{jk}\,a^j u^k .

The first term is the usual advection of the component :math:`u^i` as if it
were a scalar. The second term is the "curvature" term, which couples the
components of :math:`\mathbf{u}` to each other. For example, transporting a
vector that is constant in a Cartesian sense appears to rotate when expressed
in spherical polar components.

Diffusion term
--------------

The gradient of a vector is a rank-2 tensor. With the index :math:`j`
denoting the gradient direction, its components are obtained by raising the
index on the covariant derivative,

.. math::

    \nabla^j u^i = \sum_m g^{jm}\nabla_m u^i
    = \sum_m g^{jm}\left(\partial_m u^i + \sum_k\Gamma^i_{mk}u^k\right).

The metric tensor is needed here because the gradient involves the *dual*
basis, :math:`\mathbf{e}^j = \sum_m g^{jm}\mathbf{e}_m`. We define the
diffusive flux tensor

.. math::
    :label: flux_definition

    F^{ij} := \rho\nu\,\nabla^j u^i
    = \rho\nu\sum_m g^{jm}\left(\partial_m u^i
    + \sum_k\Gamma^i_{mk}u^k\right).

The divergence over the gradient index of a rank-2 tensor :math:`T^{ij}` is

.. math::

    \left[\nabla\cdot T\right]^i = \sum_j\nabla_j T^{ij}
    = \sum_j\partial_j T^{ij} + \sum_{j,k}\Gamma^j_{jk}T^{ik}
    + \sum_{j,k}\Gamma^i_{jk}T^{kj},

where there is one Christoffel term for each index. Using
equation :eq:`christoffel_trace` the first two terms combine to give

.. math::

    \sum_j\partial_j T^{ij} + \sum_{k}\left(\partial_k\ln\sqrt{g}\right)T^{ik}
    = \frac{1}{\sqrt{g}}\sum_j\partial_j\left(\sqrt{g}\,T^{ij}\right),

hence the :math:`i`-th component of the diffusion term is

.. math::
    :label: diffusion_term

    \frac{1}{\rho}\left[\nabla\cdot\left(\rho\nu\nabla\mathbf{u}\right)\right]^i
    = \frac{1}{\rho\sqrt{g}}\sum_j\partial_j\left(\sqrt{g}\,F^{ij}\right)
    + \frac{1}{\rho}\sum_{j,k}\Gamma^i_{jk}F^{kj}.

The first term is the divergence of a flux in the coordinate system, as would
be used for a scalar, with :math:`\sqrt{g}=\det J` converting between the
coordinate-space and physical volumes. The second term is a further curvature
term, which arises because the divergence is taken of a flux that is itself a
vector.

Evolution equation for each component
-------------------------------------

Collecting equations :eq:`transport_term` and :eq:`diffusion_term`, the
evolution of the :math:`i`-th contravariant component of
equation :eq:`vector_equation` is

.. math::
    :label: component_equation

    \frac{\partial u^i}{\partial t}
    + \sum_j a^j\left(\partial_j u^i + \sum_k\Gamma^i_{jk}u^k\right)
    = \frac{1}{\rho\sqrt{g}}\sum_j\partial_j\left(\sqrt{g}\,F^{ij}\right)
    + \frac{1}{\rho}\sum_{j,k}\Gamma^i_{jk}F^{kj},

with :math:`F^{ij}` given by equation :eq:`flux_definition`. Written out in
full, the flux is

.. math::

    \sqrt{g}\,F^{ij} = \sqrt{g}\,\rho\nu\sum_m g^{jm}\left(\partial_m u^i
    + \sum_k\Gamma^i_{mk}u^k\right).

Several remarks:

* For a given :math:`i`, the quantities required are the component
  :math:`u^i` and its derivatives, the other components :math:`u^k` (through
  the Christoffel symbols), the inverse metric :math:`g^{jm}` and the
  determinant :math:`\sqrt{g}=\det J`. All of these follow from the Jacobian
  :math:`J` and its first derivatives, via
  equations :eq:`christoffel_jacobian` and :eq:`christoffel_metric`.
* Because the :math:`\mathbb{W}_2` degrees of freedom are
  :math:`\det J\,u^i`, the contravariant components are the natural variables
  for a model using :math:`\mathbb{W}_2` winds (see
  :ref:`coordinates_and_vectors_vector_components`).
* If the coordinates are orthogonal, :math:`g^{jm}` is diagonal and the
  flux simplifies considerably. If :math:`J` is also constant, the
  Christoffel symbols vanish and equation :eq:`component_equation` reduces
  to the scalar equation for each component.
* An approximation that is sometimes made is to treat each component as a
  scalar-like quantity, by dropping the final term in
  equation :eq:`diffusion_term`. This is exact only when
  :math:`\sum_{j,k}\Gamma^i_{jk}F^{kj}` vanishes (e.g. for Cartesian
  coordinates) and is otherwise a model simplification.

