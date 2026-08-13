.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
     Some of the content of this file has been produced with the assistance of
     Met Office GitHub Copilot Enterprise.
   ------------------------------------------------------------------------------

.. _user_guide_timestepping:

Timestepping
============

This page describes the namelist options in the ``timestepping`` namelist,
found in the ``lfric-gungho`` rose-meta metadata
(``science/gungho/rose-meta/lfric-gungho/HEAD/rose-meta.conf``), which control
how the GungHo dynamical core (and its coupling to physics schemes) is advanced
in time. A description of the underlying numerical schemes can be found in the
:ref:`Science Guide <science_guide_timestepping>`.

Choosing a timestepping method
-------------------------------

The ``method`` option selects the overall timestepping algorithm. The
available values are:

* ``semi_implicit``: the Semi-Implicit Quasi-Newton (SIQN) scheme. This is
  GungHo's default timestepping scheme, evolved from the ENDGame dynamical
  core's timestepping scheme. It combines an off-centred Crank-Nicolson scheme
  with a Semi-Lagrangian treatment of transport (although in GungHo transport
  has an Eulerian formulation). The scheme consists of a nested outer loop
  (which includes transport and fast physics), and an inner loop which calls
  the linear solver. An off-centering parameter, ``alpha``, controls the
  implicitness of the scheme.
* ``tr_bdf2``: a two-stage semi-implicit method (Trapezoidal, then Backwards
  Difference), using a Quasi-Newton iterative structure like the
  Semi-Implicit scheme. Each of the two stages takes the form of a
  Semi-Implicit step. There is no off-centering parameter, which improves the
  accuracy of the scheme relative to Semi-Implicit.
* ``rk``: an explicit multi-stage Runge-Kutta method. This has not been
  implemented with physics schemes.
* ``jules``: a timestepping scheme which calls the standalone JULES surface
  scheme without coupling to a dynamical core.
* ``no_timestepping``: a null timestepping scheme, which does not call the
  dynamics or physics schemes. Time-varying fields such as ancillaries and
  LBCs are still updated on each timestep, and diagnostics are output.

Iteration counts
-----------------

The Semi-Implicit and TR-BDF2 schemes both involve nested outer and inner
iterative loops. The outer loop includes transport and fast physics terms; the
inner loop calls the implicit forcing terms and the linear solver.

* ``outer_iterations``: the number of outer iterations to perform. This
  applies to both the Semi-Implicit and TR-BDF2 schemes (and to the ``jules``
  method).
* ``inner_iterations_si``: the number of inner Newton iterations to perform
  within the Semi-Implicit scheme.
* ``inner_iterations_tr``: the number of inner Newton iterations to perform in
  the Trapezoidal (TR) stage of the TR-BDF2 scheme.
* ``inner_iterations_bdf2``: the number of inner Newton iterations to perform
  in the Backwards Difference (BDF2) stage of the TR-BDF2 scheme.

Semi-Implicit off-centering and relaxation parameters
-------------------------------------------------------

These options are only relevant when ``method='semi_implicit'``.

* ``alpha``: the off-centering parameter for the Semi-Implicit scheme, in the
  range ``0.0`` to ``1.0``. A value of ``0.5`` corresponds to a
  Crank-Nicolson scheme, while a value of ``1.0`` corresponds to a fully
  implicit scheme. A value of ``0.5`` can be appropriate for dynamics-only
  runs, but is generally unstable when orography and physics are included.
* ``spinup_alpha``: if ``.true.``, the off-centering parameter is set to
  ``1`` (fully implicit) during the spin-up period, to improve stability.

The relaxation parameters give the implicitness weighting of semi-implicit
terms (for example transport or forcing terms) on the left-hand side of the
linear solver. See the :ref:`Science Guide <science_guide_timestepping>` for
the role that these parameters play in the linear solver formulation.

* ``tau_u``: relaxation parameter for the momentum equation. It is
  recommended to set ``tau_u=alpha``.
* ``tau_r``: relaxation parameter for density. When using more than one
  Quasi-Newton inner iteration, it is recommended to set ``tau_r=1``, as the
  density transport term is not updated between inner iterations.
* ``tau_t``: relaxation parameter for potential temperature. When using more
  than one Quasi-Newton inner iteration, it is recommended to set
  ``tau_t=1``, as the potential temperature transport term is not updated
  between inner iterations.

Runge-Kutta options
--------------------

This option is only relevant when ``method='rk'``.

* ``runge_kutta_method``: selects the choice of explicit Runge-Kutta scheme
  used either to timestep the whole dynamical core, or for use with the
  Method-of-Lines transport scheme. The Strong Stability Preserving (SSP)
  methods are explicit multi-stage methods which combine the values over
  multiple stages as a weighted linear sum to provide stability. Available
  values are ``forward_euler``, ``ssp2``, ``ssp3``, ``ssp4`` and ``ssp5``.
