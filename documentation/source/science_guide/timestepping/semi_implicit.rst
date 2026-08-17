.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
    Some of the content of this file has been produced with the assistance of
    Met Office GitHub Copilot Enterprise.
   -----------------------------------------------------------------------------

.. _science_guide_timestepping_semi_implicit:

Semi-Implicit
=============

Overview
--------

The Semi-Implicit Quasi-Newton (SIQN) scheme is GungHo's default timestepping
scheme. It belongs to the family of
Semi-Implicit Semi-Lagrangian (SISL) schemes, which have been used
successfully in operational numerical weather prediction for several decades.
The appeal of SISL schemes is that they permit a large, stable timestep while
retaining good accuracy. This is achieved by combining two ideas:

* a Crank-Nicolson-like (trapezoidal) treatment of the wave-like and forcing
  terms, which is unconditionally stable with respect to the timestep for
  wave equations; and
* a transport scheme with some Semi-Lagrangian aspect so that stability does
  not depend on the Courant number.

In a classical SISL scheme, the explicit :math:`n`-level values are evaluated
at the departure points, while the implicit :math:`(n+1)`-level values are
evaluated at the arrival points. This results in a discretisation of the
material derivative which preserves the numerical properties of the
Crank-Nicolson scheme: high frequency (acoustic) waves are damped while low
frequency waves (such as Rossby waves) are represented accurately.

GungHo builds on the classical SISL picture by providing local mass conservation
through the use of conservative transport schemes. Like the ENDGame dynamical
core of [Wood2014]_ used by the Unified Model, GungHo uses a
Quasi-Newton iteration to solve the non-linear system of equations arising from
the implicit treatment of the wave-like and forcing terms.

For a full description of the formulation of GungHo, see [Melvin2019]_
and [Melvin2024]_.

Details
-------

Basic formulation
~~~~~~~~~~~~~~~~~~

Given the off-centering parameter :math:`\alpha`, the SIQN scheme advances
the state from time level :math:`n` to time level :math:`n+1` according to:

.. math:: :label: siqn_eqn

   \boldsymbol{X}^{n+1} + \alpha\Delta t \mathcal{F}(\boldsymbol{X}^{n+1})
   = \mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1}}
   \left[\boldsymbol{X}^n - (1-\alpha)\Delta t \mathcal{F}(\boldsymbol{X}^n)\right],

where the transporting velocity is the average of the start- and
end-of-timestep winds,

.. math:: :label: siqn_transport_velocity

   \overline{\boldsymbol{u}}_n^{n+1} = \tfrac{1}{2}\left(\boldsymbol{u}^n + \boldsymbol{u}^{n+1}\right).

The forcing terms are handled with a weighting of :math:`\alpha` on the
implicit (:math:`n+1`) contribution and :math:`(1-\alpha)` on the explicit
(:math:`n`) contribution.

A value of :math:`\alpha=0.5` corresponds to a centred Crank-Nicolson scheme,
which is formally second-order accurate in time but does not damp any waves.
If orography is present, it can generate spurious small-scale waves which,
without any damping, may persist and pollute the solution, leading to
instability; centred SISL schemes can also suffer from orographic resonance.
For these reasons the scheme is run slightly off-centred
(:math:`\alpha > 0.5`), which introduces damping of high-frequency (typically
acoustic) modes at the cost of reducing the formal accuracy to first order.
A value of :math:`\alpha=1.0` corresponds to a fully implicit scheme.

Quasi-Newton iteration
~~~~~~~~~~~~~~~~~~~~~~~

Equation :eq:`siqn_eqn` is non-linear in the unknown state
:math:`\boldsymbol{X}^{n+1}`, and is solved iteratively. Let us use :math:`k`
to index the iteration, so that :math:`k`-th iteration of
:math:`\boldsymbol{X}^{n+1}` is :math:`\boldsymbol{X}^{n+1}_(k)`, where
:math:`\boldsymbol{X}^{n+1}_{(0)} = \boldsymbol{X}^n` is the first ''guess''.
The :math:`k`-th transported state :math:`\boldsymbol{X}^T_{(k)}` is introduced as:

.. math::

   \boldsymbol{X}^T_{(k)} = \mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1,(k)}}
   \left[\boldsymbol{X}^n - (1-\alpha)\Delta t \mathcal{F}(\boldsymbol{X}^n)\right],

which depends upon the latest value of the end-of-timestep transporting velocity.
Then the residual of the :math:`k`-th iterate :math:`\boldsymbol{X}^{n+1}_{(k)}` is

.. math:: :label: siqn_residual

   \boldsymbol{R}_{(k)} = -\left(\boldsymbol{X}^{n+1}_{(k)}
   + \alpha\Delta t \mathcal{F}(\boldsymbol{X}^{n+1}_{(k)}) - \boldsymbol{X}^T_{(k)}\right).

A full Newton method would drive this residual to zero by solving

.. math::

   \mathcal{J}\left[\boldsymbol{X}^{n+1}_{(k)}\right]
   \left(\boldsymbol{X}^{n+1}_{(k+1)} - \boldsymbol{X}^{n+1}_{(k)}\right)
   = \boldsymbol{R}_{(k)},

where :math:`\mathcal{J}` is the Jacobian of the residual with respect to
:math:`\boldsymbol{X}^{n+1}`. Forming and inverting the exact Jacobian is
expensive, so (like ENDGame) GungHo instead uses a *Quasi-Newton* method: the Jacobian is
approximated by a fixed linearisation :math:`\mathcal{S}` of the equation set
about a reference state. Writing the increment as
:math:`\boldsymbol{X}' = \boldsymbol{X}^{n+1}_{(k+1)} - \boldsymbol{X}^{n+1}_{(k)}`,
each iteration then solves the linear system

.. math:: :label: siqn_linear_system

   \mathcal{S}\,\boldsymbol{X}' = \boldsymbol{R}_{(k)}.

This "incremental" form is the one used in GungHo, and the scheme is not
iterated to convergence: a small, fixed number of iterations is taken. Although
this can still be accurate, it can be unstable for certain numbers of iterations.

Outer and inner iterations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In practice, the change to transport increments from one iteration to the next
is smaller than the change from the implicit forcing terms. This allows us
to split the Quasi-Newton iteration into two nested loops:

* an **outer loop**, denoted :math:`n_o`, in which the transporting velocity
  is updated with the latest estimate of :math:`\boldsymbol{u}^{n+1}` and the
  prognostic variables are transported. The transport scheme is therefore
  called :math:`n_o` times per timestep.
* an **inner loop**, denoted :math:`n_i`, in which the implicit forcing terms
  are evaluated to form the residual :eq:`siqn_residual`, and the linear
  system :eq:`linear_solver_eqn` is assembled and solved for the increment.
  The number of inner iterations can be varied between outer iterations, but
  has historically been kept the same for each outer iteration, so that the
  linear solve is performed :math:`n_o \times n_i` times per timestep.

Like ENDGame, GungHo uses a default of :math:`n_o=2` outer iterations and
:math:`n_i=2` inner iterations, for a total of 4 linear solves per timestep.

Because the transport terms are only updated in the outer loop, using more
than one inner iteration (:math:`n_i > 1`), the density and potential
temperature fields do not get updated between inner iterations. To keep the
discretisation consistent, when :math:`n_i > 1` the density and
potential temperature contributions to the residual are set to zero. This is why
it is recommended to set the density and potential-temperature relaxation
parameters to 1 when more than one inner iteration is used.

Linear solver
~~~~~~~~~~~~~

The linearised operator :math:`\mathcal{S}` in Equation
:eq:`siqn_linear_system` is obtained by linearising the residuals about a
reference state :math:`\boldsymbol{X}^* = (\boldsymbol{0}, \rho^*, \theta^*,
\Pi^*)^\mathrm{T}`. The resulting linear problem takes the general form

.. math:: :label: linear_solver_eqn

   \boldsymbol{X}' + \sigma \mathcal{I}[\boldsymbol{X}']
   + \tau_i \Delta t \mathcal{L}[\boldsymbol{X}'] = \boldsymbol{R},

where :math:`\mathcal{I}` is an implicit operator (corresponding to the
damping layer) and :math:`\mathcal{L}` is the linear operator corresponding
to the forcing and transport terms. The relaxation parameters
:math:`\sigma` and :math:`\tau_i` control the implicit weighting; the
:math:`\tau_i` may differ between the momentum, density and
potential-temperature variables :math:`i`. The implicit terms in the
continuity and thermodynamic equations correspond to linearisations of the
transport terms, :math:`\nabla\cdot(\rho^*\boldsymbol{u}')` and
:math:`\boldsymbol{u}'\cdot\nabla\theta^*` respectively.

For the SIQN scheme, the relaxation parameters are typically chosen as:

.. math::

   \tau_u = \alpha, \qquad \tau_\rho=\tau_\theta=1, \qquad \sigma=1.

Right-Hand Side and Predictors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As GungHo uses a finite element discretisation, the right-hand side of the
linear system :eq:`linear_solver_eqn` is actually calculated by assembling the
residual in weak form. When field values are needed, for instance for the
transport operator or physics parametrisations, it may be necessary to
solve a mass matrix system to obtain the field. This is known as calculating
the "predictors".

Dynamical Core Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~

We now bring together the outer loop (indexed by :math:`j=1,\dots,n_o`) and
inner loop (indexed by :math:`k=1,\dots,n^j_i`) to summarise the whole
Semi-Implicit Quasi-Newton timestep as a sequence of equations. The state
:math:`\boldsymbol{X}^{n+1}_{(j,k)}` denotes the :math:`k`-th inner-loop
estimate of the end-of-timestep state within the :math:`j`-th outer loop, with
:math:`\boldsymbol{X}^{n+1}_{(j,0)} = \boldsymbol{X}^{n+1}_{(j-1)}`.
The state :math:`\boldsymbol{X}^{n+1}_{(j)}` is the final inner-loop estimate of the
end-of-timestep state within the :math:`j`-th outer loop, and is used for
carrying the final state of the previous outer loop forward.
Similarly, :math:`\boldsymbol{X}^T_{(j)}` is the transported state formed
within outer loop :math:`j`.

The first guess for the end-of-timestep state is simply the start-of-timestep
state,

.. math:: :label: dyn_core_first_guess

   \boldsymbol{X}^{n+1}_{(0)} = \boldsymbol{X}^n.

The explicit forcing contribution only depends on the start-of-timestep state,
so it is evaluated once (in weak form), before the outer loop begins,

.. math:: :label: dyn_core_explicit_forcing

   \boldsymbol{X}^n_{FE} = \boldsymbol{X}^n - (1-\alpha)\Delta t\,\mathcal{F}(\boldsymbol{X}^n).

As :math:`\boldsymbol{X}^n_{FE}` is generally calculated in weak form, a
mass matrix system is solved to obtain the "predictor" values that are
subsequently transported.

Each outer loop :math:`j` begins by updating the transporting velocity using
the latest outer-loop estimate of the end-of-timestep wind,
:math:`\boldsymbol{u}^{n+1}_{(j-1)}`,

.. math:: :label: dyn_core_transport_velocity

   \overline{\boldsymbol{u}}_n^{n+1,(j)} = \tfrac{1}{2}\left(\boldsymbol{u}^n + \boldsymbol{u}^{n+1}_{(j-1)}\right).

The transport step is then performed using this updated velocity,

.. math:: :label: dyn_core_transport

   \boldsymbol{X}^T_{(j)} = \mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1,(j)}}
   \left[\boldsymbol{X}^n_{FE}\right].

Note that :math:`\boldsymbol{X}^T_{(j)}` is also generally calculated in weak
form.

Within outer loop :math:`j`, each inner iteration :math:`k` first calculates
the implicit forcing contribution in weak form,

.. math:: :label: dyn_core_implicit_forcing

   \boldsymbol{X}^{n+1}_{(j,k-1),FI} = \boldsymbol{X}^{n+1}_{(j,k-1)}
   + \alpha\Delta t\,\mathcal{F}\!\left(\boldsymbol{X}^{n+1}_{(j,k-1)}\right).

The residual of the current iterate against this fixed transported state is then

.. math:: :label: dyn_core_residual

   \boldsymbol{R}_{(j,k)} = -\left(\boldsymbol{X}^{n+1}_{(j,k-1),FI}
   - \boldsymbol{X}^T_{(j)}\right),

solves the linearised system for the increment :math:`\boldsymbol{X}'`,

.. math:: :label: dyn_core_linear_solve

   \mathcal{S}\left[\boldsymbol{X}'\right] = \boldsymbol{R}_{(j,k)},

and updates the state accordingly,

.. math:: :label: dyn_core_increment

   \boldsymbol{X}^{n+1}_{(j,k)} = \boldsymbol{X}^{n+1}_{(j,k-1)} + \boldsymbol{X}'.

Once all :math:`n_o` outer loops have completed, the timestep is advanced by
setting

.. math:: :label: dyn_core_advance

   \boldsymbol{X}^n = \boldsymbol{X}^{n+1}_{(n_o)}.

This process is summarised by the pseudo-code below:

.. code-block:: text

   Set the first guess:
      X^{n+1}_(0) = X^n
   Explicit forcing:
      X_FE = X^n - (1 - alpha) * dt * F(X^n)
   Calculate transport predictor from X_FE

   # Outer Loop ----------------------------------------------------------------
   for j = 1 to n_o:
       Update transporting velocity:
           u_bar^(j) = 0.5 * (u_n + u_{n+1}^(j-1))
       Transport:
           X^T_(j) = Transport_{u_bar^(j)}(X_FE)

       X^{n+1}_(j,0) = X^{n+1}_(j-1)

       # Inner Loop ------------------------------------------------------------
       for k = 1 to n_i^j:
           Implicit forcing:
               X_FI^(j,k-1) = X^{n+1}_(j,k-1) + alpha * dt * F(X^{n+1}_(j,k-1))
           Form the residual:
               R_(j,k) = -( X_FI^(j,k-1) - X^T_(j) )
           Solve the linear system for the increment:
               S[X'] = R_(j,k)
           Increment the state:
               X^{n+1}_(j,k) = X^{n+1}_(j,k-1) + X'

       Set final estimate:
         X^{n+1}_(j) = X^{n+1}_(j,n_i^j)

   Advance the timestep:
      X^n = X^{n+1}_(n_o)


Inclusion of Physics Parametrisations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

So far the description of the timestep has been limited to the dynamical core.
We now describe how physics parametrisation schemes are included.
Physics parametrisations can be divided into three categories:

1. **Slow Physics**: denoted :math:`\mathcal{P}_S`. This is evaluated on the
   start-of-timestep state :math:`\boldsymbol{X}^n`. It represents explicit
   physics processes which are not strongly coupled to the dynamics. Note that
   explicit forcing terms are not evaluated on the slow physics increments.
2. **Fast Physics**: denoted :math:`\mathcal{P}_F`. This is evaluated on the
   "physics predictor", which is the latest estimate of the end-of-timestep
   state following the evaluation of transport terms. It describes implicit
   physics processes that are strongly coupled to the dynamics, and is called
   within the outer loop of the Quasi-Newton iteration.
3. **End-of-step Physics**: denoted :math:`\mathcal{P}_E`. Evaluated after the
   Quasi-Newton iteration has completed, it also represents explicit physics
   processes.

The whole system of equations then becomes:

.. math:: :label: siqn_with_physics_eqn

   \begin{aligned}
   \boldsymbol{X}^\dagger &= \boldsymbol{X}^n - \Delta t\,\mathcal{P}_S(\boldsymbol{X}^n), \\
   \boldsymbol{X}^{n+1} + \alpha\Delta t \mathcal{F}(\boldsymbol{X}^{n+1})
   + \mathcal{P}_F(\boldsymbol{X}^{n+1})
   &= \mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1}}
   \left[\boldsymbol{X}^\dagger - (1-\alpha)\Delta t \mathcal{F}(\boldsymbol{X}^n)
   \right], \\
   \boldsymbol{X}^{n+1} &= \boldsymbol{X}^{n+1} - \Delta t\,\mathcal{P}_E(\boldsymbol{X}^{n+1}).
   \end{aligned}

With the inclusion of these physics parametrisation schemes, the full timestep
can be summarised by the following pseudocode.

.. code-block:: text

   Set the first guess:
      X^{n+1}_(0) = X^n
   Slow physics:
      X^dag = X^n + dt * P_S(X^n)
   Explicit forcing:
      X_FE = X^dag - (1 - alpha) * dt * F(X^n)
   Calculate transport predictor (X_FE)

   # Outer Loop ----------------------------------------------------------------
   for j = 1 to n_o:
       Update transporting velocity:
           u_bar^(j) = 0.5 * (u_n + u_{n+1}^(j-1))
       Transport:
           X^T_(j) = T_{u_bar^(j)}(X_FE)
       Implicit forcing:
           X_FI^(j,1) = X^{n+1}_(j-1) + alpha * dt * F(X^{n+1}_(j-1))
       Fast physics predictor:
           X_star^(j) = X^T_(j) + (X_FI^(j,0) - X^{n+1}_(j-1))
       Fast physics increment:
           dX_F^(j) = dt * P_F(X_star^(j))

       X^{n+1}_(j,0) = X^{n+1}_(j-1)

       # Inner Loop ------------------------------------------------------------
       for k = 1 to n_i^j:
           if k > 1:
               Implicit forcing:
                   X_FI^(j,k) = X^{n+1}_(j,k) + alpha * dt * F(X^{n+1}_(j,k-1))
           Form the residual:
               R_(j,k) = -( X_FI^(j,k) - X^T_(j) - dX_F^(j) )
           Solve the linear system for the increment:
               S[X'] = R_(j,k)
           Increment the state:
               X^{n+1}_(j,k) = X^{n+1}_(j,k-1) + X'

       Set final estimate:
         X^{n+1}_(j) = X^{n+1}_(j,n_i^j)

   End-of-step physics:
      X^{n+1} = X^{n+1}_(n_o) + dt * P_E(X^{n+1}_(n_o))
   Advance the timestep:
      X^n = X^{n+1}

Tracers
~~~~~~~

Some tracers (for instance chemical or aerosol species) do not feed back directly
onto the dynamical core and therefore do not need including in the Quasi-Newton
iteration. These tracers are transported in only the final outer iteration of
the timestep.

These tracers may also be updated by physics parametrisations (or for example a
chemical or aerosol scheme, which for now can be considered as a physics
parametrisation). Let tracers be represented by the state vector
:math:`\boldsymbol{Y}`. The evolution of tracers through the timestep is given by:

.. math:: :label: tracer_eqn

   \begin{aligned}
   \boldsymbol{Y}^\dagger &= \boldsymbol{Y}^n
      - \Delta t\,\mathcal{P}_S(\boldsymbol{X}^n; \boldsymbol{Y}^n), \\
   \boldsymbol{Y}^T &= \mathcal{T}^{\Delta t}_{\overline{\boldsymbol{u}}_n^{n+1}}
   \left[\boldsymbol{Y}^\dagger \right] + \mathcal{P}_F(\boldsymbol{X}^{n+1}; \boldsymbol{Y}^\dagger), \\
   \boldsymbol{Y}^{n+1} &= \boldsymbol{Y}^T - \Delta t\,\mathcal{P}_E(\boldsymbol{X}^{n+1};\boldsymbol{Y}^T).
   \end{aligned}

References
----------

.. [Wood2014] Wood, N., Staniforth, A., White, A., Allen, T., Diamantakis, M.,
   Gross, M., Melvin, T., Smith, C., Vosper, S., Zerroukat, M., and Thuburn, J.
   (2014). An inherently mass-conserving semi-implicit semi-Lagrangian
   discretization of the deep-atmosphere global non-hydrostatic equations.
   *Quarterly Journal of the Royal Meteorological Society*, 140(682), 1505-1520.
   https://doi.org/10.1002/qj.2235

.. [Melvin2019] Melvin, T., Benacchio, T., Shipway, B., Wood, N., Thuburn, J.,
   and Cotter, C. (2019). A mixed finite-element, finite-volume, semi-implicit
   discretization for atmospheric dynamics: Cartesian geometry.
   *Quarterly Journal of the Royal Meteorological Society*, 145(724), 2835-2853.
   https://doi.org/10.1002/qj.3501

.. [Melvin2024] Melvin, T., Shipway, B., Wood, N., Benacchio, T., Bendall, T.,
   Boutle, I., Brown, A., Johnson, C., Kent, J., Pring, S., et al. (2024). A
   mixed finite-element, finite-volume, semi-implicit discretisation for
   atmospheric dynamics: spherical geometry.
   *Quarterly Journal of the Royal Meteorological Society*, 150(764), 4252-4269.
   https://doi.org/10.1002/qj.4887