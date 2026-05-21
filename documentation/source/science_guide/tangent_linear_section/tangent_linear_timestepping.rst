.. _tangent_linear_timestepping:

Iterated (Semi-) Implicit Timestepping
====================

Nonlinear Equations
-------------------

The nonlinear equations are

.. math::

   \begin{aligned}
     \frac{\partial u}{\partial t} &= - u . \nabla u - 2 \Omega \times u - \nabla \Phi - c_P \theta \nabla \Pi \\
     \frac{\partial \rho}{\partial t} & = - \nabla . (\rho u) \\
     \frac{\partial \theta}{\partial t} & = - u. \nabla \theta
     \end{aligned}

together with the nonlinear equation of state:

.. math:: \Pi ^{\left( \frac{1-\kappa}{\kappa} \right)} = \frac{R}{p_0} \rho \theta

Linear Equations
----------------

The linearized equations are:

.. math::

   \begin{aligned}
     \frac{\partial \delta u}{\partial t} &= - \bar{u} . \nabla \delta u - \delta {u} . \nabla \bar{ u} - 2 \Omega \times \delta u - c_P \delta \theta \nabla \bar{\Pi} - c_P \bar{\theta} \nabla \delta \Pi \\
     \frac{\partial \delta \rho}{\partial t} & = - \nabla . (\bar{\rho} \delta u) - \nabla . (\delta \rho \bar{u}) \\
     \frac{\partial \delta \theta}{\partial t} & = - \bar{u}. \nabla \delta \theta - \delta u. \nabla \bar{\theta}
     \end{aligned}

together with the linear equation of state:

.. math:: \frac{ \delta \Pi}{\bar{\Pi}}  = \frac{\kappa}{1-\kappa} \left( \frac{\delta \rho}{\bar{\rho}} + \frac{\delta \theta}{\bar{\theta}} \right)


Iterated (Semi-) implicit Timestepping
--------------------------

In the iterated-implicit timestepping, the equations are solved iteratively
using a Newton method. If the full nonlinear discrete equations are
written as

.. math:: \mathcal{R} (\mathbf{x}^{n+1}) =0

then this can be solved iteratively as

.. math::

   \begin{aligned}
     \mathbf{x}^{n+1}_{(k+1)} & =  \mathbf{x}^{n+1}_{(k)} +  \mathbf{x}'_{(k)} \\
     \mathcal{L} (\mathbf{x}^*) \mathbf{x}'_{(k)} & = - \mathcal{R} ( \mathbf{x}^{n+1}_{(k)})
     \end{aligned}

where :math:`\mathcal{L}` is an approximation to the Jacobian and
:math:`\mathbf{x}^*` is the reference state.

The reference state :math:`u^*` is zero whilst the other reference state
variables are set to :math:`\mathbf{x}^n`. Relaxation parameters
:math:`\tau <1` are also introduced to give extra weight to the diagonal
terms (and hence similar to solving the resulting linear equation using
the Jacobi method or a damped Newton method).

This gives the approximation to the Jacobian as

.. math::

   \mathcal{L}(\mathbf{x}^*) \mathbf{x}' = \left\{
     \begin{aligned}
       u' & + \tau_u \Delta t c_P (\theta' \nabla \Pi^* + \theta^* \nabla \Pi ') \\
       \rho' & + \tau_\rho \Delta t \nabla.(\rho^* u') \\
       \theta' & + \tau_\theta \Delta t u' \nabla. \nabla \theta^* \\
   & \frac{1-\kappa}{\kappa}  \frac{\Pi'}{\Pi^*}  - \frac{1}{E} \frac{\rho'}{\rho^*} - \frac{1}{E} \frac{\theta'}{\theta^*}
     \end{aligned}
   \right.

where

.. math:: E = \frac{p_0}{R} \frac{\Pi^{*(1-\kappa)/\kappa}}{\theta^* \rho^*}

We can proceed in a similar way with the linear model, so that

.. math::

   \begin{aligned}
     \delta \mathbf{x}^{n+1}_{(k+1)} & =  \delta \mathbf{x}^{n+1}_{(k)} +  \delta \mathbf{x}'_{(k)} \\
     \mathcal{L} (\mathbf{x}^*) \delta \mathbf{x}'_{(k)} & = - R ( \mathbf{x}_{(k)}) \delta \mathbf{x}_{(k)}
     \end{aligned}

where :math:`R` is the linearized version of :math:`\mathcal{R}`.

As the operator :math:`\mathcal{L}` is already linear, this does not
need to be linearized.

Discretisation
--------------------

The nonlinear model evolves the linearisation state
:math:`\bar{\mathbf{x}}`. The aim is to find the state
:math:`\bar{\mathbf{x}}^n` at time :math:`t_{n+1}`. This is solved
iteratively, by letting
:math:`{\bar{\mathbf{x}}}^{n+1}_{(k+1)} = \bar{\mathbf{x}}^{n+1}_{(k)} + \bar{\mathbf{x}}'_{(k+1)}`
where :math:`k` is the iteration count.

The increment :math:`\bar{\mathbf{x}}'_{(k+1)}` solves the linear
equation

.. math:: \mathcal{L}( \bar{\mathbf{x}}*) \bar{\mathbf{x}}'_{(k+1)} = \mathcal{R}^n(\bar{\mathbf{x}}^n) + \mathcal{R}^{n+1}(\bar{\mathbf{x}}^{n+1}_{(k)}) + \mathcal{R}^{adv}(\bar{\mathbf{x}}^n, \bar{\mathbf{x}}^{n+1}_{(k)})

where :math:`\mathcal{L}( \bar{\mathbf{x}}*)` is a linear operator based
on the reference state :math:`\bar{\mathbf{x}}*`.

The procedure is:

- | Let :math:`\bar{\mathbf{x}}^{n+1}_{(0)} = \bar{\mathbf{x}}^n`
  | The first iteration is given by the state at the previous timestep.

- | Set
    :math:`\bar{\mathbf{x}}* \equiv (\bar{\rho}*, \bar{\theta}*, \bar{\Pi}* ) = (\bar{\rho}^n, \bar{\theta}^n, \bar{\Pi}^n)`
  | The reference state is given by the state at the previous timestep
    (for density, potential temperature and exner pressure only).

- | Compute :math:`\mathcal{R}^n(\bar{\mathbf{x}}^n)`
  | This is the RHS based on the previous timestep.

- Loop k=1, K

  - | Compute :math:`\mathcal{R}^{n+1}(\bar{\mathbf{x}}^{n+1}_{(k)})`
      and
      :math:`\mathcal{R}^{adv}(\bar{\mathbf{x}}^n, \bar{\mathbf{x}}^{n+1}_{(k)})`
    | This is the RHS and advective terms based on the current timestep.

  - Solve the linear equation for :math:`\bar{\mathbf{x}}'_{(k+1)}` and
    hence for :math:`\bar{\mathbf{x}}^{n+1}_{(k+1)}`.

The linear model evolves the perturbations :math:`\mathbf{x}`. To find
the state at time :math:`t_{n+1}`, we follow the same procedure as for
the nonlinear model, but using the tangent linear versions of the right
hand side in the equation.

Let
:math:`{\mathbf{x}}^{n+1}_{(k+1)} = \mathbf{x}^{n+1}_{(k)} + \mathbf{x}'_{(k+1)}`
where :math:`\mathbf{x}'_{(k+1)}` solves

.. math:: \mathcal{L} ( \bar{\mathbf{x}}*) \mathbf{x}'_{(k+1)} = R^n(\mathbf{x}^n, \bar{\mathbf{x}}^n) + R^{n+1}(\mathbf{x}^{n+1}_{(k)},\bar{\mathbf{x}}^{n+1}_{(k)}) + R^{adv}(\mathbf{x}^n, \mathbf{x}^{n+1}_{(k)},\bar{\mathbf{x}}^n, \bar{\mathbf{x}}^{n+1}_{(k)})

where :math:`R^n` is the tangent linear of :math:`\mathcal{R}^n` etc.
The tangent linear operators are a function of both the perturbation
:math:`\mathbf{x}` and the linearisation state :math:`\bar{\mathbf{x}}`.

The procedure is:

- | Let :math:`\mathbf{x}^{n+1}_{(0)} = \mathbf{x}^n`
  | The first iteration is given by the perturbation at the previous
    timestep.

- | Set
    :math:`\bar{\mathbf{x}}* \equiv (\bar{\rho}*, \bar{\theta}*, \bar{\Pi}* ) = (\bar{\rho}^n, \bar{\theta}^n, \bar{\Pi}^n)`
  | The reference state is given by the linearisation state at the
    previous timestep - so that the linear operator is identical to the
    that in the nonlinear model.

- Compute :math:`R^n(\mathbf{x}^n, \bar{\mathbf{x}}^n)`

- Loop k=1, K

  - Compute
    :math:`R^{n+1}(\mathbf{x}^{n+1}_{(k)},\bar{\mathbf{x}}^{n+1}_{(k)})`
    and
    :math:`R^{adv}(\mathbf{x}^n, \mathbf{x}^{n+1}_{(k)},\bar{\mathbf{x}}^n, \bar{\mathbf{x}}^{n+1}_{(k)})`

  - Solve for :math:`\mathbf{x}'_{(k+1)}` and hence for
    :math:`\mathbf{x}^{n+1}_{(k+1)}`.

The linear solver :math:`\mathcal{L} ( \bar{\mathbf{x}}*)` is exactly
the same as for the nonlinear model.

The right hand side terms :math:`R^n`, :math:`R^{n+1}` and
:math:`R^{adv}` require the linearisation state at the same iteration.
Therefore, for an exact tangent linear semi-implicit algorithm, the
nonlinear semi-implicit step needs to be run first, and the intermediate
``ls_state`` values need to be stored and used in the following tangent
linear semi-implicit step.

