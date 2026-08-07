.. _tangent_linear_tests:


Tangent Linear Tests
====================

Taylor remainder convergence test
---------------------------------

By definition, the tangent linear model is the Jacobian of the nonlinear
model, and therefore the Taylor series expansion can be used to form a
test for the code.

The Taylor series expansion of :math:`N(\mathbf{x})` is

.. math:: N(\mathbf{x} + \gamma \mathbf{\delta x} ) = N(\mathbf{x}) + \gamma \mathbf{L \delta x} + \mathcal{O}(\gamma^2)

This is used to calculate the norm of the difference between the
nonlinearly evolved perturbation
:math:`N(\mathbf{x} + \gamma \mathbf{\delta x} ) - N(\mathbf{x})` and
the linearly evolved perturbation :math:`\gamma \mathbf{L \delta x}`

.. math:: E_n(\gamma) = || N(\mathbf{x} + \gamma \mathbf{\delta x}) - N(\mathbf{x}) - \gamma \mathbf{L \delta x} || = \mathcal{O}(\gamma^2)

where the norm could be the L2-norm,
:math:`||\mathbf{x}|| = (\mathbf{x}^T\mathbf{x})^{1/2}`.

The norm :math:`E_n(\gamma)` is known as the linearization error, and it
tends to zero as the size of the perturbation :math:`\gamma` tends to
zero, at order :math:`\mathcal(O)(\gamma^2)`.

If we examine the series of :math:`E_n` for values of :math:`\gamma`
that are successively halved, we can evaluate the convergence rate. For
example, with the values of :math:`E_n(\gamma)` and :math:`E_n(2\gamma)`
we have

.. math:: E_n(2\gamma) / E_n(\gamma) = \mathcal{O}(4 \gamma^2) / \mathcal{O}(\gamma^2) = 4

To create a test which gives a simple pass or fail, two values of the
sequence :math:`E_n` are chosen, and the ratio is compared to the value
:math:`4`. If the value of the ratio is :math:`4`, to within a specified
tolerance, then the code passes the test.

Implementation
--------------

The Taylor remainder convergence test is implemented in the LFRic code
through the use of integration tests. In the linear app integration
tests, the convergence test is applied to small sections of code. The
process for each test is

- Initialise the nonlinear model, linear model and linearisation state.

- Run the nonlinear model using the linearisation state as initial
  conditiosn.

- Run the nonlinear model with a perturbed linearisation state, and run
  the linear model with the perturbation, and calculate
  :math:`E_n(\gamma)`.

- Halve the size of the perturbation, and repeat.

- Compare the ratio of the values of :math:`E_n` with :math:`4`, to a
  specified tolerance.
