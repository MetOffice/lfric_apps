.. _tangent_linear_overview:

Tangent Linear Model Overview
=======================

In variational data assimilation, the aim is to compute the optimal
value of the initial state :math:`\mathbf{x}_0` based on the
minimization of a cost function:

.. math:: J(\mathbf{x}_0) = (\mathbf{x}_0 - \mathbf{x}_b)^T \mathbf{B}^{-1} (\mathbf{x}_0 - \mathbf{x}_b) + \sum_t ( \mathbf{y}_t - \mathbf{H}_t N_t(\mathbf{x}_0))^T \mathbf{R}^{-1} ( \mathbf{y}_t - \mathbf{H}_t N_t(\mathbf{x}_0))

where :math:`t` is the time, :math:`\mathbf{x}_b` is the background
state, :math:`\mathbf{y}_t` are the observations, :math:`\mathbf{H}_t`
is the observation operator (which we can assume is linear here), and
:math:`N_t` is the nonlinear model that evolves the state
:math:`\mathbf{x}_0` from time :math:`0` to time :math:`t`.
:math:`\mathbf{B}` and :math:`\mathbf{R}` are the error covariance
matrices.

As :math:`N` is nonlinear, there is no unique solution. Therefore, the
problem is rewritten in incremental form, by rewriting the initial state
in terms of a linearisation state plus a perturbation,
:math:`\mathbf{x}_0 = \tilde{\mathbf{x}} + \mathbf{\delta x}`, so that
on the first iteration, when :math:`\tilde{\mathbf{x}} = \mathbf{x}_b`,

.. math:: J(\mathbf{\delta x}_0) = (\mathbf{\delta x}_0)^T B^{-1} (\mathbf{\delta x}_0) + \sum_t ( \mathbf{\delta y}_t - \mathbf{H}_t \mathbf{L}_t\mathbf{\delta x}_0)^T \mathbf{R}^{-1} ( \mathbf{\delta y}_t - \mathbf{H}_t \mathbf{L}_t \mathbf{\delta x}_0)

where :math:`\mathbf{\delta x}_0 = \mathbf{x}_0 - \mathbf{x}_b`, and
:math:`\mathbf{\delta y}_t = \mathbf{y}_t - \mathbf{H}_t N_t (x_b)`, and
we have used the Taylor series expansion to approximate the nonlinearly
evolved perturbation,

.. math:: N( \mathbf{x}_b + \mathbf{\delta x} ) - N(\mathbf{x}_b) \approx \mathbf{L\delta x}

where the linear model
:math:`\mathbf{L}= \mathbf{L}(\tilde{\mathbf{x}})`, known as the tangent
linear model, is a function of the linearisation state
:math:`\tilde{\mathbf{x}}`, and is the Jacobian or first derivative of
the nonlinear model.
