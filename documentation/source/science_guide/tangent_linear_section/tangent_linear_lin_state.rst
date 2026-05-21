.. _tangent_linear_lin_state:

Linearisation state
===================

If the linear model is run using idealized initial conditions, then it
is possible to define an analytical linearisation state. Otherwise, the
linearisation state needs to be read and updated from file. The
linearisation state is continually updated over the linear model run
because it is not fixed in time; the linearisation state is the
nonlinear model trajectory.

The linearisation state is read and updated from file in the same way as
the ancillaries, and the lateral-boundary-conditions. This means that
the data is read with XIOS using a time-axis and on every timestep, the
data is updated and interpolated in time.
