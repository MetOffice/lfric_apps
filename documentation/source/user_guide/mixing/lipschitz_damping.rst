.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _lipschitz_damping_user_guide:

Lipschitz Damping
=================

Lipschitz damping targets wind values that are locally unstable, as measured
by their Lipschitz number (see the :ref:`science guide <lipschitz_damping_doc>`
for details of the scheme). It is controlled by two independent namelist
options, both called ``lipschitz_damping``:

``namelist:mixing=lipschitz_damping``
   Applies the scheme to the wind increment from fast physics, before it is
   added to the dynamics wind field. This is applied every timestep.

``namelist:initialization=lipschitz_damping``
   Applies the scheme once, to the initial wind field, when the model is
   initialised from a start dump (not from a checkpoint file).

Both options default to ``.false.`` and can be enabled independently of one
another.
