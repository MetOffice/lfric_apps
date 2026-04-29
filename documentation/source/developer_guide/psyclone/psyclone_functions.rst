.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _psyclone_functions:

============================================
Available Functions for Optimisation Scripts
============================================

There are a number of transformations provided by Psyclone that can be used for optimisation. These are listed in the Psyclone documentation `here <https://psyclone.readthedocs.io/en/latest/user_guide/transformations.html>`__. To reduce the amount of boiler-plate code, a few useful tranformations have been aggregated and applied in convenient functions specific to LFRic in a ``psyclone_tools.py`` file that resides in LFRic Core. These functions are described below.

---------------------------------
Profiling using ``profile_loops``
---------------------------------

Using the ``profile_loops`` function allows the user to wrap the outermost loop around a coded kernel with a set of timing calipers. This function utilises the ``ProfileTrans`` transformation that allows psyclone to inject a set of timing calipers and employ which ever profiling software is being used in conjunction with LFRic. The function takes a psyir node as an essential argument. As an optional argument, a ``colours_only`` bool can be used to profile only kernels that use coloured ordering. This argument is True by default but when set to False, every coded kernel instance will be profiled. An example that can be included in either the ``global.py`` or a specific algorithm file is below.


.. code-block:: python

   
   from psyclone_tools import profile_loops
   
   def trans(psyir):
       profile_loops(psyir,colours_only=False)


An important thing to note when including this function in an optimisation script is its position relative to other transformation functions.
The ``profile_loops`` function must come after the ``colour_loops`` function but before the ``openmp_parallelise_loops`` function.
This ensures that the callipers are injected into the correct position in the Psyclone generated code.
The build step will fail gracefully if this function is used out of order.
