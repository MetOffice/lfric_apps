.. -----------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric_atm_checkpoint:

LFRic atmosphere checkpoint/restart system
==========================================

The LFRic atmosphere ``lfric_atm`` application can be configured to generate a
checkpoint dump at the end of each model run. The checkpoint dump can be read in
by a new integration of the model allowing further timesteps to be run. The dump
is written using XIOS.

Requesting checkpoint restart
-----------------------------

Set ``checkpoint_write=.true.`` in the ``io`` namelist of the model
configuration to generate a checkpoint dump at the end of a model run. The
checkpoint dump will be named after the ``checkpoint_stem_name`` string in the
``files`` namelist appended with the number of the last timestep of the run.

Set ``checkpoint_read=.true.`` in the ``io`` namelist to restart a run from an
existing checkpoint dump. The expected start timestep will be defined by
``timestep_start`` in the ``time`` namelist. There must exist a checkpoint file
named according to the ``checkpoint_stem_name``, with the correct timestep
(which would be the ``timestep_start`` setting minus one) otherwise the run will
fail. Typically, such a file will have been created by a preceding run.

Some key ``lfric_atm`` configurations are regularly tested to check that a long
single run produces exactly the same results as a short initial run plus a
restarted run of the same total length. Maintaining such equivalence is
important for many scientific requirements.

As some configurations of ``lfric_atm`` do not run all the same science every
timestep, runs that use checkpoint and restart may need to restart at a
particular frequency. In particular, aligning restart times with radiation
timesteps is important.


The checkpoint dump name and format
-----------------------------------

The name of the checkpoint dump is computed while the XIOS context is being
defined. Its name is computed from a configurable input string (that can define
a full path name plus file name prefix, or just a file name prefix in the
current working directory), and a suffix comprising the timestep number with
leading zeros allowing timesteps up to 9,999,999,999. For example, a checkpoint
written after completion of 4 timesteps may be named ``checkpoint_0000000004``.

Field data can be written to checkpoint files in one of two formats:

* "Legacy" or "stream" format. The field data in this form are written as a 1D
  stream of data, representing the whole 3D field, just as the data is held in
  memory in an LFRic field object.
* "Non-legacy" or "layer" format. The data in this form are written out much
  like it is in diagnostic output: as a series of 2D horizontal layers that are
  stacked to form the full 3D field.

In a single checkpoint file, some fields will be written in stream format, and
some will be written in layers. There are some restrictions to how fields can
be checkpointed:

1. Fields that have data points on both full- and half-levels (so anything on a
   W1 or W2 function space) can only be written directly in stream form. When
   the data is written like this, it cannot be successfully read if the parallel
   decomposition of the reading code is different between the writing and
   reading process. To ensure correct behaviour, such fields could be split into
   the full-level and half-level components which are written and read
   separately.
2. XIOS can write discontinuous fields in layer form only when the convention is
   set to `UGRID`_. Other field outputs can be written in the `CF`_ convention.
3. XIOS does not currently support the layer output of data on the edges of
   cells in discontinuous fields on a bi-periodic mesh.

.. _CF: https://cfconventions.org/
.. _UGRID: https://ugrid-conventions.github.io/ugrid-conventions/

Currently, ``lfric_atm`` outputs at least one W2 field which means
checkpoint-restart suffers the restriction described in the first point
above.

Configuring the XIOS context
----------------------------

This section does not consider LBC prognostics or ancillary fields.

The term "prognostic fields" broadly refers to the set of model fields that must
remain in scope throughout a model timestep. The LFRic atmosphere model stores
the set of prognostic fields in a prognostic field collection held in the main
``modeldb`` data structure. This means that a field can be passed more easily
from one section of the science to another.

Not all prognostic fields need to be written to the checkpoint dump. For some
fields the need to write a field to the dump depends on the configuration. For
other fields, the data stored in a field is passed between science sections
within the same timestep and it does not need to be passed to the next timestep.

* Prognostic fields that are always or sometimes written to the dump must be
  defined in the ``lfric_dictionary.xml`` file which is included in the
  ``iodef.xml`` file that configures XIOS.
* Prognostic fields that are not written to the dump should not be defined in
  the ``lfric_dictionary.xml`` file.

A module called ``create_physics_prognostics`` holds code that uses model
configuration to determine several things about the physics prognostics:

* Identifying which fields need to be written to the checkpoint dump at the end
  of the run.
* Identifying which of many possible prognostic fields are required for a given
  model configuration.
* Initialising those fields that are required.
* Adding initialised fields to the correct science-specific field collection.

Currently, the underpinning infrastructure cannot support doing all of these
tasks at the same time. A solution has been created that involves calling
``create_physics_prognostics`` twice. The procedure takes a procedure pointer as
an argument, enabling the routine to do different tasks on each call.

.. attention::

  There is also a ``create_gungho_prognostics`` file that manages the subset of
  prognostics used by Gungho model configurations. In part, the gungho
  prognostics are handled differently because the ``gungho_model`` application
  supports both lowest and higher order finite element order, and checkpointing
  of the latter requires a different file format. When written with XIOS, all
  gungho prognostics are output using the legacy format even though only one of
  them needs to be on a 3D mesh because it is a W2 field which cannot be split
  into levels.

Within the procedure, there is code like the following for each possible
prognostic field. The first prognostic field ``lw_down_surf_rtsi`` may be
written to the checkpoint file. The second prognostic field ``qcl_at_inv_top``
is never written to the checkpoint file.

.. code-block:: fortran

  call processor%apply(make_spec('lw_down_surf_rtsi', main%radiation,    &
                       ckp=checkpoint_flag))
  call processor%apply(make_spec('qcl_at_inv_top', main%turbulence, W3,  &
                       twod=.true.))

.. attention::

   Note that one of the arguments in the calls above is a function call to
   ``make_spec``. For the field that may be checkpointed the ``make_spec`` call
   does not include the function space type or other information that describes
   the shape of the field, whereas the field that is not checkpointed provides
   the function space type ``W3`` and ``twod=.true`` (the latter indicates that
   the field is on a single layer).

   Any field that `may` be checkpointed `must` have an XIOS definition (as
   pointed out above, these definitions are stored in the
   ``lfric_dictionary.xml`` file). The XML definition is used to infer the shape
   of the field (its function space and so forth). Such informaton is required
   to initialise the field within the model. As the information can be extracted
   from the XML there is no requirement to provide the function space in the
   ``make_spec`` argument list. Indeed, doing so would result in the definition
   appearing in two places and raise the possibility that the two definitions
   may diverge.

Key parts of the above call are:

* ``make_spec`` returns a ``field spec`` data structure that summarises the
  requirements for the field. Optionally, the arguments to ``make_spec`` can
  include stating which field collection should reference the field.

  * If the ``make_spec`` call sets the ``ckp`` flag then it can be assumed that
    XIOS may need to read or write the field and therefore an XML record for the
    field will be in the ``lfric_dictionary.xml`` file that is included in the
    application's ``iodef.xml`` file. This means that the shape of the field
    (its function space, number of levels and so forth) can be inferred from the
    XIOS metadata and does not need to be supplied in the ``make_spec`` argument
    list.
  * If the ``make_spec`` call does not set the ``ckp`` flag, then the field is
    not checkpointed and it must be assumed that there is no XML record for the
    field in the ``iodef.xml`` file. Therefore, additional arguments are
    required to the ``make_spec`` call to define the function space type and
    other information about the field.

* ``apply`` is a procedure in the ``processor`` object that does work on the
  ``field spec``.
* ``processor`` can be one of two different objects depending on whether this is
  the first or second call to the ``create_physics_prognostics`` code.


The first call
^^^^^^^^^^^^^^

The first call to the code in ``create_physics_prognostics`` is done while the
XIOS context is being defined. It must be done during this period because the
aim of the call is to provide a definition of the checkpoint fields to XIOS.

During this call, the ``processor`` object is the ``persistor_type`` defined in
``gungho_model_mod``. Therefore procedure ``gungho_model_mod:persistor_apply``
is applied to the ``field spec`` for each field.

If the field is to be written to or read from the checkpoint dump the ``apply``
method calls ``lfric_xios_metafile_mod:add_field``. The ``add_field`` procedure
registers an ``xios_field``, for writing to or reading from the checkpoint dump.
The ID is the name of the field prefixed with the hard-coded string literal,
``checkpoint_``. The procedure defines the checkpoint field for XIOS using XIOS
attributes copied from the definition of the original field ID (without the
``checkpoint_`` prefix) in ``lfric_dictionary.xml``.

If the field is not to be checkpointed, then nothing is done during the first
call to ``processor%apply``. In fact, the whole first call to the
``create_physics_prognostics`` code has no purpose and does nothing if
checkpoint reading or writing is not requested.

The second call
^^^^^^^^^^^^^^^
In the second call, the ``processor`` passed to ``create_prognostics_field`` is
the ``field_maker_type`` in ``field_maker_mod``.

Briefly, the ``apply`` method in ``field_maker_mod`` initialises the model field
based on the ``field_spec`` definition returned by the ``make_spec``
call. Information in the ``field_spec`` is drawn from the
``lfric_dictionary.xml`` file (if the field has a record in this file) and
information from the arguments to the ``make_spec`` call. For multidata fields
(often referred to as "tiled fields") some hardwired or configured settings are
obtained by the ``apply`` method by calling
``multidata_field_dimensions_mod:get_multidata_field_dimension``.

Problematically, the ``checkpoint_`` prefix is also used here (hard coded as a
string literal) to check if the field is valid before creating anything.

For fields that may be checkpointed, the information from
``lfric_dictionary.xml`` is extracted using the XIOS API. This can only be done
after the context definition has been closed. That is why this call is done
separately from the first call.


Simplified call tree for setting up I/O in LFRic_atm
----------------------------------------------------


.. code-block:: rst

  lfric_atm							(lfric_atm/lfric_atm.f90)
    │
    └──initialise						(gungho/driver/gungho_driver_mod.F90)
        │
        ├──initialise_infrastructure				(gungho/driver/gungho_model_mod.F90)
        │   │
        │   └──init_io						(components/driver/driver_io_mod.F90)
        │       │
        │       └──init_xios_io_context				(components/driver/driver_io_mod.F90)
        │           │
        │           ├──populate_filelist (=> init_gungho_files)	(gungho/driver/gungho_setup_io_mod.F90)
        │           │
        │           └──io_context%initialise_xios_context	(components/lfric_xios/lfric_xios_context_mod.f90)
        │               │
        │               ├──xios_context_initialise		xios
        │               │
        │               ├──xios_get_handle			xios
        │               │
        │               ├──xios_set_current_context		xios
        │               │
        │               ├──init_xios_calendar			(components/lfric_xios/lfric_xios_setup_mod.x90)
        │               │
        │               ├──init_xios_dimensions			(components/lfric_xios/lfric_xios_setup_mod.x90)
        │               │
        │               ├──setup_xios_files			(components/lfric_xios/lfric_xios_setup_mod.x90)
        │               │
        │               ├──before_close (=> before_context_close) 	(gungho/driver/gungho_model_mod.F90)
        │               │   │
        │               │   ├──persistor%init           		(gungho/driver/gungho_model_mod.F90)
        │               │   │
        │               │   ├──process_gungho_prognostics(persistor)
        │               │   │   │					(gungho/driver/create_gungho_prognostics_mod.F90)
        │               │   │   └──persistor%apply(makespec())		(gungho/driver/gungho_model_mod.F90)
        │               │   │       │
        │               │   │       └──add_field			(components/lfric-xios/lfric_xios_metafile_mod.F90)
        │               │   │           │
        │               │   │           ├──(various xios calls…) 	xios
        │               │   │           │
        │               │   │           └──handle_legacy_field 		(components/lfric-xios/lfric_xios_metafile_mod.F90)
        │               │   │
        │               │   └──process_physics_prognostics(persistor)
        │               │       │					(gungho/driver/create_physics_prognostics_mod.F90)
        │               │       └──(…)
        │               │
        │               └──xios_close_context_definitions 	xios
        │
        │
        └──create_model_data					(gungho/driver/gungho_init_fields_mod.X90)
            │
            ├──field_mapper%init				(gungho/driver/field_mapper_mod.F90)          │
            └──create_gungho_prognostics			(gungho/driver/create_gungho_prognostics_mod.F90)
                │
                ├──creator%init					(gungho/driver/field_maker_mod.F90)
                │
                └──process_gungho_prognostics(creator)		(gungho/driver/create_gungho_prognostics_mod.F90)
                    │
                    └──creator%apply(makespec())		(gungho/driver/field_maker_mod.F90)
                        │
                        └──add_real_field			(gungho/driver/field_maker_mod.F90)
                            │
                            ├──field%initialise			(infrastructure/field/field_mod.t90)
                            │
                            ├──field%set_read/write_behaviour	(infrastructure/field/field_mod.t90)
                            │
                            ├──depository%add_field		(infrastructure/field/field_collection_mod.F90)
                            │
                            └──prognostic_fields%add_field	(infrastructure/field/field_collection_mod.F90)

(For brevity, the paths are shortened. Some prefix directories and the “source”
directory have been omitted.)
