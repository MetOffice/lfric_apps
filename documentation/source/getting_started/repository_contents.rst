.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _repository_contents:

Contents of the Repository
--------------------------

The repository contains the following directories:

- The ``Science`` directory contains several libraries of science code each
  of which can be used to in the development of applications. For
  example, the **gungho-model** library provides the dynamical core used in
  the Momentum Atmosphere application.
- The ``applications`` directory contains the different applications that have
  been developed using the LFRic Core and the science libraries, this includes
  the main Momentum Atmosphere application as well as smaller applications that
  test subsections of different models (e.g. transport, solver, etc.)
- The ``rose-stem`` directory contains the system test suites for development
  testing.
- The ``build`` directory contains additional build scripts required by the
  modelling system.

Many of the directories contain directories of unit and integration
tests, directories of Rose metadata, and Makefiles for building
applications.