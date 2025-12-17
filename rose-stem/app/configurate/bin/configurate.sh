#!/usr/bin/env bash
##############################################################################
# (c) Crown copyright 2017-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
set -e

# Sub-directory containing the "Configurator" tool
TOOLS_DIR=$LFRIC_CORE_DIR/infrastructure/build/tools

# Run rose_picker to generate rose-meta.json and config_namelists.txt
rose_picker $SOURCE_ROOT/lfric_apps/rose-meta/lfric-gungho/HEAD/rose-meta.conf -include_dirs $SOURCE_ROOT/lfric_apps/rose-meta -include_dirs $SOURCE_ROOT/lfric_core/rose-meta

# Run "GenerateNamelist". This takes a rose-meta.json. This produces the *_config_mod.F90 files
echo $TOOLS_DIR/GenerateNamelist rose-meta.json -directory LFRIC_CORE_DIR
$TOOLS_DIR/GenerateNamelist rose-meta.json -directory $LFRIC_CORE_DIR

echo Namelists:
cat config_namelists.txt

# Run "GenerateLoader". This automatically generates the top-level "read_configuration"
# subroutine which processes an input file and reads all namelists LFRic recognises.
# It (rather unhelpfully) aborts if it find an alien namelist.
echo $TOOLS_DIR/GenerateLoader $LFRIC_CORE_DIR/configuration_mod.f90 $(cat config_namelists.txt | xargs)
$TOOLS_DIR/GenerateLoader $LFRIC_CORE_DIR/configuration_mod.f90 $(cat config_namelists.txt | xargs)
