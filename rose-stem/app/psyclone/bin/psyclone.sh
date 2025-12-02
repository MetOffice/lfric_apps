#!/usr/bin/env bash
##############################################################################
# (c) Crown copyright 2017-2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

set -e

# Path to PSyclone
export PSYCLONE="$(which psyclone)"

# Version of PSyclone API
PSYCLONE_API=lfric

# Declare project array, which indicates the project source PSyclone is supposed to be run on.
declare -a project
project[0]="lfric_apps"
project[1]="lfric_core"
project[2]="scintelapi"


# Base of extracted source directory
BASE_SRC_DIR=$CYLC_SUITE_SHARE_DIR/$BUILD_NAME/extract/

# Psyclone.cfg path
PSYCLONE_CFG_PATH=$BASE_SRC_DIR/lfric_core/etc/psyclone.cfg

# Declare project source directories
declare -a project_src_dir
project_src_dir[0]="${BASE_SRC_DIR}lfric_apps"
project_src_dir[1]="${BASE_SRC_DIR}lfric_core"
project_src_dir[2]="${BASE_SRC_DIR}lfric_apps/applications/lfricinputs/source/scintelapi/generators"

# Declare project kernel directories
declare -a kernel_src_flag
kernel_src_flag[0]="-d ${BASE_SRC_DIR}lfric_apps/science/gungho/source/kernel -d ${BASE_SRC_DIR}lfric_core/components/science/source/kernel"
kernel_src_flag[1]="-d ${BASE_SRC_DIR}lfric_core/components/science/source/kernel"
kernel_src_flag[2]=


# Declare project algorithm directories
declare -a alg_src_dir
alg_src_dir[0]="${BASE_SRC_DIR}lfric_apps/science/gungho/source/algorithm"
alg_src_dir[1]="${BASE_SRC_DIR}lfric_core/infrastructure/source/field ${BASE_SRC_DIR}lfric_core/components/lfric-xios/source ${BASE_SRC_DIR}lfric_core/components/science/source/algorithm"
alg_src_dir[2]="${BASE_SRC_DIR}lfric_apps/applications/lfricinputs/source/scintelapi/generators/toolset ${BASE_SRC_DIR}lfric_apps/applications/lfricinputs/source/scintelapi/generators/analytic"

PRE_PROCESS_MACROS="RDEF_PRECISION=64"

# PSyclone input files are labelled ".x90"; for each algorithm file we find
# which matches that naming convention, PSyclone will generate two output files,
# one containing "psy" code and one containing the transformed algorithm.
# We generate appropriate file names for each, and then invoke PSyclone.
for i in "${!project[@]}"; do

  echo
  echo 'Running PSyclone on '"${project[$i]}"' source'
  echo

  # First preprocess the .x90 and .X90 files to remove ifdef's etc.
  processed_source="${BASE_SRC_DIR}../preprocess-${project[$i]}-x90"
  mkdir -p $processed_source

  DIR_LIST="${alg_src_dir[$i]}"
  for x90file in $(find $DIR_LIST -name '*.[xX]90'); do
    basename=`basename $x90file`
    processed_name="${processed_source}/$basename"
    cpp -traditional-cpp -P -D ${PRE_PROCESS_MACROS} $x90file $processed_name

  done

  # Now setup and invoke PSyclone for each processed source file
  DIR_LIST="${processed_source}"
  for x90file in $(find $DIR_LIST -name '*.[xX]90'); do

    basename=`basename $x90file`
    algname=`echo $basename | sed -e 's/.[xX]90/_alg.f90/g'`
    psyname=`echo $basename | sed -e 's/.[xX]90/_psy.f90/g'`

    if [ -z "${kernel_src_flag[$i]}" ]
    then
      FLAG_KERNEL_DIR=
    else
      FLAG_KERNEL_DIR="${kernel_src_flag[$i]}"
    fi
    PROJ_DIR="${project_src_dir[$i]}"

    echo $PSYCLONE -api $PSYCLONE_API -l all $FLAG_KERNEL_DIR --config $PSYCLONE_CFG_PATH -opsy $PROJ_DIR/$psyname -oalg $PROJ_DIR/$algname $x90file
    $PSYCLONE -api $PSYCLONE_API -l all $FLAG_KERNEL_DIR --config $PSYCLONE_CFG_PATH -opsy $PROJ_DIR/$psyname -oalg $PROJ_DIR/$algname $x90file

  done

done
