##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

#----------------------------------------------------------------------
# Plots the convergence rate of the tangent linear model.
#----------------------------------------------------------------------

# INSTRUCTIONS TO RUN
# 1. Specify CONFIG
# 2. Run using . plot_convergence.sh, from the plot_convergence directory 

# SCIENCE DETAILS
# The relative linearisation error is
# E =  || N(x+ gamma x') - N(x) - L(x) gamma x' || / || L(x) gamma x' ||
# where N=nonlinear model, L=linear model, x=linearisation state
# x'=perturbation, gamma=scalar.
# From the Taylor series expansion, E(gamma) = O(gamma) i.e. of the order gamma
# So the relative error should be a linear function of gamma

# SCRIPT STEPS
# 1. Produce the data: The integration test tl_test_timesteps is extended by
#    running over 10 values of gamma, rather than 2 values of gamma.
# 2. Plot the data: The data is plotted for each prognostic variable.

# EXTENSION
# The plot_configuration.nml can also be extended to other configurations e.g
# * increase the number of timesteps (timesteps_end)
# * increase the number of timesteps between updating the linearisation state
#   (update_ls_frequency)

#--------------------------------------------------------------------------

# CONFIG can be specified as either runge_kutta or semi_implicit
CONFIG=semi_implicit

# Define directories using the current working directory
Working_dir=$PWD
Linear_dir="$(dirname "$PWD")"

# Integration tests executable name
exe=$Linear_dir/test/$CONFIG

# Build the integration tests, unless that has already been completed
if [ -f $exe ] ; then
  echo "Do not need to build the executable as $exe exists"
else
  echo "$exe does not exist, so now building the executable"
  cd $Linear_dir
  make integration-tests

  if [$? -ne 0 ]; then
    echo "Error building the executable"
    return
  fi
fi

# Setup the configuration - to test with 10 values of gamma
cd $Linear_dir/test/test_files/$CONFIG
cp ${CONFIG}_configuration.nml plot_configuration.nml
sed -i 's/number_gamma_values=2/number_gamma_values=10/g' plot_configuration.nml
if [ $? -ne 0 ]; then
  echo "Error in creating plot_configuration.nml"
  return
fi

# Run the tl_test_timesteps integration test
echo "Running the integration test"
../../$CONFIG plot_configuration.nml test_timesteps > outfile
if [ $? -ne 0 ]; then
  echo "Error in creating outfile data"
  return
else
  echo "Data created successfully"
fi

# Plot the data, together with the expected gradient
echo "Plotting the data"
python $Working_dir/plot_convergence.py


