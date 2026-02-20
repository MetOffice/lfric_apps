##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

#----------------------------------------------------------------------
# Plots the convergence rate of the tangent linear model.
#----------------------------------------------------------------------

# INSTRUCTIONS TO RUN
# 1. Specify the config_list
# 2. Run using . plot_convergence.sh, from the plot_convergence directory 

# SCIENCE DETAILS
# The relative linearisation error is
# E =  || N(x+ gamma x') - N(x) - L(x) gamma x' || / || L(x) gamma x' ||
# where N=nonlinear model, L=linear model, x=linearisation state
# x'=perturbation, gamma=scalar.
# From the Taylor series expansion, E(gamma) = O(gamma) i.e. of the order gamma
# So the relative error should be a linear function of gamma

# SCRIPT STEPS
# 1. Build the model
# 2. Produce the data: The integration test tl_test_timesteps is extended by
#    running over 10 values of gamma, rather than 2 values of gamma.
# 3. Plot the data: The data is plotted for each prognostic variable.

# EXTENSION
# The plot_configuration.nml can also be extended to other configurations e.g
# * increase the number of timesteps (timesteps_end)
# * increase the number of timesteps between updating the linearisation state
#   (update_ls_frequency)

#--------------------------------------------------------------------------

######################## Functions ##########################################
build(){
    echo $CONFIG " Building the executable"

    # Integration tests executable name
    exe=$Linear_dir/test/$CONFIG/$CONFIG

    # Build the integration tests, unless that has already been completed
    if [ -f $exe ] ; then
        echo "Do not need to build the executable as $exe exists"
    else
        echo "$exe does not exist, so now building the executable"
        cd $Root_dir/build
        python3 local_build.py linear -t integration-tests

      if [$? -ne 0 ]; then
          echo "Error building the executable"
          return
      fi
    fi
}

integration_test(){
    # Setup the configuration - to test with 10 values of gamma
    echo $CONFIG " Setting up the configuration"
    cd $Linear_dir/integration-test/$CONFIG/resources/
    cp ${CONFIG}_configuration.nml plot_configuration.nml
    sed -i 's/number_gamma_values=2/number_gamma_values=10/g' plot_configuration.nml
    if [ $? -ne 0 ]; then
        echo "Error in creating plot_configuration.nml"
        return
    else
        echo "plot_configuration.nml created in " $PWD
    fi

    # Run the tl_test_timesteps integration test
    echo $CONFIG " Running the integration test"
    cd $Linear_dir/integration-test/$CONFIG
    echo $PWD
    $exe resources/plot_configuration.nml test_timesteps > outfile
    if [ $? -ne 0 ]; then
        echo "Error in creating outfile data"
        return
    else
        echo "Data created successfully"
    fi
}

plot(){
    # Plot the data, together with the expected gradient
    echo $CONFIG " Plotting the data"
    cd $Linear_dir/integration-test/$CONFIG
    python $Linear_dir/plot_convergence/plot_convergence.py
}

#################### MAIN PROGRAM #################################################

# Modules
module purge
module use /home/users/lfricadmin/lmod
module load lfric/vn3.0
module load scitools

config_list=(nwp_gal9 semi_implicit runge_kutta)

# Directory of this script
SCRIPT_DIR="$(dirname "$(realpath "$BASH_SOURCE")")"
# And parent directories
export Linear_dir="$(dirname $SCRIPT_DIR)"
export Parent_dir="$(dirname $Linear_dir)"
export Root_dir="$(dirname $Parent_dir)"
echo "Linear_dir" $Linear_dir
echo "Root_dir" $Root_dir

for configuration in "${config_list[@]}"; do
    echo $configuration
    export CONFIG="$configuration"

    build
    integration_test
    plot
done
