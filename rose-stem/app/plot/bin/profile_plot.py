#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Plots profiles of horizontally averaged fields.
'''

import sys
import iris

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

if iris.__version__ < "3.0.0":
    iris.FUTURE.netcdf_promote = True

#------------------------------------------------------------------------------#
# Use a variables name to obtain its data from an iris cube
#------------------------------------------------------------------------------#
def load_cube_by_varname(filename, var):

   variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var))
   return iris.load_cube(filename, constraint=variable_constraint)

#------------------------------------------------------------------------------#
# Set heights of levels for selected extrusion
#------------------------------------------------------------------------------#
def make_extrusion(extrusion, number_of_layers, domain_top):

    if extrusion == 'linear':
        z1d = np.linspace(0.0, domain_top, nlayers+1)
    elif extrusion == 'quadratic':
        z1d = domain_top * np.array([(float(k) / float(number_of_layers)) ** 2
                                      for k in range(number_of_layers+1)])
    elif extrusion == 'geometric':
        stretching = 1.03
        z1d = np.zeros([number_of_layers+1])
        dz = 1.0
        z1d[0] = 0.0
        for k in range(number_of_layers):
            z1d[k+1] = z1d[k] + dz
            dz = dz * 1.03
        z1d = z1d * domain_top / z1d[number_of_layers]
    elif extrusion == 'umL38':
        if number_of_layers != 38:
            raise ValueError('For UM extrusion need 38 layers not %d' % number_of_layers)

        z1d = domain_top * np.array([.0000000, .0005095, .0020380,
                 .0045854, .0081519, .0127373, .0183417, .0249651, .0326074,
                 .0412688, .0509491, .0616485, .0733668, .0861040, .0998603,
                 .1146356, .1304298, .1472430, .1650752, .1839264, .2037966,
                 .2246857, .2465938, .2695209, .2934670, .3184321, .3444162,
                 .3714396, .3998142, .4298913, .4620737, .4968308, .5347160,
                 .5763897, .6230643, .6772068, .7443435, .8383348, 1.000000])

        z1dh = domain_top * np.array([.0002547,  .0012737,  .0033117,  .0063686,  .0104446,
                                      .0155395,  .0216534,  .0287863,  .0369381,  .0461090,
                                      .0562988,  .0675076,  .0797354,  .0929822,  .1072479,
                                      .1225327,  .1388364,  .1561591,  .1745008,  .1938615,
                                      .2142411,  .2356398,  .2580574,  .2814940,  .3059496,
                                      .3314242,  .3579279,  .3856269,  .4148527,  .4459825,
                                      .4794523,  .5157734,  .5555529,  .5997270,  .6501355,
                                      .7107751,  .7913392,  .9191674])

    elif extrusion == 'umL70':
        if number_of_layers != 70:
            raise ValueError('For UM extrusion need 70 layers not %d' % number_of_layers)

        z1d = domain_top * np.array([
   0.0000000E+00,   0.1250000E-03,   0.5416666E-03,   0.1125000E-02,   0.1875000E-02,
   0.2791667E-02,   0.3875000E-02,   0.5125000E-02,   0.6541667E-02,   0.8125000E-02,
   0.9875000E-02,   0.1179167E-01,   0.1387500E-01,   0.1612500E-01,   0.1854167E-01,
   0.2112500E-01,   0.2387500E-01,   0.2679167E-01,   0.2987500E-01,   0.3312500E-01,
   0.3654167E-01,   0.4012500E-01,   0.4387500E-01,   0.4779167E-01,   0.5187500E-01,
   0.5612501E-01,   0.6054167E-01,   0.6512500E-01,   0.6987500E-01,   0.7479167E-01,
   0.7987500E-01,   0.8512500E-01,   0.9054167E-01,   0.9612500E-01,   0.1018750E+00,
   0.1077917E+00,   0.1138750E+00,   0.1201250E+00,   0.1265417E+00,   0.1331250E+00,
   0.1398750E+00,   0.1467917E+00,   0.1538752E+00,   0.1611287E+00,   0.1685623E+00,
   0.1761954E+00,   0.1840590E+00,   0.1921980E+00,   0.2006732E+00,   0.2095645E+00,
   0.2189729E+00,   0.2290236E+00,   0.2398690E+00,   0.2516917E+00,   0.2647077E+00,
   0.2791699E+00,   0.2953717E+00,   0.3136506E+00,   0.3343919E+00,   0.3580330E+00,
   0.3850676E+00,   0.4160496E+00,   0.4515977E+00,   0.4924007E+00,   0.5392213E+00,
   0.5929016E+00,   0.6543679E+00,   0.7246365E+00,   0.8048183E+00,   0.8961251E+00,
   0.1000000E+01])

        z1dh = domain_top * np.array([
   0.6249999E-04,   0.3333333E-03,   0.8333333E-03,   0.1500000E-02,   0.2333333E-02,
   0.3333333E-02,   0.4500000E-02,   0.5833333E-02,   0.7333333E-02,   0.9000000E-02,
   0.1083333E-01,   0.1283333E-01,   0.1500000E-01,   0.1733333E-01,   0.1983333E-01,
   0.2250000E-01,   0.2533333E-01,   0.2833333E-01,   0.3150000E-01,   0.3483333E-01,
   0.3833333E-01,   0.4200000E-01,   0.4583333E-01,   0.4983333E-01,   0.5400000E-01,
   0.5833334E-01,   0.6283334E-01,   0.6750000E-01,   0.7233334E-01,   0.7733333E-01,
   0.8250000E-01,   0.8783333E-01,   0.9333333E-01,   0.9900000E-01,   0.1048333E+00,
   0.1108333E+00,   0.1170000E+00,   0.1233333E+00,   0.1298333E+00,   0.1365000E+00,
   0.1433333E+00,   0.1503334E+00,   0.1575020E+00,   0.1648455E+00,   0.1723789E+00,
   0.1801272E+00,   0.1881285E+00,   0.1964356E+00,   0.2051189E+00,   0.2142687E+00,
   0.2239982E+00,   0.2344463E+00,   0.2457803E+00,   0.2581997E+00,   0.2719388E+00,
   0.2872708E+00,   0.3045112E+00,   0.3240212E+00,   0.3462124E+00,   0.3715503E+00,
   0.4005586E+00,   0.4338236E+00,   0.4719992E+00,   0.5158110E+00,   0.5660614E+00,
   0.6236348E+00,   0.6895022E+00,   0.7647274E+00,   0.8504717E+00,   0.9480625E+00])

    else:
        raise NotImplementedError('Extrusion %s is not implemented' % extrusion)

    z1dh = 0.5 * ( z1d[:-1] + z1d[1:] )

    return (z1d, z1dh)

#################################################################################

def make_figures(file_name, field_list, extrusion, number_of_layers, domain_top, plot_path):

    z1d, z1dh = make_extrusion(extrusion, number_of_layers, domain_top)

    for field_name in field_list:

        field_cube = load_cube_by_varname(file_name, field_name)

        field_shape = field_cube.shape

        if len(field_shape) == 3:
           (ntims, nlevs, npts) = field_shape
        else:
           (nlevs, npts) = field_shape
           ntims = 1

# Determine whether field uses full or half levels
        if nlevs == len(z1d):
            z = z1d
        elif nlevs == len(z1dh):
            z = z1dh
        else:
            raise Exception('Number of levels in ' + field_name +
                            'not same as for specified grid.')

        for itim in range(ntims):
            if len(field_shape) == 3:
                field = field_cube.data[itim, :, :]
            else:
                field = field_cube.data[:,:]

            field_prof = np.zeros(nlevs)
            for k in range(nlevs):
                field_prof[k] = np.sum(field[k, :]) / float(npts)

            plt.plot(field_prof, z)

        plt.xlabel('Domain mean '+field_name, style='italic')
        plt.ylabel('height (km)', style='italic')
        plt.savefig(plot_path+'/'+field_name+'_profile.png')
        plt.close()

if __name__ == "__main__":

    try:
        file_name, fields, extrusion, number_of_layers, domain_top, plot_path = sys.argv[1:7]
    except ValueError:
        print(
            "Usage: {0} <file_name> <field1:field2:...> <extrusion> <number_of_layers>"
            "<domain_top(km)> <plot_path>".format(sys.argv[0]))
        exit(1)

    field_list = fields.split(':')

    make_figures(file_name, field_list, extrusion, int(number_of_layers), float(domain_top), plot_path)
