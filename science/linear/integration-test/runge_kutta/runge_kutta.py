#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Run the linear model integration tests for the runge-kutta configuration
'''

import os
import re
import sys


from testframework import LFRicLoggingTest, TestEngine, TestFailed


class TLTest(LFRicLoggingTest):
    '''
    Run the linear model integration tests
    '''

    def __init__(self, flag):
        self._flag = flag
        if 'MPIEXEC_BROKEN' in os.environ:
            TLTest.set_mpiexec_broken()
        super(TLTest, self).__init__([sys.argv[1],
                                      'resources/runge_kutta_configuration.nml',
                                      'test_' + self._flag],
                                     processes=1,
                                     name='tl_test.Log')

    def test(self, return_code, out, err):
        '''
        Error messages if the test failed to run
        '''
        if return_code != 0:
            message = 'Test program failed with exit code: {code}'
            raise TestFailed(message.format(code=return_code),
                             stdout=out, stderr=err,
                             log=self.getLFRicLoggingLog())

        # "out" becomes self.getLFRicLoggingLog() when PE>1
        if not self.test_passed(out):
            message = 'Test {} failed'
            raise TestFailed(message.format(self._flag),
                             stdout=out, stderr=err,
                             log=self.getLFRicLoggingLog())

        return 'TL test : '+self._flag

    def test_passed(self, out):
        '''
        Examine the output to see if the validity test passed
        '''
        success = False
        pattern = re.compile(r'\s+test\s+.*?:\s*PASS\s*$')
        for line in out.split("\n"):
            match = pattern.search(line)
            if match:
                success = True
        return success


class tl_test_kinetic_energy_gradient(TLTest):
    '''
    Test the kinetic energy gradient kernel
    '''
    def __init__(self):
        flag = "kinetic_energy_gradient"
        super(tl_test_kinetic_energy_gradient, self).__init__(flag)


class tl_test_rk_alg(TLTest):
    '''
    Test the runge-kutta timestepping
    '''
    def __init__(self):
        flag = "rk_alg"
        super(tl_test_rk_alg, self).__init__(flag)

class tl_test_project_eos_pressure(TLTest):
    def __init__(self):
        flag = "project_eos_pressure"
        super(tl_test_project_eos_pressure, self).__init__(flag)

class tl_test_advect_density_field(TLTest):
    '''
    Test density advection
    '''
    def __init__(self):
        flag = "advect_density_field"
        super(tl_test_advect_density_field, self).__init__(flag)


class tl_test_advect_theta_field(TLTest):
    '''
    Test theta advection
    '''
    def __init__(self):
        flag = "advect_theta_field"
        super(tl_test_advect_theta_field, self).__init__(flag)


class tl_test_vorticity_advection(TLTest):
    '''
    Test the vorticity advection kernel
    '''
    def __init__(self):
        flag = "vorticity_advection"
        super(tl_test_vorticity_advection, self).__init__(flag)


class tl_test_pressure_gradient_bd(TLTest):
    '''
    Test the pressure_gradient bd kernel
    '''
    def __init__(self):
        flag = "pressure_gradient_bd"
        super(tl_test_pressure_gradient_bd, self).__init__(flag)


class tl_test_hydrostatic(TLTest):
    '''
    Test the hydrostatic kernel
    '''
    def __init__(self):
        flag = "hydrostatic"
        super(tl_test_hydrostatic, self).__init__(flag)


class tl_test_timesteps(TLTest):
    '''
    Test running over multiple timesteps
    '''
    def __init__(self):
        flag = "timesteps"
        super(tl_test_timesteps, self).__init__(flag)


if __name__ == '__main__':
    TestEngine.run( tl_test_kinetic_energy_gradient() )
    TestEngine.run( tl_test_advect_density_field() )
    TestEngine.run( tl_test_advect_theta_field() )
    TestEngine.run( tl_test_vorticity_advection() )
    TestEngine.run( tl_test_hydrostatic() )
    TestEngine.run( tl_test_pressure_gradient_bd() )
    TestEngine.run( tl_test_project_eos_pressure() )
    TestEngine.run( tl_test_rk_alg() )
