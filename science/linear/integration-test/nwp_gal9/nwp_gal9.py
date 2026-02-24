#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Run the linear model integration tests for the default (nwp_gal9) configuration

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
                                      'resources/nwp_gal9_configuration.nml',
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

class tl_test_semi_imp_alg(TLTest):
    '''
    Test the semi implicit timestep
    '''
    def __init__(self):
        flag = "semi_imp_alg"
        super(tl_test_semi_imp_alg, self).__init__(flag)


class tl_test_rhs_alg(TLTest):
    '''
    Test the right hand side forcing for the mixed solver
    '''
    def __init__(self):
        flag = "rhs_alg"
        super(tl_test_rhs_alg, self).__init__(flag)

class tl_test_transport_control(TLTest):
    '''
    Test the transport
    '''
    def __init__(self):
        flag = "transport_control"
        super(tl_test_transport_control, self).__init__(flag)


class tl_test_timesteps(TLTest):
    '''
    Test running over multiple timesteps
    '''
    def __init__(self):
        flag = "timesteps"
        super(tl_test_timesteps, self).__init__(flag)

if __name__ == '__main__':
    TestEngine.run(tl_test_rhs_alg())
    TestEngine.run(tl_test_transport_control())
    TestEngine.run(tl_test_semi_imp_alg())
    TestEngine.run(tl_test_timesteps())
 
