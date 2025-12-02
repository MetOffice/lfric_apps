#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Run the jedi_lfric integration tests
'''

import os
import re
import sys

from testframework import LFRicLoggingTest, TestEngine, TestFailed


class JediLfricTest(LFRicLoggingTest):
    '''
    Run the jedi_lfric integration tests
    '''

    def __init__(self, flag: str) -> None:
        self._flag = flag
        if 'MPIEXEC_BROKEN' in os.environ:
            JediLfricTest.set_mpiexec_broken()
        super().__init__([sys.argv[1],
                          'resources/algorithm_test_configuration.nml',
                          'test_' + self._flag],
                         processes=1,
                         name='JediLfricTest.Log')

    def test(self, return_code: int, out: str, err: str) -> str:
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

        return 'da dev test : '+self._flag

    @staticmethod
    def test_passed(out: str) -> bool:
        '''
        Examine the output to see if the validity test passed
        '''
        success = False
        pattern = re.compile(r'test\s*PASS\s*$')
        for line in out.split("\n"):
            match = pattern.search(line)
            if match:
                success = True
        return success


class jedi_lfric_increment_alg(JediLfricTest):
    '''
    Test running the jedi_lfric_increment_alg_mod.x90 algorithm
    '''

    def __init__(self):
        flag = "jedi_lfric_increment_alg_mod"
        super().__init__(flag)


if __name__ == '__main__':
    TestEngine.run(jedi_lfric_increment_alg())
