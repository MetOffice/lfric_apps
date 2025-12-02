#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

from __future__ import print_function

from __future__ import absolute_import
from abc import ABCMeta
import os
import re
import sys

from testframework import LFRicLoggingTest, TestEngine, TestFailed
import six

class CMATest(six.with_metaclass(ABCMeta, LFRicLoggingTest)):
  def __init__(self,flag):
    self._flag = flag
    if 'MPIEXEC_BROKEN' in os.environ:
      CMATest.set_mpiexec_broken()
    super(CMATest, self).__init__( [sys.argv[1],
                                    'resources/cma_test_configuration.nml',
                                    'test_' + self._flag],
                                    processes=1,
                                    name='cma_test.Log' )


  def test( self, return_code, out, err ):
    if return_code != 0:
      message = 'Test program failed with exit code: {code}'
      raise TestFailed( message.format( code=return_code ),
                        stdout=out, stderr=err,
                        log=self.getLFRicLoggingLog() )

    if not self.test_passed( out ): # "out" becomes self.getLFRicLoggingLog() when PE>1
      message = 'Test {} failed'
      raise TestFailed( message.format( self._flag ),
                        stdout=out, stderr=err,
                        log=self.getLFRicLoggingLog() )

    return 'CMA test : '+self._flag

  def test_passed(self, out):
    success = False
    pattern = re.compile( r'\s+test\s+.*?:\s*PASS\s*$' )
    for line in out.split("\n"):
      match = pattern.search( line )
      if match:
        success = True
    return success

class cma_test_apply_mass_p(CMATest):
    def __init__(self):
        flag = "apply_mass_p"
        super(cma_test_apply_mass_p, self).__init__(flag)

class cma_test_apply_mass_v(CMATest):
    def __init__(self):
        flag = "apply_mass_v"
        super(cma_test_apply_mass_v, self).__init__(flag)

class cma_test_apply_div_v(CMATest):
    def __init__(self):
        flag = "apply_div_v"
        super(cma_test_apply_div_v, self).__init__(flag)

class cma_test_multiply_div_v_mass_v(CMATest):
    def __init__(self):
        flag = "multiply_div_v_mass_v"
        super(cma_test_multiply_div_v_mass_v, self).__init__(flag)

class cma_test_multiply_grad_v_div_v(CMATest):
    def __init__(self):
        flag = "multiply_grad_v_div_v"
        super(cma_test_multiply_grad_v_div_v, self).__init__(flag)

class cma_test_add(CMATest):
    def __init__(self):
        flag = "add"
        super(cma_test_add, self).__init__(flag)

class cma_test_apply_inv(CMATest):
    def __init__(self):
        flag = "apply_inv"
        super(cma_test_apply_inv, self).__init__(flag)

class cma_test_diag_dhmdht(CMATest):
    def __init__(self):
        flag = "diag_dhmdht"
        super(cma_test_diag_dhmdht, self).__init__(flag)


if __name__ == '__main__':
    TestEngine.run( cma_test_apply_mass_v() )
    TestEngine.run( cma_test_apply_mass_p() )
    TestEngine.run( cma_test_apply_div_v()  )
    TestEngine.run( cma_test_multiply_div_v_mass_v() )
    TestEngine.run( cma_test_multiply_grad_v_div_v() )
    TestEngine.run( cma_test_add() )
    TestEngine.run( cma_test_apply_inv() )
    TestEngine.run( cma_test_diag_dhmdht() )
