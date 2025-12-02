#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Run the jedi-lfric integration tests for datetime, duration, and clock code
'''

import datetime
import os
import re
import sys

from testframework import LFRicLoggingTest, MpiTest, TestEngine, TestFailed


class LFRicTimeTest(LFRicLoggingTest):
    '''
    Run the integration tests
    '''

    def __init__(self, flag: str, optional_arg='None') -> None:
        self._flag = flag
        self._optional_arg = optional_arg
        if 'MPIEXEC_BROKEN' in os.environ:
            LFRicTimeTest.set_mpiexec_broken()
        super().__init__([sys.argv[1],
                          'resources/jedi_lfric_time_test_configuration.nml',
                          'test_' + self._flag,
                          self._optional_arg],
                         processes=1,
                         name='TimeTest.Log')

    def test(self, return_code: int, out: str, err: str) -> str:
        '''
        Error messages if the test failed to run
        '''
        if return_code != 0:
            message = (
                '[optional arg: {}] Test program failed with exit code: {}'
            )
            raise TestFailed(message.format(self._optional_arg, return_code),
                             stdout=out, stderr=err,
                             log=self.getLFRicLoggingLog())

        # "out" becomes self.getLFRicLoggingLog() when PE>1
        if not self.test_passed(out):
            message = '[optional arg: {}] Test {} failed'
            raise TestFailed(message.format(self._optional_arg, self._flag),
                             stdout=out, stderr=err,
                             log=self.getLFRicLoggingLog())

        return 'jedi-lfric test : '+self._flag

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


###############################################################################
class ErrorSerialTest(MpiTest):  # pylint: disable=too-few-public-methods
    """
    Tests logging in serial scenarios.
    """
    def __init__(self, flag: str, optional_arg='None') -> None:
        self._flag = flag
        self._optional_arg = optional_arg

        super().__init__([sys.argv[1],
                          'jedi_lfric_time_test_configuration.nml',
                          'test_' + self._flag,
                          self._optional_arg],
                         processes=1)

        self._minimum_timestamp = datetime.datetime.utcnow()

    def test(self, returncode: int, out: str, err: str):
        """
        Tests that logging an error ends execution.
        """

        if returncode == 0:  # pylint: disable=no-else-raise
            message = (
                "[optional arg: {}] Logging an error did not cause "
                "termination to end"
            )
            raise TestFailed(message.format(self._optional_arg))
        elif returncode == 127:
            raise TestFailed('Test executable not found')
        elif returncode > 128:
            message = (
                "[optional arg: {}] Execution fault such as segmentation fault"
            )
            raise TestFailed(message.format(self._optional_arg))

        if not self.test_passed(err):
            message = (
                '[optional arg: {}] Test {} failed with output: {} and err: {}'
            )
            raise TestFailed(
                message.format(self._optional_arg, self._flag, out, err)
            )

        message = 'Logging an error caused exit as expected with code {code}'
        return message.format(code=returncode)

    @staticmethod
    def test_passed(err: str) -> bool:
        '''
        Examine the output to see if the validity test passed
        '''
        success = False
        pattern = re.compile(r'ERROR:')
        for line in err.split("\n"):
            match = pattern.search(line)
            if match:
                success = True
        return success


class init_string_err(ErrorSerialTest):
    '''
    Tests logging an error when initialising a datetime
    with a bad string
    '''

    def __init__(self):
        flag = "init_string_err"
        super().__init__(flag)


class copy_from_jedi_datetime_err(ErrorSerialTest):
    '''
    Tests logging an error when attempting to initialise
    a datetime with another uninitialised datetime
    '''

    def __init__(self):
        flag = "copy_from_jedi_datetime_err"
        super().__init__(flag)


class add_duration_to_datetime(LFRicTimeTest):
    '''
    Test adding a jedi duration instance to a
    jedi datetime instance
    '''

    def __init__(self):
        flag = "add_duration_to_datetime"
        super().__init__(flag)


class duration_from_datetimes(LFRicTimeTest):
    '''
    Test getting a jedi duration by subtracting
    one datetime from another
    '''

    def __init__(self):
        flag = "duration_from_datetimes"
        super().__init__(flag)


class YYYYMMDD_to_JDN(ErrorSerialTest):
    '''
    Test logging an error in the YYYYMMDD_to_JDN
    subroutine of jedi_datetime_functions_mod
    '''

    def __init__(self):
        flag = "YYYYMMDD_to_JDN"
        super().__init__(flag)


class JDN_to_YYYYMMDD_invalid(ErrorSerialTest):
    '''
    Test logging an error in the JDN_to_YYYYMMDD
    subroutine of jedi_datetime_functions_mod
    with an invalid JDN value
    '''

    def __init__(self):
        flag = "JDN_to_YYYYMMDD_invalid"
        super().__init__(flag)


class hhmmss_to_seconds(ErrorSerialTest):
    '''
    Test logging an error in the hhmmss_to_seconds
    subroutine of jedi_datetime_functions_mod
    '''

    def __init__(self):
        flag = "hhmmss_to_seconds"
        super().__init__(flag)


class seconds_to_hhmmss_large(ErrorSerialTest):
    '''
    Test logging an error in the seconds_to_hhmmss
    subroutine of jedi_datetime_functions_mod if
    the time is over the number of seconds in a day
    '''

    def __init__(self):
        flag = "seconds_to_hhmmss_large"
        super().__init__(flag)


class seconds_to_hhmmss_neg(ErrorSerialTest):
    '''
    Test logging an error in the seconds_to_hhmmss
    subroutine of jedi_datetime_functions_mod if
    the time is negative
    '''

    def __init__(self):
        flag = "seconds_to_hhmmss_neg"
        super().__init__(flag)


class duration_init_bad_string_err(ErrorSerialTest):
    '''
    Test logging an error initialising a jedi
    duration instance with a bad iso string
    '''

    def __init__(self, bad_string: str) -> None:
        flag = "duration_init_bad_string_err"
        super().__init__(flag, bad_string)


class duration_divide_zero_err(ErrorSerialTest):
    '''
    Test logging an error when dividing
    a duration instance by zero
    '''

    def __init__(self):
        flag = "duration_divide_zero_err"
        super().__init__(flag)


class duration_divide_remainder_err(ErrorSerialTest):
    '''
    Test logging an error when two durations are
    not evenly divisible
    '''

    def __init__(self):
        flag = "duration_divide_remainder_err"
        super().__init__(flag)


class duration_divide_int_zero_err(ErrorSerialTest):
    '''
    Test logging an error when dividing
    a duration instance by an integer zero
    '''

    def __init__(self):
        flag = "duration_divide_int_zero_err"
        super().__init__(flag)


class duration_divide_int_remainder_err(ErrorSerialTest):
    '''
    Test logging an error when a durations and an
    integer are not evenly divisible
    '''

    def __init__(self):
        flag = "duration_divide_int_remainder_err"
        super().__init__(flag)


if __name__ == '__main__':
    TestEngine.run(init_string_err())
    TestEngine.run(copy_from_jedi_datetime_err())
    TestEngine.run(add_duration_to_datetime())
    TestEngine.run(duration_from_datetimes())
    TestEngine.run(YYYYMMDD_to_JDN())
    TestEngine.run(JDN_to_YYYYMMDD_invalid())
    TestEngine.run(hhmmss_to_seconds())
    TestEngine.run(seconds_to_hhmmss_large())
    TestEngine.run(seconds_to_hhmmss_neg())
    TestEngine.run(duration_init_bad_string_err("P@D"))
    TestEngine.run(duration_init_bad_string_err("PT300"))
    TestEngine.run(duration_init_bad_string_err("T0S"))
    TestEngine.run(duration_init_bad_string_err("P"))
    TestEngine.run(duration_init_bad_string_err("-P"))
    TestEngine.run(duration_divide_zero_err())
    TestEngine.run(duration_divide_remainder_err())
    TestEngine.run(duration_divide_int_zero_err())
    TestEngine.run(duration_divide_int_remainder_err())
