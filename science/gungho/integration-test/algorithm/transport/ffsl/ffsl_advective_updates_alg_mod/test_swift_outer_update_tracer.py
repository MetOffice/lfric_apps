#!/usr/bin/env python3
###############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
from re import compile as re_compile
from sys import argv as sys_argv
from testframework import MpiTest, TestEngine, TestFailed
from typing import Dict

class SwiftOuterUpdateTracerTest(MpiTest):
    def __init__(self):
        super().__init__([sys_argv[1], '64'], processes=6)
        self.__pattern = re_compile(r'tracer \((\d+)\) = (.*)')

    def test(self, return_code: int, out: str, err: str) -> str:
        expected = {
            0: [-0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036],
            1: [-0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036],
            2: [-0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036],
            3: [-0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036],
            4: [-0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036],
            5: [-0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036, -0.036,
                -0.036, -0.036, -0.036]
        }

        if return_code != 0:
            message = f"Test program failed with exit code: {return_code}"
            raise TestFailed(message,
                             stdout=out,
                             stderr=err)

        result: Dict[float] = {}
        for line in out.splitlines():
            match = self.__pattern.match(line)
            if match:
                result[int(match.group(1))] = [float(value)
                                               for value
                                               in match.group(2).split()]
        for rank in range(0, 5):
            if result[rank] != expected[rank]:
                message = ("Tracer field does not match expectation for rank "
                           + str(rank))
                raise TestFailed(message, stdout=out, stderr=err)

        return "Swift outer update tracer"


if __name__ == '__main__':
    TestEngine.run(SwiftOuterUpdateTracerTest())
