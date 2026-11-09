import pytest
import warnings
with warnings.catch_warnings():
  from pyparsing.exceptions import PyparsingDeprecationWarning
  warnings.filterwarnings("ignore", category=PyparsingDeprecationWarning)
  from psyclone.psyir.backend.fortran import FortranWriter
  from psyclone.psyir.frontend.fortran import FortranReader
  from psyclone.psyir.nodes import (
    Loop)
  from psyclone.transformations import OMPParallelLoopTrans


def test_omplooptrans_force_private():
    ''' 
    Test to assist with catching common errors with applying Transmute to
    physics source. Passing 'force_private' and 'ignore_dependencies_for' to
    a transformation is a common pattern and fail point.
    Apply OMPParallelLoopTrans and OMPLoopTrans on a nested loop structure,
    which contains a LHS RHS assignment and a single dimension array on the
    inner loop.
    The LHS RHS assignment, PSyclone will fail safely on assuming it is unsafe
    in a parallel section. 
    The single dimension array, psyclone will assume it's shared currently,
    however as its only for the inner loop, and parallelised around the outer,
    it needs to be private.
    '''
    fread = FortranReader()
    # Example of looping structure with LHS RHS which PSyclone may perceive as
    # a false dependency, and an array which has a single dimension on the,
    # inner loop
    psyir = fread.psyir_from_source('''
        module my_mod
            contains
            subroutine my_subroutine()
                integer :: ji, jj, jk, jpkm1, jpjm1, jpim1
                real, dimension(10, 10, 10) :: array1, array2
                do jk = 2, jpkm1, 1
                  do jj = 2, jpjm1, 1
                    do ji = 2, jpim1, 1
                       array2(ji) = array2(ji) + 1
                       array1(ji) = array2(ji)
                    enddo
                  enddo
                enddo
            end subroutine
        end module my_mod''')
    omplooptrans = OMPParallelLoopTrans()
    loop = psyir.walk(Loop)[0]
    omplooptrans.apply(loop, force_private=['array1','array2'], ignore_dependencies_for=['array1','array2'])
    expected = '''\
    !$omp parallel do default(shared) private(array1,array2,ji,jj,jk) schedule(auto)
    do jk = 2, jpkm1, 1
      do jj = 2, jpjm1, 1
        do ji = 2, jpim1, 1
          array2(ji) = array2(ji) + 1
          array1(ji) = array2(ji)
        enddo
      enddo
    enddo
    !$omp end parallel do\n'''

    fwrite = FortranWriter()
    gen = fwrite(psyir)
    assert expected in gen