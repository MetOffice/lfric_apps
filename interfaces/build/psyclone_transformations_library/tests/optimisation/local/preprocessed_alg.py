"""
PSyclone transformation script.
"""

from psyclone.psyGen import InvokeSchedule


def trans(psyir):
    """
    Duplicates the body of the first loop leaving two copies.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    # Loop over all of the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        print(f"Transforming invoke '{subroutine.name}' ...")

        loop = subroutine.loops()[0]
        loop_schedule = loop.loop_body
        new_node = loop_schedule[0].copy()
        loop_schedule.addchild(new_node)

        # Take a look at what we've done
        print(subroutine.view())

    return psyir
