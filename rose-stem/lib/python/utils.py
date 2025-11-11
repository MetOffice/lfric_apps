##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

from subprocess import Popen, PIPE


def get_site():
    proc = Popen(
        ["rose", "config", "rose-stem", "automatic-options"], stdout=PIPE, text=True
    )
    out, _ = proc.communicate()
    if proc.returncode or "SITE" not in out:
        raise Exception('Could not determine the rose-stem "SITE"')
    # At some sites there may be many variables that are returned by rose config rose-stem
    # Try to just grab the thing after SITE= and then ignore anything afterwards
    site = out.split("SITE=")[1].split(' ')[0].strip()
    return site
