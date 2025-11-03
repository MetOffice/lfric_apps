from subprocess import Popen, PIPE


def get_site():
    proc = Popen(
        ["rose", "config", "rose-stem", "automatic-options"], stdout=PIPE, text=True
    )
    out, _ = proc.communicate()
    if proc.returncode or "SITE" not in out:
        raise Exception('Could not determine the rose-stem "SITE"')
    return out.replace("SITE=", "").strip()
