from subprocess import Popen, PIPE


def get_site():
    proc = Popen(['rose', 'config', 'rose-stem', 'automatic-options'], stdout=PIPE, text=True)
    out, _ = proc.communicate()
    if proc.returncode or 'SITE' not in out:
        raise Exception('Could not determine the rose-stem "SITE"')
    return out.replace('SITE=', '').strip()


def get_git_ref(working_copy_path):
    proc = Popen(['git', '-C', working_copy_path, 'rev-parse', 'HEAD'], stdout=PIPE, text=True)
    if proc.wait():
        raise Exception('Could not determine the HEAD commit.')
    return proc.communicate()[0].strip()
