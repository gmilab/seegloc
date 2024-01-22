''' Use fsleyes CLI instead of Python package when not available. '''

import os
import os.path
import subprocess


def render(*args):
    ''' Run fsleyes CLI. '''
    fsleyes_path = os.environ.get('FSLDIR', None)
    if fsleyes_path is None:
        raise RuntimeError('FSLDIR environment variable not set.')
    fsleyes_path = os.path.join(fsleyes_path, 'bin', 'fsleyes')

    if len(args) == 1 and isinstance(args[0], list):
        args = args[0]

    cmd = [fsleyes_path, 'render']
    cmd.extend(args)
    subprocess.run(cmd)
