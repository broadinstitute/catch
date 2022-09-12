#!/usr/bin/python
"""This determines the catch package version, primarily based on git.

If git is available, we use the version specified by git. This can indicate
commits on top of the last numbered version and can also indicate if the
working directory is dirty (i.e., has local modifications). If git is not
available but some version from git was stored in a file, we use that. Finally,
if none of these are available, we resort to the numbered release version
manually specified in the variable RELEASE_VERSION.
"""

# Manually specify a numbered release version to use as a fallback
RELEASE_VERSION = 'v1.5.0'


import subprocess
import os

__author__ = ['Danny Park <dpark@broadinstitute.org>',
        'Hayden Metsky <hayden@mit.edu>']

# Set __version__ lazily below
__version__ = None


def get_project_path():
    """Determine absolute path to the top-level of the catch project.

    This is assumed to be the parent of the directory containing this script.

    Returns:
        path (string) to top-level of the catch project
    """
    # abspath converts relative to absolute path; expanduser interprets ~
    path = __file__  # path to this script
    path = os.path.expanduser(path)  # interpret ~
    path = os.path.abspath(path)  # convert to absolute path
    path = os.path.dirname(path)  # containing directory: utils
    path = os.path.dirname(path)  # containing directory: catch project dir
    return path


def get_version_from_git_describe():
    """Determine a version according to git.

    This calls `git describe`, if git is available.

    Returns:
        version from `git describe --tags --always --dirty` if git is
        available; otherwise, None
    """
    cwd = os.getcwd()
    try:
        os.chdir(get_project_path())
        cmd = ['git', 'describe', '--tags', '--always', '--dirty']
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
        if not isinstance(out, str):
            out = out.decode('utf-8')
        ver = out.strip()
    except Exception:
        ver = None
    os.chdir(cwd)
    return ver


def release_file():
    """Obtain path to file storing version, according to git.

    Returns:
        path to VERSION file
    """
    return os.path.join(get_project_path(), 'VERSION')


def read_release_version():
    """Read VERSION file, containing git version.

    Returns:
        if VERSION file exists, version stored in it; otherwise, None
    """
    try:
        with open(release_file(), 'rt') as inf:
            version = inf.readlines()[0].strip()
    except Exception:
        version = None
    return version


def write_release_version(version):
    """Save version, according to git, into VERSION file.

    Args:
        version: version to save
    """
    with open(release_file(), 'wt') as outf:
        outf.write(version + '\n')


def get_version():
    """Determine version from git, and save if available.
    """
    # Allow modifying the global __version__ variable
    global __version__

    if __version__ is None:
        from_git = get_version_from_git_describe()
        from_file = read_release_version()

        if from_git:
            # A version is available from git; use this
            if from_file != from_git:
                # Update the version stored in the VERSION file
                write_release_version(from_git)
            __version__ = from_git
        else:
            if from_file:
                # No version is available from git but one is in the
                # VERSION file; use this
                __version__ = from_file

        if __version__ is None:
            # No version is available from git and there is no VERSION
            # file; use the manually set release version
            __version__ = RELEASE_VERSION

    return __version__


if __name__ == "__main__":
    # Determine and print the package version
    print(get_version())
