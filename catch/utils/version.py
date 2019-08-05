#!/usr/bin/python
"""This determines the catch package version, primarily based on git.

If git is available, we use the version specified by git. This can indicate
commits on top of the last numbered version and can also indicate if the
working directory is dirty (i.e., has local modifications). If git is not
available but some version from git was stored in a file, we use that. Finally,
if none of these are available, we resort to the numbered version manually
specified in the variable VERSION.
"""


import subprocess
import os
import re
import time, datetime

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


def call_git_describe():
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


def approx_version_number():
    """
        In the event that git is unavailable and the VERSION file is not present
        this returns a "version number" in the following precedence:
            - version number from path
                downloads from GitHub tagged releases
                might be extracted into directories containing
                the version number. If they contain a version number
                in the form d.d.d, we can use it
            - modification time of this file (unix timestamp)
                file modification time for github releases corresponds to
                when the release archives were created, a rough way to ballpark
                the release date. If we can't get the version number from the path
                we can at least use the modification time of this file as a proxy
                for the true version number
            - the current time (unix timestamp)
                the current time is better than not having any version number
    """
    version_re = re.compile(r"(?:(\d+)\.)?(?:(\d+)\.)?(?:(\d+))")
    # path relative to version.py
    relative_path = os.path.basename(get_project_path())

    # for tagged releases, the version number might be part of
    # the root directory name
    matches = version_re.search(relative_path)

    if matches and len([n for n in matches.groups() if n]) == 3:
        version = ".".join(map(str, matches.groups()))
    else:
        try:
            # Try to use modification time of the current file
            version = str(int(os.path.getmtime(__file__)))
        except OSError:
            # Just use the current time
            version = str(int(time.time()))

    return version


def get_version():
    """Determine version from git, and save if available.
    """
    global __version__
    if __version__ is None:
        from_git = call_git_describe()
        from_file = read_release_version()

        if from_git:
            if from_file != from_git:
                write_release_version(from_git)
            __version__ = from_git
        else:
            __version__ = from_file

        if __version__ is None:
            __version__ = approx_version_number()

    return __version__


if __name__ == "__main__":
    # Determine and print the package version
    print(get_version())
