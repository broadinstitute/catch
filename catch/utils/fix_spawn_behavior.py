"""Function for fixing a macOS multiprocessing bug.
"""

import logging
import multiprocessing
import os

__author__ = 'Hayden Metsky <hmetsky@broadinstitute.org>'

logger = logging.getLogger(__name__)


already_checked = False
def fix_spawn_behavior():
    global already_checked
    if already_checked:
        return

    # On macOS, starting with Python 3.8, new processes begin following the
    #   spawn behavior rather than fork; apparently, forking processes in
    #   macOS can cause crashes, but CATCH with older versions of
    #   Python has not experienced those issues. Parts of CATCH --
    #   especially the pools in the `probe` module -- are written to follow
    #   the fork behavior, where child processes inherit memory from their
    #   parent. When spawning processes, those child processes do not inherit
    #   memory, and CATCH crashes because the global variables are not
    #   accessible.
    # See https://github.com/python/cpython/issues/84112 for another user
    #   reporting the change in Python 3.8.
    # A quick fix is to simply change the behavior, on macOS, to fork
    #   processes. This is implemented below. A longer term, more appropriate
    #   fix might be to change the global variables accessed by children to
    #   use multiprocessing.shared_memory.SharedMemory objects.
    if os.uname().sysname == 'Darwin':
        logger.debug(("Setting multiprocessing start method to 'fork'"))
        multiprocessing.set_start_method('fork')

    already_checked = True

