"""Functions for limiting the execution time of function calls.

This is based on the answers in the StackOverflow post:
  http://stackoverflow.com/questions/366682/how-to-limit-execution-time-of-a-function-call-in-python
"""

from contextlib import contextmanager
import signal

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TimeoutException(Exception): pass


@contextmanager
def time_limit(sec):
    """Limit the execution time of a function.

    Use like:
    `
    try:
        with time_limit(sec):
            # Python code
    except TimeoutException:
        # Python code timed out
    `

    Args:
        sec: number of seconds before raising a timeout

    Raises:
        TimeoutException if the wrapped code executes for more than
        sec seconds
    """
    def signal_handler(signum, frame):
        raise TimeoutException

    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(sec)

    try:
        yield
    finally:
        signal.alarm(0)
