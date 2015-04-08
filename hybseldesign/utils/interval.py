"""Functions for working with intervals.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def merge_overlapping(intervals):
  """Merge a list of possibly overlapping intervals.

  Args:
      intervals: list of intervals, each of which is a tuple of the
          form (start, end)

  Returns: 
      list of intervals in which overlapping ones are merged. For
      example, the input [(1,5), (3,7), (9,12)] yields
      [(1,7), (9,12)]. Intervals that are touching (e.g., (1,3) and
      (3,1)) are merged into one.
  """
  if len(intervals) == 0:
    return []

  intervals = sorted(intervals)
  intervals_merged = []
  curr_start, curr_end = intervals[0][0], intervals[0][1]
  for start, end in intervals:
    if start <= curr_end:
      curr_end = max(curr_end, end)
    else:
      intervals_merged += [(curr_start, curr_end)]
      curr_start = start
      curr_end = end
  intervals_merged += [(curr_start, curr_end)]

  return intervals_merged


def schedule(intervals):
  """Schedule the maximum number of compatible intervals.

  This uses to greedy interval scheduling algorithm to find (schedule)
  the maximum number of compatible (non-overlapping) intervals.

  Args:
      intervals: list of intervals, each of which is a tuple (x,y)
      in which x is a tuple of the form (start, end) and y is a
      reference to an object represented by the interval.

  Returns:
      list of the objects corresponding to the chosen intervals
      (i.e., the 'y' for each chosen element)
  """
  # Sort all intervals by their endpoint (the "finishing time")
  # x[0] gives the interval and x[0][1] gives the endpoint
  intervals = sorted(intervals, key=lambda x: x[0][1])

  # Scan through the intervals in sorted order and choose
  # compatible ones with the earliest endpoint
  last_chosen_interval = None
  chosen_objects = []
  for interval, obj in intervals:
    is_compatible = False
    if last_chosen_interval == None:
      # No intervals have been chosen yet, so interval is
      # of course compatible
      is_compatible = True
    else:
      # interval is compatible with the chosen intervals iff
      # its start comes after the finish of the last chosen
      # interval
      if interval[0] >= last_chosen_interval[1]:
        is_compatible = True
    if is_compatible:
      last_chosen_interval = interval
      chosen_objects += [obj]

  return chosen_objects

