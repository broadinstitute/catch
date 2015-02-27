"""Functions for working with intervals.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


"""Merges a list of possibly overlapping intervals.

Each interval in 'intervals' is a tuple of the form (start, end).
This returns a list of intervals in which overlapping ones are
merged. For example, the input [(1,5), (3,7), (9,12)] yields
[(1,7), (9,12)].

Intervals that are touching (e.g., (1,3) and (3,1)) are merged into
one.
"""
def merge_overlapping(intervals):
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
