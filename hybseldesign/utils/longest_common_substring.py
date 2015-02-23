"""Functions for computing the longest common substring between two
sequences, a problem that is potentially relevant to probe design.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from collections import deque


"""Computes the longest common substring with k mismatches.

The two sequences, a and b, can be numpy arrays or Python strings.
The algorithm runs in O(|a|*|b|) time and requires O(k) space.
It returns, in order, the length of the longest common substring,
the starting position of the substring in a, and the starting
position of the substring in b.

The implementation is of an algorithm published in the 2014 paper
"Longest common substrings with k mismatches"
  by Flouri, Giaquinta, Kobert, and Ukkonen
  (http://arxiv.org/pdf/1409.1694v1.pdf)
It is an implementation of the pseudocode in Figure 1 of the paper.
Note that the pseudocode is meant for k>0; the code below
includes a modification to support k=0.
"""
def k_lcf(a, b, k):
  n = len(a)
  m = len(b)
  ell, r_a, r_b = 0, 0, 0
  for d in xrange(-m+1, n):
    i = max(-d, 0) + d
    j = max(-d, 0)
    Q = deque([])
    s, l = 0, 0
    while l <= min(n-i, m-j) - 1:
      if a[i+l] != b[j+l]:
        if k == 0:
          s = l+1
        else:
          if len(Q) == k:
            s = min(Q) + 1
            Q.popleft()
          Q.append(l)
      l = l+1
      if l-s > ell:
        ell = l-s
        r_a = i+s
        r_b = j+s
  return ell, r_a, r_b

