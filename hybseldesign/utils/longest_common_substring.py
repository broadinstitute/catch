"""Functions for computing the longest common substring between two sequences.
"""

from collections import deque

import numpy as np

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def k_lcf(a, b, k):
    """Computes the longest common substring with k mismatches.

    The algorithm runs in O(|a|*|b|) time and requires O(k) space.
    The implementation is of an algorithm published in the 2014 paper
    "Longest common substrings with k mismatches"
      by Flouri, Giaquinta, Kobert, and Ukkonen
      (http://arxiv.org/pdf/1409.1694v1.pdf)
    It is an implementation of the pseudocode in Figure 1 of the paper.
    Note that the pseudocode is meant for k>0; the code below
    includes a modification to support k=0.

    Args:
        a: sequence; either numpy array or Python string
        b: sequence; either numpy array or Python string
        k: find the longest common substring with this number of
            mismatches

    Returns:
        a tuple (l, s_a, s_b) where l is the length of the longest
        common substring, s_a is the starting position of the substring
        in a, and s_b is the starting position of the substring in b
    """
    n = len(a)
    m = len(b)
    ell, r_a, r_b = 0, 0, 0
    for d in range(-m + 1, n):
        i = max(-d, 0) + d
        j = max(-d, 0)
        Q = deque([])
        s, l = 0, 0
        while l <= min(n - i, m - j) - 1:
            if a[i + l] != b[j + l]:
                if k == 0:
                    s = l + 1
                else:
                    if len(Q) == k:
                        s = min(Q) + 1
                        Q.popleft()
                    Q.append(l)
            l = l + 1
            if l - s > ell:
                ell = l - s
                r_a = i + s
                r_b = j + s
    return ell, r_a, r_b


def k_lcf_around_anchor(a, b, anchor_start, anchor_end, k):
    """Compute longest common substring around a shared anchor substring.

    The longest common substring found will always contain the anchor
    substring (i.e., a[anchor_start:anchor:end]) as a substring of
    itself. Thus, what is returned may not in fact be the longest common
    substring with k mismatches between a and b; rather, it is the
    longest common substring with k mismatches that contains the given
    anchor.

    The algorithm works by expanding outward from the shared anchor. It
    runs in O(|a|+|b|+k) time and requires O(|a|+|b|) space.

    Args:
        a: sequence; either numpy array or Python string
        b: sequence; either numpy array or Python string
        anchor_start/anchor_end: a and b must be "anchored" around
            the substring between anchor_start and anchor_end (exclusve);
            that is, it must be true that
                a[anchor_start:anchor_end] == b[anchor_start:anchor_end]
        k: find the longest common substring with this number of
            mismatches

    Returns:
        a tuple (l, s) where l is the length of the longest common
        substring found and s is the starting index of the longest
        common substring. Because this substring is found around a
        shared anchor of a and b, the starting index is the same in
        a and b. (For example, when k=0, the returned tuple (l,s)
        indicates that a[s:(s+l)] == b[s:(s+l)].)
    """
    # Make a and b the same length
    if len(a) > len(b):
        # Ignore the end of a
        a = a[:len(b)]
    elif len(b) > len(a):
        # Ignore the end of b
        b = b[:len(a)]

    # Convert a and b to NumPy arrays, which is useful for quickly
    # computing mismatch positions
    if not isinstance(a, np.ndarray):
        a = np.fromiter(a, dtype='U1')
    if not isinstance(b, np.ndarray):
        b = np.fromiter(b, dtype='U1')

    # Check that the anchor is the same in a and b
    if np.any(a[anchor_start:anchor_end] != b[anchor_start:anchor_end]):
        raise ValueError("anchors are different in a and b")

    # Find bools indicating the mismatches between a and b.
    mismatches = a != b
    # Find the indices of the mismatches before the anchor, but
    # reverse them so that the mismatch closest to the anchor
    # has the smallest index and the one further has the largest.
    mismatch_ind_before = np.nonzero(mismatches[:anchor_start][::-1])[0]
    # Find the indices of the mismatches after the anchor
    mismatch_ind_after = np.nonzero(mismatches[anchor_end:])[0]
    # As one moves away from the anchor toward the beginning of the
    # sequence, mismatch_ind_before gives the number of bases
    # that match between a and b at each successive mismatch. For
    # example, mismatch_ind_before[1] gives the number bases between
    # the anchor and the second mismatch encountered moving toward
    # the beginning of the sequence (not including that second
    # mismatched base). Similarly, as one moves from the end of the
    # anchor toward the end of the sequence, mismatch_ind_after gives
    # the number of bases that match between a and b at each successive
    # mismatch. For example, mismatch_ind_after[0] gives the number of
    # bases between the anchor and the first mismatch encountered after
    # the anchor (not including that first mismatched base).

    anchor_len = anchor_end - anchor_start
    max_common_substring_len = -1
    max_common_substring_start = -1
    for i in range(k + 1):
        # Consider the longest common substring that includes the anchor
        # and has i mismatches left of the anchor and k-i mismatches
        # right of the anchor. (This will have k mismatches in total,
        # as desired.) It has anchor_len bases from the anchor,
        # mismatch_ind_before[i] bases from before the anchor, and
        # mismatch_ind_after[k-i] bases from after the anchor.
        if i >= len(mismatch_ind_before):
            # There are no more mismatches to encounter between the anchor
            # and the beginning of the sequence, so all of these bases
            # count.
            before_len = anchor_start
        else:
            before_len = mismatch_ind_before[i]
        if k - i >= len(mismatch_ind_after):
            # There are no more mismatches to encounter between the anchor
            # and the end of the sequence, so all of these bases count.
            after_len = len(a) - anchor_end
        else:
            after_len = mismatch_ind_after[k - i]
        len_of_substring = before_len + anchor_len + after_len
        if len_of_substring > max_common_substring_len:
            max_common_substring_len = len_of_substring
            # Compute where this common substring starts.
            max_common_substring_start = anchor_start - before_len

    return max_common_substring_len, max_common_substring_start
