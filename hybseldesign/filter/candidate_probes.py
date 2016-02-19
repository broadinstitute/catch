"""Functions for generating lists of candidate probes from sequences.

These functions compute lists of many (likely redundant) probes, termed
candidate probes, from a sequence of list of sequences.
"""

import re
import sys

import numpy as np

from hybseldesign import probe
from hybseldesign.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def make_candidate_probes_from_sequence(seq,
                                        probe_length,
                                        probe_stride,
                                        min_n_string_length=2,
                                        insert_bugs=False):
    """Generate a list of candidate probes from a sequence.

    It is possible (especially when there are strings of N's) that
    duplicate probes are returned.

    Args:
        seq: sequence as a string or np.array from which to generate
            candidate probes
        probe_length: generate candidate probes with this number of bp
        probe_stride: generate probes from seq separated by this number
             of bp
        min_n_string_length: possible probes that would contain strings
            of this number or more N's are discarded and, instead, new
            probes flanking the string are added
        insert_bugs: add bugs to the code in order to replicate past
            software (the first version of probe design, as Matlab code)
            and its results

    Returns:
        list of candidate probes as instances of probe.Probe
    """
    if probe_length > len(seq):
        raise ValueError("Invalid probe_length " + str(probe_length))

    if isinstance(seq, np.ndarray):
        seq = ''.join(seq)
    n_string_query = re.compile('(N{' + str(min_n_string_length) + ',})')

    # Make a probe based on the subsequence seq[start:end].
    # Namely, if that subsequence contains no string of N's, then it
    # is a probe to be added and the probe is returned in a single-
    # valued list. Otherwise, an empty list is returned.
    def add_probe_from_subsequence(start, end,
                                   is_flanking_n_string=False,
                                   is_bug_location=False):
        subseq = seq[start:end]
        probes = []

        if insert_bugs and is_bug_location:
            if 'N' in subseq:
                # A bug was present that considered a probe with just a
                # single 'N' (even if min_n_string_length > 1) to be invalid
                # and otherwise to be valid
                # (line 116 of Matlab code)
                pass
            else:
                probes += [subseq]
        else:
            # Search for strings of min_n_string_length or more N's in subseq
            # and only add a probe if there is not such a string
            if not n_string_query.search(subseq):
                # There's no string of N's, so this subsequence is a valid
                # probe
                probes += [subseq]

        # Convert the probes from a Python list of Python strings to a
        # list of probe.Probe
        probes = [probe.Probe.from_str(p) for p in probes]
        for p in probes:
            p.is_flanking_n_string = is_flanking_n_string

        return probes

    # Populate a list of probes
    probes = []
    for start in np.arange(0, len(seq), probe_stride):
        if start + probe_length > len(seq):
            break
        if insert_bugs:
            if start + probe_length == len(seq) and \
                    start % probe_length == 0:
                # A bug was present that stopped creating candidate probes
                # one bp short of where it should have stopped (but only
                # for the so-called 'A adapter'; this is not a problem for
                # the 'B-adapter')
                # (line 86 of Matlab code)
                break

        # A particular bug was present only in designing the
        # so-called 'B adapter'
        # (line 116 of Matlab code)
        is_bug_location = start % probe_length == probe_stride

        probes += add_probe_from_subsequence(start, start + probe_length,
                                             is_bug_location=is_bug_location)
    if len(seq) % probe_stride != 0:
        # There are bases on the right that were never covered, so add
        # another probe for this
        if insert_bugs:
            # A bug was present that ignored the bases on the end
            pass
        else:
            probes += add_probe_from_subsequence(len(seq) - probe_length,
                                                 len(seq))

    # Add probes flanking each string of N's. Specifically, add a probe
    # to the left of a string and to the right. The called function
    # must check that the flanking probe does not contain a string of
    # N's before adding. (Don't recursively chase flanking probes.)
    for match in n_string_query.finditer(seq):
        if match.start() - probe_length >= 0:
            # Add the left flanking probe for match
            probes += add_probe_from_subsequence(match.start() - probe_length,
                                                 match.start(),
                                                 is_flanking_n_string=True)
        if match.end() + probe_length <= len(seq):
            # Add the right flanking probe for match
            probes += add_probe_from_subsequence(match.end(),
                                                 match.end() + probe_length,
                                                 is_flanking_n_string=True)

    return probes


def make_candidate_probes_from_sequences(
        seqs,
        probe_length,
        probe_stride,
        min_n_string_length=2,
        insert_bugs=False,
        move_all_n_string_flanking_probes_to_end=False):
    """Generate a list of candidate probes from a list of sequences.

    It is possible (perhaps even likely depending on where
    the sequences come from) that duplicate probes are returned.

    Args:
        seqs: list of sequences, each as a string or np.array from which
            to generate candidate probes
        probe_length: generate candidate probes with this number of bp
        probe_stride: generate probes from each sequence separated by this
             number of bp
        min_n_string_length: possible probes that would contain strings
            of this number or more N's are discarded and, instead, new
            probes flanking the string are added
        insert_bugs: add bugs to the code in order to replicate past
            software (the first version of probe design, as Matlab code)
            and its results
        move_all_n_string_flanking_probes_to_end: takes all probes that
            were generated because they flank a string of N's, and moves
            these to the end of the returned list; this is necessary to
            replicate past code and results

    Returns:
        list of candidate probes as instances of probe.Probe
    """
    if not isinstance(seqs, list):
        raise ValueError("seqs must be a list of sequences")
    if len(seqs) == 0:
        raise ValueError("seqs must have at least one sequence")
    for seq in seqs:
        if not isinstance(seq, str):
            raise ValueError("seqs must be a list of Python strings")

    probes = []
    for seq in seqs:
        probes += make_candidate_probes_from_sequence(
            seq,
            probe_length=probe_length,
            probe_stride=probe_stride,
            min_n_string_length=min_n_string_length,
            insert_bugs=insert_bugs)

    if move_all_n_string_flanking_probes_to_end:
        # To precisely replicate past software/results, this option
        # moves all probes that flank a string of N's to the end
        # of the returned list.
        probes_non_n = []
        probes_n = []
        for p in probes:
            if p.is_flanking_n_string:
                probes_n += [p]
            else:
                probes_non_n += [p]
        probes = probes_non_n + probes_n

    return probes
