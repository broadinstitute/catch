"""Functions for generating lists of candidate probes from sequences.

These functions compute lists of many (likely redundant) probes, termed
candidate probes, from a sequence of list of sequences.
"""

import logging
import re
import sys

import numpy as np

from catch import probe
from catch.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def make_candidate_probes_from_sequence(seq,
                                        probe_length,
                                        probe_stride,
                                        min_n_string_length=2,
                                        allow_small_seqs=None):
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
        allow_small_seqs: if set, allow sequences that are smaller than the
            probe length by creating candidate probes equal to the sequence;
            the value gives the minimum allowed probe (sequence) length

    Returns:
        list of candidate probes as instances of probe.Probe
    """
    n_string_query = re.compile('(N{' + str(min_n_string_length) + ',})')

    if len(seq) < probe_length:
        if allow_small_seqs:
            if len(seq) < allow_small_seqs:
                raise ValueError(("Allowing sequences smaller than the probe "
                                  "length (" + str(probe_length) + "), but "
                                  "input sequence is smaller than minimum "
                                  "allowed length"))
            else:
                if n_string_query.search(seq):
                    raise Exception(("Only possible probe from input "
                                     "sequence has too long a stretch of N's"))
                else:
                    # Make a probe equal to this sequence
                    return [probe.Probe.from_str(seq)]
        else:
            raise ValueError(("An input sequence is smaller than the probe "
                              "length (" + str(probe_length) + ")"))

    if isinstance(seq, np.ndarray):
        seq = ''.join(seq)

    # Make a probe based on the subsequence seq[start:end].
    # Namely, if that subsequence contains no string of N's, then it
    # is a probe to be added and the probe is returned in a single-
    # valued list. Otherwise, an empty list is returned.
    def add_probe_from_subsequence(start, end,
                                   is_flanking_n_string=False):
        subseq = seq[start:end]
        probes = []

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
        probes += add_probe_from_subsequence(start, start + probe_length)

    if len(seq) % probe_stride != 0:
        # There are bases on the right that were never covered, so add
        # another probe for this
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
        allow_small_seqs=None,
        seq_length_to_skip=None):
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
        allow_small_seqs: if set, allow sequences that are smaller than the
            probe length by creating candidate probes equal to the sequence;
            the value gives the minimum allowed probe (sequence) length
        seq_length_to_skip: if set, skip sequences whose length is <=
            the given value (i.e., do not design candidate probes for
            them)

    Returns:
        list of candidate probes as instances of probe.Probe
    """
    if not isinstance(seqs, list):
        raise TypeError("seqs must be a list of sequences")
    if len(seqs) == 0:
        raise ValueError("seqs must have at least one sequence")
    for seq in seqs:
        if not isinstance(seq, str):
            raise TypeError("seqs must be a list of Python strings")

    probes = []
    for seq in seqs:
        if seq_length_to_skip is not None:
            if len(seq) <= seq_length_to_skip:
                logger.info(("Not designing candidate probes for a "
                    "sequence with length %d, since it is <= %d"),
                    len(seq), seq_length_to_skip)
                continue

        probes += make_candidate_probes_from_sequence(
            seq,
            probe_length=probe_length,
            probe_stride=probe_stride,
            min_n_string_length=min_n_string_length,
            allow_small_seqs=allow_small_seqs)

    return probes
