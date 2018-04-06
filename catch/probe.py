"""Structure(s) and functions for directly working with probes.
"""

import atexit
import bisect
import ctypes
from collections import defaultdict
from functools import partial
import gc
import hashlib
import logging
import multiprocessing
from multiprocessing import sharedctypes

import numpy as np

from catch.utils import interval
from catch.utils import longest_common_substring
from catch.utils import timeout

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class Probe:
    """Immutable sequence representing a probe/bait.
    """

    def __init__(self, seq):
        """
        Args:
            seq: np.array representing the sequence of a probe
        """
        self.seq = seq
        self.seq_str = ''.join(seq)
        self.is_flanking_n_string = False
        self.header = None

        self.kmers = defaultdict(set)
        self.kmers_rand_choices = defaultdict(lambda: defaultdict(set))

    def mismatches(self, other):
        """Count number of mismatches with other.

        Args:
            other: another Probe, which must be of the same length as self

        Returns:
            number of mismatches between self and other
        """
        return self.mismatches_at_offset(other, 0)

    def mismatches_at_offset(self, other, offset):
        """Count number of mismatches with other given shift.

        Args:
            other: another Probe, which must be of the same length as self
            offset: number of bp by which to shift 'other'; can be negative
                (corresponding to 'other' being shifted left) or
                positive (corresponding to 'other' being shifted right)

        Returns:
            number of mismatches between self and 'other' after 'other' is
            shifted by 'offset' bp
        """
        if len(self.seq) != len(other.seq):
            raise ValueError("Sequences must be of same length")
        if abs(offset) >= len(other.seq):
            raise ValueError("Invalid offset value " + str(offset))
        if offset == 0:
            return np.sum(self.seq != other.seq)
        elif offset < 0:
            return np.sum(self.seq[:offset] != other.seq[-offset:])
        else:
            return np.sum(self.seq[offset:] != other.seq[:-offset])

    def min_mismatches_within_shift(self, other, max_shift):
        """Compute minimum number of mismatches while shifting.

        Args:
            other: another Probe, which must be of the same length as self
            max_shift: number of bp by which to shift 'other' (in both
                directions)

        Returns:
            the minimum number of mismatches between self and 'other' as
            'other' is shifted with an offset between -max_shift and
            +max_shift relative to self
        """
        return min(self.mismatches_at_offset(other, offset)
                   for offset in range(-max_shift, max_shift + 1))

    def longest_common_substring_length(self, other, k):
        """Compute length of longest common substring with other.

        Args:
            other: another Probe
            k: maximum number of mismatches to tolerate in a common
                substring

        Returns:
            length of the longest common substring with at most k
            mismatches between self and other
        """
        l, _, _ = longest_common_substring.k_lcf(self.seq, other.seq, k)
        return l

    def reverse_complement(self):
        """Create reverse complement of this probe.

        Returns:
            a Probe that is the reverse complement of this probe
        """
        rc_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        # Using rc_map.get(b, b) ensures that this can process bases
        # like 'N'. It returns the base itself (e.g., 'N') if it is
        # not either 'A', 'T', 'C', or 'G'.
        rc_seq = np.array([rc_map.get(b, b) for b in self.seq[::-1]],
                          dtype='U1')
        return Probe(rc_seq)

    def with_prepended_str(self, s):
        """Create a probe with 's' prepended to this probe.

        Args:
            s: string to prepend

        Returns:
            a Probe with 's' prepended to the sequence of this probe
        """
        s_seq = np.fromiter(s, dtype='U1')
        new_seq = np.concatenate([s_seq, self.seq])
        return Probe(new_seq)

    def with_appended_str(self, s):
        """Create a probe with 's' appended to this probe.

        Args:
            s: string to append

        Returns:
            a Probe with 's' appended to the sequence of this probe
        """
        s_seq = np.fromiter(s, dtype='U1')
        new_seq = np.concatenate([self.seq, s_seq])
        return Probe(new_seq)

    def construct_kmers(self, k, include_positions=False):
        """Return a list of k-mers in this probe.

        Args:
            k: the number of bp in a k-mer
            include_positions: when True, the set consists of tuples in
                which the first element is a k-mer and the second is its
                position in the probe

        Returns:
            list of all k-mers of length k in this probe, where the list
            is ordered according to the positions of the k-mers in the
            probe
        """
        kmers = []
        for i in range(len(self.seq) - k + 1):
            kmer = self.seq_str[i:(i + k)]
            if include_positions:
                kmers += [(kmer, i)]
            else:
                kmers += [kmer]
        return kmers

    def shares_some_kmers(self, other,
                          k=20,
                          num_kmers_to_test=10,
                          memoize_kmers=True,
                          return_kmer=False):
        """Determine whether this probe likely shares one or more k-mers with other.

        This heuristic outputs whether it is likely that self and
        other share at least one k-mer. Note that, depending on the
        sequences being compared and the parameter values, false negatives
        are a very real possibility.

        False negatives are possible (indeed even likely for two sequences
        with very few k-mers in common) and their probability of occurring
        depends on the number of k-mers in common and on num_kmers_to_test.

        This heuristic is intended primarily for determining whether it
        is possible that two sequences are 'redundant'. If two sequences
        are 'redundant', they ought to (by definition) share many k-mers
        and therefore this heuristic should have a small probability
        of giving outputting False (a false negative). If two sequences
        are not 'redundant', they should have few or no k-mers in common
        and therefore this heuristic should have a small probability of
        outputting True (a false positive). [In this paragraph, the notion
        of a false positive/negative relates to whether the sequences are
        'redundant', whereas in the earlier paragraph it relates to whether
        they share at least one k-mer.]

        In the case where two sequences are 'redundant', assume they have
        N k-mers in common. The probability of this outputting False
        (falsely indicating that the two are likely not redundant) is the
        probability that this heuristic does not select any of those N
        k-mers. That is:
          ( 1 - N/(len(seq)-k+1) )^{num_kmers_to_test}

        In the case where two sequences are not 'redundant', assume they
        are each produced at random with uniformity and independence across
        possible k-mers. The probability of this outputting True (falsely
        indicating that the two might be redundant) is the probability that
        this heuristic finds at least one k-mer in common. That is:
          1 - ( (1-(1/4)^k )^{len(seq)-k+1} )^{num_kmers_to_test}

        Args:
            other: another Probe
            k: the number of bp in a k-mer
            num_kmers_to_test: the number of k-mers to randomly pick from
                one probe and lookup in the other (i.e., the number of
                times to sample)
            memoize_kmers: save the randomly selected k-mers for self and
                for other; randomly selected k-mers is costly, so this is
                useful when this function is called repeatedly on the same
                probe(s)
            return_kmer: when set, if this function determines that the
                probes do share k-mers, then it returns one of the
                k-mers that they share (rather than the value True)

        Returns:
            True or False depending on whether this probe likely shares
            one or more k-mers with other; if True and return_kmer is
            True, then this returns the tuple a k-mer that the two probes
            share rather than the value True
        """
        if memoize_kmers:
            # Construct the k-mers for self and other if they have
            # not yet been constructed for the given k
            if len(self.kmers[k]) == 0:
                self.kmers[k] = set(self.construct_kmers(k))
            if len(other.kmers[k]) == 0:
                other.kmers[k] = set(other.construct_kmers(k))

            if len(self.kmers_rand_choices[k][num_kmers_to_test]) == 0:
                # Instead of using self.kmers[k], redetermine a list
                # of kmers. self.kmers[k] is a set, and so if this
                # probe is repetitive (i.e., k-mers appear multiple
                # times) those k-mers would only be weighted once
                # in the random selection; we want them to be picked
                # more often if they appear more often
                kmers_list = self.construct_kmers(k)
                rand_kmers = np.random.choice(kmers_list,
                                              size=num_kmers_to_test,
                                              replace=True)
                # rand_kmers can be treated as a set because if
                # a k-mer appears multiple times in rand_kmers, it
                # only needs to be compared against other k-mers
                # once
                rand_kmers = set(rand_kmers)
                # Memoize the random choices too because the calls to
                # list(self.kmers[k]) and to np.random.choice are
                # slow
                self.kmers_rand_choices[k][num_kmers_to_test] = \
                    rand_kmers
            else:
                rand_kmers = \
                    self.kmers_rand_choices[k][num_kmers_to_test]

            kmers_intrst = rand_kmers & other.kmers[k]
            if kmers_intrst:
                # Pick out an arbitrary k-mer from the intersection
                # of rand_kmers and other.kmers[k] (i.e., one that
                # is shared between self and other)
                shared_kmer = next(iter(kmers_intrst))
                return shared_kmer if return_kmer else True
            else:
                return False
        else:
            rand_kmer_positions = np.random.randint(0,
                                                    len(self.seq) - k + 1,
                                                    num_kmers_to_test)
            for n in range(num_kmers_to_test):
                # Read a random k-mer from self and explicitly test for
                # its presence in other
                rand_kmer_pos = rand_kmer_positions[n]
                rand_kmer = self.seq_str[rand_kmer_pos:(rand_kmer_pos + k)]
                if rand_kmer in other.seq_str:
                    return rand_kmer if return_kmer else True
            return False

    def identifier(self, length=10):
        """Return an identifier for this probe, based on its sequence.

        The identifier is probably unique among all the probes being
        considered.

        The identifier is computed from a hash of this probe's sequence
        (self.seq_str); it is the final 'length' hex digits of the
        hash. Python's hash(..) function could be used, but the size of
        the hashes it produces depends on the size of the input (longer
        input yield larger hashes); using the SHA-224 hash function
        should produce more uniform hash values.

        For example, when length=10, this is equivalent to taking the final
        40 bits of the SHA-224 digest since each hex digit is 4 bits.
        Thus, it is the SHA-224 digest modulo 2^40. There are 2^40 (roughly
        one trillion) possible identifiers for a probe.

        Returns:
            a (probably) unique identifier for this probe, as a string
        """
        return hashlib.sha224(self.seq_str.encode()).hexdigest()[-length:]

    def __hash__(self):
        return hash(self.seq_str)

    def __eq__(self, other):
        return isinstance(other, Probe) and \
            np.array_equal(self.seq, other.seq)

    def __cmp__(self, other):
        c = np.where(self.seq != other.seq)[0]
        if len(c) == 0:
            return 0
        else:
            # c[0] holds the first index where a char in self.seq does
            # not equal the corresponding char in other.seq
            return cmp(self.seq[c[0]], other.seq[c[0]])

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, i):
        return self.seq[i]

    def __str__(self):
        return self.seq_str

    def __repr__(self):
        return self.seq_str

    @staticmethod
    def from_str(seq_str):
        """Construct a Probe from a string.

        Args:
            seq_str: sequence as a Python string

        Returns:
            instance of Probe, whose sequence is seq_str
        """
        return Probe(np.fromiter(seq_str, dtype='U1'))


def _construct_rand_kmer_probe_map(probes,
                                   k=20,
                                   num_kmers_per_probe=20,
                                   include_positions=False):
    """Construct k-mer/probe map by choosing k-mers randomly.

    Given a collection of probes, this finds the k-mers (of length k)
    in each probe and randomly selects num_kmers_per_probe from each.
    Then, it builds a map from the randomly selected k-mers to a set
    of probes from which the k-mer is located. If more than one probe
    share a k-mer that is randomly selected from those probes, then
    each of those probes are in the set mapped to by the shared k-mer.

    Args:
        probes: list of probes from which to construct the map
        k: the number of bp in a k-mer
        num_kmers_per_probe: the number of k-mers to add to the map
            (as keys) for each probe
        include_positions: when True, the set mapped to by each k-mer
            key consists of tuples in which the first element is a probe
            containing the k-mer and the second is the k-mer's position
            in the probe; if a k-mer appears more than once in a probe
            and it is randomly selected more than once, that probe may
            appear in more than one tuple mapped to by that k-mer

    Returns:
        dict mapping k-mers to sets of probes that contains those
        k-mers
    """
    kmer_probe_map = defaultdict(set)
    for probe in probes:
        if k > len(probe.seq):
            raise ValueError("k is larger than the length of a probe")
        kmers = probe.construct_kmers(k, include_positions)
        if include_positions:
            # np.random.choice won't directly pick tuples from a list,
            # so instead randomly select indices
            rand_kmers = [kmers[i]
                          for i in np.random.choice(len(kmers),
                                                    size=num_kmers_per_probe,
                                                    replace=True)]
            for kmer, pos in rand_kmers:
                kmer_probe_map[kmer].add((probe, pos))
        else:
            rand_kmers = np.random.choice(kmers,
                                          size=num_kmers_per_probe,
                                          replace=True)
            for kmer in rand_kmers:
                kmer_probe_map[kmer].add(probe)
    return dict(kmer_probe_map)


class PigeonholeRequiresTooSmallKmerSizeError(Exception):
    """The pigeonhole approach requires a k-mer length that is too small
    """
    pass


def _construct_pigeonholed_kmer_probe_map(probes,
                                          mismatches,
                                          min_k=20,
                                          include_positions=False):
    """Construct k-mer/probe map by pigeonholing mismatches into k-mers.

    Given a collection of probes (all of the same length) and some
    number of mismatches, this selects a k-mer length k so that k-mers
    can be selected (non-overlapping) to ensure that all mismatches are
    pigeonholed into the k-mers such that there is at least one k-mer
    without a mismatch. Consider some probe p that is some number of
    mismatches away from another probe q; if all the k-mers in p were
    looked up in the k-mer probe map generated from q as an input,
    this method guarantees that at least one match will be found between
    a k-mer in p and a k-mer in q.

    It builds a map from the k-mers to a set of probes from which
    the k-mer is located. If more than one probe share a k-mer that is
    randomly selected from those probes, then each of those probes are
    in the set mapped to by the shared k-mer.

    Note that this function is intended for use when finding probe
    covers based on longest common substring. For other approaches, it
    may yield too few k-mers per probe.

    Args:
        probes: list of probes from which to construct the map; these
            must all have the same length
        mismatches: number of mismatches that will be tolerated when
            this k-mer probe map is used to search for probe coverage;
            chooses a k-mer length k such that there is always at least
            one k-mer that does not have a mismatch with a sequence that
            a probe covers
        min_k: the smallest k-mer length allowed (if the k-mer length is
            too small, a map from k-mers to probes can become useless
            as its use in find_probe_covers_in_sequence will yield too
            many false positives); an exception will be thrown if the
            chosen value of k is less than this value, which may happen
            if mismatches is too large compared to the probe length or if
            the probe length is not easily divisible
        include_positions: when True, the set mapped to by each k-mer
            key consists of tuples in which the first element is a probe
            containing the k-mer and the second is the k-mer's position
            in the probe; if a k-mer appears more than once in a probe
            and it is randomly selected more than once, that probe may
            appear in more than one tuple mapped to by that k-mer

    Returns:
        dict mapping k-mers to sets of probes that contains those
        k-mers
    """
    # Find the probe length
    if len(probes) == 0:
        return {}
    probe_length = len(probes[0].seq)
    for p in probes:
        if len(p.seq) != probe_length:
            raise ValueError("All probes must have the same length")

    if mismatches == 0:
        # Just one k-mer of length probe_length suffices
        k = probe_length
    else:
        # For some k, let us have probe_length/k k-mers. In the worst-case,
        # we place one mismatch in each k-mer. We want to have at least one
        # k-mer with no mismatches. So we need probe_length/k > mismatches.
        # That is, k < probe_length/mismatches.
        k = int(probe_length / mismatches)
        if k == float(probe_length) / mismatches:
            # mismatches divides probe_length, so decrement k
            k -= 1
        # We need k to divide probe_length, so keep decrementing k until
        # this is true
        while probe_length % k != 0:
            k -= 1

    if k < min_k:
        raise PigeonholeRequiresTooSmallKmerSizeError()

    # Construct a map with k-mers from each probe separated by k bp
    kmer_probe_map = defaultdict(set)
    for p in probes:
        kmers = p.construct_kmers(k, include_positions)
        for i in range(0, len(p.seq), k):
            if include_positions:
                kmer, pos = kmers[i]
                kmer_probe_map[kmer].add((p, pos))
            else:
                kmer = kmers[i]
                kmer_probe_map[kmer].add(p)
    return dict(kmer_probe_map)


def construct_kmer_probe_map_to_find_probe_covers(probes,
                                                  mismatches,
                                                  lcf_thres,
                                                  min_k=20,
                                                  k=20,
                                                  include_positions=True):
    """Construct map from k-mers to probes that contain these k-mers.

    This wraps around two other functions for constructing k-mer probe
    maps: _construct_rand_kmer_probe_map and
    _construct_pigeonholed_kmer_probe_map. It only calls the "pigeonhole"
    function if all probes have the same length and this length is
    equal to lcf_thres. If this length is greater than lcf_thres, then
    the "pigeonhole" method would not suffice for finding probe covers
    because it selects k-mers by spreading across the full probe length.
    If the "pigeonhole" function fails because it requires too small a
    value for k, then this resorts to calling the "random" function.

    Args:
        probes: list of probes from which to construct the map
        mismatches: number of mismatches that will be tolerated when
            this k-mer probe map is used to search for probe coverage
        lcf_thres: stretch of aligned bp required when this k-mer probe
            map is used to search for probe coverage
        min_k: the smallest k-mer length allowed (if the k-mer length is
            too small, a map from k-mers to probes can become useless
            as its use in find_probe_covers_in_sequence will yield too
            many false positives) when using the pigeonhole approach
        k: when using the random approach, use this k-mer length
        include_positions: when True, the set mapped to by each k-mer
            key consists of tuples in which the first element is a probe
            containing the k-mer and the second is the k-mer's position
            in the probe; if a k-mer appears more than once in a probe
            and it is randomly selected more than once, that probe may
            appear in more than one tuple mapped to by that k-mer

    Returns:
        dict mapping k-mers to sets of probes that contains those
        k-mers
    """
    # Find the probe length
    if len(probes) == 0:
        return {}
    probe_length = len(probes[0].seq)
    probe_lengths_differ = False
    for p in probes:
        if len(p.seq) != probe_length:
            probe_lengths_differ = True
            break

    if probe_lengths_differ or lcf_thres < probe_length:
        # Use the random method with its default values for k and
        # num_kmers_per_probe
        return _construct_rand_kmer_probe_map(
            probes, k=k, include_positions=include_positions)

    # Try the pigeonhole approach
    try:
        return _construct_pigeonholed_kmer_probe_map(
            probes, mismatches, min_k=min_k,
            include_positions=include_positions)
    except PigeonholeRequiresTooSmallKmerSizeError:
        # Resort to the random approach
        return _construct_rand_kmer_probe_map(
            probes, k=k, include_positions=include_positions)


class SharedKmerProbeMap:
    """A read-only kmer_probe_map that can be shared by processes.

    Suppose we are using multiple processes (via Python's multiprocessing
    module) and want these processes to read from the dict kmer_probe_map.
    Because this dict may be large, we do not want to copy it to each of
    the processes; it would be preferable it they could access the dict
    in shared memory. This should be simplified by the fact that the dict
    is intended to be read-only. This is, on some level, possible using
    a kmer_probe_map as a global variable. Since multiprocessing calls
    os.fork, and on Linux a fork gives copy-on-write memory, as long as
    we do not modify the dict it should not be copied to the processes
    that do the reading. Unfortunately, when we access a key from a dict
    in Python, Python's dict implementation increments a reference counter
    for the value of that key. (Python wraps each value in the dict with
    its own object.) When this reference counter is incremented, that is
    a write and therefore the memory page containing the reference is
    copied to the memory space of the process doing the read. Since a page
    is about 4 KB and we have a lot of keys/values in this dict, a
    substantial amount of memory is copied to each process that reads from
    the dict.

    This class avoids that problem by using lower-level structures (from
    the multiprocessing.sharedctypes module) that can be shared by
    processes. It implements the same functionality as kmer_probe_map. 
    While a dict is not an available shared ctype, arrays of primitive types
    are. Thus, this implementation is based on arrays of shared primitive
    types. For example, it does not store instances of probe.Probe, but does
    store just the sequence strings of the probes (as these are all that
    are needed).
    """

    def __init__(self, keys, probe_seqs_ind, probe_pos, probe_seqs, k,
            probe_seqs_to_probe, native_dict):
        """Accepts arrays containing the information of a kmer_probe_map.

        These arrays are allocated using multiprocessing.sharedctypes.RawArray
        because the read/write synchronization (as offered by Array) is not
        needed.

        Args:
            keys: Contains the k-mers (keys) from kmer_probe_map, as strings,
                in sorted order. This is used for lookup. The same k-mer may
                appear multiple times if there are multiple probes that contain
                that k-mer (e.g., ['abc', 'def', 'def', 'ghi'] means that the
                k-mer 'def' appears at positions in two probes.
            probe_seqs_ind: For a k-mer keys[i], the value probe_seqs_ind[i]
                gives an index in the array probe_seqs whose value contains the
                sequence of a probe that contains the k-mer keys[i]. That is,
                the sequence probe_seqs[probe_seqs_ind[i]] contains the k-mer
                keys[i].
            probe_pos: Contains the position of k-mers in probes. A k-mer
                keys[i] appears at the position probe_pos[i] in the probe whose
                sequence is given by probe_seqs[probe_seqs_ind[i]].
            probe_seqs: The sequences of all the probes that appear in values
                in the kmer_probe_map. Note that there may be many k-mers/keys
                that map to the same probe; that probe's sequence appears just
                once in the probe_seqs array. If we were to make a direct
                mapping between indices in keys and indices in probe_seqs, we
                would have to store the same probe sequence many times.
            k: length of the k-mers (as an int) in keys
            probe_seqs_to_probe: dict mapping probe sequences (as strings)
                to the instances of probe.Probe from which these sequences
                came
            native_dict: kmer_probe_map as a native Python dict
        """
        self.keys = keys
        self.probe_seqs_ind = probe_seqs_ind
        self.probe_pos = probe_pos
        self.probe_seqs = probe_seqs
        self.k = k
        self.probe_seqs_to_probe = probe_seqs_to_probe
        self.native_dict = native_dict

    def get(self, kmer):
        """Get the value in kmer_probe_map for the given kmer.

        Args:
            kmer: k-mer (string) to lookup

        Returns:
            list of tuples (seq, pos) where seq is the sequence (string) of
            a probe that contains kmer and pos is the position of kmer in
            the sequence; returns None if kmer is not found as a key
        """
        # The kmers in self.keys are sorted, so do a binary search to
        # find kmer
        kmer_bytes = kmer.encode()
        i = bisect.bisect_left(self.keys, kmer_bytes)
        if i == len(self.keys) or self.keys[i] != kmer_bytes:
            # The key kmer is not present
            return None

        # There may be more than one match for the key kmer, so keep
        # scanning while there is a match
        matches = []
        while i < len(self.keys) and self.keys[i] == kmer_bytes:
            seq = self.probe_seqs[self.probe_seqs_ind[i]].decode()
            pos = self.probe_pos[i]
            matches += [(seq, pos)]
            i += 1
        return matches

    @staticmethod
    def construct(kmer_probe_map):
        """Construct a SharedKmerProbeMap instance from a kmer_probe_map dict.

        Args:
            kmer_probe_map: dict as output by the function
                probe.construct_kmer_probe_map_to_find_probe_covers

        Returns:
            instance of SharedKmerProbeMap that offers the same functionality
            and stores the same information as kmer_probe_map

        Raises:
            ValueError if k-mers have different lengths or the kmer_probe_map
            does not include positions
        """
        # Find the k-mer length k, check that all k-mers in the map are
        # of length k, and check that the k-mers in the map come with
        # positions
        k = None
        for kmer in kmer_probe_map.keys():
            if k is None:
                k = len(kmer)
            if len(kmer) != k:
                raise ValueError("Inconsistent kmer lengths in kmer_probe_map")
            for v in kmer_probe_map[kmer]:
                if not isinstance(v, tuple):
                    raise ValueError(("Given kmer_probe_map must include kmer "
                                      "positions"))

        # First copy all the (unique) probe sequences to an array probe_seqs
        # While filling in probe_seqs, save the position of a seq in probe_seqs
        # as unique_probe_seqs[seq]
        # Also, save a mapping of all the probe sequences back to the instances
        # of Probe
        unique_probe_seqs = {}
        probe_seqs_to_probe = {}
        for kmer, kmer_alignments in kmer_probe_map.items():
            for probe, pos in kmer_alignments:
                unique_probe_seqs[probe.seq_str] = True
                probe_seqs_to_probe[probe.seq_str] = probe
        probe_seqs = multiprocessing.sharedctypes.RawArray(
            ctypes.c_char_p, len(unique_probe_seqs))
        for i, (seq, _) in enumerate(unique_probe_seqs.items()):
            probe_seqs[i] = seq.encode()
            unique_probe_seqs[seq] = i

        # Find the number of keys that are needed and allocate the keys array
        num_keys = sum(len(kmer_alignments)
                       for kmer, kmer_alignments in kmer_probe_map.items())
        keys = multiprocessing.sharedctypes.RawArray(ctypes.c_char_p, num_keys)

        # Allocate the probe_seqs_ind and probe_pos arrays; their indices map
        # to indices in key, so their length is also num_keys
        probe_seqs_ind = multiprocessing.sharedctypes.RawArray(
            ctypes.c_uint, num_keys)
        probe_pos = multiprocessing.sharedctypes.RawArray(
            ctypes.c_uint, num_keys)

        # Fill in keys, probe_seqs_ind, and probe_pos
        i = 0
        for kmer in sorted(kmer_probe_map.keys()):
            num_alignments = len(kmer_probe_map[kmer])
            for probe, pos in kmer_probe_map[kmer]:
                keys[i] = kmer.encode()
                probe_seqs_ind[i] = unique_probe_seqs[probe.seq_str]
                probe_pos[i] = pos
                i += 1

        # Fill in native_dict
        native_dict = defaultdict(list)
        for kmer in kmer_probe_map.keys():
            for probe, pos in kmer_probe_map[kmer]:
                native_dict[kmer].append((probe.seq_str, pos))
        native_dict = dict(native_dict)

        return SharedKmerProbeMap(keys, probe_seqs_ind, probe_pos, probe_seqs,
                                  k, probe_seqs_to_probe, native_dict)


def set_max_num_processes_for_probe_finding_pools(max_num_processes=8):
    """Set the maximum number of processes to use in a probe finding pool.

    Args:
        max_num_processes: an int (>= 1) specifying the maximum number of
            processes to use in a multiprocessing.Pool when a num_processes
            argument is not provided to probe.open_probe_finding_pool; uses
            min(the number of CPUs in the system, max_num_processes) processes
            when num_processes is not provided to
            probe.open_probe_finding_pool
    """
    global _pfp_max_num_processes
    _pfp_max_num_processes = max_num_processes
set_max_num_processes_for_probe_finding_pools()


def open_probe_finding_pool(kmer_probe_map,
                            cover_range_for_probe_in_subsequence_fn,
                            num_processes=None,
                            use_native_dict=False):
    """Open a pool for calling find_probe_covers_in_sequence().

    The variables to share with the processes (e.g., kmer_probe_map.keys)
    cannot be pickled and are not intended to be. But all variables that are
    global in a module (prior to making the pool) are accessible to processes
    in the pool that are executing a top-level function in the module
    (like _find_probe_covers_in_subsequence()). Thus, this function -- along
    with opening a multiprocessing pool -- also makes global the variables
    that should be shared with the worker processes.

    All the global variables that are part of this probe finding pool are
    prefixed with '_pfp'.

    Args:
        kmer_probe_map: instance of SharedKmerProbeMap
        cover_range_for_probe_in_subsequence_fn: function that
            determines whether a probe "covers" a part of a subsequence
            of sequence; if it returns None, there is no coverage;
            otherwise it returns the range of the subsequence covered
            by the probe
        num_processes: number of processes/workers to have in the pool;
            if None, uses min(the number of CPUs in the system,
            _pfp_max_num_processes)
        use_native_dict: use the native Python dict in kmer_probe_map
            rather than the primitive types that are more suited to
            sharing across processes; depending on the input, this can
            result in considerably more memory use (see SharedKmerProbeMap
            for an explanation of why) but may provide an improvement
            in runtime

    Raises:
        RuntimeError if the pool is already open; only one pool may be
        open at a time
    """
    global _pfp_is_open
    global _pfp_max_num_processes
    global _pfp_pool
    global _pfp_work_was_submitted
    global _pfp_cover_range_for_probe_in_subsequence_fn
    global _pfp_kmer_probe_map_keys
    global _pfp_kmer_probe_map_probe_seqs_ind
    global _pfp_kmer_probe_map_probe_pos
    global _pfp_kmer_probe_map_probe_seqs
    global _pfp_kmer_probe_map_probe_seqs_to_probe
    global _pfp_kmer_probe_map_k
    global _pfp_kmer_probe_map_native
    global _pfp_kmer_probe_map_use_native

    try:
        if _pfp_is_open:
            raise RuntimeError("Probe finding pool is already open")
    except NameError:
        pass

    if num_processes is None:
        num_processes = min(multiprocessing.cpu_count(),
                            _pfp_max_num_processes)

    logger.debug("Opening a probe finding pool with %d processes",
                 num_processes)

    _pfp_is_open = True

    _pfp_cover_range_for_probe_in_subsequence_fn = \
        cover_range_for_probe_in_subsequence_fn

    # Rather than saving kmer_probe_map directly, save pointers to individual
    # variables and have the function a process executes reconstruct an
    # instance of SharedKmerProbeMap from these variables
    # This way, we can be careful to only share in memory with other processes
    # variables that do not have to be copied -- i.e., those processes explicitly
    # access certain variables and we can ensure that variables that would
    # need to be copied (like kmer_probe_map.probe_seqs_to_probe) are not
    # accidentally accessed by a process
    _pfp_kmer_probe_map_keys = kmer_probe_map.keys
    _pfp_kmer_probe_map_probe_seqs_ind = kmer_probe_map.probe_seqs_ind
    _pfp_kmer_probe_map_probe_pos = kmer_probe_map.probe_pos
    _pfp_kmer_probe_map_probe_seqs = kmer_probe_map.probe_seqs
    _pfp_kmer_probe_map_probe_seqs_to_probe = \
        kmer_probe_map.probe_seqs_to_probe
    _pfp_kmer_probe_map_k = kmer_probe_map.k
    _pfp_kmer_probe_map_native = kmer_probe_map.native_dict
    _pfp_kmer_probe_map_use_native = use_native_dict

    # Note that the pool must be created at the very end of this function
    # because the only global variables shared with processes in this
    # pool are those that are created prior to creating the pool

    # Sometimes opening a pool (via multiprocessing.Pool) hangs indefinitely,
    # particularly when many pools are opened/closed repeatedly by a master
    # process; this likely stems from issues in multiprocessing.Pool. So set
    # a timeout on opening the pool, and try again if it times out. It
    # appears, from testing, that opening a pool may timeout a few times in
    # a row, but eventually succeeds.
    time_limit = 60
    while True:
        try:
            with timeout.time_limit(time_limit):
                _pfp_pool = multiprocessing.Pool(num_processes)
            break
        except timeout.TimeoutException:
            # Try again
            logger.debug("Pool initialization timed out; trying again")
            time_limit *= 2
            continue

    _pfp_work_was_submitted = False
    logger.debug("Successfully opened a probe finding pool")


def close_probe_finding_pool():
    """Close the pool for calling find_probe_covers_in_sequence().

    This closes the multiprocessing pool and also deletes pointers to the
    variables that were made global in this module in order to be shared
    with worker processes.

    Raises:
        RuntimeError if the pool is not open
    """
    global _pfp_is_open
    global _pfp_pool
    global _pfp_work_was_submitted
    global _pfp_cover_range_for_probe_in_subsequence_fn
    global _pfp_kmer_probe_map_keys
    global _pfp_kmer_probe_map_probe_seqs_ind
    global _pfp_kmer_probe_map_probe_pos
    global _pfp_kmer_probe_map_probe_seqs
    global _pfp_kmer_probe_map_probe_seqs_to_probe
    global _pfp_kmer_probe_map_k
    global _pfp_kmer_probe_map_native
    global _pfp_kmer_probe_map_use_native

    pfp_is_open = False
    try:
        if _pfp_is_open:
            pfp_is_open = True
    except NameError:
        pass
    if not pfp_is_open:
        raise RuntimeError("Probe finding pool is not open")

    logger.debug("Closing the probe finding pool of processes")

    del _pfp_cover_range_for_probe_in_subsequence_fn

    del _pfp_kmer_probe_map_keys
    del _pfp_kmer_probe_map_probe_seqs_ind
    del _pfp_kmer_probe_map_probe_pos
    del _pfp_kmer_probe_map_probe_seqs
    del _pfp_kmer_probe_map_probe_seqs_to_probe
    del _pfp_kmer_probe_map_k
    del _pfp_kmer_probe_map_native
    del _pfp_kmer_probe_map_use_native

    # In Python versions earlier than 2.7.3 there is a bug (see
    # http://bugs.python.org/issue12157) that occurs if a pool p is
    # created and p.join() is called, but p.map() is never called (i.e.,
    # no work is submitted to the processes in the pool); the bug
    # causes p.join() to sometimes hang indefinitely.
    # That could happen here if a probe finding pool is opened/closed
    # but find_probe_covers_in_sequence() is never called; the variable
    # _pfp_work_was_submitted ensures that join() is only called on
    # the pool if work was indeed submitted.
    # Similarly, when no work is submitted, a call to p.close() may yield
    # a RuntimeError that is printed but ignored; so only call close()
    # when work was indeed submitted.
    if _pfp_work_was_submitted:
        _pfp_pool.close()
        # Due to issues that likely stem from bugs in the multiprocessing
        # module, calls to _pfp_pool.terminate() and _pfp_pool.join()
        # sometimes hang indefinitely (even when work was indeed submitted
        # to the processes). So make a best effort in calling these functions
        # -- i.e., use a timeout around calls to these functions
        try:
            with timeout.time_limit(60):
                _pfp_pool.terminate()
        except timeout.TimeoutException:
            # Ignore the timeout
            # If _pfp_pool.terminate() or _pfp_pool.join() fails this will
            # not affect correctness and will not necessarily prevent
            # additional pools from being created, so let the program continue
            # to execute because it will generally be able to keep making
            # progress
            logger.debug(("Terminating the probe finding pool timed out; "
                          "ignoring"))
            pass
        except:
            # _pfp_pool.terminate() occassionally raises another exception
            # (NoneType) if it tries to terminate a process that has already
            # been terminated; ignoring that exception should not affect
            # correctness or prevent additional pools from being created, so
            # is better to ignore it than to let the exception crash the
            # program
            pass

        try:
            with timeout.time_limit(60):
                _pfp_pool.join()
        except timeout.TimeoutException:
            # Ignore the timeout
            # If _pfp_pool.terminate() or _pfp_pool.join() fails this will
            # not affect correctness and will not necessarily prevent
            # additional pools from being created, so let the program continue
            # to execute because it will generally be able to keep making
            # progress
            logger.debug(("Joining the probe finding pool timed out; "
                          "ignoring"))
            pass
        except:
            # Ignore any additional exception from _pfp_pool.join() rather
            # than letting it crash the program
            pass

    del _pfp_pool
    _pfp_is_open = False
    del _pfp_work_was_submitted

    gc.collect()
    logger.debug("Successfully closed the probe finding pool")


def _find_probe_covers_in_subsequence(bounds,
                                      sequence,
                                      merge_overlapping=True):
    """Helper function for find_probe_covers_in_sequence().

    Scans through a subsequence of sequence, as specified by bounds, and
    looks for probes that cover a range of the subsequence.

    Args:
        bounds: tuple of the form (start, end); scan through each k-mer
            in sequence beginning with the k-mer whose first base is
            at start and ending with the k-mer whose first base is at
            end-1
        sequence: sequence (as a string) in which to find ranges that
            probes cover
        merge_overlapping: when True, merges overlapping ranges into
            a single range and returns the ranges in sorted order; when
            False, intervals returned may be overlapping (e.g., if a
            probe covers two regions that overlap)

    Returns:
        dict mapping probe sequences (as strings) to the set of ranges
        (each range is a tuple of the form (start, end)) that each probe
        "covers" in the scanned subsequence
    """
    if bounds is None:
        return {}

    global _pfp_cover_range_for_probe_in_subsequence_fn
    global _pfp_kmer_probe_map_keys
    global _pfp_kmer_probe_map_probe_seqs_ind
    global _pfp_kmer_probe_map_probe_pos
    global _pfp_kmer_probe_map_probe_seqs
    global _pfp_kmer_probe_map_k
    global _pfp_kmer_probe_map_use_native

    if _pfp_kmer_probe_map_use_native:
        global _pfp_kmer_probe_map_native
        shared_kmer_probe_map = _pfp_kmer_probe_map_native
    else:
        shared_kmer_probe_map = SharedKmerProbeMap(
            _pfp_kmer_probe_map_keys,
            _pfp_kmer_probe_map_probe_seqs_ind,
            _pfp_kmer_probe_map_probe_pos,
            _pfp_kmer_probe_map_probe_seqs,
            _pfp_kmer_probe_map_k,
            None,
            None)
    k = _pfp_kmer_probe_map_k
    # Each time a probe is found to cover a range of sequence,
    # add that range, as a tuple, to the probe's entry in
    # subseq_probe_cover_ranges
    start, end = bounds
    subseq_probe_cover_ranges = defaultdict(list)
    for i in range(start, end):
        kmer = sequence[i:(i + k)]
        # Find the probes with this kmer (with the potential to miss
        # some probes due to false negatives)
        probes_to_align = shared_kmer_probe_map.get(kmer)
        if probes_to_align is None:
            # No probes (from kmer_probe_map) share this kmer
            continue
        for probe_seq_str, pos in probes_to_align:
            # kmer appears in probe at position pos. So align probe
            # to sequence at i-pos and see how much of the subsequence
            # starting here the probe covers.
            probe_seq_full = np.fromiter(probe_seq_str, dtype='U1')
            subseq_left = max(0, i - pos)
            subseq_right = min(len(sequence), i - pos + len(probe_seq_full))
            subsequence = sequence[subseq_left:subseq_right]
            if i - pos < 0:
                # An edge case where probe is cutoff on left end because it
                # extends further left than where sequence begins
                probe_seq = probe_seq_full[-(i - pos):]
                # Shift kmer_start left from pos to determine its new
                # position in probe_seq (equivalently its position in
                # subsequence, which is i)
                kmer_start = pos + (i - pos)
            elif i - pos + len(probe_seq_full) > len(sequence):
                # An edge case where probe is cutoff on right end because it
                # extends further right than where sequence ends
                probe_seq = probe_seq_full[:-(i - pos + len(probe_seq_full) -
                                            len(sequence))]
                kmer_start = pos
            else:
                probe_seq = probe_seq_full
                kmer_start = pos
            cover_range = \
                _pfp_cover_range_for_probe_in_subsequence_fn(
                    probe_seq, subsequence, kmer_start, kmer_start + k,
                    len(probe_seq_full), len(sequence))
            if cover_range is None:
                # probe does not meet the threshold for covering this
                # subsequence
                continue
            cover_start, cover_end = cover_range
            # cover_start and cover_end are relative to subsequence, so
            # adjust these to be relative to sequence
            cover_start += subseq_left
            cover_end += subseq_left
            subseq_probe_cover_ranges[probe_seq_str].append(
                (cover_start, cover_end))
            if merge_overlapping:
                # Save some memory in each process by merging cover ranges,
                # since many found by this method will overlap
                # (This is not necessary because all the cover ranges for
                # each probe will be merged across processes at the end of
                # find_probe_covers_in_sequence(), but it can save
                # considerable memory before that final merge.)
                subseq_probe_cover_ranges[probe_seq_str] = interval.\
                    merge_overlapping(subseq_probe_cover_ranges[probe_seq_str])
    return dict(subseq_probe_cover_ranges)


def find_probe_covers_in_sequence(sequence,
                                  merge_overlapping=True):
    """Find ranges in sequence that a collection of probes cover.

    This uses multiple processes to scan through sequence in parallel.
    Prior to calling this function, a pool of processes must have been
    created by calling open_probe_finding_pool(). That function takes
    arguments (like kmer_probe_map and cover_range_for_probe_in_sequence_fn)
    that the worker processes use in finding ranges that the probes cover.
    Those variables are made global in this module so that the worker
    processes can access them without having to copy the memory.

    Probes are from the values of kmer_probe_map. A probe is said
    to "cover" (i.e., hybridize to) a region as determined by the
    function cover_range_for_probe_in_subsequence_fn.

    This works by scanning through sequence, reading the k-mer at each
    position, and looking this up in kmer_probe_map to retrieve a set of
    probes that are "candidates" for covering sequence around the current
    position. This aligns each candidate probe to sequence around the
    shared k-mer and then calls cover_range_for_probe_in_sequence_fn to
    determine whether the probe "covers" sequence in this region, where
    coverage is determined by that function.

    This determines the k-mer length k from the kmer_probe_map that was
    given to open_probe_finding_pool().

    Note that kmer_probe_map could be generated in a manner that randomly
    selects a subset of the k-mers from each probe (i.e., by
    _construct_rand_kmer_probe_map). In such a case, this algorithm
    is a Monte Carlo algorithm and may yield false negatives. That is,
    it is possible that when scanning through the sequence we encounter
    a region that a probe ought to cover, but none of the k-mers in
    that region map to the probe in kmer_probe_map. This should be
    unlikely as long as kmer_probe_map is generated with a small enough
    k and large enough num_kmers_per_probe. (There is a balance: as k
    decreases and num_kmers_per_probe increases there is a larger set
    of "candidate" probes for covering a region, and the runtime of this
    function increases.) Assume that there is a region that should be
    covered by some probe of length L and that the region and the probe
    share N k-mers. (Note that as k decreases, N should increase.) The
    probability that this function does not detect the coverage is the
    probability that none of those N k-mers in the probe are selected
    when generating kmer_probe_map. That is:
        ( 1 - N/(L-k+1) )^{num_kmers_per_probe}
    where num_kmers_per_probe is a parameter used when constructing
    kmer_probe_map.

    Args:
        sequence: sequence (as a string) in which to find ranges that
            probes cover
        merge_overlapping: when True, merges overlapping ranges into
            a single range and returns the ranges in sorted order; when
            False, intervals returned may be overlapping (e.g., if a
            probe covers two regions that overlap)

    Returns:
        dict mapping probes to the set of ranges (each range is a tuple
        of the form (start, end)) that each probe "covers"

    Raises:
        RuntimeError if a pool for finding probes is not open; a pool
        must be opened prior to calling this function by calling
        open_probe_finding_pool()
    """
    global _pfp_is_open
    global _pfp_pool
    global _pfp_work_was_submitted
    global _pfp_kmer_probe_map_probe_seqs_to_probe
    global _pfp_kmer_probe_map_k

    pfp_is_open = False
    try:
        if _pfp_is_open:
            pfp_is_open = True
    except NameError:
        pass
    if not pfp_is_open:
        raise RuntimeError("Probe finding pool is not open")

    k = _pfp_kmer_probe_map_k

    # Setup a function that the processes can execute; do this using
    # functools.partial so that the created function (scan_subsequence)
    # takes just the argument 'bounds' and all the other arguments to
    # _find_probe_covers_in_subsequence are filled in
    scan_subsequence = partial(_find_probe_covers_in_subsequence,
                               sequence=sequence,
                               merge_overlapping=merge_overlapping)

    # Create bounds for each process
    # The first num_processes-1 processes should be given bounds
    # with size bound_size, and the final process may have a smaller
    # range to scan
    # (Rather than having processes that are never sent any work -- which
    # seems to sometimes cause trouble for a multiprocessing Pool -- send
    # 'None' as the bounds for this process; the helper function
    # _find_probe_covers_in_subsequence treats bounds='None' as a NO-OP)
    num_processes = _pfp_pool._processes
    bounds_size = int((len(sequence) - k + 1) / num_processes + 1)
    bounds_by_process = []
    for start in range(0, len(sequence) - k + 1, bounds_size):
        end = min(len(sequence) - k + 1, start + bounds_size)
        bounds_by_process += [(start, end)]
    while len(bounds_by_process) < num_processes:
        bounds_by_process += [None]

    # Run the processes
    try:
        _pfp_work_was_submitted = True
        all_subseq_probe_cover_ranges = _pfp_pool.map(scan_subsequence,
                                                      bounds_by_process)
    except KeyboardInterrupt:
        _pfp_pool.terminate()
        _pfp_pool.join()

    # Merge the outputs from the different processes. Namely:
    # all_subseq_probe_cover_ranges is a list of dicts, where each
    # dict is keyed on probe sequences and has values that are lists.
    # Merge these to create one dict, keyed on probes, by concatenating
    # all the lists (across the dicts) for each probe.
    probe_cover_ranges = defaultdict(list)
    for subseq_probe_cover_ranges in all_subseq_probe_cover_ranges:
        for probe_seq, cover_ranges in subseq_probe_cover_ranges.items():
            probe = _pfp_kmer_probe_map_probe_seqs_to_probe[probe_seq]
            probe_cover_ranges[probe].extend(cover_ranges)

    # It's possible that the list of cover ranges for a probe has
    # overlapping ranges. Clean the list of cover ranges by "merging"
    # overlapping ones, if desired. Also, convert the defaultdict to
    # a regular dict.
    probe_cover_ranges_cleaned = {}
    for probe, cover_ranges in probe_cover_ranges.items():
        if merge_overlapping:
            probe_cover_ranges_cleaned[probe] = interval.\
                merge_overlapping(cover_ranges)
        else:
            # Remove duplicate cover ranges
            probe_cover_ranges_cleaned[probe] = sorted(list(set(cover_ranges)))
    return probe_cover_ranges_cleaned


def probe_covers_sequence_by_longest_common_substring(mismatches,
                                                      lcf_thres,
                                                      island_of_exact_match=0):
    """Return a function that determines coverage of a probe in a sequence.

    The returned function lcf takes a probe sequence (probe_seq) and a
    sequence (intended to be the same length), as well as the indices
    of a shared k-mer around which both are anchored/aligned. That is,
    it should be true that
      probe_seq[kmer_start:kmer_end] == sequence[kmer_start:kmer_end].
    lcf computes the longest common substring based around this anchor
    that has at most 'mismatches' mismatches. If lcf is below the
    specified length (lcf_thres) for the probe to be "covering" a
    portion of sequence, then lcf returns None. Otherwise, we say that
    a portion (namely, the common substring) of the probe covers
    sequence, and lcf returns the range (shared by both the probe
    sequence and sequence) of sequence that the probe covers, where
    the range is the bounds of the longest common substring. lcf also
    accepts the length of the probe (probe_seq can actually be a
    subsequence of the probe if the probe hangs off the end of a
    sequence, so len(probe_seq) is not helpful) and the length of the
    sequence (the argument sequence is generally a subsequence of the
    full sequence, so len(sequence) is not helpful); if either is less than
    lcf_thres, then lcf_thres is taken to be the length of the full
    probe or the full sequence, whichever is shorter.

    Furthermore, lcf requires that there be a longest common substring
    with 0 mismatches, based around the anchor, of length at least
    'island_of_exact_match'. If this is not met, lcf returns None.
    This helps simulate hybridization, which often requires an island
    of 100% identity between a probe and a fragment in order for the
    probe to hybridize to (and capture) the fragment; the length of
    this island is a parameter, often about 20-30 bp. When
    'island_of_exact_match' is unset and given a default value of 0,
    this requirement not effectively not applied.

    Args:
        mismatches/lcf_thres: if the length of the longest common
            substring with at most 'mismatches' mismatches is >=
            'lcf_thres', then the returned function (lcf) outputs the
            bounds of the longest common substring; otherwise it outputs
            None, indicating that the provided probe that does cover
            the provided sequence
        island_of_exact_match: in order for the returned function (lcf)
            to output bounds, require that there be an exact match of
            at least length 'island_of_exact_match' between the probe
            and the provided sequence (i.e., a longest common substring
            with 0 mismatches that has this length)

    Returns:
        function that, given a probe and sequence anchored at a shared
        k-mer, returns whether the probe covers part of the sequence and,
        if so, which part
    """
    def lcf(probe_seq, sequence, kmer_start, kmer_end,
            full_probe_len, full_sequence_len):
        l, start = longest_common_substring.k_lcf_around_anchor(
            probe_seq, sequence, kmer_start, kmer_end, mismatches)
        if l < min(lcf_thres, full_probe_len, full_sequence_len):
            return None

        if island_of_exact_match > 0:
            if mismatches == 0:
                exact_match_l = l
            else:
                exact_match_l, _ = longest_common_substring.k_lcf_around_anchor(
                    probe_seq, sequence, kmer_start, kmer_end, 0)
            if exact_match_l < island_of_exact_match:
                return None

        return (start, start + l)

    return lcf
