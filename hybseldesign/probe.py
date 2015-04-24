"""Structure(s) and functions for directly working with probes.
"""

from collections import defaultdict
import hashlib

import numpy as np

from hybseldesign.utils import interval
from hybseldesign.utils import longest_common_substring

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class Probe:
    """Immutable sequence representing a probe/bait.
    """

    def __init__(self, seq):
        """
        Args:
            seq: np.array representing the sequence of a probe
        """
        self.seq = seq
        self.seq_str = seq.tostring()
        self.is_flanking_n_string = False
        self.header = None

        self.kmers = defaultdict(set)
        self.kmers_rand_choices = defaultdict(lambda: defaultdict(list))

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
                   for offset in xrange(-max_shift, max_shift + 1))

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
                          dtype='S1')
        return Probe(rc_seq)

    def with_prepended_str(self, s):
        """Create a probe with 's' prepended to this probe.

        Args:
            s: string to prepend

        Returns:
            a Probe with 's' prepended to the sequence of this probe
        """
        s_seq = np.fromstring(s, dtype='S1')
        new_seq = np.concatenate([s_seq, self.seq])
        return Probe(new_seq)

    def with_appended_str(self, s):
        """Create a probe with 's' appended to this probe.

        Args:
            s: string to append

        Returns:
            a Probe with 's' appended to the sequence of this probe
        """
        s_seq = np.fromstring(s, dtype='S1')
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
        for i in xrange(len(self.seq) - k + 1):
            kmer = self.seq_str[i:(i + k)]
            if include_positions:
                kmers += [(kmer, i)]
            else:
                kmers += [kmer]
        return kmers

    def shares_some_kmers(self, other,
                          k=10,
                          num_kmers_to_test=10,
                          memoize_kmers=True):
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

        Returns:
            True or False depending on whether this probe likely shares
            one or more k-mers with other
        """
        if memoize_kmers:
            # Construct the k-mers for self and other if they have
            # not yet been constructed for the given k
            if len(self.kmers[k]) == 0:
                self.kmers[k] = set(self.construct_kmers(k))
            if len(other.kmers[k]) == 0:
                other.kmers[k] = set(other.construct_kmers(k))

            if len(self.kmers_rand_choices[k][num_kmers_to_test]) == 0:
                rand_kmers = np.random.choice(list(self.kmers[k]),
                                              size=num_kmers_to_test,
                                              replace=True)
                # Memoize the random choices too because the calls to
                # list(self.kmers[k]) and to np.random.choice are
                # slow
                self.kmers_rand_choices[k][num_kmers_to_test] = \
                    rand_kmers
            else:
                rand_kmers = \
                    self.kmers_rand_choices[k][num_kmers_to_test]

            other_kmers = other.kmers[k]
            for rand_kmer in rand_kmers:
                if rand_kmer in other_kmers:
                    return True
            return False
        else:
            rand_kmer_positions = np.random.random_integers(0,
                                                            len(self.seq) - k,
                                                            num_kmers_to_test)
            for n in xrange(num_kmers_to_test):
                # Read a random k-mer from self and explicitly test for
                # its presence in other
                rand_kmer_pos = rand_kmer_positions[n]
                rand_kmer = self.seq_str[rand_kmer_pos:(rand_kmer_pos + k)]
                if rand_kmer in other.seq_str:
                    return True
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
        return hashlib.sha224(self.seq_str).hexdigest()[-length:]

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
        return Probe(np.fromstring(seq_str, dtype='S1'))


def _construct_rand_kmer_probe_map(probes,
                                   k=10,
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
                                          min_k=10,
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
        k = probe_length / mismatches
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
        for i in xrange(0, len(p.seq), k):
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
                                                  min_k=10,
                                                  k=10,
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


def find_probe_covers_in_sequence(
        sequence,
        kmer_probe_map,
        cover_range_for_probe_in_subsequence_fn=None,
        merge_overlapping=True):
    """Find ranges in sequence that a collection of probes cover.

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

    This determines the k-mer length k from the given kmer_probe_map.

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
        kmer_probe_map: dict mapping kmers to probes; must map k-mers
            (of length k) to a set of tuples, in which each tuple
            contains a probe whose sequence contains the k-mer as well
            as the position of the k-mer in the probe.
        cover_range_for_probe_in_subsequence_fn: function that
            determines whether a probe "covers" a part of a subsequence
            of sequence; if it returns None, there is no coverage;
            otherwise it returns the range of the subsequence covered
            by the probe
        merge_overlapping: when True, merges overlapping ranges into
            a single range and returns the ranges in sorted order; when
            False, intervals returned may be overlapping (e.g., if a probe
            covers two regions that overlap)

    Returns:
        dict mapping probes to the set of ranges (each range is a tuple
        of the form (start, end)) that each probe "covers"
    """
    if cover_range_for_probe_in_subsequence_fn is None:
        # By default, determine a cover range using a longest common
        # substring with its default parameters
        cover_range_for_probe_in_subsequence_fn = \
            probe_covers_sequence_by_longest_common_substring()

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

    # Iterate through every kmer in sequence
    # Each time a probe is found to cover a range of sequence,
    # add that range, as a tuple, to the probe's entry in
    # probe_cover_ranges
    probe_cover_ranges = defaultdict(list)
    for i in xrange(len(sequence) - k + 1):
        kmer = sequence[i:(i + k)]
        # Find the probes with this kmer (with the potential to miss
        # some probes due to false negatives)
        if kmer not in kmer_probe_map:
            # In an earlier version, kmer_probe_map was a defaultdict
            # and this 'if' condition was left out -- if the kmer was not
            # in kmer_probe_map, then probes_to_align (below) would simply
            # be empty. However, this leads to a *huge* memory leak. Each
            # time kmer is looked up in kmer_probe_map, a new set is
            # created by the defaultdict when kmer does not exist as a key.
            # This new set is entered into kmer_probe_map, with kmer as
            # the key, and is never removed. Because of the many possible
            # kmers and the large overhead of a Python set, the size of
            # kmer_probe_map grows enormously (into hundreds of GB very
            # quickly) when many kmers are looked up that do not exist.
            # Thus, a better solution is for kmer_probe_map to be a
            # regular dict, and to have this check instead.
            continue
        probes_to_align = kmer_probe_map[kmer]
        for probe, pos in probes_to_align:
            # kmer appears in probe at position pos. So align probe
            # to sequence at i-pos and see how much of the subsequence
            # starting here the probe covers.
            subseq_left = max(0, i - pos)
            subseq_right = min(len(sequence), i - pos + len(probe.seq))
            subsequence = sequence[subseq_left:subseq_right]
            if i - pos < 0:
                # An edge case where probe is cutoff on left end because it
                # extends further left than where sequence begins
                probe_seq = probe.seq[-(i - pos):]
                # Shift kmer_start left from pos to determine its new
                # position in probe_seq (equivalently its position in
                # subsequence, which is i)
                kmer_start = pos + (i - pos)
            elif i - pos + len(probe.seq) > len(sequence):
                # An edge case where probe is cutoff on right end because it
                # extends further right than where sequence ends
                probe_seq = probe.seq[:-(i - pos + len(probe.seq) -
                                         len(sequence))]
                kmer_start = pos
            else:
                probe_seq = probe.seq
                kmer_start = pos
            cover_range = \
                cover_range_for_probe_in_subsequence_fn(
                    probe_seq, subsequence, kmer_start, kmer_start + k)
            if cover_range is None:
                # probe does not meet the threshold for covering this
                # subsequence
                continue
            cover_start, cover_end = cover_range
            # cover_start and cover_end are relative to subsequence, so
            # adjust these to be relative to sequence
            cover_start += subseq_left
            cover_end += subseq_left
            probe_cover_ranges[probe].append((cover_start, cover_end))

    # It's possible that the list of cover ranges for a probe has
    # overlapping ranges. Clean the list of cover ranges by "merging"
    # overlapping ones, if desired. Also, convert the defaultdict to
    # a regular dict.
    probe_cover_ranges_cleaned = {}
    for probe, cover_ranges in probe_cover_ranges.iteritems():
        if merge_overlapping:
            probe_cover_ranges_cleaned[probe] = interval.\
                merge_overlapping(cover_ranges)
        else:
            # The method above will generally give many of the same cover
            # ranges for a probe, so remove duplicates
            probe_cover_ranges_cleaned[probe] = sorted(list(set(cover_ranges)))
    return probe_cover_ranges_cleaned


def probe_covers_sequence_by_longest_common_substring(mismatches=0,
                                                      lcf_thres=100):
    """Return a function that determines coverage of a probe in a sequence.

    The returned function lcf takes a probe sequence (probe.seq) and a
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
    the range is the bounds of the longest common substring.

    Args:
        mismatches/lcf_thres: if the length of the longest common
            substring with at most 'mismatches' mismatches is >=
            'lcf_thres', then the returned function (lcf) outputs the
            bounds of the longest common substring; otherwise it outputs
            None, indicating that the provided probe that does cover
            the provided sequence

    Returns:
        function that, given a probe and sequence anchored at a shared
        k-mer, returns whether the probe covers part of the sequence and,
        if so, which part
    """
    def lcf(probe_seq, sequence, kmer_start, kmer_end):
        l, start = longest_common_substring.k_lcf_around_anchor(
            probe_seq, sequence, kmer_start, kmer_end, mismatches)
        if l >= lcf_thres:
            return (start, start + l)
        else:
            return None

    return lcf
