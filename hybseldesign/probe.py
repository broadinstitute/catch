"""Structure(s) storing probes, as well as methods and functions
for directly working with them.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from collections import defaultdict
import numpy as np

from hybseldesign.utils import longest_common_substring
from hybseldesign.utils import interval


"""Immutable sequence representing a probe/bait.
"""
class Probe:
  def __init__(self, seq):
    self.seq = seq
    self.seq_str = seq.tostring()
    self.is_flanking_n_string = False

    self.kmers = defaultdict(set)
    self.kmers_rand_choices = defaultdict(lambda : defaultdict(list))

  """Counts the number of mismatches between self and other.

  They must be the same length.
  """
  def mismatches(self, other):
    return self.mismatches_at_offset(other, 0)
  
  """Counts the number of mismatches between self and other after
  other is shifted by offset bp.

  offset can be negative (corresponding to other being shifted left)
  or positive (corresponding to other being shifted right).
  """
  def mismatches_at_offset(self, other, offset):
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

  """Counts the minimum number of mismatches between self and other
  as other is shifted with an offset between -max_shift and
  +max_shift relative to self.
  """
  def min_mismatches_within_shift(self, other, max_shift):
    return min(self.mismatches_at_offset(other, offset) \
                for offset in xrange(-max_shift, max_shift+1))

  """Returns the length of the longest common substring with at most
  k mismatches between self and other.
  """
  def longest_common_substring_length(self, other, k):
    l, _, _ = longest_common_substring.k_lcf(self.seq, other.seq, k)
    return l

  """Returns a probe that is the reverse-complement of this
  probe.

  Using rc_map.get(b, b) ensures that this can process bases
  like 'N'. It returns the base itself (e.g., 'N') if it is
  not either 'A', 'T', 'C', or 'G'.
  """
  def reverse_complement(self):
    rc_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rc_seq = np.array([rc_map.get(b, b) for b in self.seq[::-1]],
                dtype='S1')
    return Probe(rc_seq)

  """Returns the set of k-mers in this probe.

  When include_positions is True, the set consists of tuples in
  which the first element is a kmer and the second is its position
  in the probe.
  """
  def construct_kmers(self, k, include_positions=False):
    kmers = set()
    for i in xrange(len(self.seq)-k+1):
      kmer = self.seq_str[i:(i+k)]
      if include_positions:
        kmers.add((kmer, i))
      else:
        kmers.add(kmer)
    return kmers

  """A heuristic that outputs whether it is likely that self and
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
  """
  def shares_some_kmers(self, other, k=10, num_kmers_to_test=10,
      memoize_kmers=True):
    if memoize_kmers:
      # Construct the k-mers for self and other if they have
      # not yet been constructed for the given k
      if len(self.kmers[k]) == 0:
        self.kmers[k] = self.construct_kmers(k)
      if len(other.kmers[k]) == 0:
        other.kmers[k] = other.construct_kmers(k)

      if len(self.kmers_rand_choices[k][num_kmers_to_test]) == 0:
        rand_kmers = np.random.choice(list(self.kmers[k]),
            size=num_kmers_to_test, replace=True)
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
      rand_kmer_positions = np.random.random_integers(
          0, len(self.seq)-k, num_kmers_to_test)
      for n in xrange(num_kmers_to_test):
        # Read a random k-mer from self and explicitly test for
        # its presence in other
        rand_kmer_pos = rand_kmer_positions[n]
        rand_kmer = self.seq_str[rand_kmer_pos:(rand_kmer_pos+k)]
        if rand_kmer in other.seq_str:
          return True
      return False

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
    return str(self.seq)

  """Make a Probe from a Python str.
  """
  @staticmethod
  def from_str(seqStr):
    return Probe(np.fromstring(seqStr, dtype='S1'))


"""Constructs a map from k-mers to the probes that contain those
k-mers.

Given a collection of probes, this finds the k-mers (of length k)
in each probe and randomly selects num_kmers_per_probe from each.
Then, it builds a map from the randomly selected k-mers to a set
of probes from which the k-mer is located. If more than one probe
share a k-mer that is randomly selected from those probes, then
each of those probes are in the set mapped to by the shared k-mer.
When include_positions is True, the set mapped to by each k-mer
consists of a 2-element tuple in which the first element is the
probe containing the k-mer and the second is the position of the
k-mer in the probe. If a k-mer appears more than once in a probe
and is randomly selected more than once, the probe may appear
in more than one tuple mapped to by that k-mer.
"""
def construct_kmer_probe_map(probes, k=15, num_kmers_per_probe=10,
    include_positions=False):
  kmer_probe_map = defaultdict(set)
  for probe in probes:
    kmers = list(probe.construct_kmers(k, include_positions))
    if include_positions:
      # np.random.choice won't directly pick tuples from a list,
      # so instead randomly select indices
      rand_kmers = [kmers[i] for i in np.random.choice(len(kmers),
          size=num_kmers_per_probe, replace=True)]
      for kmer, pos in rand_kmers:
        kmer_probe_map[kmer].add((probe, pos))
    else:
      rand_kmers = np.random.choice(kmers,
          size=num_kmers_per_probe, replace=True)
      for kmer in rand_kmers:
        kmer_probe_map[kmer].add(probe)
  return dict(kmer_probe_map)


"""Finds the ranges in sequence that each probe in a given set of
probes "covers" (i.e., should hybridize to), where coverage is
determined by a given function.

The probes in the values of kmer_probe_map make up the set of probes
for which we seek coverage. kmer_probe_map must map k-mers (of length
k) to a set of tuples, in which each tuple contains a probe
whose sequence contains the k-mer as well as the position of the
k-mer in the probe. This works by scanning through sequence, reading
the k-mer at each position, and looking this up in kmer_probe_map to
retrieve a set of probes that are "candidates" for covering sequence
around the current position. This aligns each candidate probe to
sequence around the shared k-mer and then calls
cover_range_for_probe_in_sequence_fn to determine whether the probe
"covers" sequence in this region, where coverage is determined by
that function. If that function returns None, there is no coverage,
and otherwise it returns the range covered by the probe. The
output of this function is a dict mapping probes to the set of
ranges that each "covers".

Note that kmer_probe_map is generated in a manner that randomly
selects a subset of the k-mers from each probe. Thus, this algorithm
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
"""
def find_probe_covers_in_sequence(sequence,
    kmer_probe_map, k=15,
    cover_range_for_probe_in_subsequence_fn=None):
  if cover_range_for_probe_in_subsequence_fn == None:
    # By default, determine a cover range using a longest common
    # substring with its default parameters
    cover_range_for_probe_in_subsequence_fn = \
      probe_covers_sequence_by_longest_common_substring()

  # Check that the kmers in the given map are of length k and
  # that the kmers in the given map come with positions
  for kmer in kmer_probe_map.keys():
    if len(kmer) != k:
      raise ValueError("Given kmer has length %d but expected %d" % \
                        (len(kmer), k))
    for v in kmer_probe_map[kmer]:
      if type(v) != tuple:
        raise ValueError(("Given kmer_probe_map must include kmer "
            "positions"))

  # Iterate through every kmer in sequence
  # Each time a probe is found to cover a range of sequence,
  # add that range, as a tuple, to the probe's entry in
  # probe_cover_ranges
  probe_cover_ranges = defaultdict(list)
  for i in xrange(len(sequence)-k+1):
    kmer = sequence[i:(i+k)]
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
      subseq_left = max(0, i-pos)
      subseq_right = min(len(sequence), i-pos+len(probe.seq))
      subsequence = sequence[subseq_left:subseq_right]
      if i-pos < 0:
        # An edge case where probe is cutoff on left end because it
        # extends further left than where sequence begins
        probe_seq = probe.seq[-(i-pos):]
        # Shift kmer_start left from pos to determine its new
        # position in probe_seq (equivalently its position in
        # subsequence, which is i)
        kmer_start = pos + (i-pos)
      elif i-pos+len(probe.seq) > len(sequence):
        # An edge case where probe is cutoff on right end because it
        # extends further right than where sequence ends
        probe_seq = probe.seq[:-(i-pos+len(probe.seq)-len(sequence))]
        kmer_start = pos
      else:
        probe_seq = probe.seq
        kmer_start = pos
      cover_range = \
        cover_range_for_probe_in_subsequence_fn(
            probe_seq, subsequence, kmer_start, kmer_start+k)
      if cover_range == None:
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
  # overlapping ones.
  probe_cover_ranges_merged = defaultdict(list)
  for probe, cover_ranges in probe_cover_ranges.iteritems():
    probe_cover_ranges_merged[probe] = interval.\
        merge_overlapping(cover_ranges)
  return dict(probe_cover_ranges_merged)


"""Returns a function that, given a probe and sequence anchored at
a shared k-mer, returns the portion of the sequence covered by the
probe, where coverage is determined by longest common substring.

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
"""
def probe_covers_sequence_by_longest_common_substring(mismatches=0,
    lcf_thres=100):
  def lcf(probe_seq, sequence, kmer_start, kmer_end):
    l, start = longest_common_substring.k_lcf_around_anchor(
        probe_seq, sequence, kmer_start, kmer_end, mismatches)
    if l >= lcf_thres:
      return (start, start + l)
    else:
      return None
  return lcf
