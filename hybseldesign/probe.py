"""Structure(s) storing probes and methods for directly working
with them.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'

from collections import defaultdict
import numpy as np

from hybseldesign.utils import longest_common_substring


"""Immutable sequence representing a probe/bait.
"""
class Probe:
  def __init__(self, seq):
    self.seq = seq
    self.seq_str = seq.tostring()
    self.is_flanking_n_string = False

    self.kmer_hashes = defaultdict(set)
    self.kmer_hashes_rand_choices = defaultdict(
                                       lambda : defaultdict(list))

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

  """Adds the hash of each k-mer of this probe to self.kmer_hashes.

  This could add the k-mer itself, rather than the hash, to avoid
  false positives (due to collisions) on lookup; however, storing
  the hash rather than a k-mer gives space savings and avoids having
  to compute the hash of a k-mer repeatedly on lookup.

  Use hash(kmer) & 0xffffffff so that the stored hash is 32-bit
  rather than 64-bit; 32-bit suffices here and halves the space
  requirements.
  """
  def construct_kmer_hashes(self, k):
    for i in xrange(len(self.seq)-k+1):
      kmer = self.seq_str[i:(i+k)]
      self.kmer_hashes[k].add(hash(kmer) & 0xffffffff)

  """A heuristic that outputs whether it is likely that self and
  other share at least one k-mer. Note that, depending on the
  sequences being compared and the parameter values, false negatives
  are a very real possibility.

  When memoize_kmer_hashes is True, false positives are possible
  although extremely unlikely; when it is False, there are no false
  positives. In either case, false negatives are possible (indeed
  even likely for two sequences with very few k-mers in common) and
  their probability of occurring depends on the number of k-mers in
  common and on num_kmers_to_test.

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
  def shares_some_kmers(self, other, k=10, num_kmers_to_test=5,
      memoize_kmer_hashes=True):
    if memoize_kmer_hashes:
      # Construct the k-mer hashes for self and other if they have
      # not yet been constructed for the given k
      if len(self.kmer_hashes[k]) == 0:
        self.construct_kmer_hashes(k)
      if len(other.kmer_hashes[k]) == 0:
        other.construct_kmer_hashes(k)

      if len(self.kmer_hashes_rand_choices[k][num_kmers_to_test]) == 0:
        rand_kmer_hashes = np.random.choice(list(self.kmer_hashes[k]),
            size=num_kmers_to_test, replace=True)
        # Memoize the random choices too because the calls to
        # list(self.kmer_hashes[k]) and to np.random.choice are
        # slow
        self.kmer_hashes_rand_choices[k][num_kmers_to_test] = \
            rand_kmer_hashes
      else:
        rand_kmer_hashes = \
            self.kmer_hashes_rand_choices[k][num_kmers_to_test]

      other_kmer_hashes = other.kmer_hashes[k]
      for rand_kmer_hash in rand_kmer_hashes:
        if rand_kmer_hash in other_kmer_hashes:
          # This may be a collision and it may be that in fact
          # self and other do not share this k-mer, but we allow
          # false positives
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
