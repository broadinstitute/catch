"""Structure(s) storing genomes, as well as methods and functions
for directly working with them.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


"""Immutable collection of one or more sequences (chromosomes)
representing a genome.

This is designed to accomodate both:
 - genomes that are not divided into multiple chromosomes and thus
   do not have a label for even a single sequence
 - genomes that are divided into multiple chromosomes and must
   support the ability to obtain a sequence from a chromosome label
"""
class Genome:
  """seqs is a list of sequences (chromosomes) that make up this
  genome. If len(seqs)==1, then this genome is not divided into
  multiple chromosomes.

  chrs is an OrderedDict mapping the labels of chromosomes to the
  sequence corresponding to the label. It should be true that
  chrs.values()==seqs -- i.e., they contain the same sequences in
  the same order (although this is not asserted).
  """
  def __init__(self, seqs, chrs=None):
    if len(seqs) > 1 and chrs == None:
      raise ValueError(("When there is more than one sequence, chrs "
                        "should also be specified"))
    self.seqs = seqs
    self.chrs = chrs

  """Returns True iff this genome is broken into more than one
  chromosome.
  """
  def divided_into_chrs(self):
    return len(self.seqs) > 1

  """Returns the total length of the genome across all chromosomes.
  """
  def size(self):
    return sum(len(seq) for seq in self.seqs)

  def __hash__(self):
    return hash(self.seqs)

  def __eq__(self, other):
    return isinstance(other, Genome) and \
      self.seqs == other.seqs and \
      self.chrs == other.chrs

  """Constructs and returns a Genome, where seqs_by_chr is an
  OrderedDict mapping chromosome names to sequences.
  """
  @staticmethod
  def from_chrs(seqs_by_chr):
    return Genome(seqs_by_chr.values(), seqs_by_chr)

  """Constructs and returns a Genome with just one sequence
  (chromosome), seq.
  """
  @staticmethod
  def from_one_seq(seq):
    return Genome([seq])

