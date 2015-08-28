"""Structure(s) for storing and directly working with genomes.
"""

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class Genome:
    """Immutable collection of sequences (chromosomes) representing a genome.

    This is designed to accomodate both:
     - genomes that are not divided into multiple chromosomes and thus
       do not have a label for even a single sequence
     - genomes that are divided into multiple chromosomes and must
       support the ability to obtain a sequence from a chromosome label
    """

    def __init__(self, seqs, chrs=None):
        """
        Args:
            seqs: list of sequences (chromosomes) that make up this genome;
                if len(seqs)==1, then this genome is not divided into
                multiple chromosomes
            chrs: OrderedDict mapping the labels of chromosomes to the
                sequence corresponding to the label; it should be true that
                chrs.values()==seqs -- i.e., they contain the same sequences
                in the same order (although this is not asserted)
        """
        if len(seqs) > 1 and chrs is None:
            raise ValueError(("When there is more than one sequence, chrs "
                              "should also be specified"))
        self.seqs = seqs
        self.chrs = chrs

        self.hash_cached = None
        self.size_cached = None
        self.size_unambig_cached = None

    def divided_into_chrs(self):
        """Return if the genome is broken into more than one chromosome.
        """
        return len(self.seqs) > 1

    def size(self, only_unambig=False):
        """Return the total length of the genome across all chromosomes.

        Args:
            only_unambig: when True, only counts unambiguous bases (i.e.,
                'A', 'T', 'C', and 'G'); when False, counts all bases
        """
        if only_unambig:
            if self.size_unambig_cached is None:
                count = 0
                for b in ['A', 'T', 'C', 'G']:
                    count += sum(seq.count(b) for seq in self.seqs)
                self.size_unambig_cached = count
            return self.size_unambig_cached
        else:
            if self.size_cached is None:
                self.size_cached = sum(len(seq) for seq in self.seqs)
            return self.size_cached

    def __hash__(self):
        if self.hash_cached is None:
            self.hash_cached = hash(tuple(self.seqs))
        return self.hash_cached

    def __eq__(self, other):
        return isinstance(other, Genome) and \
            self.seqs == other.seqs and \
            self.chrs == other.chrs

    @staticmethod
    def from_chrs(seqs_by_chr):
        """Construct a Genome from multiple chromosomes.

        Args:
            seqs_by_chr: OrderedDict mapping chromosomes names to
                sequences

        Returns:
            instance of Genome from the input chromosome sequences
        """
        for seq in seqs_by_chr.itervalues():
            if not isinstance(seq, str):
                raise TypeError("Sequences must be strings")
        return Genome(seqs_by_chr.values(), seqs_by_chr)

    @staticmethod
    def from_one_seq(seq):
        """Construct a Genome from a single sequence.

        The returned Genome is composed of just a single sequence (i.e.,
        is not divided into more than one chromosome).

        Args:
            seq: one sequence, as a string

        Returns:
            instance of Genome from the input sequence
        """
        if not isinstance(seq, str):
            raise TypeError("seq must be a string")
        return Genome([seq])
