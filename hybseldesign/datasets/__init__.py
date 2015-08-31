"""Convenience classes and functions for working with datasets.
"""

import logging
from os.path import dirname
from os.path import join

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class GenomesDataset:
    """Individual (sampled) genomes, usually from the same species.
    """

    def __init__(self, __name__, __file__, __spec__):
        """
        Args:
            __name__: name of this dataset (module)
            __file__: path to the module corresponding to this dataset
            __spec__: import-related information about the module
                corresponding to this dataset
        """
        self.__name__ = __name__
        self.__file__ = __file__
        self.__spec__ = __spec__

    def is_multi_chr(self):
        """Return whether genomes in this dataset have multiple chromosomes.
        """
        try:
            self.chrs
        except AttributeError:
            return False
        else:
            return True


class GenomesDatasetSingleChrom(GenomesDataset):
    """Dataset whose genomes have one chromosome.

    This is used when the genome held by this dataset can be fully
    specified by one sequence (i.e., there are not multiple chromosomes).
    All of the sample/individual genomes should be in a single FASTA file
    set by set_fasta_path().
    """

    def __init__(self, __name__, __file__, __spec__):
        """
        Args:
            __name__: name of this dataset (module)
            __file__: path to the module corresponding to this dataset
            __spec__: import-related information about the module
                corresponding to this dataset
        """
        GenomesDataset.__init__(self, __name__, __file__, __spec__)
        self.fasta_path = None

    def set_fasta_path(self, path, relative=False):
        """Set the FASTA file for the genome samples stored by this dataset.

        Args:
            path: path to FASTA file
            relative: when True, path is assumed to be relative to the
                path of the dataset file (self.__file__)
        """
        if relative:
            path = join(dirname(self.__file__), path)
        self.fasta_path = path


class GenomesDatasetMultiChrom(GenomesDataset):
    """Dataset whose genomes have multiple chromosomes.

    This is used when the genome held by this dataset has more than one
    chromosome.  Each FASTA file offers just one sample/individual genome
    in which each sequence in the file corresponds to a chromosome; there
    is typically more than one FASTA file for the dataset, and they should
    each be added with add_fasta_path().
    """

    def __init__(self, __name__, __file__, __spec__,
                 chrs, seq_header_to_chr):
        """
        Args:
            __name__: name of this dataset (module)
            __file__: path to the module corresponding to this dataset
            __spec__: import-related information about the module
                corresponding to this dataset
            chrs: list of chromosomes for the genome held by this dataset;
                these should be the header of each sequence in the FASTA
                files
            seq_header_to_chr: function that takes a sequence header from
                a FASTA file and returns the chr (as a string) from that
                header
        """
        if len(chrs) == 1:
            logger.critical(("Typically when there is just one chromosome, "
                             "the GenomesDatasetSingleChrom class should "
                             "used. Even though there is just one "
                             "chromosome, GenomesDatasetMultiChrom will "
                             "assume that each sequence (chromosome) is in "
                             "a separate FASTA file."))

        GenomesDataset.__init__(self, __name__, __file__, __spec__)
        self.chrs = chrs
        self.seq_header_to_chr = seq_header_to_chr
        self.fasta_paths = []

    def add_fasta_path(self, path, relative=False):
        """Add FASTA file giving one sample/individual genome for this dataset.

        Args:
            path: path to FASTA file
            relative: when True, path is assumed to be relative to the
                path of the dataset file (self.__file__)
        """
        if relative:
            path = join(dirname(self.__file__), path)
        self.fasta_paths += [path]
