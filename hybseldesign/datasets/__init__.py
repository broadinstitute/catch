"""Convenience classes and functions for working with datasets.
"""

import importlib
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

    def add_fasta_path(self, path, relative=False):
        """Add a FASTA file giving sequences corresponding to this dataset.

        When is_multi_chr() is True, the FASTA file should give one
        sample/individual genome for this dataset (each sequence in the
        file corresponding to a chromosome). When is_multi_chr() is False,
        the FASTA file can give multiple samples/individual genomes, where
        each sequence in the file corresponds to a genome; adding more than
        one FASTA file is equivalent to specifying one FASTA file that
        is the concatenation of the multiple files.

        Args:
            path: path to FASTA file
            relative: when True, path is assumed to be relative to the
                path of the dataset file (self.__file__)
        """
        if relative:
            path = join(dirname(self.__file__), path)
        self.fasta_paths += [path]


class GenomesDatasetSingleChrom(GenomesDataset):
    """Dataset whose genomes have one chromosome.

    This is used when the genome held by this dataset can be fully
    specified by one sequence (i.e., there are not multiple chromosomes).
    All of the sample/individual genomes should be specified as individual
    sequences across one or more FASTA files added with add_fasta_path().
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
        self.fasta_paths = []


class GenomesDatasetMultiChrom(GenomesDataset):
    """Dataset whose genomes have multiple chromosomes.

    This is used when the genome held by this dataset has more than one
    chromosome.  Each FASTA file offers just one sample/individual genome
    in which each sequence in the file corresponds to a chromosome; there
    is typically more than one FASTA file for the dataset, and they should
    each be added with add_fasta_path().

    If there are many genomes, keeping each in a separate FASTA file may
    result in an unwiedly number of files. If seq_header_to_genome is
    provided, then this will assume that the provided FASTA paths (which
    may be just one) have an ID for each genome in the sequence headers.
    This allows multiple genomes to be placed in a single FASTA file,
    and for the segments of each genome to be grouped.
    """

    def __init__(self, __name__, __file__, __spec__,
                 chrs, seq_header_to_chr, seq_header_to_genome=None):
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
            seq_header_to_genome: (optional) function that takes a
                sequence header from a FASTA file and returns an ID for
                the genome from that header; if not provided, then this
                assumes that all of the sequences in a FASTA file are from
                a single genome
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
        self.seq_header_to_genome = seq_header_to_genome
        self.fasta_paths = []


class DatasetCollection:
    """Immutable collection of datasets.
    """

    def __init__(self, dataset_names):
        """
        Args:
            dataset_names: list of the names of the datasets (modules)
                in this collection
        """
        self.dataset_names = dataset_names

    def import_all(self):
        """Import and return the datasets (modules) in this collection.

        Returns:
            List of tuples (name, dataset) where 'name' is the name of a
            dataset and 'dataset' is the imported dataset object
        """
        datasets = []
        for name in self.dataset_names:
            try:
                dataset = importlib.import_module(
                                'hybseldesign.datasets.' + name)
            except ImportError:
                raise ValueError(("Unknown dataset %s; cannot import into "
                                  "collection") % name)
            datasets += [(name, dataset)]
        return datasets
