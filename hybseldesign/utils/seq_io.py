"""Utilities for working with sequence i/o.
"""

from collections import defaultdict
from collections import OrderedDict
import gzip
import logging
import re

import numpy as np

from catch import genome

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def read_dataset_genomes(dataset):
    """Read genomes of the given dataset.

    Args:
        dataset: instance of datasets.GenomesDataset

    Returns:
        list of genome.Genome
    """
    genomes = []

    if dataset.is_multi_chr():
        # The genomes in this dataset have more than one chromosome.
        logger.debug("Reading dataset %s broken up by chromosome",
                     dataset.__name__)
        if dataset.seq_header_to_genome is None:
            # The dataset should have one or more paths to FASTA files
            # (dataset.fasta_paths), where each file gives one genome
            # (sample) and the sequences in that file correspond to the
            # different chromosomes. The header of each sequence should
            # specify a chromosome in dataset.chrs.
            for fn in dataset.fasta_paths:
                seq_map = read_fasta(fn)
                seq_map_by_chr = {dataset.seq_header_to_chr(header): seq_map[header]
                                  for header in seq_map.keys()}
                seqs = OrderedDict(seq_map_by_chr)
                genomes += [genome.Genome.from_chrs(seqs)]
        else:
            # The dataset should have one or more paths to FASTA files
            # (dataset.fasta_paths), where each file gives one or more
            # genomes and their corresponding chromosomes. The header of
            # each sequence should specify a genome ID and a chromosome
            # in dataset.chrs.
            genomes_by_id = defaultdict(list)
            for fn in dataset.fasta_paths:
                seq_map = read_fasta(fn)
                for header, seq in seq_map.items():
                    genome_id = dataset.seq_header_to_genome(header)
                    chrom = dataset.seq_header_to_chr(header)
                    genomes_by_id[genome_id].append((chrom, seq))
            for genome_id, seq_tups in genomes_by_id.items():
                seqs = OrderedDict(seq_tups)
                genomes += [genome.Genome.from_chrs(seqs)]
    else:
        # There is just one sequence (chromosome) for each genome in
        # this dataset. The dataset should have a path to one or more
        # FASTA files (dataset.fasta_paths), and the sequences in each
        # file should correspond to the different genomes. The headers of
        # each sequence are ignored.
        logger.debug("Reading dataset %s with one chromosome per genome",
                     dataset.__name__)
        for fn in dataset.fasta_paths:
            seqs = list(read_fasta(fn).values())
            for seq in seqs:
                genomes += [genome.Genome.from_one_seq(seq)]

    return genomes


def read_genomes_from_fasta(fn):
    """Read genomes from a fasta file.

    This assumes that the genomes each have one chromosome -- i.e.,
    that genomes are not broken up into separate sequences within
    the fasta.

    Returns:
        list of genome.Genome
    """
    logger.debug("Reading fasta %s; assuming one sequence per genome", fn)

    genomes = []
    seqs = list(read_fasta(fn).values())
    for seq in seqs:
        genomes += [genome.Genome.from_one_seq(seq)]
    return genomes


def read_fasta(fn, data_type='str', replace_degenerate=True,
               skip_gaps=True, make_uppercase=True):
    """Read a FASTA file.

    Args:
        fn: path to FASTA file to read
        data_type: determines whether to store a sequence as
            a native Python string ('str') or as a numpy array
            ('np')
        replace_degenerate: when True, replace the degenerate
            bases ('Y','R','W','S','M','K','B','D','H','V')
            with 'N'
        skip_gaps: when True, do not read dashes ('-'), which
            represent gaps
        make_uppercase: when True, change all bases to be
            uppercase

    Returns:
        dict mapping the name of each sequence to the sequence
        itself. The mapping is ordered by the order in which
        the sequence is encountered in the FASTA file; this
        helps in particular with replicating past results,
        where the input order could affect the output.
    """
    logger.info("Reading fasta file %s", fn)

    degenerate_pattern = re.compile('[YRWSMKBDHV]')

    def process(f):
        m = OrderedDict()
        curr_seq_name = ""
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                # Reset the sequence being read on an empty line
                curr_seq_name = ""
                continue
            if curr_seq_name == "":
                # Must encounter a new sequence
                assert line.startswith('>')
            if line.startswith('>'):
                curr_seq_name = line[1:]
                m[curr_seq_name] = ''
            else:
                # Append the sequence
                if make_uppercase:
                    line = line.upper()
                if replace_degenerate:
                    line = degenerate_pattern.sub('N', line)
                if skip_gaps:
                    line = line.replace('-', '')
                m[curr_seq_name] += line
        return m

    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as f:
            m = process(f)
    else:
        with open(fn, 'r') as f:
            m = process(f)

    if data_type == 'str':
        # Already stored sequence as string
        m_converted = m
    elif data_type == 'np':
        m_converted = OrderedDict()
        for seq_name, seq in m.items():
            m_converted[seq_name] = np.fromiter(seq, dtype='U1')
    else:
        raise ValueError("Unknown data_type " + data_type)

    return m_converted


def iterate_fasta(fn, data_type='str', replace_degenerate=True):
    """Scan through a FASTA file and yield each sequence.

    This is a generator that scans through a given FASTA file and,
    upon completing the read of a sequence, yields that sequence.

    Args:
        fn: path to FASTA file to read
        data_type: determines whether to store a sequence as
            a native Python string ('str') or as a numpy array
            ('np')
        replace_degenerate: when True, replace the degenerate
            bases ('Y','R','W','S','M','K','B','D','H','V')
            with 'N'

    Yields:
        each sequence in the FASTA file
    """
    degenerate_pattern = re.compile('[YRWSMKBDHV]')

    def format_seq(seq):
        if data_type == 'str':
            # Already stored as str
            return seq
        elif data_type == 'np':
            return np.fromiter(seq, dtype='U1')
        else:
            raise ValueError("Unknown data_type " + data_tyoe)

    def process(f):
        curr_seq = ''
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                # Skip the blank line
                continue
            if line.startswith('>'):
                # Yield the current sequence (if there is one) and reset the
                # sequence being read
                if len(curr_seq) > 0:
                    yield format_seq(curr_seq)
                curr_seq = ''
            else:
                # Append the sequence
                if replace_degenerate:
                    line = degenerate_pattern.sub('N', line)
                curr_seq += line
        if len(curr_seq) > 0:
            yield format_seq(curr_seq)

    if fn.endswith('.gz'):
        with gzip.open(fn, 'rt') as f:
            yield from process(f)
    else:
        with open(fn, 'r') as f:
            yield from process(f)

def write_probe_fasta(probes, out_fn):
    """Write probe sequences to a FASTA file.

    This writes one probe sequence per line, with a header immediately
    preceding the sequence. If set, the header written is the one in
    probe.Probe.header. If not set, the probe.Probe.identifier() is used.

    Args:
        probes: list of instances of probe.Probe
        out_fn: path to FASTA file to write
    """
    with open(out_fn, 'w') as f:
        for p in probes:
            if p.header:
                f.write('>' + p.header + '\n')
            else:
                f.write('>probe_%s\n' % p.identifier())
            f.write(p.seq_str + '\n')
