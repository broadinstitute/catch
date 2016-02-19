"""Provides analysis of a probe set for targeted genomes.

This computes the number of bp across each target genome that the
probes cover, as well as the percentage of each target genome that
the probes cover. It computes a percentage against the full
length of the target genome, as well as a percentage against the
length of the target genome when only counting unambiguous bases
(non-'N' bases). It also computes the average coverage that the
probes yield across each target genome: again, one value against
the full target genome and one against just the unambiguous bases.

A few notes about the percentage of bases covered (in each note,
assume that we specified a desired coverage of 100% when designing
the probes, though this could be generalized to any desired coverage):
 - The percentage of bases covered in the full target genome
   (including ambiguous bases) may be less than 100%. This is expected
   if the target genome contains ambiguous bases ('N') because these
   are not covered. However, if this coverage analysis is run on probes
   to which adapters have been added, it may also be less than 100%
   even if the target genome contains no ambiguous bases. The reason is
   that, when this finds the ranges covered by each probe, it includes
   the adapters; some number of mismatches may be "used up" by the
   adapter, leaving fewer mismatches at the end of the true probe for
   aligning to the sequence, and meaning that the bases at the end of
   the true probe fail to cover bases of the target genome that they
   were originally intended to cover. For example, suppose we have
   a probe "ATCG" which was intended to cover the sequence "ATCC" in
   a target genome. Now suppose that (2-base) adapters are added to the
   probe, making it "GGATCGCC" (left adapter is 'GG' and right adapter
   is 'CC'). Suppose that the region in the target genome around "ATCC"
   is "GTATCCCC". And suppose we design the probes allowing for one
   mismatch. The one mismatch is "used up" by the second base on the
   left adapter ('G' in the adapter mismatches with 'T' in the target
   genome), so the base 'C' at the end of "ATCC" is not covered. Running
   the coverage analysis before adding adapters to the probes should
   alleviate this, making the percentage 100% as desired.
 - The percentage of bases covered in the target genome, when only
   counting unambiguous bases, may be less than 100% when adapters have
   been added to probes. The reason is the same as above.
 - The percentage of bases covered in the target genome, when only
   counting unambiguous bases, may be less than 100% even when adapters
   have not been added to probes. This can occur if there are bases
   in between two stretches of ambiguous bases ('N'). It may not be
   possible to cover the (unambiguous) bases in between these two
   stretches using a probe, so these bases then go uncovered. For example,
   consider the region 'NNNNNNATCGNNNNNN'; the four unambiguous bases
   ('ATCG') may not be covered by any probe.
 - The percentage of bases covered in the target genome, when only
   counting unambiguous bases, may be more than 100%. This could happen
   when allowing mismatches. A probe could cover an ambiguous base
   (i.e., 'N') of a target genome by simply using a mismatch when
   aligning to that base. Then, the probes collectively cover more bases
   than there are unambiguous bases in the target genome, making the
   percentage more than 100%.
"""

from collections import defaultdict
import logging

import numpy as np

from hybseldesign import probe
from hybseldesign.utils import interval
from hybseldesign.utils import pretty_print

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


class Analyzer:
    """Methods for testing quality control of a probe set.
    """

    def __init__(self,
                 probes,
                 mismatches,
                 lcf_thres,
                 target_genomes,
                 target_genomes_names=None,
                 cover_extension=0,
                 kmer_probe_map_k=10):
        """
        Args:
            probes: collection of instances of probe.Probe that form a
                complete probe set
            mismatches/lcf_thres: consider a probe to hybridize to a sequence
                if a stretch of 'lcf_thres' or more bp aligns with
                'mismatches' or fewer mismatched bp; used to compute whether
                a probe "covers" a portion of a sequence
            target_genomes: list [g_1, g_2, ..., g_m] of m groupings of
                genomes, where each g_i is a list of genome.Genomes belonging
                to group i. For example, a group may be a species and each g_i
                would be a list of the target genomes of species i.
            target_genomes_names: list [s_1, s_2, ..., s_m] of strings where
                the name of the i'th genome grouping (from target_genomes) is
                s_i. When None, the name of the i'th grouping is "Group i".
            cover_extension: number of bp by which to extend the coverage on
                each side of a probe; a probe "covers" the portion of the
                sequence that it hybridizes to, as well as 'cover_extension'
                bp on each side of that portion
            kmer_probe_map_k: in calls to probe.construct_kmer_probe_map...,
                uses this value as min_k and k
        """
        self.probes = probes
        self.target_genomes = target_genomes
        if target_genomes_names:
            if len(target_genomes_names) != len(target_genomes):
                raise ValueError(("Number of target genome names must be same "
                                  "as the number of target genomes"))
            self.target_genomes_names = target_genomes_names
        else:
            self.target_genomes_names = ["Group %d" % i
                                         for i in range(len(target_genomes))]

        self.mismatches = mismatches
        self.lcf_thres = lcf_thres
        self.cover_range_fn = \
            probe.probe_covers_sequence_by_longest_common_substring(
                mismatches, lcf_thres)
        self.cover_extension = cover_extension
        self.kmer_probe_map_k = kmer_probe_map_k

    def _iter_target_genomes(self, rc_too=False):
        """Yield target genomes across groupings to iterate over.

        Args:
            rc_too: when True, also yields False and True for each
                target genome

        Yields:
            i, j, gnm [, rc]
                - i is the index of a target genome grouping
                - j is the index of a genome in grouping i
                - gnm is an instance of genome.Genome corresponding to
                  to genome j in grouping i
                - if rc_too is True, rc cycles through False and True
                  so that an iterator can take the reverse complement
                  of gnm's sequences
        """
        for i, genomes_from_group in enumerate(self.target_genomes):
            for j, gnm in enumerate(genomes_from_group):
                if rc_too:
                    yield i, j, gnm, False
                    yield i, j, gnm, True
                else:
                    yield i, j, gnm

    def _find_covers_in_target_genomes(self):
        """Find intervals across the target genomes covered by the probe set.

        This considers the given probe set (self.probes) and determines the
        intervals, in each genome of the target genomes (as well as their
        reverse complements), that are covered by the probes. This saves a
        dict, self.target_covers, as follows: self.target_covers[i][j][b]
        is a list of all the intervals covered by the probes in the target
        genome j of grouping i (in the reverse complement of the genome if
        b is True, and provided sequence if b is False).

        The endpoints of the intervals are offset so as to give unique integer
        positions in the genome (e.g., endpoints in the second chromosome
        are offset based on the length of the first chromosome). There may
        be duplicate intervals if two probes cover the same region of a
        sequence.
        """
        logger.info("Finding probe covers across target genomes")
        logger.info("Building map from k-mers to probes")
        # Note that if adapters are added to the probes before this filter
        # is run (which would be typical), then self.lcf_thres will likely
        # be less than the probe length. So the k-mer to probe map will
        # be constructed using the random approach (yielding many k-mers
        # and thus a slower runtime in finding probe covers) rather than
        # the pigeonhole approach.
        kmer_probe_map = probe.SharedKmerProbeMap.construct(
            probe.construct_kmer_probe_map_to_find_probe_covers(
                self.probes, self.mismatches, self.lcf_thres,
                min_k=self.kmer_probe_map_k, k=self.kmer_probe_map_k)
        )
        probe.open_probe_finding_pool(kmer_probe_map,
                                      self.cover_range_fn)

        self.target_covers = {}
        for i, j, gnm, rc in self._iter_target_genomes(True):
            if not rc:
                logger.info(("Computing coverage in grouping %d (of %d), "
                             "with target genome %d (of %d)"), i + 1,
                            len(self.target_genomes), j + 1,
                            len(self.target_genomes[i]))
            if i not in self.target_covers:
                self.target_covers[i] = {}
            if j not in self.target_covers[i]:
                self.target_covers[i][j] = {False: None, True: None}

            gnm_covers = []
            length_so_far = 0
            for sequence in gnm.seqs:
                if rc:
                    # Take the reverse complement of sequence
                    rc_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                    sequence = ''.join([rc_map.get(b, b)
                                       for b in sequence[::-1]])

                # Find cover ranges of the probes, while allowing the ranges
                # to overlap (e.g., if one probe covers two regions that
                # overlap)
                probe_cover_ranges = probe.find_probe_covers_in_sequence(
                    sequence,
                    merge_overlapping=False)
                for p, cover_ranges in probe_cover_ranges.items():
                    for cover_range in cover_ranges:
                        # Extend the range covered by probe p on both sides
                        # by self.cover_extension
                        cover_start = max(0,
                            cover_range[0] - self.cover_extension)
                        cover_end = min(len(sequence),
                            cover_range[1] + self.cover_extension)
                        # The endpoints of the cover give positions in just
                        # this sequence (chromosome), so adjust them (according
                        # to length_so_far) to give a unique integer position
                        # in the genome gnm
                        adjusted_cover = (cover_start + length_so_far,
                                          cover_end + length_so_far)
                        gnm_covers += [adjusted_cover]
                length_so_far += len(sequence)
            self.target_covers[i][j][rc] = gnm_covers

        probe.close_probe_finding_pool()

    def _compute_bp_covered_in_target_genomes(self):
        """Count number of bp covered by probes in each target genome.

        self._find_covers_in_target_genomes() must be called prior to this,
        so that self.target_covers can be accessed.

        This saves a dict, self.bp_covered, as follows:
        self.bp_covered[i][j][b] gives the number of bp covered by the
        probes in genome j of target genome grouping i (in the reverse
        complement of j if b is True, and in the provided sequence if b
        if False).
        """
        logger.info("Computing bases covered across target genomes")
        self.bp_covered = {}
        for i, j, gnm, rc in self._iter_target_genomes(True):
            if i not in self.bp_covered:
                self.bp_covered[i] = {}
            if j not in self.bp_covered[i]:
                self.bp_covered[i][j] = {False: None, True: None}
            covers = self.target_covers[i][j][rc]

            # Make an IntervalSet out of all covers to merge overlapping
            # ones and make it easy to count the number of bp covered
            covers_set = interval.IntervalSet(covers)
            self.bp_covered[i][j][rc] = len(covers_set)

    def _compute_average_coverage_in_target_genomes(self):
        """Calculate the average coverage/depth in each target genome.

        self._find_covers_in_target_genomes() must be called prior to this,
        so that self.target_covers can be accessed.

        This saves a dict, self.average_coverage, as follows:
        self.average_coverage[i][j][b] gives the average coverage/depth
        provided by the probes in genome j of target genome grouping i
        (in the reverse complement of j if b is True, and in the provided
        sequence if b is False).

        Specifically, the value is the average, taken across all bases,
        of the number of probes that hybridize to a region that includes
        a given base.
        """
        logger.info("Computing average coverage across target genomes")
        self.average_coverage = {}
        for i, j, gnm, rc in self._iter_target_genomes(True):
            if i not in self.average_coverage:
                self.average_coverage[i] = {}
            if j not in self.average_coverage[i]:
                self.average_coverage[i][j] = {False: None, True: None}
            covers = self.target_covers[i][j][rc]

            # Count the total number of bases covered by all the probe
            # hybridizations
            # (covers may include duplicates if two probes hybridize to
            # the same region, so it is important not to convert probes
            # to an IntervalSet or merge its intervals)
            total_covered = sum(c[1] - c[0] for c in covers)

            # Divide by the genome length to average across bases
            # (do this including ambiguous bases ('N') and not including them)
            avg_covg_over_all = float(total_covered) / gnm.size(False)
            avg_covg_over_unambig = float(total_covered) / gnm.size(True)
            self.average_coverage[i][j][rc] = (avg_covg_over_all,
                                               avg_covg_over_unambig)

    def _compute_sliding_coverage_in_target_genomes(self,
                                                    window_length,
                                                    window_stride):
        """Calculate the coverage over a sliding window in each target genome.

        self._find_covers_in_target_genomes() must be called prior to this,
        so that self.target_covers can be accessed.

        This saves a dict, self.sliding_coverage, as follows:
        self.sliding_coverage[i][j][b] gives a dict D provided by the probes
        in genome j of target genome grouping i (in the reverse complement
        of j if b is True, and in the provided sequence if b is False).
        D is a dict mapping a position p of a target genome to the average
        coverage/depth provided by the probes in a window centered at p.

        The value is the average, taken across the bases in the window, of
        the number of probes that hybridize to a region that includes each base
        in the window.

        This includes ambiguous bases in the sliding windows.
        """
        logger.info("Computing sliding coverage across target genomes")
        self.sliding_coverage = {}
        for i, j, gnm, rc in self._iter_target_genomes(True):
            if i not in self.sliding_coverage:
                self.sliding_coverage[i] = {}
            if j not in self.sliding_coverage[i]:
                self.sliding_coverage[i][j] = {False: None, True: None}
            covers = self.target_covers[i][j][rc]

            # Build a dict mapping endpoints to a list of False/True
            # values, where False indicates the endpoint is an end of
            # a probe and True indicates the endpoint is a start of a
            # probe
            endpoints = defaultdict(list)
            for c in covers:
                endpoints[c[0]].append(True)
                endpoints[c[1]].append(False)

            # Store, at each base in gnm, the number of probes that cover
            # a region including the base
            probe_counts = np.zeros(gnm.size(False), dtype='uint16')
            curr_count = 0
            prev_endpoint = None
            for endpoint in sorted(endpoints.keys()):
                if prev_endpoint != None:
                    # Update number of probes between prev_endpoint and
                    # endpoint
                    for pos in range(prev_endpoint, endpoint):
                        probe_counts[pos] = curr_count
                # Compute the net change in number of probe covers at
                # this endpoint, which will last until the next endpoint
                endpoint_types = endpoints[endpoint]
                net_change = sum([1 if t == True else -1
                                  for t in endpoint_types])
                curr_count += net_change
                prev_endpoint = endpoint
            # Don't bother updating counts between prev_endpoint and the end of
            # the genome because all probes should have been ended, so
            # curr_count should equal 0
            # (also, if probes cover the very end of the genome, prev_endpoint
            # might equal the end of the genome)

            # Slide over windows and fill in the dict gnm_sliding_coverage
            gnm_sliding_coverage = {}
            for window_start in np.arange(0, gnm.size(False), window_stride):
                window_end = window_start + window_length
                if window_end > gnm.size(False):
                    # This window stretches past the end of the genome, so
                    # make a new one that goes up to the very end
                    window_end = gnm.size(False)
                    window_start = window_end - window_length
                middle = window_start + (window_length / 2)
                window_average_count = np.average(
                    probe_counts[window_start:window_end])
                gnm_sliding_coverage[middle] = window_average_count

            self.sliding_coverage[i][j][rc] = gnm_sliding_coverage

    def run(self, window_length=50, window_stride=25):
        """Run all analysis methods.

        The methods called save their output to self.

        Args:
            window_length: number of bp in a window (for use in computing
                coverage over sliding windows)
            window_stride: number of bp by which to step (for use in computing
                coverage over sliding windows)
        """
        self._find_covers_in_target_genomes()
        self._compute_bp_covered_in_target_genomes()
        self._compute_average_coverage_in_target_genomes()
        self._compute_sliding_coverage_in_target_genomes(
            window_length, window_stride)

    def write_data_matrix_as_tsv(self, fn):
        """Write 2D array representing results as a TSV file.

        Args:
            fn: path to file to write to
        """
        # Make row headers
        data = [["Genome",
                 "Num bases covered",
                 "Frac bases covered",
                 "Frac bases covered over unambig",
                 "Average coverage/depth",
                 "Average coverage/depth over unambig"]]

        # Create a row for every genome, including reverse complements
        for i, j, gnm, rc in self._iter_target_genomes(True):
            col_header = "%s, genome %d" % (self.target_genomes_names[i], j)
            if rc:
                col_header += " (rc)"

            bp_covered = self.bp_covered[i][j][rc]
            frac_covered_all = float(bp_covered) / gnm.size(False)
            frac_covered_unambig = float(bp_covered) / gnm.size(True)

            avg_covg_all, avg_covg_unambig = self.average_coverage[i][j][rc]

            row = [col_header,
                   bp_covered,
                   frac_covered_all,
                   frac_covered_unambig,
                   avg_covg_all,
                   avg_covg_unambig]
            data += [row]

        # Write to fn as a TSV
        with open(fn, 'w') as f:
            for row in data:
                line = '\t'.join([str(entry) for entry in row])
                f.write(line + '\n')

    def _make_data_matrix_string(self):
        """Return 2D array representing results (as strings) to output.

        Returns:
            2D array, with row and column headers, containing data to
            output as a table
        """
        # Make row headers
        data = [["Genome",
                 "Num bases covered\n[over unambig]",
                 "Average coverage/depth\n[over unambig]"]]

        # Create a row for every genome, including reverse complements
        for i, j, gnm, rc in self._iter_target_genomes(True):
            col_header = "%s, genome %d" % (self.target_genomes_names[i], j)
            if rc:
                col_header += " (rc)"

            # Format bp covered
            bp_covered = self.bp_covered[i][j][rc]
            frac_covered_all = float(bp_covered) / gnm.size(False)
            frac_covered_unambig = float(bp_covered) / gnm.size(True)
            if frac_covered_all < 0.0001:
                prct_covered_all_str = "<0.01%"
            else:
                prct_covered_all_str = "{0:.2%}".format(frac_covered_all)
            if frac_covered_unambig < 0.0001:
                prct_covered_unambig_str = "<0.01%"
            else:
                prct_covered_unambig_str = "{0:.2%}".format(frac_covered_unambig)
            bp_covered_str = "%d (%s) [%s]" % (bp_covered,
                                                prct_covered_all_str,
                                                prct_covered_unambig_str)

            # Format average covered
            avg_covg_all, avg_covg_unambig = self.average_coverage[i][j][rc]
            if avg_covg_all < 0.01:
                avg_covg_all_str = "<0.01"
            else:
                avg_covg_all_str = "{0:.2f}".format(avg_covg_all)
            if avg_covg_unambig < 0.01:
                avg_covg_unambig_str = "<0.01"
            else:
                avg_covg_unambig_str = "{0:.2f}".format(avg_covg_unambig)
            avg_covg_str = "%s [%s]" % (avg_covg_all_str,
                                         avg_covg_unambig_str)

            row = [col_header,
                   bp_covered_str,
                   avg_covg_str]
            data += [row]

        return data

    def print_analysis(self):
        """Print the number of probes and a table of results of the analysis.
        """
        print("NUMBER OF PROBES: %d" % len(self.probes))
        print()
        print(pretty_print.table(self._make_data_matrix_string(),
                                 ["left", "right", "right"],
                                 header_underline=True))

    def write_sliding_window_coverage(self, fn):
        """Write coverage in sliding windows to a file.

        Args:
            fn: path to file to write to
        """
        with open(fn, 'w') as f:
            # Create an entry for every genome, including reverse complements
            for i, j, gnm, rc in self._iter_target_genomes(True):
                header = "%s, genome %d" % (self.target_genomes_names[i], j)
                if rc:
                    header += " (rc)"
                gnm_sliding_coverage = self.sliding_coverage[i][j][rc]
                for pos in sorted(gnm_sliding_coverage.keys()):
                    covg = gnm_sliding_coverage[pos]
                    line = '\t'.join([str(x) for x in [header, pos, covg]])
                    f.write(line + '\n')
