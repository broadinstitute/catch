"""Tests for ncbi_neighbors module.
"""

import logging
import unittest

from catch.utils import ncbi_neighbors as nn
from catch.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'

TAXIDS = {'ZIKV': 64320, 'IAV': 11320}
ZIKV_ACCS = ['KY785462', 'KY785463']


class TestURLConstruction(unittest.TestCase):
    """Tests constructing URLs.
    """

    def encoded_url_check(self, url, expected_start, expected_fields):
        url_split = url.split('?')
        self.assertEqual(url_split[0], expected_start)
        fields = url_split[1].split('&')
        self.assertCountEqual(fields, expected_fields)

    def test_ncbi_neighbors_url(self):
        url = nn.ncbi_neighbors_url(123)
        self.encoded_url_check(url,
            'https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi',
            ['taxid=123', 'cmd=download2'])

    def test_ncbi_fasta_download_url(self):
        url = nn.ncbi_fasta_download_url(['A123', 'A456', 'B789'])
        self.encoded_url_check(url,
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
            ['id=A123,A456,B789', 'db=nuccore', 'rettype=fasta',
                'retmode=text'])


class TestConstructNeighborsUnit(unittest.TestCase):
    """Unit tests for the construct_neighbors() function.

    The function construct_neighbors() calls fetch_neighbors_table(),
    which makes a request to NCBI. To avoid the request, thsi overrides
    fetch_neighbors_table() to return a known neighbors table.
    """

    def setUp(self):
        self.expected_table_contents = \
            ("## Comment line 1\n"
             "## Comment line 2\n"
             "## Columns:\tRepresentative\tNeighbor\tHost\tSelected lineage\tTaxonomy name\tSegment name\n"
             "NC_0123\tKY456\tvertebrate,human\tFamilyA,GenusA,SpeciesA\tSpeciesA\tsegment \n"
             "NC_0456\tAB123\tinvertebrate\tFamilyB,GenusB,SpeciesB\tSpeciesB\tsegment 1\n"
             "NC_0456\tAB456\tinvertebrate\tFamilyB,GenusB,SpeciesB\tSpeciesB\tsegment 2\n")
        self.expected_neighbors = [
            nn.Neighbor('KY456', 'NC_0123', ['vertebrate', 'human'],
                ('FamilyA', 'GenusA', 'SpeciesA'), 'SpeciesA', ''),
            nn.Neighbor('AB123', 'NC_0456', ['invertebrate'],
                ('FamilyB', 'GenusB', 'SpeciesB'), 'SpeciesB', '1'),
            nn.Neighbor('AB456', 'NC_0456', ['invertebrate'],
                ('FamilyB', 'GenusB', 'SpeciesB'), 'SpeciesB', '2')
        ]

        # Override nn.fetch_neighbors_table() to return expected_table_lines,
        # but keep the real function
        self.fetch_neighbors_table_real = nn.fetch_neighbors_table
        nn.fetch_neighbors_table = lambda taxid: self.expected_table_contents.split('\n')

    def test_construct_neighbors(self):
        neighbors = nn.construct_neighbors(123)
        self.assertEqual(neighbors, self.expected_neighbors)

    def tearDown(self):
        # Reset nn.fetch_neighbors_table()
        nn.fetch_neighbors_table = self.fetch_neighbors_table_real


class TestConstructNeighborsIntegration(unittest.TestCase):
    """Integration tests for the construct_neighbors() function.

    This tests for non-influenza and influenza species. Note that it
    makes remote calls to NCBI's API; if that is not working, these
    tests will fail.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_construct_neighbors_for_zika(self):
        # Download Zika virus neighbors
        neighbors = nn.construct_neighbors(TAXIDS['ZIKV'])

        # Check that there are at least 100 neighbors (there should
        # be many more)
        self.assertGreaterEqual(len(neighbors), 100)

        # Check that most sequences have 'Zika virus' in their lineage (all
        # should but there may be outliers so only check
        # if at least half are)
        num_with_zika = 0
        for n in neighbors:
            if 'Zika virus' in n.lineage:
                num_with_zika += 1
        self.assertGreaterEqual(num_with_zika, len(neighbors) / 2)

    def test_construct_neighbors_for_influenza(self):
        # Download Influenza A virus neighbors
        neighbors = nn.construct_influenza_genome_neighbors(TAXIDS['IAV'])

        # Check that there are at least 100 neighbors (there should
        # be many more)
        self.assertGreaterEqual(len(neighbors), 100)

        # Check that most sequences have 'Influenza A virus' in their lineage
        # (all should but there may be outliers so only check
        # if at least half are)
        num_with_influenza = 0
        for n in neighbors:
            if 'Influenza A virus' in n.lineage:
                num_with_influenza += 1
        self.assertGreaterEqual(num_with_influenza, len(neighbors) / 2)

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestFetchFastasFromAccessionsIntegration(unittest.TestCase):
    """Integration tests for the fetch_fastas() function.

    Note that this makes remote calls to NCBI's API; if that is not working,
    these tests will fail.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_fetch_fastas(self):
        # Download Zika virus accessions
        fasta_tf = nn.fetch_fastas(ZIKV_ACCS)

        # Read the fasta
        seqs = seq_io.read_fasta(fasta_tf.name)

        # Verify that the right number of sequences were fetched
        self.assertEqual(len(seqs), len(ZIKV_ACCS))

        # Verify that each accession appears in a sequence header (it may
        # not match exactly because the sequence header is
        # [accession].[verison], but the accession should be a substring)
        for acc in ZIKV_ACCS:
            num_with_acc = sum(1 for seq_header in seqs.keys() if acc in
                    seq_header)
            self.assertEqual(num_with_acc, 1)

        fasta_tf.close()

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)


class TestConstructFastaFromTaxIDIntegration(unittest.TestCase):
    """Integration tests for the construct_fasta_for_taxid() function.

    Note that this makes remote calls to NCBI's API; if that is not working,
    these tests will fail.
    """

    def setUp(self):
        # Disable logging
        logging.disable(logging.INFO)

    def test_construct_fasta_for_taxid(self):
        # Download Zika virus sequences
        fasta_tf = nn.construct_fasta_for_taxid(TAXIDS['ZIKV'])

        # Read the fasta
        seqs = seq_io.read_fasta(fasta_tf.name)

        # Check that there are at least 100 sequences (there should be
        # many more)
        self.assertGreaterEqual(len(seqs), 100)

        fasta_tf.close()

    def tearDown(self):
        # Re-enable logging
        logging.disable(logging.NOTSET)
