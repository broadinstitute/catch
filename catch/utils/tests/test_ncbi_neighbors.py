"""Tests for ncbi_neighbors module.
"""

import unittest

from catch.utils import ncbi_neighbors as nn

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class TestURLConstruction(unittest.TestCase):
    """Tests constructing URLs.
    """

    def test_ncbi_neighbors_url(self):
        url = nn.ncbi_neighbors_url(123)
        expected_url = ('https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup'
            '.cgi?taxid=123&cmd=download2')
        self.assertEqual(url, expected_url)

    def test_ncbi_fasta_download_url(self):
        url = nn.ncbi_fasta_download_url(['A123', 'A456', 'B789'])
        expected_url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch'
            '.fcgi?id=A123,A456,B789&db=nuccore&rettype=fasta&retmode=text')
        self.assertEqual(url, expected_url)


class TestConstructNeighbors(unittest.TestCase):
    """Tests the construct_neighbors() function.

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
