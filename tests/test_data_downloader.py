import os
import unittest
import tempfile

from Bio.PDB.Structure import Structure

from pdb_superimposer.data_downloader import pdb_search_query, ftp_download
from pdb_superimposer.superimposition import Helper


class TestNoClassMethods(unittest.TestCase):
    def test_pdb_search_query(self):
        parent = os.path.dirname(__file__)
        with open(os.path.join(parent, "testdata", "test_search.json"), "r") as r:
            query = r.read()

        result = pdb_search_query(query)

        # Successful searches will:
        #     1: return a json not None
        #     2: have a result set
        self.assertIn("result_set", result)

    def testftp_download(self):
        pdb = "2VTA"
        out = tempfile.mkdtemp()
        response = ftp_download(pdb, out)

        # does the file exist
        fpath = os.path.join(out, f"{pdb}.pdb")
        self.assertTrue(os.path.exists(fpath))

        # can it be read?
        a = Helper.protein_from_file(pdb, fpath)
        self.assertIsInstance(a, Structure)


if __name__ == '__main__':
    unittest.main()
