import unittest
from ../src/data_downloader import pdb_search_query, ftp_download
import os
from src.superimposition import Helper
from Bio.PDB.Structure import Structure


class TestNoClassMethods(unittest.TestCase):
    def test_pdb_search_query(self):
        with open("testdata/test_search.json", "r") as r:
            query = r.read()

        result = pdb_search_query(query)

        # Successful searches will:
        #     1: return a json not None
        #     2: have a result set
        self.assertIn("result_set", result)

    def testftp_download(self):

        pdb = "2VTA"
        out = "testdata"
        response = ftp_download(pdb, out)

        # does the file exist
        fpath = os.path.join(out, f"{pdb}.pdb")
        self.assertTrue(os.path.exists(fpath))

        # can it be read?
        a = Helper.protein_from_file(pdb, fpath)
        self.assertIsInstance(a, Structure)


if __name__ == '__main__':
    unittest.main()
