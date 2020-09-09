import os
import unittest

import numpy as np

from pdb_superimposer.superimposition import ChainSuperimposer, Helper


class TestChainSuperimposer(unittest.TestCase):
    def setUp(self):
        parent = os.path.dirname(__file__)
        # ref
        self.ref = Helper.protein_from_file("2VTA",
                                            os.path.join(parent, "testdata", "2VTA.pdb"))
        self.chain_ref = [c for c in self.ref[0]][0]
        # other
        pdb = "6YLK"
        self.a = Helper.protein_from_file(pdb,
                                          os.path.join(parent, "testdata", "{pdb}.pdb"))
        self.coords_a = [atm.coord for atm in self.a.get_atoms()]
        self.chain_a = [c for c in self.a[0]][0]

        self.cs = ChainSuperimposer(reference=self.chain_ref, other=self.chain_a, other_struc=self.a)

    def testChainSuperimposerConstruction(self):
        # artifact removal?
        self.assertIn(-1, self.cs.other_seq)
        self.assertNotIn(-1, self.cs.selected_index)

        # mutation
        self.assertIn(177, self.cs.other_seq)
        self.assertNotIn(177, self.cs.selected_index)

        # missing atoms
        self.assertIn(178, self.cs.other_seq)
        self.assertNotIn(178, self.cs.selected_index)

        # not modelled res
        self.assertNotIn(45, self.cs.other_seq)

    def testsuperimpose(self):
        self.cs.superimpose()

        # check alignment is sensible
        self.assertLess(self.cs.rms, 5)

        # did all atoms move?
        moved = [atm.coord for atm in self.a.get_atoms()]
        self.assertTrue(all([np.linalg.norm(a - b) != 0 for a, b in zip(self.coords_a, moved)]))


if __name__ == '__main__':
    unittest.main()
