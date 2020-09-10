"""

A Module to superimpose two protein chains.

"""

import numpy as np
from Bio import PDB
from scipy.spatial import distance


class ChainSuperimposer:
    """ A class to organise the superimposition of two chains. First the atoms involved in the operation
    are selected. Then :class:`Bio.PDB.Superimposer.Superimposer` is employed to do the superimposition.
    """

    def __init__(self, reference, other, other_struc):
        """

        Initialise the class.

        :param reference: The chain used as a template in the superimposition.
        :type reference: :class:`Bio.PDB.Chain.Chain`
        :param other: The chain to be transformed.
        :type other: :class:`Bio.PDB.Chain.Chain`
        :param other_struc: The structure object of the other chain. This enables the transformation matrix to be applied across the whole structure.
        :type other_struc: :class:`Bio.PDB.Structure.Structure`

        """
        # Input data
        self.reference = reference
        self.other = other
        self.other_struc = other_struc  # Needed to translate the whole structure

        # ensure only single chains are supplied
        assert (reference.get_level() == 'C')
        assert (other.get_level() == 'C')

        # Setup a dictionary matching uniprot res index with res. Discard ligands and water.
        # res.get_id()[0] ---> hetid flag
        # res.get_id()[1] ---> res index
        self.reference_seq = {res.get_id()[1]: res for res in self.reference if res.get_id()[0] == " "}
        self.other_seq = {res.get_id()[1]: res for res in self.other if res.get_id()[0] == " "}
        self.selected_index = self._residue_selection()

    def _residue_selection(self):
        """
        A class method to do the essential residue filtering using the UniProt sequence positions.

            - Find intersecting sequence indices between the reference and other.
            - Remove index if at a given UniProt index the residue type in reference and other is not the same
            - Remove index if at a given UniProt index either the reference or other residue is incomplete

        :return: A set of sequence indices
        :rtype: set
        """
        sele = set(self.reference_seq.keys()).intersection(set(self.other_seq.keys()))

        # for testing
        # for i in sele:
            # print(self.reference_seq[i].resname != self.other_seq[i].resname)
            # print(self.is_residue_incomplete(self.reference_seq[i]))
            # print(self.is_residue_incomplete(self.other_seq[i]))

        rm = {i for i in sele if
              any((self.reference_seq[i].resname != self.other_seq[i].resname,
                   self.is_residue_incomplete(self.reference_seq[i]),
                   self.is_residue_incomplete(self.other_seq[i])
                   ))
              }

        for r in rm:
            sele.remove(r)
        return sele

    def superimpose(self, hetid=None, within=8.0):
        """

        The superimposition method. If optional `hetid` is supplied, only binding site
        residues will used for the superimposition.

        :param hetid: Ligand identifier
        :type hetid: str
        :param binding_site: flag
        :type binding_site: bool
        :param within: binding site cutoff distance
        :type within: float

        .. code-block:: python

            # Example useage
            from pdb_superimposer import ChainSuperimposer, Helper

            ref = Helper.protein_from_file("2VTA", "2VTA.pdb")
            ref_chain = [c for c in ref[0]][0]

            other = Helper.protein_from_file("6YLK", "6YLK.pdb")
            ref_chain = [c for c in other[0]][0]

            cs = ChainSuperimposer(reference=ref_chain, other=other_chain, other_struc=other)
            cs.superimpose()


        """
        if hetid:
            # detect binding site
            bs = self.binding_site(chain=self.other, hetid=hetid, within=within)
            self.selected_index = self.selected_index.intersection(bs)

        super_imposer = PDB.Superimposer()

        # set the active atoms
        reference_atms = [atm for resi in self.reference_seq.values() for atm in resi
                          if resi.get_id()[1] in self.selected_index and atm.element != "H"]
        other_atms = [atm for resi in self.other_seq.values() for atm in resi
                      if resi.get_id()[1] in self.selected_index and atm.element != "H"]
        super_imposer.set_atoms(reference_atms, other_atms)

        # apply the transformation matrix to the whole chain
        super_imposer.apply(self.other_struc.get_atoms())
        self.rms = super_imposer.rms

    @staticmethod
    def is_residue_incomplete(res):
        """
        Determine whether the residue is incomplete (i.e. whether there are missing residue atoms)

        :param res: A PDB residue
        :type: `Bio.PDB.Residue.Residue`
        :return: A boolean expressing whether the residue is incomplete
        :rtype: bool

        >>> from pdb_superimposer import ChainSuperimposer, Helper
        >>> ref = Helper.protein_from_file("2VTA", "2VTA.pdb")
        >>> ref_chain = [c for c in ref[0]][0]
        >>> first_res = [r for r in ref_chain][0]
        >>> bs = ChainSuperimposer.is_residue_incomplete(first_res)
        False

        """
        aa_atm_dic = {'MET': 8, 'GLU': 9, 'ASN': 8, 'PHE': 11, 'GLN': 9, 'LYS': 9,
                      'VAL': 7, 'ILE': 8, 'GLY': 4, 'THR': 7, 'TYR': 12, 'ALA': 5,
                      'ARG': 11, 'LEU': 8, 'SER': 6, 'HIS': 10, 'PRO': 7, 'ASP': 8,
                      'CYS': 6, 'TRP': 14}

        return len([a for a in res if a.element != "H"]) != aa_atm_dic[res.resname]

    @staticmethod
    def binding_site(chain, hetid, within=8.0):
        """
        Detect binding site residues for a supplied hetid

        :param chain: A protein chain
        :param hetid: A ligand identifier
        :param within: A cutoff distance
        :type chain: `Bio.PDB.Chain`
        :type hetid: str
        :type within: float
        :return: A boolean list denoting inclusion in the superimposition operation
        :rtype: list

        >>> from pdb_superimposer import ChainSuperimposer, Helper
        >>> ref = Helper.protein_from_file("2VTA", "2VTA.pdb")
        >>> ref_chain = [c for c in ref[0]][0]
        >>> ans = ChainSuperimposer.binding_site(ref_chain, "LZ1", 10.0)
        {1, 2, 3, ...}

        """
        # array of ligand coords
        ligand = [resi for resi in chain if hetid in resi.get_id()[0]][0]
        ligand_atms = np.array([atm.get_coord() for atm in ligand])
        # list of array's containing residue coords
        resi_atms = [np.array([atm.get_coord() for atm in resi]) for resi in chain if resi.get_id()[0] == " "]
        # loop through the residues,
        # calculate and all-by-all atom distances
        # if any atom is within the cutoff distance, return True
        return {r.get_id()[1] for r in resi_atms if np.amin(distance.cdist(r, ligand_atms)) <= within}


class Helper:
    @staticmethod
    def protein_from_file(name, fpath):
        """
        Reads a PDB file.

        :param name: name of object
        :param fpath: path to pdb file
        :type name: str
        :type fpath: str
        :return: a protein structure
        :rtype: `Bio.PDB.Structure.Structure`
        """
        return PDB.PDBParser().get_structure(name, fpath)

    @staticmethod
    def protein_to_file(structure, fpath):
        """
        Writes a PDB to file.

        :param structure: a protein structure
        :param fpath: path for output pdb file
        :type structure: `Bio.PDB.Structure.Structure`
        :type fpath: str
        """
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(fpath)
