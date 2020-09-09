"""

Functionality to superimpose protein chains

"""
import numpy as np
from Bio import PDB
from scipy.spatial import distance


class ChainSuperimposer:
    def __init__(self, reference, other, other_struc):
        """
        Class to managed the data and methods for protein chain superimposition.
        Initialise with a reference and query chain.

        :param reference: Superimposition template
        :param other: Chain to be transformed
        :type reference: `Bio.PDB.Chain.Chain`
        :type other: `Bio.PDB.Chain.Chain`
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
            - Find intersecting sequence indices
            - Remove index if the residue type at a given position is not the same
            - Remove index if either the reference or other residue is incomplete

        :return: A set of sequence indices
        :rtype: set
        """
        sele = set(self.reference_seq.keys()).intersection(set(self.other_seq.keys()))

        rm = {i for i in sele if
              any((self.reference_seq[i].resname != self.other_seq[i].resname,
                   self.is_incomplete(self.reference_seq[i]),
                   self.is_incomplete(self.other_seq[i])
                   ))
              }
        for r in rm:
            sele.remove(r)
        return sele

    def superimpose(self, hetid=None, binding_site=False, within=8.0):
        """
        Superimpose other to reference, if hetid is and supplied and binding_site set to True
        then only the binding site residues will be used for the superimposition.

        :param hetid: Ligand identifier
        :param binding_site: flag
        :param within: binding site cutoff distance
        :type hetid: str
        :type binding_site: bool
        :type within: float
        """
        if binding_site and hetid:
            # detect binding site
            bs = self.binding_site(chain=self.other, hetid=hetid, within=within)
            self.selected_index = self.selected_index.intersection(bs)

        super_imposer = PDB.Superimposer()

        # set the active atoms
        reference_atms = [atm for resi in self.reference_seq.values() for atm in resi
                          if resi.get_id()[1] in self.selected_index]
        other_atms = [atm for resi in self.other_seq.values() for atm in resi
                      if resi.get_id()[1] in self.selected_index]
        super_imposer.set_atoms(reference_atms, other_atms)

        # apply the transformation matrix to the whole chain
        super_imposer.apply(self.other_struc.get_atoms())
        self.rms = super_imposer.rms

    @staticmethod
    def is_incomplete(res):
        """
        Determine whether the residue is incomplete (i.e. whether there are missing residues)

        :param res: A PDB residue
        :type: `Bio.PDB.Residue.Residue`
        :return: A boolean expressing whether the residue is incomplete
        :rtype: bool
        """
        aa_atm_dic = {'MET': 8, 'GLU': 9, 'ASN': 8, 'PHE': 11, 'GLN': 9, 'LYS': 9,
                      'VAL': 7, 'ILE': 8, 'GLY': 4, 'THR': 7, 'TYR': 12, 'ALA': 5,
                      'ARG': 11, 'LEU': 8, 'SER': 6, 'HIS': 10, 'PRO': 7, 'ASP': 8,
                      'CYS': 6, 'TRP': 14}
        return len(res) != aa_atm_dic[res.resname]

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
