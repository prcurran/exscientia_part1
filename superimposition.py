"""
Used for reference:

https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_superposition/crystallography_pdb_alignment.py

"""
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Superimposer, PPBuilder, Polypeptide
from scipy.spatial import distance
import numpy as np


class ModifiedResidue(PDB.Residue.Residue):
    def __init__(self, id, resname, segid):
        super().__init__(id, resname, segid)

    def __str__(self):
        return Polypeptide.three_to_one(self.resname)


PDB.Residue.Residue = ModifiedResidue


class ChainSuperimposer:
    def __init__(self, reference, other):
        """
        Class to managed the data and methods for protein chain superimposition. Initialise with
        a reference and query chain.

        :param reference: Superimposition template
        :param other: Chain to be transformed
        :type reference: `Bio.PDB.Chain`
        :type other: `Bio.PDB.Chain`
        """
        self.spacer = "-"
        self.reference = reference
        assert (reference.get_level() == 'C')
        self.reference_seq = self.sequence_with_spacers(self.reference)

        self.other = other
        assert (other.get_level() == 'C')
        self.other_seq = self.sequence_with_spacers(self.other)

        # the sequence with spacer must be equal (and match the UniProt sequence)
        assert(len(self.reference_seq) == len(self.other_seq))

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
            spacer_indices = [i for i, e in enumerate(self.other_seq) if e == self.spacer]
            # insert missing residue spacers
            for i in spacer_indices:
                bs.insert(i, self.spacer)

            # use the residue if:
            #     - same as reference residue
            #     - not a missing residue
            #     - class as within the binding site
            use = [all((str(res_a) == str(res_b), str(res_a) != "-", bs_res))
                   for res_a, res_b, bs_res in zip(self.reference_seq, self.other_seq, bs)]
        else:
            use = [all((str(res_a) == str(res_b), str(res_a) != "-"))
                   for res_a, res_b in zip(self.reference_seq, self.other_seq)]

        super_imposer = Superimposer()
        # set the active atoms
        # TODO: Could have some atom-based inclusion crtiera here e.g. backbone? etc.
        reference_atms = [atm for resi, flag in zip(self.reference_seq, use) for atm in resi if flag]
        other_atms = [atm for resi, flag in zip(self.other_seq, use) for atm in resi if flag]
        super_imposer.set_atoms(reference_atms, other_atms)
        # apply the transformation matrix to the whole chain
        super_imposer.apply(self.other.get_atoms())

        self.rms = super_imposer.rms

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
        # TODO: A residue is included if any residue atom is within the cutoff.
        #  (whole_residue=True). Would be nice to also provide the option to only include atoms
        #  within the cutoff. (whole_residue=False)
        # array of ligand coords
        ligand = [resi for resi in chain if hetid in resi.get_id()[0]][0]
        ligand_atms = np.array([atm.get_coord() for atm in ligand])
        # list of array's containing residue coords
        resi_atms = [np.array([atm.get_coord() for atm in resi]) for resi in chain if resi.get_id()[0] == " "]
        # loop through the residues,
        # calculate and all-by-all atom distances
        # if any atom is within the cutoff distance, return True
        return [np.amin(distance.cdist(r, ligand_atms)) <= within for r in resi_atms]

    @staticmethod
    def sequence_with_spacers(chain, chara="-"):
        """
        Method to add spacer character for missing residues

        :param chain: A protein chain
        :param chara: The character used to signify a missing residue
        :type chain: `Bio.PDB.Chain`
        :type chara: str
        :return: PDB sequence as a list
        :rtype: list
        """
        sequence = []
        prev_pos = 0
        for i, residue in enumerate(chain):
            # Only look at amino acids (residues have no identifier, HET have HETID)
            if residue.get_id()[0] == " ":
                # Find missing residues
                current_pos = int(residue.get_id()[1])
                diff = current_pos - (prev_pos + 1)
                if diff > 0:
                    sequence.extend([chara] * diff)

                sequence.append(residue)
                prev_pos = current_pos
        return sequence


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
    return PDBParser().get_structure(name, fpath)


def protein_to_file(structure, fpath):
    """
    Writes a PDB to file.

    :param structure: a protein structure
    :param fpath: path for output pdb file
    :type structure: `Bio.PDB.Structure.Structure`
    :type fpath: str
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(fpath)


if __name__ == "__main__":

    a = protein_from_file("1aq1", "testdata/1aq1.pdb")
    b = protein_from_file("2vta", "testdata/2vta.pdb")

    chain_a = [c for c in a[0]][0]
    chain_b = [c for c in b[0]][0]

    cs = ChainSuperimposer(reference=chain_b, other=chain_a)
    cs.superimpose("STU", binding_site=True)

    protein_to_file(a, "testdata/aligned_1aq1.pdb")

