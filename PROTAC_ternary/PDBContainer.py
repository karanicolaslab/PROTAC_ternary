from rdkit import Chem
from rdkit.Chem import rdMolAlign


class PDBContainer:
    """Container for storage PDB information

    Several operations with input structure as add/remove atoms, align against
    reference structure and save to file using RDKit package.

    Attributes:
        atom_ids (dict): Relationship between PDB <residue name>/<atom name>
                         and atom ids of RDKit molecule representation

        structure (rdkit.Chem.rdchem.Mol): RDKit container for structure
    """

    def __init__(self, file):
        """Inititializing of container

        Note: we disabled sanitizing step as well as removing of hydrogen atoms
        to escape structural errors during reading PDB file. Also assumptions
        about bonding of closely located also disabled.

        Args:
            file (str): File name of input structure
        """
        self.atom_ids = {}
        self.structure = Chem.rdmolfiles.MolFromPDBFile(file,
                                                        sanitize=False,
                                                        removeHs=False,
                                                        flavor=0,
                                                        proximityBonding=False)
        self.update_atom_cache()

    def update_atom_cache(self):
        """
        Update state of atom_ids dictionary
        """
        self.atom_ids = {}

        for i, atom in enumerate(self.structure.GetAtoms()):
            pdb_info = atom.GetPDBResidueInfo()

            if pdb_info is not None:
                atom_name = pdb_info.GetName().strip()
                resi_name = pdb_info.GetResidueName().strip()
                key = f"{resi_name} {atom_name}"
                self.atom_ids[key] = i

    def add_hetatm(self,
                   symbol,
                   coords=(0, 0, 0),
                   atom_name="D1",
                   res_name="LG1",
                   res_numb=1):
        """Add hetatm into structure

        Args:
            symbol (str): Any chemical element
            coords (tuple, optional): 3D coordinate of added atom
            atom_name (str, optional): PDB name of added atom
            res_name (str, optional): PDB residue name of added atom
            res_numb (int, optional): PDB residue number of added atom

        Returns:
            str: PDB <residue name>/<atom name> of added atom
        """
        # prepare pdb information
        pdb_info = Chem.AtomPDBResidueInfo()
        pdb_info.SetName(atom_name)
        pdb_info.SetResidueName(res_name)
        pdb_info.SetResidueNumber(res_numb)
        pdb_info.SetIsHeteroAtom(True)

        # add pdb information into atom
        atom = Chem.rdchem.Atom(symbol)
        atom.SetMonomerInfo(pdb_info)

        # add atom to structure
        structure = Chem.rdchem.EditableMol(self.structure)
        idx = structure.AddAtom(atom)
        self.structure = structure.GetMol()
        self.structure.GetConformer(0).SetAtomPosition(idx, coords)

        a = self.structure.GetAtomWithIdx(idx)
        pdb_info = a.GetPDBResidueInfo()
        atom_name = pdb_info.GetName().strip()
        resi_name = pdb_info.GetResidueName().strip()
        key = f"{resi_name} {atom_name}"

        self.update_atom_cache()

        return key

    def get_coordinates(self, atom_labels):
        """Return coordinates of given atoms

        Args:
            atom_labels (list): PDB labels of interested atoms

        Returns:
            list: coordinates of given atoms as nested list of X,Y,Z
            coordinates
        """
        x_coords, y_coords, z_coords = [], [], []

        for label in atom_labels:
            atom_id = self.atom_ids[label]
            coord = self.structure.GetConformer(0).GetAtomPosition(atom_id)
            x_coords.append(coord.x)
            y_coords.append(coord.y)
            z_coords.append(coord.z)

        return x_coords, y_coords, z_coords

    def align(self, reference, targ_atoms, ref_atoms, max_iters=50):
        """Align structure against reference

        Args:
            reference (PDBContainer): Container of reference structure
            targ_atoms (list): PDB atom labels of target using for alignment
            ref_atoms (list): PDB atom labels of reference using for alignment
            max_iters (int, optional): Max number of iterations for alignment

        Returns:
            float: RMSD between aligned targ_atoms and ref_atoms
        """
        targ_dummy_atom_ids = []
        targ_atom_ids = []

        for i, atom in enumerate(targ_atoms):
            if len(atom) > 1:
                coords = self.get_coordinates(atom)
                x_avg = sum(coords[0]) / len(coords[0])
                y_avg = sum(coords[1]) / len(coords[1])
                z_avg = sum(coords[2]) / len(coords[2])
                label = self.add_hetatm("H",
                                        coords=(x_avg, y_avg, z_avg),
                                        atom_name=f"D{i}",
                                        res_name="LG1",
                                        res_numb=1)
                atom_id = self.atom_ids[label]
                targ_dummy_atom_ids.append(label)
            else:
                atom_id = self.atom_ids[atom[0]]

            targ_atom_ids.append(atom_id)

        ref_dummy_atom_ids = []
        ref_atom_ids = []

        for i, atom in enumerate(ref_atoms):
            if len(atom) > 1:
                coords = reference.get_coordinates(atom)
                x_avg = sum(coords[0]) / len(coords[0])
                y_avg = sum(coords[1]) / len(coords[1])
                z_avg = sum(coords[2]) / len(coords[2])
                label = reference.add_hetatm("H",
                                             coords=(x_avg, y_avg, z_avg),
                                             atom_name=f"D{i}",
                                             res_name="LG1",
                                             res_numb=1)
                atom_id = reference.atom_ids[label]
                ref_dummy_atom_ids.append(label)
            else:
                atom_id = reference.atom_ids[atom[0]]

            ref_atom_ids.append(atom_id)

        atom_map = list(zip(targ_atom_ids, ref_atom_ids))

        rmsd = Chem.rdMolAlign.AlignMol(self.structure,
                                        reference.structure,
                                        atomMap=atom_map,
                                        maxIters=max_iters)

        self.delete_atoms(targ_dummy_atom_ids)
        reference.delete_atoms(ref_dummy_atom_ids)

        return rmsd

    def delete_atoms(self, atoms):
        """Delete set of atoms

        Args:
            atoms (list): PDB labels of deleting atoms
        """
        structure = Chem.rdchem.EditableMol(self.structure)

        ids = set([self.atom_ids[a] for a in atoms if a in self.atom_ids])

        for idx in sorted(ids, reverse=True):
            structure.RemoveAtom(idx)

        for a in atoms:
            if a in self.atom_ids:
                del self.atom_ids[a]

        self.structure = structure.GetMol()
        self.update_atom_cache()

    def merge(self, compound):
        """Copy atoms and bonds of some structure to self

        Args:
            compound (PDBContainer): Structure that will be pasted into given
            structure
        """
        structure = Chem.rdchem.EditableMol(self.structure)
        atoms_ids = {}
        coords = {}

        for i, a in enumerate(compound.structure.GetAtoms()):
            idx = structure.AddAtom(a)
            atoms_ids[i] = idx
            coords[idx] = compound.structure.GetConformer(0).GetAtomPosition(i)

        for bond in compound.structure.GetBonds():
            bgn = atoms_ids[bond.GetBeginAtomIdx()]
            end = atoms_ids[bond.GetEndAtomIdx()]
            print(bgn, end)
            structure.AddBond(bgn, end, bond.GetBondType())

        self.structure = structure.GetMol()

        for idx, coord in coords.items():
            self.structure.GetConformer(0).SetAtomPosition(idx, coord)

        self.update_atom_cache()

    def save(self, file):
        """Summary

        Args:
            file (TYPE): Description
        """
        Chem.rdmolfiles.MolToPDBFile(self.structure, file, flavor=2)

    def delete_hetHs(self):
        """Remove HETATM hydrogens from the structure
        """
        hetHs = []

        for i, a in enumerate(self.structure.GetAtoms()):
            pdb_info = a.GetPDBResidueInfo()

            is_hetero = pdb_info.GetIsHeteroAtom()
            atom_name = pdb_info.GetName().strip()
            resi_name = pdb_info.GetResidueName().strip()

            if is_hetero and a.GetSymbol() == "H":
                key = resi_name + " " + atom_name
                hetHs.append(key)

        self.delete_atoms(hetHs)
        self.update_atom_cache()

    def rename(self, old_resn, new_resn, new_chain=None, new_resi=None):
        """Replace residue name and number as well as chain label

        Args:
            old_resn (str): Old residue name
            new_resn (str): New residue name
            new_chain (str, optional): New chain name
            new_resi (str, optional): New residue number
        """
        for i, a in enumerate(self.structure.GetAtoms()):
            pdb_info = a.GetPDBResidueInfo()

            if pdb_info.GetResidueName().strip() == old_resn:
                pdb_info.SetResidueName(new_resn)
                if new_chain:
                    pdb_info.SetChainId(new_chain)
                if new_resi:
                    pdb_info.SetResidueNumber(new_resi)

        self.update_atom_cache()
