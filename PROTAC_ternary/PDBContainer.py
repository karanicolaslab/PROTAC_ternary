from rdkit import Chem
from rdkit.Chem import rdMolAlign

class PDBContainer:
    def __init__(self, file):
        self.filename = file
        self.atom_ids = {}
        self.structure = Chem.rdmolfiles.MolFromPDBFile(file,
                                                        sanitize=False,
                                                        removeHs=False,
                                                        flavor=0,
                                                        proximityBonding=False)
        self.set_atom_ids()
                    
    def set_atom_ids(self):
        self.atom_ids = {}
        for i, a in enumerate(self.structure.GetAtoms()):
            pdb_info = a.GetPDBResidueInfo()
            atom_name = pdb_info.GetName().strip()
            resi_name = pdb_info.GetResidueName().strip()
            key = resi_name + " " + atom_name
            self.atom_ids[key] = i
        
    def align(self, reference, atom_map, max_iters=50):
        rmsd = Chem.rdMolAlign.AlignMol(self.structure, 
                                        reference.structure,
                                        atomMap=atom_map,
                                        maxIters=max_iters)
        return rmsd
    
    def save(self, file):
        Chem.rdmolfiles.MolToPDBFile(self.structure, file, flavor=2)
        
    def delete_hetHs(self):
        hetHs = []
        
        for i, a in enumerate(self.structure.GetAtoms()):
            pdb_info = a.GetPDBResidueInfo()

            is_hetero = pdb_info.GetIsHeteroAtom()
            atom_name = pdb_info.GetName().strip()
            resi_name = pdb_info.GetResidueName().strip()
            
            if is_hetero and a.GetSymbol() == "H":
                key = resi_name + " " + atom_name
                hetHs.append(key)    
                        
        self.delete(hetHs)
        self.set_atom_ids()
        
    def delete(self, atoms):
        structure = Chem.rdchem.EditableMol(self.structure)
        ids = set([self.atom_ids[a] for a in atoms if a in self.atom_ids])
        
        for idx in sorted(ids, reverse=True):
            structure.RemoveAtom(idx)
            
        for a in atoms:
            if a in self.atom_ids:
                del self.atom_ids[a]
            
        self.structure = structure.GetMol()
        
        self.set_atom_ids()
        
    def merge(self,compound):
        
        structure = Chem.rdchem.EditableMol(self.structure)        
        atoms_ids = {}
        coords = {}
        
        for i, a in enumerate(compound.structure.GetAtoms()):
            idx = structure.AddAtom(a)
            atoms_ids[i] = idx
            coords[idx] = compound.structure.GetConformer(0).GetAtomPosition(i)
                    
        for bond in compound.structure.GetBonds():
            bgn = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            structure.AddBond(bgn, end, bond.GetBondType())             
                    
        self.structure = structure.GetMol()
        
        for idx, coord in coords.items():
            self.structure.GetConformer(0).SetAtomPosition(idx, coord)
        
        self.set_atom_ids()
            
    def rename(self, old_resn, new_resn, new_chain=None, new_resi=None):
        for i, a in enumerate(self.structure.GetAtoms()):
            pdb_info = a.GetPDBResidueInfo()
            
            if pdb_info.GetResidueName().strip() == old_resn:
                pdb_info.SetResidueName(new_resn)
                if new_chain:
                    pdb_info.SetChainId(new_chain)
                if new_resi:
                    pdb_info.SetResidueNumber(new_resi)
        
        self.set_atom_ids()
