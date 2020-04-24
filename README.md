# PROTAC_ternary

## Abstract
PROteolysis TArgeting Chimeras (PROTACs) are heterobifunctional small molecules designed to recruit an E3 ubiquitin ligase to degrade some protein of interest (POI). 
  
We would like to develop a computational approach for modeling the ensemble of ternary complexes that can be formed by a given POI / PROTAC / E3 ligase. This method involves protein-protein docking to identify complementary binding modes, followed by screening of low-energy linker conformations to determine which candidate binding modes are accessible to a given PROTAC linker.
  
This "ternary_model_prediction.py" script can take protein-protein docking decoy(s) and linker comformer(s) as input, align linker conformer to the decoy, and output ternary structure(s) if the alignment rmsd is less than the cutoff value.

## Dependency
Rosetta Software Suite; OpenEye Software Suite

**The docking decoys should be genearted with two ligands along with the two proteins using Rosetta.

**The linker conformers should have overlap part at each end with the two ligands in the decoy using OMEGA.

## Input files
1) protein-protein docking decoy(s), can be eitehr a pdb file or a list of pdb structures (see docking_decoy.pdb as example);

2) linker conformers, can be either a pdb file or a list of pdb structures (see linker_conformer.pdb as example);

3) decoy_atom_list.txt & linker_atom_list.txt: atoms would be used to do the alignment (see decoy_atom_list.txt & linker_atom_list.txt as examples);

4) decoy_atom_delete.txt & linker_atom_delte.txt: atoms which are repeated in decoy / conformer and need to be deleted (see decoy_atom_delete.txt & linker_atom_delte.txt as examples).

## Example command
python2.7 ternary_model_prediction.py -linker_atoms_alignment linker_atom_list.txt -decoy_atoms_alignment decoy_atom_list.txt -decoy_file docking_decoy.pdb -linker_conformer_file linker_conformer.pdb -linker_atoms_delete linker_atom_delete.txt -warheads_atoms_delete decoy_atom_delete.txt -ternary_model default -rmsd_file rmsd.txt

**python2.7 ternary_model_prediction.py -help # More information of flags
