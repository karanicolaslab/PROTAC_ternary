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
1) protein-protein docking decoy(s), can be eitehr a pdb file or a list of pdb structures (see docking_decoy.pdb and decoy_list.txt as example);

2) linker conformers, can be either a pdb file or a list of pdb structures (see linker_conformer.pdb and linker_list.txt as example);

3) decoy_atom_list.txt & linker_atom_list.txt: atoms would be used to do the alignment (see decoy_atom_list.txt & linker_atom_list.txt as examples);

4) decoy_atom_delete.txt & linker_atom_delte.txt: atoms which are repeated in decoy / conformer and need to be deleted (see decoy_atom_delete.txt & linker_atom_delte.txt as examples).

## Flags information
-da/--decoy_aligment

## Example command
python2.7 ternary_model_prediction.py -la linker_atom_list.txt -da decoy_atom_list.txt -d docking_decoy.pdb -l linker_conformer.pdb -ld linker_atom_delete.txt -wd decoy_atom_delete.txt -t default -r rmsd.txt

## More information
python2.7 ternary_model_prediction.py -h 
