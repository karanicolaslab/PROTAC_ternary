# PROTAC_ternary

## Abstract
PROteolysis TArgeting Chimeras (PROTACs) are heterobifunctional small molecules designed to recruit an E3 ubiquitin ligase to degrade the protein of interest (POI). 
  
We would like to develop a computational approach for modeling the ensemble of ternary complexes that can be formed by a given POI / PROTAC / E3 ligase. This method involves protein-protein docking to identify complementary binding modes, followed by screening of low-energy linker conformations to determine which candidate binding modes are accessible to a given PROTAC linker.
  
This `ternary_model_prediction.py` script can take protein-protein docking decoy(s) (please see example in Example folder) and linker comformer(s) as input (please see example in Example folder), align linker conformer to the decoy with the atoms you choose, and output ternary structure(s) if the alignment rmsd is less than the cutoff value.

## Dependency
[Rosetta](https://www.rosettacommons.org/software/license-and-download) Software Suite; [OpenEye](https://www.eyesopen.com/) Software Suite

* The docking decoys should be genearted with two ligands along with the two proteins using Rosetta, example command:

```
$ ./Rosetta/main/source/bin/docking_protocol.linuxgccelease –database path/to/Rosetta/main/database \
                                                            –s POI_ligand1_E3ligase_ligand2_prepacked.pdb \
                                                            –nstruct 50000 \
                                                            –s POI_ligand1_E3ligase_ligand2_prepacked.pdb \
                                                            –use_input_sc \
                                                            –spin \
                                                            –dock_pert 5 20 \
                                                            –partners XY_MN \
                                                            –ex1 \
                                                            –ex2aro \
                                                            –extra_res_fa ligand1.params ligand2_params \
                                                            –out:file:scorefile score.sc –score:docking_interface_score 1
```

(# X and Y are the chain IDs of POI and its ligand1 and M and N are the chain IDs of E3 ligase and its ligand2. The `–partners XY_MN` flag is used to make the ligands only move together with their paired proteins.)

* The linker conformers should have overlap part at each end with the two ligands in the decoy using OMEGA (`linker_conformer.pdb` as example)

## Input files
1) protein-protein docking decoy(s), can be either a one pdb file or multiple (see `docking_decoy.pdb` as example in Example folder);

2) linker conformers, can be either a one pdb file or multiple (see `linker_conformer.pdb` as example in Example folder);

3) `decoy_atom_list.txt` & `linker_atom_list.txt`: atoms would be used to do the alignment (see `decoy_atom_list.txt` & `linker_atom_list.txt` as examples in Example folder);

4) `decoy_atom_delete.txt` & `linker_atom_delte.txt`: atoms which are repeated in decoy / conformer and need to be deleted (see `decoy_atom_delete.txt` & `linker_atom_delte.txt` as examples in Example folder).

## Flags information
```
-da/--decoy_aligment          # take the decoy_atom_list.txt
-la/--linker_aligment         # take the linker_atom_list.txt
-wd/--warheads_delete         # take the decoy_atom_delete.txt
-ld/--linkerd_delete          # take the linker_atom_delete.txt
-d/--decoy                    # take the decoy pdbs
-l/--linker                   # take the linker pdbs
-c/--cutoff                   # take a float as the cutoff of alignment rmsd, if not applied, the default value (0.4) will be used
-r/--rmsd                     # output file with alignemnet rmsd value(s), if not applied, the default name (rmsd.txt) will be used
-t/--ternary                  # output ternary structure, if not applied, no ternary structure will be generated, if applied, choose either default (the output would be ternary0.pdb, etc.) or specify (the output would be decoy_linker.pdb,, etc)
-ai/--alignment_iterations    # number of rounds for alignment
```

## Example command
```
$ python ternary_model_prediction.py -la linker_atom_list.txt \
                                     -da decoy_atom_list.txt \
                                     -d docking_decoy.pdb \
                                     -l linker_conformer.pdb \
                                     -ld linker_atom_delete.txt \
                                     -wd decoy_atom_delete.txt \
                                     -t default \
                                     -r rmsd.txt
                                     -ai 50
```

## More information
`$ python ternary_model_prediction.py -h` 
