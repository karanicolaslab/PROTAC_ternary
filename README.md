# PROTAC_ternary

## Abstract
PROteolysis TArgeting Chimeras (PROTACs) are heterobifunctional small molecules designed to recruit an E3 ubiquitin ligase to degrade the protein of interest (POI). 
  
We would like to develop a computational approach for modeling the ensemble of ternary complexes that can be formed by a given POI / PROTAC / E3 ligase. This method involves protein-protein docking to identify complementary binding modes, followed by screening of low-energy linker conformations to determine which candidate binding modes are accessible to a given PROTAC linker.
  
This `ternary_model_prediction.py` script can take protein-protein docking decoy(s) (please see example in Example folder) and linker comformer(s) as input (please see example in Example folder), align linker conformer to the decoy with the atoms you choose, and output ternary structure(s) if the alignment rmsd is less than the cutoff value.

## Dependency
[Rosetta](https://www.rosettacommons.org/software/license-and-download) Software Suite; [OpenEye](https://www.eyesopen.com/) Software Suite; [RDKit](www.rdkit.org) Software Suite

* The docking decoys should be generated with two ligands along with the two proteins using Rosetta, example command:

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

## Input files for ternary model prediction
1) protein-protein docking decoy(s), can be either a one pdb file or multiple (see `docking_decoy.pdb` as example in Example folder);

2) linker conformers, can be either a one pdb file or multiple (see `linker_conformer.pdb` as example in Example folder);

3) `decoy_atom_list.txt` & `linker_atom_list.txt`: atoms would be used to do the alignment (see `decoy_atom_list.txt` & `linker_atom_list.txt` as examples in Example folder);

4) `decoy_atom_delete.txt` & `linker_atom_delte.txt`: atoms which are repeated in decoy / conformer and need to be deleted (see `decoy_atom_delete.txt` & `linker_atom_delte.txt` as examples in Example folder).

5) the flags `-d`/`-dl` and `-l`/`-ll` cannot not both be used in a single command for decoy or linker files; a single decoy can be modeling against a list of linker conformers (`-d`/`-ll`) or a list of decoys can be modeled against single linker conformer (`-dl`/`-l`) or single decoy and linker (`-d`/`-l`) or lists of both (`-dl`/`-ll`) can be modeled

## Flags information
```
-da/--decoy_aligment          # take the decoy_atom_list.txt
-la/--linker_aligment         # take the linker_atom_list.txt
-wd/--warheads_delete         # take the decoy_atom_delete.txt
-ld/--linkerd_delete          # take the linker_atom_delete.txt
-d/--decoy                    # take the decoy pdbs
-d/--decoy_list		      # take the list of decoy pdbs in decoy_list.txt
-l/--linker                   # take the linker pdbs
-ll/--linker_list             # take the list of linker conformations in linker_list.txt
-c/--cutoff                   # take a float as the cutoff of alignment rmsd, if not applied, the default value (0.4) will be used
-r/--rmsd                     # output file with alignemnet rmsd value(s), if not applied, the default name (rmsd.txt) will be used
-t/--ternary                  # output ternary structure, if not applied, no ternary structure will be generated, if applied, choose either default (the output would be ternary0.pdb, etc.) or specify (the output would be decoy_linker.pdb,, etc)
-ai/--alignment_iterations    # number of rounds for alignment (default is set to 5 iterations, shows high convergence in most cases)
```

For a visual representation of the linker/decoy atom and delete lists, see: 

![alt text](https://github.com/karanicolaslab/PROTAC_ternary/raw/master/Atom_list_Figure.png)
 
## Example command
```
$ python ternary_model_prediction.py -la linker_atom_list.txt \
                                     -da decoy_atom_list.txt \
                                     -dl decoy_list.xt \ 
                                     -ll linker_list.txt \
				     -d  decoy_model.pdb \
				     -l  linker_conf.pdb \ 
                                     -ld linker_atom_delete.txt \
                                     -wd decoy_atom_delete.txt \
                                     -t default \
                                     -r rmsd.txt
```

## More information on ternary model predition

```python ternary_model_prediction.py -h```

## Modify generated ternary models for Rosetta minimization
This section will ensure consistent atom labeling of ternary complexes for use of a single descriptive params file for minimization with Rosetta

The script requires:
1) The path to obabel be defined: 

```export BABEL=/path/to/openbabel/2.4.1/bin/obabel```

2) The path to Rosetta molfile_to_params.py be defined:

```export MOL2PARAMS=/path/to/Rosetta/main/source/scripts/python/public/molfile_to_params.py```

3) An input list of ternary models 

Example Command 

```python ternary_modify.py -s pdb_list.txt```

## Modified Ternary model outputs
Note: Only ternary PROTAC models with acceptable bond connectivity properties (angle, length) defined by Rosetta will pass this step

1) Successfully converted ternary PROTAC models will be appended with _mod.pdb
2) Ternary atom labels will be converted to TRN
3) Generate a single params file (from one of the top scoring ternary models, see RMSD output file) to be compatible with this new atom labeling:

```/path/to/Rosetta/main/source/scripts/python/public/molfile_to_params.py TRN.mol2 -n TRN```

The ternary PROTAC pdb must first be converted to a mol2 file with obabel
The flag -n TRN is critical during params file generation to maintain consistency with the output from ternary_modify.py

## Minimize the ternary models

Example Command

```
/path/to/Rosetta/main/source/bin/minimize_ppi.linuxgccrelease \
						-database /path/to/Rosetta/main/database \
						-s ternary_model_mod.pdb \
						-extra_res_fa TRN.params \
						-jump_all \
						-out:file:scorefile score.sc
```

## Calculating the FFC
First, you will need gather the interface scores from the ternary minimiztion output, run the below command:

```python ppi_ternary_scores.py```

This will ask you for your the name of your output file name (i.e., `score.sc`  - from above)

Next, you will need to minimize the skeleton (docked decoys prior to ternary model prediction) using:
```
/path/to/Rosetta/main/source/bin/minimize_ppi.linuxgccrelease \
						-database /path/to/Rosetta/main/database \
						-s ternary_model_mod.pdb \
						-extra_res_fa ligand1.params ligand2_params \
						-jump_all \
						-out:file:scorefile score.sc
```

The `-extra_res_fa ligand1.params ligand2_params` will be the same as the initial docking submission at the beginning of the pipeline


Now, calculate the median interface score the minimized skeletons using:

```python ppi_skeleton_median.py```

The output will be `skeleton_median.txt`, which contains a summary of minimization scores for the decoy skeletons

Finally, calculate the ffc using;
```python ffc_calculator.py```
This script requires the inputs of
  1) the csv output file from ppi_ternary_scores.py (ppi_ternary_scores.csv)
  2) number of docked decoys: `5,000` (benchmark from paper)
  3) number of ligand conformations: `1,000` (benchmark from paper)
  4) the median decoy skeleton interface energy (for example, `-61.71495`)
		
This script will print the ffc value, the number of ternary models that passed the energy filter, and txt file with the names/scores of each ternary model that passed this filter 




